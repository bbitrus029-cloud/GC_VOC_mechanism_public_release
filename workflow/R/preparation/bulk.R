# Source-specific raw-array/count preparation; never runs on import.
# Reader supplies exact acquired inputs. Computational formulas retain accepted order.
bulk_preparation <- function(context) {
  s02_root <- context$root
  s02_path <- function(p) if(startsWith(p,'tmp/') || startsWith(p,'cache/')) context$allowed(p) else context$input(p)
  s02_raw <- function(p) context$input(paste0('data/raw/s02_bulk_ingress/',p))
  s02_qc <- function(p) context$output('qc',p)
  s02_out <- function(p) context$output('processed',p)
  s02_write <- context$write_tsv; s02_json <- context$write_json; s02_rds <- context$write_rds
  s02_firewall <- function(...) context$guard()
  s02_resource <- function(label,heavy=FALSE) context$resource_check(if(heavy)8 else 6,label)
s02_pair_check <- function(d) {
  s01_assert_unique_key(d,'sample_id');s01_assert_unique_key(d,c('patient_id','tissue'))
  s01_require(all(d$tissue %in% c('NAT','Tumor')),'invalid tissue')
  tab<-table(d$patient_id,factor(d$tissue,levels=c('NAT','Tumor')))
  s01_require(all(tab==1L),'incomplete paired identities')
  design<-model.matrix(~tissue+patient_id,data=transform(d,tissue=factor(tissue,levels=c('NAT','Tumor')),patient_id=factor(patient_id)))
  s01_require(qr(design)$rank==ncol(design),'paired design rank failure')
  invisible(design)
}

s02_soft <- function(path) {
  con<-gzfile(path,'rt');on.exit(close(con));res<-list();current<-NULL
  repeat {lines<-readLines(con,n=10000,warn=FALSE);if(!length(lines))break
    lines<-lines[startsWith(lines,'^SAMPLE = ') | grepl('^!Sample_(title|platform_id|characteristics_ch1|supplementary_file) = ',lines)]
    for(l in lines) {
      if(startsWith(l,'^SAMPLE = ')) {if(!is.null(current))res[[current$sample_id]]<-current;current<-list(sample_id=sub('^\\^SAMPLE = ','',l),title='',platform='',characteristics=character(),cel_url='')}
      if(is.null(current))next
      if(startsWith(l,'!Sample_title = ')) current$title<-sub('^!Sample_title = ','',l)
      if(startsWith(l,'!Sample_platform_id = ')) current$platform<-sub('^!Sample_platform_id = ','',l)
      if(startsWith(l,'!Sample_characteristics_ch1 = '))current$characteristics<-c(current$characteristics,sub('^!Sample_characteristics_ch1 = ','',l))
      if(startsWith(l,'!Sample_supplementary_file = ') && grepl('\\.[Cc][Ee][Ll]\\.gz$',l))current$cel_url<-sub('^ftp:','https:',sub('^!Sample_supplementary_file = ','',l))
    }
  }
  if(!is.null(current))res[[current$sample_id]]<-current
  do.call(rbind,lapply(res,function(x)data.frame(sample_id=x$sample_id,title=x$title,platform=x$platform,characteristics=paste(x$characteristics[grepl('^tissue:',x$characteristics,ignore.case=TRUE)],collapse=' | '),cel_url=x$cel_url)))
}

s02_annotations <- function() {
  s02_firewall('annotation')
  p<-s02_raw('frozen_metadata/GPL96.annot.gz');con<-gzfile(p,'rt');lines<-readLines(con,warn=FALSE);close(con)
  a<-match('!platform_table_begin',lines);b<-match('!platform_table_end',lines)
  s01_require(!is.na(a)&&!is.na(b)&&b>a,'GPL96 table delimiters')
  tab<-read.delim(text=paste(lines[(a+1):(b-1)],collapse='\n'),sep='\t',quote='',comment.char='',check.names=FALSE,colClasses='character')
  gpl<-data.frame(feature_id=tab$ID,source_symbol=tab[['Gene symbol']],entrez_id=tab[['Gene ID']])
  archive<-s02_raw('frozen_metadata/huex10sttranscriptcluster.db_8.8.0.tar.gz')
  f<-'huex10sttranscriptcluster.db/inst/extdata/huex10sttranscriptcluster.sqlite'
  dir<-s02_path('cache/s02_bulk_ingress/annotation');dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  entries<-utils::untar(archive,list=TRUE);s01_require(sum(entries==f)==1L,'Annotation archive identity')
  utils::untar(archive,files=f,exdir=dir)
  s01_require(!nzchar(Sys.readlink(file.path(dir,f))),'Annotation SQLite must be regular')
  db<-DBI::dbConnect(RSQLite::SQLite(),file.path(dir,f),flags=RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(db))
  tables<-DBI::dbListTables(db)
  s02_write(data.frame(table=tables),s02_qc('exon_annotation_tables.tsv'))
  ex<-DBI::dbReadTable(db,'probes')
  s01_require(all(c('probe_id','gene_id','is_multiple') %in% names(ex)),'exon annotation schema mismatch')
  ex<-data.frame(feature_id=as.character(ex$probe_id),source_symbol='',entrez_id=as.character(ex$gene_id))
  canonical<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype='ENTREZID'),keytype='ENTREZID',columns='SYMBOL')
  convert<-function(x,source) {
    x[is.na(x)]<-'';x<-unique(x)
    splitx<-split(seq_len(nrow(x)),x$feature_id)
    splitx<-splitx[order(names(splitx),method='radix')]
    ids<-lapply(splitx,function(ii){v<-unique(trimws(unlist(strsplit(x$entrez_id[ii],'///',fixed=TRUE))));v[nzchar(v)]})
    single<-lengths(ids)==1L
    idstring<-vapply(ids,paste,character(1),collapse=';')
    symbol<-canonical$SYMBOL[match(idstring,canonical$ENTREZID)]
    ambiguous_entrez<-canonical$ENTREZID[duplicated(canonical$ENTREZID)|duplicated(canonical$ENTREZID,fromLast=TRUE)]
    reliable<-single & !is.na(symbol) & nzchar(symbol) & !(idstring %in% ambiguous_entrez)
    reason<-ifelse(!single,'EXCLUDED_NONUNIQUE_OR_MISSING_ENTREZ',ifelse(!reliable,'EXCLUDED_NONUNIQUE_OR_MISSING_CANONICAL_SYMBOL','RELIABLE'))
    symbol[!reliable]<-'';symbol[symbol=='FSP1']<-'AIFM2'
    bridge<-unique(data.frame(symbol=symbol[reliable],entrez=idstring[reliable]))
    ambiguous_symbol<-bridge$symbol[duplicated(bridge$symbol)|duplicated(bridge$symbol,fromLast=TRUE)]
    collision<-symbol %in% ambiguous_symbol
    reason[collision]<-'EXCLUDED_CANONICAL_SYMBOL_MULTIPLE_ENTREZ';symbol[collision]<-''
    y<-data.frame(feature_id=names(splitx),source_symbol=vapply(splitx,function(ii)paste(sort(unique(x$source_symbol[ii])),collapse=';'),character(1)),entrez_id=idstring,canonical_gene=symbol,status=reason)
    rownames(y)<-NULL;s01_assert_unique_key(y,'feature_id');s02_write(y,s02_qc(paste0(source,'_mapping.tsv')))
  }
  c(convert(gpl,'GSE29272'),convert(ex,'GSE27342'))
}

s02_qc_flags <- function(metrics,domains,tails=rep('both',ncol(metrics)),expected_unavailable=character()) {
  s01_require(ncol(metrics)==length(domains)&&length(tails)==ncol(metrics),'QC metric contract mismatch')
  s01_require(!('median_NUSE' %in% expected_unavailable),'NUSE is an intended available metric')
  flags<-matrix(FALSE,nrow(metrics),ncol(metrics),dimnames=dimnames(metrics));thresholds<-list()
  for(j in seq_len(ncol(metrics))) {
    x<-metrics[,j];available<-is.finite(x)
    s01_require(all(available)||!any(available),'partially nonfinite QC metric requires diagnosis')
    s01_require(any(available)||colnames(metrics)[j] %in% expected_unavailable,'unexpected missing QC metric requires diagnosis')
    if(colnames(metrics)[j] %in% expected_unavailable)s01_require(all(is.na(x)),'declared unavailable metric must contain explicit NA only')
    med<-if(all(available))median(x)else NA_real_;mad<-if(all(available))stats::mad(x,constant=1)else NA_real_
    if(is.finite(mad)&&mad>0) flags[,j]<-if(tails[j]=='lower')x<med-3*mad else x<med-3*mad|x>med+3*mad
    thresholds[[j]]<-data.frame(metric=colnames(metrics)[j],domain=domains[j],median=med,unscaled_MAD=mad,lower=med-3*mad,upper=med+3*mad,tail=tails[j],status=if(all(available))'AVAILABLE'else 'NA_UNAVAILABLE')
  }
  dm<-sapply(unique(domains),function(d)rowSums(flags[,domains==d,drop=FALSE])>0)
  list(metric_flags=flags,domain_flags=dm,abnormal_domains=rowSums(dm),exclude=rowSums(dm)>=2,thresholds=do.call(rbind,thresholds))
}

s02_propagate <- function(d,excluded) {
  s02_pair_check(d);s01_require(length(excluded)==nrow(d)&&!anyNA(excluded),'exclusion vector mismatch')
  bad<-unique(d$patient_id[excluded]);d$sample_QC_excluded<-excluded;d$pair_excluded<-d$patient_id %in% bad
  d$exclusion_reason<-ifelse(excluded,'AT_LEAST_TWO_AVAILABLE_ABNORMAL_QC_DOMAINS',ifelse(d$pair_excluded,'PAIRED_PARTNER_QC_EXCLUDED','NONE'))
  if(any(!d$pair_excluded))s02_pair_check(d[!d$pair_excluded,]);d
}

s02_collapse <- function(x,mapping,prepare_alternative=TRUE) {
  s01_validate_identity_vector(rownames(x),'expression feature');s01_validate_identity_vector(colnames(x),'expression sample')
  s01_assert_unique_key(mapping,'feature_id')
  mi<-match(rownames(x),mapping$feature_id)
  audit<-data.frame(feature_id=rownames(x),canonical_gene=mapping$canonical_gene[mi],mapping_status=mapping$status[mi],finite_in_all_formal_samples=rowSums(!is.finite(x))==0)
  audit$canonical_gene[is.na(mi)]<-'';audit$mapping_status[is.na(mi)]<-'EXCLUDED_FEATURE_ABSENT_FROM_ANNOTATION'
  audit<-audit[order(audit$feature_id,method='radix'),]
  mapping<-mapping[mapping$status=='RELIABLE' & mapping$feature_id %in% rownames(x),]
  s01_require(all(!is.na(mapping$canonical_gene)&nzchar(mapping$canonical_gene)),'incomplete canonical mapping')
  mapping<-mapping[order(mapping$canonical_gene,mapping$feature_id,method='radix'),]
  x<-x[match(mapping$feature_id,rownames(x)),,drop=FALSE]
  finite<-rowSums(!is.finite(x))==0
  mapping$finite_in_all_formal_samples<-finite
  excluded<-mapping[!finite,,drop=FALSE];mapping<-mapping[finite,,drop=FALSE];x<-x[finite,,drop=FALSE]
  s01_require(nrow(x)>0,'no reliable finite features')
  idx<-split(seq_len(nrow(x)),mapping$canonical_gene);genes<-sort(names(idx),method='radix')
  primary<-t(vapply(genes,function(g)apply(x[idx[[g]],,drop=FALSE],2,median),numeric(ncol(x))))
  dimnames(primary)<-list(genes,colnames(x))
  alternative<-NULL;choices<-NULL
  # Means are used internally only; freeze choices without value-bearing displays.
  if(prepare_alternative) {
    means<-rowMeans(x);chosen<-vapply(genes,function(g){ii<-idx[[g]];ii[which.max(means[ii])]},integer(1))
    alternative<-x[chosen,,drop=FALSE];rownames(alternative)<-genes
    choices<-data.frame(canonical_gene=genes,feature_id=mapping$feature_id[chosen],tie_rule='exact_mean_tie_ASCII_radix_smallest_feature')
  }
  list(primary=primary,alternative=alternative,choices=choices,mapping=mapping,nonfinite_exclusions=excluded,feature_audit=audit)
}

s02_coverage <- function(genes,source,core_only=FALSE) {
  cfg<-yaml::read_yaml(s02_path('config/analysis_config.yml'))
  programs<-c(list(G_core=cfg$programs$G_core$members,H_core=cfg$programs$H_core$members),cfg$coverage_contracts$context_programs,cfg$coverage_contracts$G_submodules[c('substrate_assembly','oxidation_support')])
  if(core_only)programs<-programs[c('G_core','H_core')]
  do.call(rbind,lapply(names(programs),function(p) {
    required<-programs[[p]];present<-intersect(required,genes)
    # Empty overlap is a coverage failure, never substitute a member.
    pass<-identical(present,required)
    if(length(present))pass<-s01_coverage_decision(present,required,'exact',genes)$estimable
    data.frame(source=source,program=p,required_count=length(required),observed_count=length(present),required_members=paste(required,collapse=';'),observed_members=paste(present,collapse=';'),missing_members=paste(setdiff(required,genes),collapse=';'),state=if(pass)'PASS'else 'NA_NONESTIMABLE')
  }))
}

s02_array_headers <- function(d,source) {
  headers<-lapply(d$raw_path,function(p)affyio::read.celfile.header(p,info='full'))
  # Retain source header structure locally; display only verified structural fields.
  s02_rds(headers,s02_out(paste0(source,'_CEL_headers.rds')))
  rows<-do.call(rbind,lapply(seq_along(headers),function(i) {
    h<-headers[[i]];text<-paste(unlist(h),collapse=' ')
    date<-regmatches(text,regexpr('[0-9]{2}/[0-9]{2}/[0-9]{2}',text))
    s01_require(length(date)==1&&nzchar(date),'CEL scan date unresolved')
    expected<-if(source=='GSE29272')'HG-U133A'else 'HuEx-1_0-st-v2'
    s01_require(identical(h$cdfName,expected),'CEL platform mismatch')
    s01_require(identical(h[['CEL dimensions']],headers[[1]][['CEL dimensions']]),'CEL dimensions mismatch')
    scan_time<-sub('^[^ ]+ +','',h$ScanDate)
    scanner<-strsplit(sub(paste0('^.*',h$ScanDate,' +'),'',h$DatHeader),' +')[[1]][1]
    data.frame(sample_id=d$sample_id[i],neutral_id=d$neutral_id[i],scan_date=date,scan_time=scan_time,scanner=scanner,columns=h[['CEL dimensions']][1],rows=h[['CEL dimensions']][2],platform=expected,header_status='PASS')
  }))
  if(source=='GSE27342') {
    frozen<-s01_read_tsv(s02_raw('frozen_metadata/gse27342_scan_header_audit.tsv'))
    s01_validate_identity_vector(rows$sample_id,'CEL header identity',frozen$gsm)
    s01_require(all(rows$scan_date==frozen$scan_date[match(rows$sample_id,frozen$gsm)]),'frozen scan header mismatch')
    s01_require(all(rows$scan_time==frozen$scan_time[match(rows$sample_id,frozen$gsm)])&&all(rows$scanner==frozen$scanner[match(rows$sample_id,frozen$gsm)]),'frozen scanner/time mismatch')
    frozen_bytes<-as.numeric(sub('^.*/','',frozen$content_range))
    s01_require(all(d$bytes==frozen_bytes[match(d$sample_id,frozen$gsm)]),'frozen CEL compressed byte count mismatch')
  }
  rows
}

s02_array <- function(source,inventory,annotation) {
  s02_firewall('rma');before<-s02_resource(paste0(source,'_start'),TRUE)
  settings<-s02_path('config/implementation/s02_processing_settings.yml');s01_require(file.exists(settings),'settings absent')
  d<-s01_read_tsv(inventory);s02_pair_check(d)
  d$raw_path<-context$input(d$raw_path)
  for(i in seq_len(nrow(d)))s01_validate_fixed_authority_file(d$raw_path[i],d$sha256[i],'raw CEL')
  h<-s02_array_headers(d,source);s02_write(h,s02_qc(paste0(source,'_header_audit.tsv')))
  d$scan_date<-h$scan_date
  if(source=='GSE27342') {
    dates<-split(d$scan_date,d$patient_id);same<-names(dates)[vapply(dates,function(x)length(unique(x))==1,logical(1))]
    s01_require(length(same)==64&&length(dates)==80,'scan populations mismatch');d$core_population<-d$patient_id %in% same
  }else d$core_population<-TRUE
  z<-s02_prepared_metrics(source,inventory);x<-z$expression;metrics<-z$metrics
  domains<-c('distribution','distribution','RLE_NUSE','RLE_NUSE','RLE_NUSE','RNA_degradation','sample_correlation')
  qc<-s02_qc_flags(metrics,domains,c(rep('both',6),'lower'),if(source=='GSE27342')'RNA_degradation_slope'else character())
  s02_write(cbind(data.frame(neutral_id=d$neutral_id),as.data.frame(metrics),as.data.frame(qc$domain_flags),abnormal_domains=qc$abnormal_domains),s02_qc(paste0(source,'_technical_QC.tsv')))
  s02_write(qc$thresholds,s02_qc(paste0(source,'_QC_thresholds.tsv')))
  d<-s02_propagate(d,qc$exclude);s02_write(d,s02_qc(paste0(source,'_pair_ledger.tsv')))
  # Diagnostic PCA uses neutral identifiers and order only.
  pc<-prcomp(t(x),center=TRUE,scale.=FALSE,rank.=2)
  s02_write(data.frame(neutral_id=d$neutral_id,PC1=pc$x[,1],PC2=pc$x[,2]),s02_qc(paste0(source,'_neutral_PCA.tsv')))
  pdf(s02_qc(paste0(source,'_neutral_QC.pdf')),width=8,height=6)
  plot(pc$x[,1:2],pch=16,col='grey40',main=paste(source,'technical PCA'),xlab='PC1',ylab='PC2');text(pc$x[,1:2],labels=d$neutral_id,cex=.35,pos=3)
  dev.off()
  map<-s01_read_tsv(annotation[grepl(source,annotation)])
  core<-!d$pair_excluded & d$core_population;s01_require(any(core),'no core pairs retained')
  collapsed<-s02_collapse(x[,core,drop=FALSE],map)
  cov<-s02_coverage(rownames(collapsed$primary),source)
  files<-c(s02_rds(collapsed$primary,s02_out(paste0(source,'_primary.rds'))),s02_rds(collapsed$alternative,s02_out(paste0(source,'_alternative.rds'))),s02_write(collapsed$choices,s02_qc(paste0(source,'_alternative_probe_choices.tsv'))),s02_write(collapsed$mapping,s02_qc(paste0(source,'_reliable_mapping.tsv'))),s02_write(collapsed$nonfinite_exclusions,s02_qc(paste0(source,'_nonfinite_features.tsv'))),s02_write(data.frame(gene=rownames(collapsed$primary)),s02_qc(paste0(source,'_background.tsv'))))
  files<-c(files,s02_write(collapsed$feature_audit,s02_qc(paste0(source,'_processed_feature_audit.tsv'))))
  if(source=='GSE27342') {
    allpop<-s02_collapse(x[,!d$pair_excluded,drop=FALSE],map,prepare_alternative=FALSE)$primary
    files<-c(files,s02_rds(allpop,s02_out('GSE27342_80_primary.rds')),s02_write(data.frame(gene=rownames(allpop)),s02_qc('GSE27342_80_background.tsv')))
    cov<-rbind(cov,s02_coverage(rownames(allpop),'GSE27342_80',TRUE))
  }
  files<-c(files,s02_write(cov,s02_qc(paste0(source,'_coverage.tsv'))))
  eligible<-all(cov$state[cov$program %in% c('G_core','H_core')]=='PASS')
  end<-s02_resource(paste0(source,'_end'))
  status<-list(source=source,expected_pairs=nrow(d)/2,retained_core_pairs=sum(core)/2,retained_all_pairs=sum(!d$pair_excluded)/2,QC_excluded_samples=sum(d$sample_QC_excluded),QC_excluded_pairs=sum(d$pair_excluded)/2,core_eligible=eligible,coverage=cov,settings_hash=s01_sha256(settings),artifacts=setNames(lapply(files,s01_sha256),files),resource_stop_before_next_heavy=end$available_GiB<6||(end$swap_GiB>2&&end$swap_GiB>before$swap_GiB))
  c(files,s02_json(status,s02_qc(paste0(source,'_processing_state.json'))))
}

s02_parse_star <- function(path) {
  cols<-c('gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')
  hdr<-readLines(path,n=2,warn=FALSE)
  s01_require(grepl('GENCODE v36',hdr[1],fixed=TRUE),'GENCODE v36 header mismatch')
  t<-read.delim(path,skip=1,sep='\t',quote='',comment.char='',check.names=FALSE,colClasses='character')
  s01_require(identical(names(t),cols),'STAR schema mismatch')
  s01_require(identical(t$gene_id[1:4],c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous')),'STAR summary row mismatch')
  t<-t[-(1:4),];s01_validate_identity_vector(t$gene_id,'STAR gene')
  s01_require(all(grepl('^ENSG[0-9]+[.][0-9]+(_PAR_Y)?$',t$gene_id)),'STAR versioned Ensembl identity mismatch')
  s01_require(all(grepl('^[0-9]+$',t$unstranded)),'invalid unstranded integer token')
  count<-as.numeric(t$unstranded);s01_require(all(is.finite(count)&count>=0&count==floor(count)),'invalid STAR count')
  list(annotation=t[,1:3],counts=count)
}

s02_tcga_prepare <- function(ledger) {
  s02_firewall('counts_filter_TMM');s02_resource('TCGA_prepare_start',TRUE)
  d<-s01_read_tsv(ledger);s02_pair_check(d);s01_require(nrow(d)==66,'TCGA count universe mismatch')
  cols<-c('gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded')
  matrices<-vector('list',nrow(d));annotation<-NULL
  for(i in seq_len(nrow(d))) {
    s01_validate_fixed_authority_file(d$raw_path[i],d$sha256[i],'TCGA raw file')
    parsed<-s02_parse_star(d$raw_path[i]);ann<-parsed$annotation
    if(is.null(annotation))annotation<-ann else s01_require(identical(ann,annotation),'STAR gene order or annotation mismatch')
    matrices[[i]]<-parsed$counts
  }
  x<-do.call(cbind,matrices);dimnames(x)<-list(annotation$gene_id,d$sample_id)
  s01_require(all(colSums(x)>0),'zero count library')
  detected<-rowSums(edgeR::cpm(x,log=FALSE)>=1)>=ceiling(.2*66)
  symbol<-annotation$gene_name;symbol[symbol=='FSP1']<-'AIFM2'
  # Ambiguous gene symbols are excluded as mappings, never collapsed as counts.
  unique_symbol<-nzchar(symbol)&!is.na(symbol)&!duplicated(symbol)&!duplicated(symbol,fromLast=TRUE)
  annotation$canonical_gene<-symbol;annotation$CPM_filter_pass<-detected;annotation$unique_symbol<-unique_symbol
  annotation$status<-ifelse(!unique_symbol,'EXCLUDED_AMBIGUOUS_OR_MISSING_SYMBOL',ifelse(!detected,'EXCLUDED_FROZEN_CPM_RULE','RELIABLE'))
  keep<-detected&unique_symbol
  counts<-x[keep,,drop=FALSE];rownames(counts)<-symbol[keep];counts<-counts[order(rownames(counts),method='radix'),,drop=FALSE]
  dge<-edgeR::normLibSizes(edgeR::DGEList(counts=counts),method='TMM')
  design<-s02_pair_check(d);rownames(design)<-colnames(counts)
  s01_require(all(is.finite(dge$samples$norm.factors)&dge$samples$norm.factors>0),'invalid TMM')
  coverage<-s02_coverage(rownames(counts),'TCGA_STAD')
  files<-c(s02_rds(list(dge=dge,design=design,sample_data=d,voom_state='APPROVED_DEFERRED_TO_S05_FIRST_EFFECT_TRANSACTION'),s02_out('TCGA_TMM_preparation.rds')),s02_write(annotation,s02_qc('TCGA_mapping_filter.tsv')),s02_write(data.frame(gene=rownames(counts)),s02_qc('TCGA_background.tsv')),s02_write(coverage,s02_qc('TCGA_coverage.tsv')))
  s02_resource('TCGA_prepare_end')
  c(files,s02_json(list(pairs=33,files=66,genes=nrow(counts),design_columns=ncol(design),design_rank=qr(design)$rank,filter='CPM >= 1 in >= 14/66 files',normalization='TMM',voom_completed=FALSE,voom_state='APPROVED_DEFERRED_TO_S05_FIRST_EFFECT_TRANSACTION',artifacts=setNames(lapply(files,s01_sha256),files)),s02_qc('TCGA_processing_state.json')))
}

s02_decision <- function(first,replication,backup) {
  a<-jsonlite::read_json(s02_qc('GSE29272_processing_state.json'));b<-jsonlite::read_json(s02_qc('GSE27342_processing_state.json'))
  branch<-if(!a$core_eligible&&!b$core_eligible)list(branch='BOTH_CORE_FAILED_OWNER_REVIEW',first=NA_character_,replication=NA_character_)else s01_resolve_bulk_roles(c(GSE29272=a$core_eligible,GSE27342=b$core_eligible),TRUE)
  ga<-s01_read_tsv(s02_qc('GSE29272_background.tsv'))$gene;gb<-s01_read_tsv(s02_qc('GSE27342_background.tsv'))$gene
  common<-sort(intersect(ga,gb),method='radix');p<-s02_write(data.frame(gene=common),s02_qc('common_microarray_background.tsv'))
  cov<-s02_coverage(common,'COMMON_MICROARRAY_BACKGROUND')
  c(p,s02_write(cov,s02_qc('common_background_coverage.tsv')),s02_json(c(branch,list(backup_activated=FALSE,activation_authorization=if(branch$branch=='normal')'NOT_REQUIRED_STAYS_INACTIVE'else 'REQUIRED_PENDING',core_eligibility=c(GSE29272=a$core_eligible,GSE27342=b$core_eligible),common_background_genes=length(common),common_background_sha256=s01_sha256(p))),s02_qc('backup_decision.json')))
}

s02_nuse_rows <- function(se,probe_counts) {
  s01_require(is.matrix(se)&&is.numeric(se)&&ncol(se)>1,'NUSE SE matrix invalid')
  s01_validate_identity_vector(rownames(se),'NUSE feature');s01_validate_identity_vector(colnames(se),'NUSE sample')
  s01_require(identical(names(probe_counts),rownames(se))&&all(probe_counts>=1),'NUSE probe identity mismatch')
  finite<-rowSums(is.finite(se));positive<-rowSums(is.finite(se)&se>0)
  denominator<-apply(se,1,median)
  reason<-ifelse(probe_counts==1,'STRUCTURAL_SINGLETON_NO_RESIDUAL_DF',
    ifelse(finite!=ncol(se),'UNEXPECTED_NONFINITE_SE',ifelse(!is.finite(denominator)|denominator<=0,'INVALID_DENOMINATOR',ifelse(positive!=ncol(se),'UNEXPECTED_NONPOSITIVE_SE','CONTRIBUTING'))))
  data.frame(feature_id=rownames(se),probe_count=as.integer(probe_counts),residual_degrees_of_freedom=(probe_counts-1)*(ncol(se)-1),finite_SE_samples=finite,positive_SE_samples=positive,denominator_finite=is.finite(denominator),denominator_positive=is.finite(denominator)&denominator>0,reason=reason,row.names=NULL)
}

s02_nuse_summary <- function(se,probe_counts) {
  rows<-s02_nuse_rows(se,probe_counts)
  s01_require(!any(!rows$reason %in% c('CONTRIBUTING','STRUCTURAL_SINGLETON_NO_RESIDUAL_DF')),'UNRESOLVED_NUSE_ROWS: unexpected SE or denominator failure')
  use<-rows$reason=='CONTRIBUTING';s01_require(any(use),'ALL_INVALID_NUSE_ROWS')
  normalized<-sweep(se[use,,drop=FALSE],1,apply(se[use,,drop=FALSE],1,median),'/')
  s01_require(all(is.finite(normalized)&normalized>0),'nonfinite normalized SE')
  values<-apply(normalized,2,median)
  s01_require(all(is.finite(values)&values>0),'unexpected missing sample NUSE')
  list(values=values,rows=rows,contributing_features=rownames(se)[use])
}

s02_fit_technical_plm <- function(d) {
  library(oligo);library(pd.huex.1.0.st.v2);library(ff);foreach::registerDoSEQ()
  ffdir<-s02_path('tmp/s02_bulk_ingress/ff');dir.create(ffdir,recursive=TRUE,showWarnings=FALSE)
  s01_require(length(list.files(ffdir,all.files=TRUE,no..=TRUE))==0,'existing ff files require ownership diagnosis')
  options(fftempdir=ffdir);oligoClasses::ldSetOptions(nsamples=10,nprobesets=1000,path=ffdir,verbose=TRUE)
  s01_require(foreach::getDoParWorkers()==1,'parallel backend prohibited')
  raw<-oligo::read.celfiles(d$raw_path,pkgname='pd.huex.1.0.st.v2',sampleNames=d$neutral_id)
  s01_require(inherits(Biobase::exprs(raw),'ff_matrix')&&startsWith(ff::filename(Biobase::exprs(raw)),paste0(ffdir,'/')),'ff ownership/path failure')
  s01_require(oligoClasses::ocSamples()==10&&oligoClasses::ocProbesets()==1000&&oligoClasses::ldPath()==ffdir,'oligo options reset')
  pmi<-get('stArrayPmInfo',asNamespace('oligo'))(raw,target='core')
  pms<-get('ffSubset',asNamespace('oligo'))(rows=pmi$fid,object=Biobase::exprs(raw),prefix='s02-repair-plm-pm-')
  for(j in seq_len(ncol(pms)))s01_require(all(is.finite(pms[,j])&pms[,j]>0),'invalid core PM')
  fit<-oligo::basicPLM(pms,as.character(pmi$fsetid),background=TRUE,normalize=TRUE,method='plm')
  se<-fit$StdErrors[,];dimnames(se)<-dimnames(fit$StdErrors)
  s01_require(identical(colnames(se),d$neutral_id),'PLM sample order mismatch')
  counts<-table(as.character(pmi$fsetid));counts<-setNames(as.integer(counts[rownames(se)]),rownames(se))
  s01_require(sum(counts)==nrow(pmi)&&!anyNA(counts)&&sum(counts==1)==174,'frozen core probe cardinality mismatch')
  rows<-s02_nuse_rows(se,counts);bad<-rows$feature_id[rows$reason!='CONTRIBUTING']
  # Compact row-level residual diagnostics identify actual failures without any tissue split.
  residual_audit<-do.call(rbind,lapply(bad,function(g){ii<-which(as.character(pmi$fsetid)==g);r<-fit$Residuals[ii,,drop=FALSE]
    data.frame(feature_id=g,residual_values=length(r),finite_residual_values=sum(is.finite(r)),nonzero_finite_residual_values=sum(is.finite(r)&r!=0))}))
  unexpected<-rows$feature_id[!rows$reason %in% c('CONTRIBUTING','STRUCTURAL_SINGLETON_NO_RESIDUAL_DF')]
  if(length(unexpected)){
    ii<-which(as.character(pmi$fsetid) %in% unexpected)
    s02_rds(list(probe_ids=pmi[ii,],raw_PM=pms[ii,,drop=FALSE],residuals=fit$Residuals[ii,,drop=FALSE]),s02_out('GSE27342_unresolved_PLM_rows.rds'))
  }
  result<-list(SE=se,probe_counts=counts,residual_audit=residual_audit,core_probe_identity_sha256=digest::digest(pmi[,c('fid','fsetid')],algo='sha256'))
  rm(fit,pms,raw);gc();result
}

s02_prepared_metrics <- function(source,inventory) {
 d<-s01_read_tsv(inventory);d$raw_path<-context$input(d$raw_path)
 stopifnot(source %in% c('GSE29272','GSE27342'))
    if(source=='GSE29272') {
      raw<-affy::ReadAffy(filenames=d$raw_path,cdfname='hgu133acdf',sampleNames=d$neutral_id)
      pm<-affy::pm(raw);s01_require(all(is.finite(pm))&&all(pm>0),'invalid raw PM intensity')
      distribution<-t(apply(pm,2,function(v)c(median_log2=median(log2(v)),IQR_log2=IQR(log2(v)))))
      rm(pm);gc()
      deg<-affy::AffyRNAdeg(raw,log.it=TRUE)$slope
      x<-Biobase::exprs(affy::rma(raw))
      plm<-affyPLM::fitPLM(raw)
      nuse<-apply(affyPLM::NUSE(plm,type='values'),2,median)
      rm(plm,raw);gc()
    } else {
      # Attach platform dependencies first: oligoClasses .onAttach and its ff
      # attach hook reset ldSetOptions() to defaults and the working directory.
      library(oligo);library(pd.huex.1.0.st.v2);library(ff);foreach::registerDoSEQ()
      ffdir<-s02_path('tmp/s02_bulk_ingress/ff');dir.create(ffdir,recursive=TRUE,showWarnings=FALSE)
      options(fftempdir=ffdir);s01_require(foreach::getDoParWorkers()==1,'parallel backend prohibited')
      oligoClasses::ldSetOptions(nsamples=10,nprobesets=1000,path=ffdir,verbose=TRUE)
      raw<-oligo::read.celfiles(filenames=d$raw_path,pkgname='pd.huex.1.0.st.v2',sampleNames=d$neutral_id)
      s01_require(identical(oligoClasses::ldPath(),ffdir)&&oligoClasses::ocSamples()==10&&oligoClasses::ocProbesets()==1000,'oligo resource options reset')
      s01_require(inherits(Biobase::exprs(raw),'ff_matrix'),'oligo raw object is not file backed')
      s01_require(startsWith(ff::filename(Biobase::exprs(raw)),paste0(ffdir,'/')),'raw ff file escaped temporary directory')
      info<-oligo::getProbeInfo(raw,target='core',field=c('fid','fsetid'))
      idx<-sort(unique(info$fid));distribution<-t(vapply(seq_len(nrow(d)),function(j){v<-Biobase::exprs(raw)[idx,j];s01_require(all(is.finite(v))&&all(v>0),'invalid exon PM');c(median_log2=median(log2(v)),IQR_log2=IQR(log2(v)))},numeric(2)))
      eset<-oligo::rma(raw,target='core',background=TRUE,normalize=TRUE)
      e<-Biobase::exprs(eset);x<-e[,];dimnames(x)<-dimnames(e)
      # Native basicPLM supports ff directly, avoiding fitProbeLevelModel's
      # full-platform in-memory weights/residual arrays. Same core probes.
      pmi<-get('stArrayPmInfo',asNamespace('oligo'))(raw,target='core')
      pms<-get('ffSubset',asNamespace('oligo'))(rows=pmi$fid,object=Biobase::exprs(raw),prefix='s02-plm-pm-')
      fit<-oligo::basicPLM(pms,as.character(pmi$fsetid),background=TRUE,normalize=TRUE,method='plm')
      se<-fit$StdErrors[,];dimnames(se)<-dimnames(fit$StdErrors)
      counts<-table(as.character(pmi$fsetid));counts<-setNames(as.integer(counts[rownames(se)]),rownames(se))
      nuse<-s02_nuse_summary(se,counts)$values
      deg<-rep(NA_real_,nrow(d))
      rm(fit,pms,se);gc()
      rm(eset,e,raw);gc()
    }
    s01_require(ncol(x)==nrow(d)&&all(is.finite(x)),'RMA finite/dimension failure')
    s01_require(identical(colnames(x),d$neutral_id),'RMA sample order failure')
    rle<-sweep(x,1,apply(x,1,median),'-')
    cc<-cor(x,method='pearson');diag(cc)<-NA_real_
    metrics<-cbind(distribution,absolute_median_RLE=abs(apply(rle,2,median)),RLE_IQR=apply(rle,2,IQR),median_NUSE=nuse,RNA_degradation_slope=deg,median_correlation=apply(cc,2,median,na.rm=TRUE))
 list(expression=x,metrics=metrics)
}
list(s02_soft=s02_soft,s02_annotations=s02_annotations,s02_qc_flags=s02_qc_flags,s02_propagate=s02_propagate,s02_collapse=s02_collapse,s02_coverage=s02_coverage,s02_array_headers=s02_array_headers,s02_prepared_metrics=s02_prepared_metrics,s02_array=s02_array,s02_nuse_rows=s02_nuse_rows,s02_nuse_summary=s02_nuse_summary,s02_fit_technical_plm=s02_fit_technical_plm,s02_parse_star=s02_parse_star,s02_tcga_prepare=s02_tcga_prepare,s02_decision=s02_decision)
}
