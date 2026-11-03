# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Public adaptation: explicit inputs and outputs; original scientific rules retained.

original_bulk_module <- function(context) {
  historical_namespace <- 's05_bulk_effects'
  stopifnot(is.list(context), is.function(context$input), is.function(context$output))
  resolve_path <- function(path) {
    if (startsWith(path, context$output_root)) return(path)
    prefix <- paste0('data/processed/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('processed',substring(path,nchar(prefix)+1L)))
    prefix <- paste0('results/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('results',substring(path,nchar(prefix)+1L)))
    prefix <- paste0('reports/qc/', historical_namespace, '/')
    if (startsWith(path,prefix)) return(context$output('qc',substring(path,nchar(prefix)+1L)))
    context$input(path)
  }
  read_rds <- function(path,...) base::readRDS(context$input(resolve_path(path)),...)
  read_yaml <- function(path) yaml::read_yaml(context$input(resolve_path(path)))
  read_lines <- function(path,...) readLines(context$input(resolve_path(path)),...)
  write_rds <- function(x,path) {context$write_rds(x,resolve_path(path));path}
  utc <- function() format(Sys.time(),'%Y-%m-%dT%H:%M:%OS6Z',tz='UTC')
  s05_path <- resolve_path
  s05_read <- function(path) utils::read.delim(context$input(resolve_path(path)),quote='',comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
  s05_json_read <- function(path) jsonlite::fromJSON(context$input(resolve_path(path)),simplifyVector=TRUE)
  s05_utc <- utc
  s05_entry_path <- function() 'data/manifests/s05_bulk_effects/ENTRY_RELEASE.json'
  s05_guard <- function() context$guard()
  s05_rds <- write_rds
  implementation_paths <- c('workflow/R/common/contracts.R','workflow/R/common/statistics.R','workflow/R/analysis/original_bulk.R')
  s05_hash_rows <- function(paths) data.frame(path=paths,bytes=as.numeric(file.info(vapply(paths,resolve_path,''))$size),sha256=vapply(paths,function(p)s01_sha256(resolve_path(p)),''),row.names=NULL)
  s05_verify_rows <- function(rows) {
    s01_require(!anyNA(rows$path)&&!anyDuplicated(rows$path),'Invalid frozen input binding')
    actual<-s05_hash_rows(rows$path)
    s01_require(identical(as.numeric(rows$bytes),actual$bytes)&&identical(as.character(rows$sha256),actual$sha256),'Frozen input or current implementation binding mismatch')
    invisible(TRUE)
  }
s05_pairs <- function(metadata, sample_ids) {
 s01_require(identical(as.character(metadata$sample_id),sample_ids),'sample order mismatch')
 s01_assert_unique_key(metadata,c('patient_id','tissue'),'paired sample')
 s01_require(all(metadata$tissue%in%c('Tumor','NAT')),'invalid tissue')
 id<-sort(unique(metadata$patient_id),method='radix')
 n<-match(paste(id,'NAT'),paste(metadata$patient_id,metadata$tissue))
 t<-match(paste(id,'Tumor'),paste(metadata$patient_id,metadata$tissue))
 s01_require(!anyNA(c(n,t))&&length(c(n,t))==nrow(metadata),'incomplete pairing')
 data.frame(patient_id=id,NAT_index=n,Tumor_index=t)
}

s05_bulk <- function() {e<-s05_read('data/manifests/s04_pre_effect/global_hypothesis_eligibility.tsv');e[e$contract=='C1',]}

s05_members <- function(construct,config) {
 if(construct%in%c('G_core','H_core'))return(config$programs[[construct]]$members)
 config$coverage_contracts$context_programs[[construct]]
}

s05_registry <- function(bulk=s05_bulk(),config=yaml::read_yaml(s05_path('config/analysis_config.yml'))) {
 rows<-list();add<-function(e,kind,spec='primary_rank',member='NA_NOT_APPLICABLE',sub='NA_NOT_APPLICABLE'){
  rows[[length(rows)+1L]]<<-data.frame(endpoint_id=e$endpoint_id,source=e$source,construct=e$construct,
   category=kind,specification=spec,deleted_member=member,submodule=sub,
   diagnostic_id=paste(kind,spec,member,sub,sep='__'),stringsAsFactors=FALSE)
 }
 for(i in seq_len(nrow(bulk))) {
  e<-bulk[i,]; if(e$primary_method!='paired_delta_t_test_Student_t_CI_plus_g_av'||e$source=='GSE27342_80')next
  for(k in c('paired_Wilcoxon','paired_percentile_bootstrap','influence'))add(e,k)
  add(e,'standardized_member_mean_direction','standardized_member_mean')
  if(e$source%in%c('GSE29272','GSE27342_64')) {
   add(e,'common_background','common_background');add(e,'alternative_probe_collapse','alternative_probe_collapse')
  }
  if(e$construct%in%c('G_core','H_core')) {
   add(e,'member_direction_concordance')
   for(g in s05_members(e$construct,config)) add(e,'leave_one_gene_out','primary_rank',g)
   if(e$source=='GSE29272')for(site in c('cardia','non_cardia'))add(e,'site_stratum_estimate',site)
   if(e$construct=='G_core')for(sub in c('substrate_assembly','oxidation_support'))
    for(k in c('G_submodule','paired_Wilcoxon','paired_percentile_bootstrap','influence'))add(e,k,'primary_rank',sub=sub)
  }
 }
 out<-do.call(rbind,rows);rownames(out)<-NULL
 s01_assert_unique_key(out,c('endpoint_id','diagnostic_id'),'diagnostic registry')
 s01_require(!any(out$source=='GSE27342_80'),'80 population extensions prohibited')
 out
}

s05_array_metadata <- function(d,source) {
 keep<-!d$pair_excluded
 if(source!='GSE27342_80')keep<-keep&d$core_population
 fields<-c('neutral_id','patient_id','tissue',if(source=='GSE29272')'site')
 md<-d[keep,fields];names(md)[1]<-'sample_id';md
}

s05_input <- function(source) {
 m<-context$bindings() # Current reader-supplied identities, not the historical pre-effect manifest
 if(source=='TCGA_STAD') {
  paths<-c('data/processed/s02_bulk_ingress/TCGA_TMM_preparation.rds')
  paths<-unique(paths);s05_verify_rows(m[match(paths,m$path),])
  tc<-read_rds(s05_path(paths[1])); md<-tc$sample_data
  p<-s01_prepare_paired_model(tc$dge$counts,md,'patient_id','patient',TRUE)
  s01_require(isTRUE(all.equal(unname(p$design),unname(tc$design),check.attributes=FALSE)),'TCGA frozen paired design mismatch')
  s01_require(identical(rownames(tc$design),md$sample_id)&&qr(tc$design)$rank==34L,'TCGA design identity/rank mismatch')
  s01_require(all(colSums(tc$dge$counts)==tc$dge$samples$lib.size)&&all(is.finite(tc$dge$samples$norm.factors))&&all(tc$dge$samples$norm.factors>0),'TCGA normalization integrity mismatch')
  mapping<-s05_read('reports/qc/s02_bulk_ingress/TCGA_mapping_filter.tsv')
  bg<-s05_read('reports/qc/s02_bulk_ingress/TCGA_background.tsv')$gene
  s01_require(identical(bg,rownames(tc$dge$counts))&&identical(bg,sort(mapping$canonical_gene[mapping$status=='RELIABLE'],method='radix')),'TCGA mapping/filter/background mismatch')
  contract<-s05_json_read('data/manifests/s04_pre_effect/D01_TCGA_execution_contract.json')
  s01_require(digest::digest(limma::voom,algo='sha256')==contract$implementation_sha256&&packageDescription('limma')$Version==contract$package,'voom implementation changed')
  s01_require(identical(lapply(formals(limma::voom),function(v)paste(deparse(v),collapse=' ')),contract$remaining_formals),'voom defaults changed')
  return(list(tc=tc,metadata=md,pairs=s05_pairs(md,colnames(tc$dge$counts)),background=bg))
 }
 name<-if(source=='GSE29272')source else 'GSE27342'
 stems<-if(source=='GSE27342_80')'GSE27342_80_primary' else paste0(name,c('_primary','_alternative'))
 paths<-paste0('data/processed/s02_bulk_ingress/',stems,'.rds')
 ledger<-paste0('reports/qc/s02_bulk_ingress/',name,'_pair_ledger.tsv')
 bgpath<-paste0('reports/qc/s02_bulk_ingress/',if(source=='GSE27342_80')source else name,'_background.tsv')
 s05_verify_rows(m[match(c(paths,ledger,bgpath),m$path),])
 d<-s05_read(ledger);md<-s05_array_metadata(d,source)
 x<-read_rds(s05_path(paths[1]));bg<-s05_read(bgpath)$gene
 s01_validate_numeric_matrix(x,'array expression');s01_require(identical(rownames(x),bg),'array background order mismatch')
 pairs<-s05_pairs(md,colnames(x)); alt<-NULL
 if(length(paths)>1L) {alt<-read_rds(s05_path(paths[2]));s01_require(identical(dimnames(alt),dimnames(x))&&all(is.finite(alt)),'alternative input identity mismatch')}
 list(expression=x,alternative=alt,metadata=md,pairs=pairs,background=bg)
}

s05_voom <- function(tc) {
 # Continue directly from the frozen DGEList; never recalculate filtering or TMM.
 p<-s01_prepare_paired_model(tc$dge$counts,tc$sample_data,'patient_id','patient',TRUE)
 p$design<-tc$design
 v<-limma::voom(tc$dge,design=tc$design,plot=FALSE,normalize.method='none')
 s01_require(identical(dimnames(v$E),dimnames(tc$dge$counts))&&all(is.finite(v$E))&&all(is.finite(v$weights))&&all(v$weights>0),'invalid voom output')
 fit<-limma::eBayes(limma::lmFit(v,tc$design))
 list(voom=v,fit=fit,member=s01_extract_limma_contract(fit,p,'expression~tissue+patient','frozen_TMM_then_voom'))
}

s05_score <- function(input,members,spec='primary_rank') {
 if(spec=='standardized_member_mean')return(s01_standardized_member_mean(input$expression,members))
 x<-if(spec=='alternative_probe_collapse')input$alternative else input$expression
 bg<-if(spec=='common_background')s05_read('reports/qc/s02_bulk_ingress/common_microarray_background.tsv')$gene else input$background
 s01_rank_program_score(x,members,bg)
}

s05_estimate <- function(score,pairs,subset=seq_len(nrow(pairs))) {
 if(!isTRUE(score$estimable))return(list(estimable=FALSE,state='NA_NONESTIMABLE',reason=score$reason,biological_n=as.integer(length(subset)),effect=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,p_value=NA_real_,g_av=NA_real_,g_av_ci_lower=NA_real_,g_av_ci_upper=NA_real_))
 pp<-pairs[subset,]; t<-score$score[pp$Tumor_index]; n<-score$score[pp$NAT_index]
 a<-s01_paired_delta_t(t,n,pp$patient_id);g<-s01_paired_g_av(t,n,pp$patient_id)
 a$g_av<-g$g_av;a$g_av_ci_lower<-g$ci_lower;a$g_av_ci_upper<-g$ci_upper
 a$g_av_state<-g$state;a$g_av_reason<-if(is.null(g$reason))'NA_NOT_APPLICABLE'else g$reason
 a$score_data<-data.frame(patient_id=pp$patient_id,NAT_score=n,Tumor_score=t,delta=t-n)
 a
}

s05_empty_row <- function(e) {
 e$effect<-e$ci_lower<-e$ci_upper<-e$p_value<-e$q_value<-NA_real_
 e$g_av<-e$g_av_ci_lower<-e$g_av_ci_upper<-NA_real_
 e$standard_error<-e$statistic<-e$numerator_df<-e$denominator_df<-e$model_degrees_of_freedom<-e$raw_member_paired_mean_difference<-NA_real_
 e$standard_error_status<-e$confidence_interval_status<-'NA_NONESTIMABLE'
 e$uncertainty_reason<-e$reason;e$ql_df_status<-'NA_NOT_APPLICABLE'
 e$estimability<-if(e$structural_eligibility=='INACTIVE_CONDITIONAL_BACKUP')'NA_NOT_APPLICABLE'else 'NA_NONESTIMABLE'
 e$missing_reason<-e$reason;e$multiplicity_status<-'FROZEN_FAMILY_BH'
 e$timing<-'PRESPECIFIED_BEFORE_DISEASE_EFFECT_ACCESS'
 e$dependence<-switch(e$source,GSE29272='independent_population_first_estimate',GSE27342_64='independent_population_replication',GSE27342_80='same_data_robustness',TCGA_STAD='cross_platform_sensitivity','inactive_conditional_backup')
 e$provenance_id<-s05_entry_path();e$native_direction<-'NA_NONESTIMABLE'
 e$member<-if(grepl('limma',e$primary_method))e$family_hypothesis else 'NA_NOT_APPLICABLE'
 e$g_av_state<-'NA_NONESTIMABLE';e$g_av_reason<-e$reason
 e
}

s05_apply_estimate <- function(row,a) {
 for(k in setdiff(intersect(names(a),names(row)),'reason')) if(length(a[[k]])==1L)row[[k]]<-a[[k]]
 row$estimability<-a$state;row$missing_reason<-if(isTRUE(a$estimable))'NA_NOT_APPLICABLE'else a$reason
 if(isTRUE(a$estimable)) {
  row$confidence_interval_status<-'AVAILABLE'
  row$standard_error_status<-if(is.finite(row$standard_error))'AVAILABLE'else 'NA_NOT_APPLICABLE'
  row$uncertainty_reason<-'NA_NOT_APPLICABLE'
 }
 row$native_direction<-s01_direction_label(row$effect);row
}

s05_diagnostic <- function(r,a=NULL,state=NULL,reason=NULL) {
 r$diagnostic_state<-if(is.null(state))a$state else state
 r$missing_reason<-if(is.null(reason))if(is.null(a$reason))'NA_NOT_APPLICABLE'else a$reason else reason
 for(k in c('effect','ci_lower','ci_upper','p_value','g_av','g_av_ci_lower','g_av_ci_upper','statistic')) r[[k]]<-if(!is.null(a[[k]]))a[[k]] else NA_real_
 r$q_value<-NA_real_;r$multiplicity_status<-'NO_FORMAL_FAMILY';r$role<-'same_data_robustness'
 r$biological_n<-if(!is.null(a$biological_n))a$biological_n else NA_integer_
 r$valid_bootstrap_replicates<-if(!is.null(a$valid_replicates))a$valid_replicates else NA_integer_
 r$concordant_members<-r$fixed_denominator<-NA_integer_
 r$operation_name<-'NA_NOT_APPLICABLE';r$seed<-NA_integer_
 r$primary_effect<-NA_real_;r$direction_status<-'NA_NOT_APPLICABLE';r$provenance_id<-s05_entry_path();r
}

s05_analyze_source_internal <- function(source,release,registry) {
 s01_require(isTRUE(release),'release dependency missing');s05_guard()
 start<-s05_utc();bulk<-s05_bulk();ee<-bulk[bulk$source==source,]; config<-yaml::read_yaml(s05_path('config/analysis_config.yml'))
 input<-s05_input(source); n<-nrow(input$pairs)
 expected<-c(GSE29272=112L,GSE27342_64=49L,GSE27342_80=60L,TCGA_STAD=33L)
 s01_require(n==expected[[source]],'retained patient n mismatch')
 model<-NULL
 if(source=='TCGA_STAD') {model<-s05_voom(input$tc);input$expression<-model$voom$E}
 else if(source!='GSE27342_80')model<-list(member=s01_fit_limma_paired(input$expression,input$metadata))
 scores<-list();estimates<-list();formal<-list();diagnostics<-list();influence<-list()
 for(i in seq_len(nrow(ee))) {
  e<-ee[i,];row<-s05_empty_row(e)
  if(grepl('limma',e$primary_method)) {
   if(e$structural_eligibility=='ELIGIBLE') {
    j<-match(e$family_hypothesis,model$member$feature_id);s01_require(!is.na(j),'eligible member absent')
    s01_require(model$member$state[j]=='ESTIMABLE','Unexpected registered member model failure; owner review required')
    for(k in c('effect','standard_error','statistic','p_value','ci_lower','ci_upper'))row[[k]]<-model$member[[k]][j]
    row$model_degrees_of_freedom<-model$member$degrees_of_freedom[j]
    row$raw_member_paired_mean_difference<-mean(input$expression[e$family_hypothesis,input$pairs$Tumor_index]-input$expression[e$family_hypothesis,input$pairs$NAT_index])
    row$estimability<-'ESTIMABLE';row$missing_reason<-'NA_NOT_APPLICABLE'
    row$standard_error_status<-row$confidence_interval_status<-'AVAILABLE';row$uncertainty_reason<-'NA_NOT_APPLICABLE'
    row$native_direction<-s01_direction_label(row$effect)
   }
  } else if(grepl('linear_model',e$primary_method)) {
   parent<-ee$endpoint_id[ee$construct==e$construct&ee$primary_method=='paired_delta_t_test_Student_t_CI_plus_g_av']
   a<-estimates[[parent]];s01_require(!is.null(a$score_data),'site parent unavailable')
   sites<-input$metadata$site[input$pairs$NAT_index]
   s01_require(identical(sites,input$metadata$site[input$pairs$Tumor_index]),'within-pair site disagreement')
   z<-s01_fit_site_interaction(a$score_data$delta,input$pairs$patient_id,sites,'non_cardia','cardia')
   row<-s05_apply_estimate(row,z);row$contrast_direction<-'cardia_minus_non_cardia'
  } else {
   sc<-s05_score(input,s05_members(e$construct,config));scores[[e$endpoint_id]]<-sc
   a<-s05_estimate(sc,input$pairs);estimates[[e$endpoint_id]]<-a
   s01_require(isTRUE(sc$estimable)==(e$structural_eligibility=='ELIGIBLE'),'frozen coverage and scoring disagree')
   row<-s05_apply_estimate(row,a)
   if(e$structural_eligibility!='ELIGIBLE')row$missing_reason<-e$reason
  }
  if(is.null(row$contrast_direction))row$contrast_direction<-'Tumor_minus_NAT'
  formal[[length(formal)+1L]]<-row
 }
 rr<-registry[registry$source==source,];seedmap<-s05_read('config/implementation/s01_target_seed_map.tsv')
 for(i in seq_len(nrow(rr))) {
  r<-rr[i,];members<-s05_members(r$construct,config);base<-estimates[[r$endpoint_id]]
  sub<-r$submodule!='NA_NOT_APPLICABLE';loo<-r$category=='leave_one_gene_out'
  if(sub)members<-config$coverage_contracts$G_submodules[[r$submodule]]
  if(loo)members<-setdiff(members,r$deleted_member)
  # A failed full program cannot be replaced by a reduced-member LOO score.
  if(loo&&!isTRUE(scores[[r$endpoint_id]]$estimable)) {
   d<-s05_diagnostic(r,state='NA_NONESTIMABLE',reason='FULL_PROGRAM_FROZEN_COVERAGE_FAILED_NO_REDUCED_REPLACEMENT')
   diagnostics[[length(diagnostics)+1L]]<-d;next
  }
  spec<-if(r$specification%in%c('standardized_member_mean','common_background','alternative_probe_collapse'))r$specification else 'primary_rank'
  key<-paste(r$endpoint_id,if(r$category=='site_stratum_estimate')r$specification else spec,r$deleted_member,r$submodule,sep='__')
  if(!is.null(scores[[key]]))sc<-scores[[key]] else {
   sc<-if(!sub&&!loo&&spec=='primary_rank')scores[[r$endpoint_id]]else s05_score(input,members,spec)
   scores[[key]]<-sc
  }
  subset<-seq_len(n)
  if(r$category=='site_stratum_estimate')subset<-which(input$metadata$site[input$pairs$NAT_index]==r$specification)
  a<-s05_estimate(sc,input$pairs,subset);estimates[[key]]<-a
  d<-s05_diagnostic(r,a);d$biological_n<-as.integer(length(subset));d$primary_effect<-base$effect
  if(r$category=='paired_percentile_bootstrap') {
   op<-paste0(r$endpoint_id,if(sub)paste0('__G_submodule_',r$submodule)else '', '__paired_percentile_bootstrap')
   seed<-seedmap$seed[match(op,seedmap$operation_name)];s01_require(length(seed)==1L&&!is.na(seed),'unregistered bootstrap operation')
   if(sc$estimable) {
    dat<-a$score_data;b<-s01_pair_bootstrap(dat$Tumor_score,dat$NAT_score,dat$patient_id,seed)
    d<-s05_diagnostic(r,b);d$biological_n<-n;d$effect<-a$effect
   }
   d$operation_name<-op;d$seed<-as.integer(seed)
  } else if(r$category=='paired_Wilcoxon'&&sc$estimable) {
   dat<-a$score_data;w<-s01_paired_wilcoxon(dat$Tumor_score,dat$NAT_score,dat$patient_id)
   d<-s05_diagnostic(r,w);d$biological_n<-n;d$effect<-a$effect
  } else if(r$category=='influence'&&sc$estimable) {
   for(j in seq_len(n)) {
    inf<-s05_estimate(sc,input$pairs,setdiff(seq_len(n),j))
    influence[[length(influence)+1L]]<-data.frame(endpoint_id=r$endpoint_id,diagnostic_id=r$diagnostic_id,
     deleted_patient_id=input$pairs$patient_id[j],biological_n=n-1L,effect=inf$effect,ci_lower=inf$ci_lower,ci_upper=inf$ci_upper,
     p_value=inf$p_value,g_av=inf$g_av,g_av_ci_lower=inf$g_av_ci_lower,g_av_ci_upper=inf$g_av_ci_upper,
     estimability=inf$state,missing_reason=if(is.null(inf$reason))'NA_NOT_APPLICABLE'else inf$reason,
     primary_effect=a$effect,direction_status=s01_loo_direction_status(a$effect,inf$effect,inf$estimable),q_value=NA_real_,multiplicity_status='NO_FORMAL_FAMILY')
   }
  } else if(r$category=='member_direction_concordance') {
   member_delta<-vapply(members,function(g)if(g%in%rownames(input$expression))mean(input$expression[g,input$pairs$Tumor_index]-input$expression[g,input$pairs$NAT_index])else NA_real_,numeric(1))
   con<-s01_member_direction_concordance(base$effect,member_delta,length(members))
   d$diagnostic_state<-con$state;d$effect<-con$proportion;d$p_value<-NA_real_;d$ci_lower<-d$ci_upper<-NA_real_
   d$concordant_members<-if(is.null(con$concordant_members))NA_integer_ else con$concordant_members;d$fixed_denominator<-length(members)
   d$g_av<-d$g_av_ci_lower<-d$g_av_ci_upper<-d$statistic<-NA_real_
  }
  if(loo)d$direction_status<-s01_loo_direction_status(base$effect,a$effect,a$estimable)
  diagnostics[[length(diagnostics)+1L]]<-d
 }
 result<-list(source=source,formal=do.call(rbind,formal),diagnostics=if(length(diagnostics))do.call(rbind,diagnostics)else NULL,
  influence=if(length(influence))do.call(rbind,influence)else NULL,scores=scores,estimates=estimates,
  pairs=input$pairs,metadata=input$metadata,model=model,background=input$background,
  start_utc=start,completed_utc=s05_utc(),status='COMPLETE_UNINSPECTED',
  implementation=s05_hash_rows(implementation_paths))
 revision<-substr(digest::digest(result$implementation,algo='sha256'),1,12)
 s05_rds(result,paste0('data/processed/s05_bulk_effects/',source,'_complete_',revision,'.rds'))
}

s05_analyze_source <- function(source,release,registry) {
 withCallingHandlers(s05_analyze_source_internal(source,release,registry),
  warning=function(w)stop(paste0('UNEXPECTED_S05_WARNING_REQUIRES_REVIEW: ',conditionMessage(w)),call.=FALSE))
}

s05_adjust <- function(rows,config) {
 for(f in unique(rows$multiplicity_family)) {
  spec<-config$multiplicity$families[[f]];ii<-which(rows$multiplicity_family==f)
  s01_require(length(ii)==spec$frozen_size,'incomplete multiplicity family')
  order<-match(spec$hypotheses,rows$family_hypothesis[ii]);s01_require(!anyNA(order),'family identity mismatch')
  ii<-ii[order];rows$q_value[ii]<-unname(s01_frozen_bh(setNames(rows$p_value[ii],spec$hypotheses),f,config))
 }
 rows
}

s05_synthesis <- function(formal,diagnostics,config) {
 states<-list();loo<-list();scoremethod<-list()
 for(axis in c('G','H')) {
  ids<-c(paste0('GSE29272_',axis),paste0('GSE27342_primary_',axis));p<-formal[match(ids,formal$endpoint_id),]
  a<-s01_dimension_a(p$effect[1],p$effect[2],p$q_value[1],p$q_value[2],p$estimability[1]=='ESTIMABLE',p$estimability[2]=='ESTIMABLE')
  b<-s01_dimension_b(p$g_av,p$g_av_ci_lower,p$g_av_ci_upper)
  alt<-diagnostics[match(paste(ids,'standardized_member_mean_direction'),paste(diagnostics$endpoint_id,diagnostics$category)),]
  sm<-s01_score_method_sensitivity(p$effect,alt$effect,alt$diagnostic_state=='ESTIMABLE')
  ld<-diagnostics[diagnostics$endpoint_id%in%ids&diagnostics$category=='leave_one_gene_out',]
  labels<-list()
  for(g in config$programs[[paste0(axis,'_core')]]$members) {
   z<-ld[match(paste(ids,g),paste(ld$endpoint_id,ld$deleted_member)),]
   q<-vapply(seq_along(ids),function(j) {
    other<-formal$p_value[formal$source==p$source[j]&formal$construct==(if(axis=='G')'H_core'else 'G_core')&formal$primary_method=='paired_delta_t_test_Student_t_CI_plus_g_av']
    pv<-if(axis=='G')c(z$p_value[j],other)else c(other,z$p_value[j])
    unname(s01_frozen_bh_values(setNames(as.numeric(pv),c('G_core','H_core')),2L,c('G_core','H_core'))[if(axis=='G')1 else 2])
   },numeric(1))
   la<-s01_dimension_a(z$effect[1],z$effect[2],q[1],q[2],z$diagnostic_state[1]=='ESTIMABLE',z$diagnostic_state[2]=='ESTIMABLE')
   labels[[g]]<-data.frame(construct=axis,deleted_member=g,source_1= p$source[1],source_2=p$source[2],
    effect_1=z$effect[1],effect_2=z$effect[2],diagnostic_BH2_q_1=q[1],diagnostic_BH2_q_2=q[2],
    dimension_a_state=la$state,primary_dimension_a_state=a$state,label_changed=la$state!='non_estimable'&&la$state!=a$state,
    estimable=la$state!='non_estimable',multiplicity_status='DIAGNOSTIC_BH2_RECALCULATION_NOT_DISCOVERY')
  }
  l<-do.call(rbind,labels);loo[[axis]]<-l
  state<-data.frame(construct=axis,active_source_1_endpoint_id=ids[1],active_source_2_endpoint_id=ids[2],
   dimension_a_state=a$state,zero_point_estimate_flag=a$zero_point_estimate_flag,dimension_b_state=b$magnitude_class,
   both_core_CI_bounded_within_small_effect_region=b$both_core_CI_bounded_within_small_effect_region,precision_uncertain=b$precision_uncertain,
   primary_rank_unrounded_raw_paired_mean_difference=p$effect[1],standardized_member_mean_unrounded_raw_paired_mean_difference=alt$effect[1],
   primary_rank_direction=sm$primary_rank_direction[1],standardized_member_mean_direction=sm$standardized_member_mean_direction[1],
   source_2_primary_rank_unrounded_raw_paired_mean_difference=p$effect[2],source_2_standardized_member_mean_unrounded_raw_paired_mean_difference=alt$effect[2],
   source_2_primary_rank_direction=sm$primary_rank_direction[2],source_2_standardized_member_mean_direction=sm$standardized_member_mean_direction[2],
   score_method_sensitive=sm$score_method_sensitive,score_method_nonestimable=sm$score_method_nonestimable,exact_zero_direction_flag=sm$exact_zero_direction_flag,
   loo_direction_stable_both_core=s01_aggregate_loo_status(ld$direction_status),
   loo_replication_label_sensitive=s01_loo_replication_sensitivity(l$label_changed,l$estimable),
   claim_ceiling=p$claim_ceiling[1],provenance_id=s05_entry_path())
  states[[axis]]<-state
 }
 # Complete source-level method diagnostics, including frozen unavailable axes.
 d<-diagnostics[diagnostics$category=='standardized_member_mean_direction',]
 for(i in seq_len(nrow(d))) {
  p<-formal[match(d$endpoint_id[i],formal$endpoint_id),]
  sm<-s01_score_method_sensitivity(p$effect,d$effect[i],d$diagnostic_state[i]=='ESTIMABLE')
  scoremethod[[i]]<-cbind(data.frame(endpoint_id=p$endpoint_id,source=p$source,construct=p$construct),as.data.frame(sm))
 }
 list(states=do.call(rbind,states),loo=do.call(rbind,loo),scoremethod=do.call(rbind,scoremethod))
}

s05_synthesis_identity <- function() {
 vapply(list(adjust=s05_adjust,synthesis=s05_synthesis),function(fn)
  digest::digest(list(formals(fn),body(fn)),algo='sha256'),character(1))
}

s05_source_paths <- function() {
 implementation<-s05_hash_rows(implementation_paths)
 revision<-substr(digest::digest(implementation,algo='sha256'),1,12)
 paste0('data/processed/s05_bulk_effects/',c('GSE29272','GSE27342_64','GSE27342_80','TCGA_STAD'),'_complete_',revision,'.rds')
}

s05_source_projection <- function(paths,registry) {
 s01_require(identical(unname(paths),s05_source_paths()),'source identities/paths missing, duplicated, substituted or reordered')
 objects<-lapply(paths,function(p)read_rds(s05_path(p)))
 s01_require(identical(unname(vapply(objects,`[[`,character(1),'source')),c('GSE29272','GSE27342_64','GSE27342_80','TCGA_STAD')),'complete graph sources missing or reordered')
 s01_require(all(vapply(objects,function(x)x$status=='COMPLETE_UNINSPECTED',logical(1))),'incomplete source cannot pass barrier')
 for(x in objects)s05_verify_rows(x$implementation)
 bulk<-s05_bulk()
 for(x in objects) {
  expected<-bulk$endpoint_id[bulk$source==x$source]
  s01_require(identical(x$formal$endpoint_id,expected),'source formal identities/order mismatch')
  dk<-function(z)paste(z$endpoint_id,z$diagnostic_id,sep='\t')
  expected_d<-registry[registry$source==x$source,]
  s01_require(setequal(dk(x$diagnostics),dk(expected_d))&&!anyDuplicated(dk(x$diagnostics)), 'source diagnostic identities mismatch')
 }
 inactive<-bulk[bulk$structural_eligibility=='INACTIVE_CONDITIONAL_BACKUP',]
 na<-lapply(seq_len(nrow(inactive)),function(i){z<-s05_empty_row(inactive[i,]);z$contrast_direction<-'Tumor_minus_NAT';z})
 formal<-do.call(rbind,c(lapply(objects,`[[`,'formal'),na));formal<-formal[match(bulk$endpoint_id,formal$endpoint_id),]
 config<-yaml::read_yaml(s05_path('config/analysis_config.yml'));formal<-s05_adjust(formal,config)
 dr<-registry[registry$source=='GSE122401_backup',]
 # Use actual registered source spelling without activating or loading the backup matrix.
 if(!nrow(dr))dr<-registry[grepl('GSE122401',registry$source),]
 dd<-lapply(seq_len(nrow(dr)),function(i)s05_diagnostic(dr[i,],state='NA_NOT_APPLICABLE',reason='NORMAL_BRANCH_NO_ACTIVATION'))
 diagnostics<-do.call(rbind,c(lapply(objects,`[[`,'diagnostics'),dd))
 key<-function(x)paste(x$endpoint_id,x$diagnostic_id)
 diagnostics<-diagnostics[match(key(registry),key(diagnostics)),]
 synthesis<-s05_synthesis(formal,diagnostics,config)
 influence<-do.call(rbind,lapply(objects,`[[`,'influence'))
 sources<-list()
 for(x in objects)for(id in names(x$estimates)) {
  a<-x$estimates[[id]];if(is.null(a$score_data))next
  z<-a$score_data;z$analysis_id<-id;z$source<-x$source;sources[[length(sources)+1L]]<-z
 }
 source_data<-do.call(rbind,sources)
 list(formal=formal,diagnostics=diagnostics,registry=registry,synthesis=synthesis,
  influence=influence,source_data=source_data)
}

s05_compare_rows <- function(actual,expected,keys,label) {
 key<-function(z)do.call(paste,c(z[keys],sep='\t'))
 s01_require(all(keys%in%names(actual))&&all(keys%in%names(expected)),paste(label,'missing key'))
 s01_require(!anyDuplicated(key(actual))&&!anyDuplicated(key(expected))&&
  identical(key(actual),key(expected)),paste(label,'identities/order mismatch'))
 rownames(actual)<-NULL;rownames(expected)<-NULL
 s01_require(identical(actual,expected),paste(label,'source content mismatch'))
 invisible(TRUE)
}

s05_reconcile_sources <- function(result) {
 expected<-s05_source_projection(result$source_paths,result$registry)
 # Compare every formal field to source records, then independently restore frozen BH.
 s01_require(identical(names(result$formal),names(expected$formal)),'formal schema mismatch')
 fields<-setdiff(names(expected$formal),'q_value')
 s05_compare_rows(result$formal[fields],expected$formal[fields],'endpoint_id','formal')
 s01_require(identical(result$formal$q_value,expected$formal$q_value),'source-based frozen BH mismatch')
 s05_compare_rows(result$diagnostics,expected$diagnostics,c('endpoint_id','diagnostic_id'),'diagnostic')
 s05_compare_rows(result$influence,expected$influence,c('endpoint_id','diagnostic_id','deleted_patient_id'),'influence')
 z<-result$source_data
 # RDS retains the original IEEE doubles: the identical subtraction must be exact.
 s01_require(identical(z$delta,z$Tumor_score-z$NAT_score),'patient source-data delta mismatch')
 s05_compare_rows(z,expected$source_data,c('source','analysis_id','patient_id'),'patient source-data')
 invisible(TRUE)
}

s05_assemble <- function(paths,registry,release) {
 s01_require(isTRUE(release),'release dependency missing');s05_guard()
 projection<-s05_source_projection(paths,registry)
 path<-'data/processed/s05_bulk_effects/assembled.rds'
 # Preserve the reviewed RDS byte-for-byte when all content and provenance still agree.
 if(all(dirname(paths)=='data/processed/s05_bulk_effects')&&file.exists(s05_path(path))) {
  prior<-read_rds(s05_path(path))
  if(identical(prior[names(projection)],projection)&&identical(prior$source_paths,paths)&&
   identical(prior$synthesis_implementation,s05_synthesis_identity())&&
   identical(prior$status,'ALL_ANALYSES_COMPLETE_UNINSPECTED'))return(path)
 }
 result<-c(projection,list(source_paths=paths,synthesis_implementation=s05_synthesis_identity(),
  completed_utc=s05_utc(),status='ALL_ANALYSES_COMPLETE_UNINSPECTED'))
 s05_rds(result,path)
}

s05_validate_result <- function(result,bulk=s05_bulk()) {
 f<-result$formal;d<-result$diagnostics;r<-result$registry
 s01_require(result$status=='ALL_ANALYSES_COMPLETE_UNINSPECTED','incomplete graph cannot finalize')
 s01_require(identical(result$synthesis_implementation,s05_synthesis_identity()),'stale synthesis implementation')
 s01_require(identical(r,s05_registry()),'stale diagnostic execution registry')
 s05_reconcile_sources(result)
 s01_require(nrow(f)==72L&&identical(f$endpoint_id,bulk$endpoint_id)&&!anyDuplicated(f$endpoint_id),'bulk crosswalk reconciliation failed')
 for(k in names(bulk))s01_require(identical(f[[k]],bulk[[k]]),paste('frozen row metadata changed',k))
 s01_require(identical(paste(d$endpoint_id,d$diagnostic_id),paste(r$endpoint_id,r$diagnostic_id))&&!anyDuplicated(paste(d$endpoint_id,d$diagnostic_id)),'diagnostic coverage reconciliation failed')
 s01_require(!any(d$source=='GSE27342_80')&&sum(f$source=='GSE27342_80')==2L,'80 scope violation')
 ok<-f$estimability=='ESTIMABLE';inactive<-f$structural_eligibility=='INACTIVE_CONDITIONAL_BACKUP'
 s01_require(sum(inactive)==17L&&all(is.na(f$effect[inactive]))&&all(is.na(f$p_value[!ok]))&&all(is.na(f$q_value[!ok])),'inactive/nonestimable effect or family leakage')
 s01_require(all(f$estimability[f$structural_eligibility=='NA_NONESTIMABLE']=='NA_NONESTIMABLE'),'frozen unavailable endpoint replaced')
 s01_require(all(is.finite(f$effect[ok]))&&all(is.finite(f$ci_lower[ok]))&&all(is.finite(f$ci_upper[ok]))&&all(is.finite(f$p_value[ok]))&&all(is.finite(f$q_value[ok])),'invalid formal numerical output')
 s01_require(all(f$ci_lower[ok]<=f$effect[ok]&f$effect[ok]<=f$ci_upper[ok])&&all(f$p_value[ok]>=0&f$p_value[ok]<=1),'invalid formal interval or probability')
 expected<-s05_adjust(f,yaml::read_yaml(s05_path('config/analysis_config.yml')))
 s01_require(identical(expected$q_value,f$q_value),'BH numerical mismatch')
 s01_require(all(!is.na(f$missing_reason))&&all(nzchar(f$missing_reason)),'missing terminal reasons')
 numeric<-names(f)[vapply(f,is.numeric,logical(1))]
 s01_require(!any(vapply(f[numeric],function(v)any(is.nan(v)|is.infinite(v)),logical(1))),'nonfinite numerical table values')
 # Reconstruct native t and g_av from persisted patient scores; no cache-skip claim.
 checked<-0L
 for(p in result$source_paths) {
  x<-read_rds(s05_path(p));s05_verify_rows(x$implementation)
  s05_pairs(x$metadata,x$metadata$sample_id)
  for(a in x$estimates) {
   if(is.null(a$score_data))next
   z<-a$score_data;s01_require(identical(z$delta,z$Tumor_score-z$NAT_score),'paired delta mismatch')
   native<-s01_paired_delta_t(z$Tumor_score,z$NAT_score,z$patient_id)
   s01_require(identical(a$effect,native$effect)&&identical(a$p_value,native$p_value)&&identical(a$ci_lower,native$ci_lower)&&identical(a$ci_upper,native$ci_upper),'persisted native estimate mismatch')
   if(native$estimable) {
    anchor<-stats::t.test(z$delta);s01_require(isTRUE(all.equal(native$p_value,anchor$p.value,tolerance=1e-12))&&isTRUE(all.equal(c(native$ci_lower,native$ci_upper),as.numeric(anchor$conf.int),tolerance=1e-12)),'native Student t anchor mismatch')
   }
   g<-s01_paired_g_av(z$Tumor_score,z$NAT_score,z$patient_id)
   s01_require(identical(a$g_av,g$g_av)&&identical(a$g_av_ci_lower,g$ci_lower)&&identical(a$g_av_ci_upper,g$ci_upper),'persisted magnitude mismatch')
   checked<-checked+1L
  }
  if(x$source!='GSE27342_80') {
   bd<-x$diagnostics[x$diagnostics$category=='paired_percentile_bootstrap'&x$diagnostics$diagnostic_state=='ESTIMABLE',]
   for(i in seq_len(nrow(bd))) {
    key<-paste(bd$endpoint_id[i],'primary_rank','NA_NOT_APPLICABLE',bd$submodule[i],sep='__')
    dat<-x$estimates[[key]]$score_data
    boot<-s01_pair_bootstrap(dat$Tumor_score,dat$NAT_score,dat$patient_id,bd$seed[i])
    s01_require(identical(boot$ci_lower,bd$ci_lower[i])&&identical(boot$ci_upper,bd$ci_upper[i])&&
     identical(boot$valid_replicates,bd$valid_bootstrap_replicates[i]),'persisted bootstrap reconstruction mismatch')
   }
   mm<-x$formal[grepl('limma',x$formal$primary_method)&x$formal$estimability=='ESTIMABLE',]
   s01_require(all(mm$standard_error>0)&&all(mm$confidence_interval_status=='AVAILABLE'),'limma uncertainty exception wrongly applied')
   j<-match(mm$member,x$model$member$feature_id)
   crit<-stats::qt(.975,x$model$member$degrees_of_freedom[j])
   s01_require(isTRUE(all.equal(mm$ci_lower,mm$effect-crit*mm$standard_error,tolerance=1e-14)),'limma model interval mismatch')
  }
 }
 bs<-d$category=='paired_percentile_bootstrap'&d$diagnostic_state=='ESTIMABLE'
 s01_require(all(d$valid_bootstrap_replicates[bs]>=9500L)&all(d$valid_bootstrap_replicates[bs]<=10000L),'bootstrap finite count mismatch')
 s01_require(all(is.na(d$q_value)),'nominal diagnostic acquired q value')
 s01_require(nrow(result$synthesis$states)==2L&&nrow(result$synthesis$loo)==11L,'bulk synthesis coverage failed')
 reconstructed<-s05_synthesis(f,d,yaml::read_yaml(s05_path('config/analysis_config.yml')))
 s01_require(identical(result$synthesis,reconstructed),'synthesis numerical/state reconstruction mismatch')
 list(status='PASS',verified_utc=s05_utc(),crosswalk_hypotheses=72L,active_population_hypotheses=55L,
  inactive_backup_hypotheses=17L,frozen_nonestimable_hypotheses=5L,formal_estimable=sum(ok),
  diagnostic_records=nrow(d),patient_score_estimate_reconstructions=checked,influence_rows=nrow(result$influence),
  assembled_source_reconciliation=TRUE,patient_source_data_rows=nrow(result$source_data),
  exact_patient_delta_verified=TRUE,biological_n_verified=TRUE,complete_families_verified=TRUE,model_intervals_verified=TRUE,inspection_barrier_complete=TRUE)
}

s05_auxiliary_tables <- function(result) {
 d<-result$diagnostics
 exclusions<-s05_read('reports/qc/s02_bulk_ingress/exclusions_and_deviations.tsv')
 exclusions$provenance_reference<-'reports/qc/s02_bulk_ingress/exclusions_and_deviations.tsv'
 list(loo_direction_stable_by_source=d[d$category=='leave_one_gene_out',],
  submodule_estimates=d[d$category=='G_submodule',],
  bootstrap_wilcoxon=d[d$category%in%c('paired_percentile_bootstrap','paired_Wilcoxon'),],
  site_stratum_estimates=d[d$category=='site_stratum_estimate',],
  nonestimability=result$formal[result$formal$estimability!='ESTIMABLE',],
  exclusions_and_deviations=exclusions)
}

s05_table_projection <- function(result) {
 f<-result$formal;member<-grepl('limma',f$primary_method);site<-grepl('linear_model',f$primary_method)
 schema<-yaml::read_yaml(s05_path('config/implementation/s01_output_schemas.yml'))
 tables<-list(all_program_estimates=f[!member,],all_member_estimates=f[member,],site_interaction_estimates=f[site,],
  bulk_reproducibility_states=result$synthesis$states,robustness_diagnostics=result$diagnostics,
  loo_replication_labels=result$synthesis$loo,score_method_diagnostics=result$synthesis$scoremethod,
  influence_diagnostics=result$influence,patient_score_deltas=result$source_data,
  hypothesis_coverage=f[,c('endpoint_id','source','biological_n','structural_eligibility','estimability','missing_reason','multiplicity_family','family_hypothesis','frozen_family_size')],
  nominal_coverage=result$diagnostics[,c('endpoint_id','diagnostic_id','source','category','diagnostic_state','missing_reason')])
 for(name in names(tables)) {
  x<-tables[[name]];rownames(x)<-NULL
  if('model_validity'%in%names(x)) {
   names(x)[names(x)=='model_validity']<-'pre_effect_model_validity'
   x$model_validity<-ifelse(x$estimability=='ESTIMABLE','VALID_REGISTERED_MODEL',x$estimability)
   names(x)[names(x)=='reason']<-'pre_effect_reason'
  }
  if(name%in%names(schema$outputs))s01_require(all(schema$outputs[[name]]$required_columns%in%names(x)),paste('output schema missing',name))
  if('endpoint_id'%in%names(x))x<-x[do.call(order,c(x[intersect(c('endpoint_id','diagnostic_id','analysis_id','patient_id','deleted_patient_id'),names(x))],list(method='radix'))),]
  tables[[name]]<-x
 }
 c(tables,s05_auxiliary_tables(result))
}

s05_compare_table <- function(x,y,label) {
 s01_require(identical(names(x),names(y))&&nrow(x)==nrow(y),paste(label,'table schema mismatch'))
 for(k in names(x)) {
  # Existing 17-digit TSV representation: allow only serialization-scale error.
  if(is.numeric(x[[k]])) {
   u<-as.numeric(x[[k]]);v<-as.numeric(y[[k]])
   s01_require(identical(is.na(u),is.na(v))&&all(abs(u[!is.na(u)]-v[!is.na(v)])<=1e-13*pmax(1,abs(u[!is.na(u)]))),paste(label,k,'table content mismatch'))
  }
  else s01_require(identical(as.character(x[[k]]),as.character(y[[k]])),paste(label,k,'table identity mismatch'))
 }
 if(all(c('NAT_score','Tumor_score','delta')%in%names(y))) {
  error<-abs(y$delta-(y$Tumor_score-y$NAT_score))
  bound<-1e-13*pmax(1,abs(y$Tumor_score),abs(y$NAT_score))
  s01_require(all(is.finite(error)&error<=bound),paste(label,'exported patient delta mismatch'))
 }
 invisible(TRUE)
}
  run <- function() {
    context$guard();context$resource_check()
    registry<-s05_registry()
    paths<-vapply(c('GSE29272','GSE27342_64','GSE27342_80','TCGA_STAD'),function(id)s05_analyze_source(id,TRUE,registry),'',USE.NAMES=FALSE)
    ap<-s05_assemble(paths,registry,TRUE);result<-read_rds(ap)
    validation<-s05_validate_result(result);tables<-s05_table_projection(result)
    out<-vapply(names(tables),function(n){p<-context$output('results',paste0(n,'.tsv'));context$write_tsv(tables[[n]],p);s05_compare_table(tables[[n]],s05_read(p),n);p},'')
    context$write_rds(validation,context$output('qc','numerical_validation.rds'))
    list(scope='original_bulk',tables=out,objects=c(vapply(paths,resolve_path,''),resolve_path(ap)),validation=validation)
  }
  list(run=run, input=s05_input, registry=s05_registry, project=s05_table_projection, validate=s05_validate_result,
       kernels=list(score=s05_score,estimate=s05_estimate,voom=s05_voom,synthesis=s05_synthesis), resolve=resolve_path)
}
