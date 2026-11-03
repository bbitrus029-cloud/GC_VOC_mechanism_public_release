# GSE122401 abundance and D7 reference preparation; no acquisition or disease fitting.
supplementary_preparation <- function(context) {
 s10_m<-context$output('metadata');s10_p<-context$output('processed');s10_q<-context$output('qc')
 s10_guard<-context$guard;s10_memory<-function(stage,minimum_GiB=8)context$resource_check(minimum_GiB,stage)
 s10_read<-context$read_tsv;s10_tsv<-context$write_tsv;s10_rds<-context$write_rds;s10_json<-context$write_json
 s10_sha<-function(p)digest::digest(file=context$input(p),algo='sha256')
 s10_verify_rows<-function(m,path_field='path') {
   q<-context$input(m[[path_field]]);stopifnot(identical(unname(vapply(q,s01_sha256,'')),m$sha256))
 }
s10_A <- function(ledger_path,file_inventory_path,background_path,annotation_path) {
 s10_guard();s10_memory("A_start")
 ledger<-s10_read(ledger_path); inventory<-s10_read(file_inventory_path)
 files<-context$input(inventory$logical_input[match(ledger$processed_gene_file,inventory$file_basename)])
 expected<-inventory$sha256[match(ledger$processed_gene_file,inventory$file_basename)]
 s01_require(!anyNA(expected)&&all(vapply(files,s10_sha,"")==expected),"GSE122401 gene payload identity mismatch")
 schema<-gse122401_inspect_processed_schemas(files,160L,TRUE,TRUE)
 bg<-s10_read(background_path)$ensembl_gene_id
 first<-gse122401_read_processed_gene_file(files[1],schema)
 s01_require(!anyDuplicated(first$gene_id),"GSE122401 duplicate feature")
 mat<-matrix(NA_real_,nrow(first),160L,dimnames=list(first$gene_id,ledger$gsm))
 for(i in seq_along(files)) {
  s10_memory("A_atomic_file",6)
  z<-gse122401_read_processed_gene_file(files[i],schema)
  s01_require(identical(z$gene_id,rownames(mat)),"GSE122401 feature ordering changed")
  mat[,i]<-z[[schema$expression_column]]
 }
 s01_validate_numeric_matrix(mat,"GSE122401 abundance")
 s01_require(all(mat>=0)&&any(mat!=round(mat)),"GSE122401 continuous nonnegative scale mismatch")
 observed<-sort(rownames(mat)[rowSums(mat>0)>=32L],method="radix")
 s01_require(identical(bg,observed)&&length(bg)==23818L&&nrow(mat)==34908L,"GSE122401 pooled background mismatch")
 md<-data.frame(sample_id=ledger$gsm,patient_id=ledger$subject_id,tissue=ifelse(ledger$tissue_category=="adjacent_normal","NAT","Tumor"),
  original_tissue=ledger$tissue_category,source="GSE122401",role="external_cohort_post_effect")
 s01_require(nrow(md)==160L&&length(unique(md$patient_id))==80L&&!anyDuplicated(paste(md$patient_id,md$tissue))&&
  all(table(md$patient_id)==2L)&&all(ledger$paired_gsm[match(ledger$paired_gsm,ledger$gsm)]==ledger$gsm),"GSE122401 pair reconciliation failed")
 native_design<-s01_prepare_paired_model(mat[bg,,drop=FALSE],md,"patient_id","patient")
 s01_require(qr(native_design$design)$rank==81L,"GSE122401 design rank mismatch")
 cfg<-context$supplementary_config
 symbols<-c(unlist(cfg$programs,use.names=FALSE),cfg$D7)
 annotation<-s10_read(annotation_path)
 annotation$ensembl_id<-sub("[.][0-9]+$","",annotation$gene_id)
 mapping<-do.call(rbind,lapply(symbols,function(g) {
  j<-which(annotation$gene_name==g & annotation$unique_symbol)
  s01_require(length(j)==1L,paste("Ambiguous canonical mapping",g))
  data.frame(gene=g,ensembl_id=annotation$ensembl_id[j],mapping_authority="frozen_GENCODE_v36_GDC_gene_identity",present=annotation$ensembl_id[j]%in%rownames(mat),
   background=annotation$ensembl_id[j]%in%bg)
 }))
 historical<-gse122401_frozen_members()
 s01_require(all(mapping$ensembl_id[match(historical$gene_symbol,mapping$gene)]==historical$ensembl_gene_id),"G7/H4 historical mapping changed")
 s01_require(all(mapping$present)&&all(mapping$background),"Required program or D7 member unavailable")
 s10_tsv(mapping,file.path(s10_m,"GSE122401_feature_bridge.tsv"));s10_tsv(md,file.path(s10_m,"GSE122401_sample_pair_ledger.tsv"))
 s10_tsv(data.frame(ensembl_gene_id=bg),file.path(s10_m,"GSE122401_background.tsv"))
 s10_tsv(data.frame(path=files,bytes=file.info(files)$size,sha256=expected,source="GSE122401",role="external_cohort_post_effect"),file.path(s10_m,"GSE122401_payload_identity.tsv"))
 coverage<-do.call(rbind,lapply(names(cfg$programs),function(p)data.frame(program=p,member=cfg$programs[[p]],mapping_state="UNIQUE",background_state="AVAILABLE")))
 s10_tsv(coverage,file.path(s10_q,"A_all_program_coverage.tsv"))
 diagnostic<-data.frame(analysis_id=c(paste0("PEA.standardized_member_mean.",c("G_core7","H_core4")),
  paste0("PEA.LOO.",unlist(cfg$programs[1:2],use.names=FALSE)),paste0("PEA.Wilcoxon.",c("G_core7","H_core4")),
  paste0("PEA.bootstrap.",c("G_core7","H_core4")),paste0("PEA.member_concordance.",c("G_core7","H_core4"))),
  source="GSE122401",role="same_data_robustness",timing="post_effect_exploratory",status="PENDING_EFFECT_RELEASE")
 s10_tsv(diagnostic,file.path(s10_m,"A_diagnostic_registry.tsv"))
 o<-list(source="GSE122401",role="external_cohort_post_effect",abundance=mat[bg,,drop=FALSE],metadata=md,background=bg,
  member_mapping=mapping,schema=schema,historical_activation="INACTIVE_ORIGINAL_CONDITIONAL_BACKUP",model_fitted=FALSE,program_scores_calculated=FALSE)
 path<-s10_rds(o,file.path(s10_p,"GSE122401_adapter.rds"))
 s10_json(list(state="READY_FOR_INDEPENDENT_REVIEW",biological_pairs=80L,samples=160L,background_genes=23818L,
  excluded_pseudoheaders=160L,G7=7L,H4=4L,context_members=15L,clinical_context="Mun2019 early_onset_diffuse_predominant; patient_clinical_fields_unknown_unless_official_metadata",
  historical_backup="INACTIVE_ORIGINAL_CONDITIONAL_BACKUP",input_manifest_sha256=context$manifest_sha256,
  adapter_sha256=s10_sha(path),scale="continuous_UQ_RSEM_derived_abundance",real_contrasts="NOT_RUN",
  overlap="same_study_protein_candidate_dependent; cross_cohort_patient_overlap_not_proven_absent",
  effect_release="NOT_RELEASED"),file.path(s10_m,"A_TECHNICAL_FREEZE.json"))
 path
}

s10_D7 <- function(A) {
 s10_guard();s10_memory("D7_start");cfg<-context$supplementary_config
 before<-context$bindings()
 references<-list();moments<-list();ledger<-list();coverage<-list();inputs<-list();mapping_audit<-list()
 for(s in cfg$C_D7$cohorts) {
  if(s=="GSE122401") {
   a<-context$read_rds(A);mapping<-a$member_mapping;genes<-mapping$ensembl_id[match(cfg$D7,mapping$gene)]
   x<-log2(a$abundance+1);rownames(x)[match(genes,rownames(x))]<-cfg$D7;md<-a$metadata;path<-A;scale<-"log2_UQ_RSEM_abundance_plus_1"
  } else if(s=="TCGA_STAD") {
   path<-"data/processed/s09_reproducibility/s05_bulk_effects/TCGA_STAD_complete_3d541cd5368a.rds"
   m<-context$bindings();s10_verify_rows(m[match(path,m$path),])
   a<-context$read_rds(path);x<-a$model$voom$E;md<-a$metadata;scale<-"frozen_TMM_voom_log2_expression_from_released_S09_representation"
  } else {
   name<-if(s=="GSE27342_core")"GSE27342"else s
   path<-paste0("data/processed/s02_bulk_ingress/",name,"_primary.rds")
   lp<-paste0("reports/qc/s02_bulk_ingress/",name,"_pair_ledger.tsv")
   j<-match(c(path,lp),before$path);s01_require(!anyNA(j),"D7 original input absent from freeze");s10_verify_rows(before[j,])
   x<-context$read_rds(path);d<-s10_read(lp);d<-d[!d$pair_excluded&d$core_population,,drop=FALSE]
   md<-data.frame(sample_id=d$neutral_id,patient_id=d$patient_id,tissue=d$tissue)
   md<-md[match(colnames(x),md$sample_id),,drop=FALSE];scale<-"frozen_RMA_log2_gene_summary"
   mp<-paste0("reports/qc/s02_bulk_ingress/",name,"_mapping.tsv")
   ap<-paste0("reports/qc/s02_bulk_ingress/",name,"_processed_feature_audit.tsv")
   j<-match(c(mp,ap),before$path);s01_require(!anyNA(j),"D7 array mapping audit absent from freeze");s10_verify_rows(before[j,])
   map<-s10_read(mp);audit<-s10_read(ap)
   mapping_audit[[s]]<-do.call(rbind,lapply(cfg$D7,function(g) {
    rows<-map[map$canonical_gene==g&map$status=="RELIABLE",,drop=FALSE]
    measured<-rows$feature_id %in% audit$feature_id
    data.frame(source=s,member=g,reliable_mapped_features=nrow(rows),measured_features=sum(measured),
     canonical_entrez=paste(sort(unique(rows$entrez_id),method="radix"),collapse=";"),
     state=if(g%in%rownames(x))"UNIQUE_CANONICAL_GENE_AVAILABLE"else if(nrow(rows)>0&&!any(measured))"MAPPED_FEATURE_ABSENT_FROM_FROZEN_RMA_UNIVERSE"else"NO_RELIABLE_CANONICAL_MEASUREMENT",
     mapping_sha256=s10_sha(mp),processed_feature_audit_sha256=s10_sha(ap))
   }))
  }
  md<-md[,c("sample_id","patient_id","tissue"),drop=FALSE]
  s01_require(!anyNA(md)&&identical(md$sample_id,colnames(x)),paste("D7 population mismatch",s))
  expected<-c(GSE29272=112L,GSE27342_core=49L,GSE122401=80L,TCGA_STAD=33L)
  s01_require(length(unique(md$patient_id))==expected[[s]],paste("D7 formal population count changed",s))
  r<-pe_d7_reference(x,md,cfg$D7);r$source<-s;r$input_path<-path;r$input_sha256<-s10_sha(path);r$scale<-scale
  references[[s]]<-r
  coverage[[s]]<-data.frame(source=s,member=cfg$D7,present=cfg$D7%in%rownames(x),state=r$state,reason=r$reason)
  inputs[[s]]<-data.frame(source=s,path=path,bytes=file.info(context$input(path))$size,sha256=s10_sha(path),scale=scale,formal_pairs=expected[[s]])
  if(!is.null(r$metadata)) {
   ledger[[s]]<-cbind(source=s,r$metadata)
   moments[[s]]<-data.frame(source=s,member=cfg$D7,pooled_mean=unname(r$mean),sample_SD=unname(r$sd),constant=unname(r$constant),
    reference_samples=r$N_samples,reference_pairs=length(r$patient_ids),sd_denominator=r$sd_denominator,reference_sha256=r$reference_sha256)
  }
 }
 s10_tsv(do.call(rbind,inputs),file.path(s10_m,"D7_input_identity.tsv"))
 s10_tsv(do.call(rbind,coverage),file.path(s10_q,"D7_coverage.tsv"))
 s10_tsv(do.call(rbind,mapping_audit),file.path(s10_m,"D7_array_mapping_audit.tsv"))
 s10_tsv(do.call(rbind,moments),file.path(s10_m,"D7_pooled_reference_moments.tsv"))
 s10_tsv(do.call(rbind,ledger),file.path(s10_m,"D7_reference_population.tsv"))
 path<-s10_rds(references,file.path(s10_p,"D7_reference_freeze.rds"))
 s10_json(list(state="READY_WITH_DOCUMENTED_COHORT_ESTIMABILITY_STATES",
  reference_sha256=s10_sha(path),ordered_id_manifest_sha256=s10_sha(file.path(s10_m,"D7_reference_population.tsv")),
  pooled_moments_sha256=s10_sha(file.path(s10_m,"D7_pooled_reference_moments.tsv")),
  actual_paired_deltas="PENDING_EFFECT_RELEASE",actual_model_rank_covariate_variation_support_HC3_influence="PENDING_EFFECT_RELEASE",
  actual_D7_scores="NOT_CALCULATED",composite_restandardization=FALSE,effect_release="NOT_RELEASED"),file.path(s10_m,"C_D7_TECHNICAL_FREEZE.json"))
 path
}

list(A=s10_A,D7=s10_D7)
}
