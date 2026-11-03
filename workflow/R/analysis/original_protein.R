# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Public adaptation: explicit inputs and outputs; original scientific rules retained.

s07_boot_spearman <- function(x,i) stats::cor(x$RNA_delta[i],x$protein_delta[i],method='spearman')
original_protein_module <- function(context) {
  historical_namespace <- 's07_protein_bridge'
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
  s07_path <- resolve_path
  s07_m <- function(p) file.path('data/manifests',historical_namespace,p)
  s07_o <- function(p) file.path('data/processed',historical_namespace,p)
  s07_read <- function(path) utils::read.delim(context$input(resolve_path(path)),quote='',comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
  s07_cfg <- function() context$config
  s07_utc <- utc
  s07_guard <- function() context$guard()
  s07_resource <- function() context$resource_check()
  s07_save <- write_rds
s07_lf_hash <- function(x) digest::digest(paste0(paste(x,collapse='\n'),'\n'),algo='sha256',serialize=FALSE)

s07_compare <- function(x,y,label) {
 s01_require(is.data.frame(x)&&is.data.frame(y)&&identical(names(x),names(y))&&nrow(x)==nrow(y),paste(label,'schema mismatch'))
 for(k in names(x)) {
  if(is.numeric(x[[k]])) {
   a<-as.numeric(x[[k]]);b<-as.numeric(y[[k]])
   s01_require(identical(is.na(a),is.na(b))&&all(a[!is.na(a)]==b[!is.na(b)]),paste(label,k,'numeric mismatch'))
  } else s01_require(identical(as.character(x[[k]]),as.character(y[[k]])),paste(label,k,'identity mismatch'))
 };invisible(TRUE)
}

s07_crosswalk <- function() {
 e<-s07_read('config/hypothesis_estimand_crosswalk.tsv');e<-e[e$contract%in%c('C3','C4'),];rownames(e)<-NULL;e
}

s07_population <- function(endpoint) {
 if(startsWith(endpoint,'CPTAC_bridge_'))return(sub('CPTAC_','',endpoint,fixed=TRUE))
 if(grepl('_member_',endpoint))return(paste0('protein_',sub('^PDC000614_([GH])_member_.*','\\1',endpoint),'_program'))
 if(startsWith(endpoint,'PDC000614_protein_'))return(paste0(sub('PDC000614_','',endpoint,fixed=TRUE),'_program'))
 sub('PDC000614_context_','',endpoint,fixed=TRUE)
}

s07_members <- function(population,cfg=s07_cfg()) {
 if(startsWith(population,'protein_'))return(cfg$programs[[paste0(substr(population,9,9),'_core')]]$protein$members)
 if(startsWith(population,'bridge_'))return(cfg$programs[[paste0(substr(population,8,8),'_core')]]$bridge$members)
 cfg$coverage_contracts$context_programs[[population]]
}

s07_pair_index <- function(identity,patients,layer) {
 required<-c('opaque_identifier','original_matrix_column','authoritative_aliquot_identifier','canonical_case_identifier','canonical_tissue_label','excluded')
 s01_require(all(required%in%names(identity)),'identity schema missing')
 s01_validate_identity_vector(identity$opaque_identifier,'opaque sample')
 s01_validate_identity_vector(identity$original_matrix_column,'source sample')
 s01_validate_identity_vector(identity$authoritative_aliquot_identifier,'aliquot')
 s01_validate_identity_vector(patients,'patient')
 s01_assert_unique_key(identity,c('canonical_case_identifier','canonical_tissue_label'),'patient tissue')
 s01_require(!anyNA(identity$canonical_case_identifier)&&all(nzchar(identity$canonical_case_identifier))&&all(identity$canonical_tissue_label%in%c('Tumor','NAT')),'invalid patient or tissue')
 key<-paste(identity$canonical_case_identifier,identity$canonical_tissue_label,sep='::')
 i<-match(paste(rep(patients,each=2),rep(c('Tumor','NAT'),length(patients)),sep='::'),key)
 s01_require(!anyNA(i),'missing complete pair identity')
 md<-identity[i,];s01_require(!anyNA(md$excluded)&&!any(md$excluded),'excluded sample re-entered')
 data.frame(layer=layer,patient_id=md$canonical_case_identifier,tissue=md$canonical_tissue_label,
  sample_id=md$opaque_identifier,source_record_id=md$original_matrix_column,aliquot_id=md$authoritative_aliquot_identifier,stringsAsFactors=FALSE)
}

s07_validate_prepared <- function(x) {
 cfg<-s07_cfg();identities<-s07_read('data/manifests/s04_pre_effect/sample_identity.tsv');features<-s07_read('data/manifests/s04_pre_effect/feature_identity.tsv')
 ledgers<-s07_read('data/manifests/s04_pre_effect/pair_quadruplet_membership.tsv');coverage<-s07_read('data/manifests/s04_pre_effect/endpoint_coverage.tsv')
 backgrounds<-s07_read('data/manifests/s04_pre_effect/ordered_backgrounds.tsv');summary<-s07_read('data/manifests/s04_pre_effect/background_identity_summary.tsv')
 for(layer in c('protein','RNA')) {
  a<-x[[if(layer=='RNA')'rna'else layer]];ids<-identities[identities$layer==layer,setdiff(names(identities),'layer')];s07_compare(a$identity,ids,paste(layer,'identity'))
  f<-features[features$layer==layer,];s01_require(identical(f$gene_id,rownames(a$values))&&identical(f$source_row,seq_len(nrow(a$values))),'feature identities/order changed')
  s01_validate_identity_vector(rownames(a$values),'feature');s01_validate_identity_vector(colnames(a$values),'matrix column')
  s01_require(identical(colnames(a$values),a$identity$opaque_identifier),'matrix column identities changed')
 }
 s01_require(identical(dim(x$protein$values),c(11294L,189L))&&identical(dim(x$rna$values),c(62696L,172L)),'matrix dimensions changed')
 s01_require(sum(x$protein$identity$canonical_tissue_label=='Tumor')==159L&&sum(x$protein$identity$canonical_tissue_label=='NAT')==30L,'protein tissue universe changed')
 s01_require(sum(!is.finite(x$protein$values))==254539L&&all(is.finite(x$rna$values)),'matrix finite masks changed')
 s01_require(!'ALOX12'%in%rownames(x$protein$values),'unexpected ALOX12 presence')
 settings<-read_yaml('config/implementation/s04_settings.yml')
 s01_require(x$protein$sha256==settings$protein_sha256&&x$rna$sha256==settings$rna_sha256,'matrix provenance hashes changed')
 ex<-s07_read('data/manifests/s04_pre_effect/publication_exclusions.tsv')
 s01_require(!any(ex$authoritative_sample_identifier%in%x$protein$identity$original_matrix_column),'publication exclusion re-entered')
 p<-list();audit<-list()
 for(id in c(names(x$ordinary),names(x$bridge))) {
  bridge<-startsWith(id,'bridge_');z<-if(bridge)x$bridge[[id]]else x$ordinary[[id]];members<-s07_members(id,cfg)
  s01_require(identical(members,z$members),'exact ordered membership mismatch')
  l<-ledgers[ledgers$endpoint_id==id,];s07_compare(z$ledger,l,paste(id,'whole eligibility ledger'))
  s01_validate_identity_vector(l$patient_id,'ledger patient')
  # Reconstruct complete eligibility across every patient in the frozen ledger.
  eligible<-vapply(l$patient_id,function(patient) {
   all(vapply(if(bridge)c('rna','protein')else 'protein',function(layer) {
    a<-x[[layer]];d<-a$identity[a$identity$canonical_case_identifier==patient,]
    nrow(d)==2L&&setequal(d$canonical_tissue_label,c('Tumor','NAT'))&&!any(d$excluded)&&all(is.finite(a$values[members,d$opaque_identifier,drop=FALSE]))
   },logical(1)))
  },logical(1))
  s01_require(identical(unname(eligible),l$complete_case_status=='complete'),'complete pair/quadruplet finite-mask eligibility mismatch')
  patients<-l$patient_id[eligible];s01_require(identical(patients,sort(patients,method='radix')),'patient order changed')
  cv<-coverage[coverage$endpoint_id==id,];s01_require(length(patients)==z$complete_n&&z$complete_n==cv$complete_n&&paste(members,collapse=';')==cv$members,'coverage reconciliation failed')
  pp<-list(population_id=id,members=members,patients=patients,metadata=list(),backgrounds=list())
  for(layer in if(bridge)c('RNA','protein')else 'protein') {
   a<-x[[if(layer=='RNA')'rna'else layer]];md<-s07_pair_index(a$identity,patients,layer)
   # Source record IDs are joined independently to each frozen tissue field.
   for(tissue in c('Tumor','NAT')) {
    field<-paste0(if(layer=='RNA')'rna'else 'protein','_',tolower(tissue),'_record_id')
    s01_require(identical(md$source_record_id[md$tissue==tissue],l[[field]][eligible]),'pair/quadruplet tissue record mismatch')
   }
   frozen_bg<-if(bridge)z[[paste0(if(layer=='RNA')'rna'else 'protein','_background')]]else z$background
   b<-backgrounds[backgrounds$endpoint_id==id&backgrounds$layer==layer,]
   direct<-sort(rownames(a$values)[rowSums(is.finite(a$values[,md$sample_id,drop=FALSE]))==nrow(md)],method='radix')
   s01_require(identical(frozen_bg,b$gene_id)&&identical(b$position,seq_along(frozen_bg))&&identical(direct,frozen_bg),'endpoint-specific background identity/order mismatch')
   ss<-summary[summary$endpoint_id==id&summary$layer==layer,]
   s01_require(nrow(ss)==1&&s07_lf_hash(frozen_bg)==ss$ordered_LF_sha256&&s07_lf_hash(patients)==ss$patient_order_sha256&&length(frozen_bg)==ss$background_n,'ordered background/patient hash mismatch')
   pp$metadata[[layer]]<-md;pp$backgrounds[[layer]]<-frozen_bg
   audit[[length(audit)+1L]]<-data.frame(population_id=id,layer=layer,members=paste(members,collapse=';'),biological_n=length(patients),background_n=length(frozen_bg),patient_order_sha256=s07_lf_hash(patients),background_order_sha256=s07_lf_hash(frozen_bg),sample_order_sha256=s07_lf_hash(md$sample_id),identity_status='EXACT_FROZEN_MATCH')
  }
  if(!bridge)s01_require(identical(z$retained_cases,patients)&&identical(z$retained_opaque,pp$metadata$protein$sample_id),'ordinary retained order mismatch')
  else {
   s01_require(identical(vapply(z$retained,`[[`,'','case'),patients),'bridge retained patient mismatch')
   for(layer in c('rna','protein'))s01_require(identical(unlist(lapply(z$retained,`[[`,layer),use.names=FALSE),pp$metadata[[if(layer=='rna')'RNA'else layer]]$sample_id),'bridge retained sample mismatch')
  }
  p[[id]]<-pp
 }
 e<-s07_crosswalk();g<-s07_read('data/manifests/s04_pre_effect/global_hypothesis_eligibility.tsv');g<-g[g$contract%in%c('C3','C4'),]
 s07_compare(e,g[,names(e)],'global crosswalk');s01_require(all(g$structural_eligibility=='ELIGIBLE'),'frozen endpoint eligibility changed')
 for(i in seq_len(nrow(e)))s01_require(length(p[[s07_population(e$endpoint_id[i])]]$patients)==g$biological_n[i],'hypothesis biological n mismatch')
 list(populations=p,audit=do.call(rbind,audit))
}

s07_inputs <- function() {
 x<-read_rds(s07_path('data/processed/s04_pre_effect/unscored_inputs.rds'));v<-s07_validate_prepared(x)
 list(raw=x,populations=v$populations,audit=v$audit)
}

s07_scope <- function() {
 e<-s07_crosswalk();cfg<-s07_cfg();rows<-list()
 add<-function(endpoint,pop,category,construct,members,parent='NA_NOT_APPLICABLE',deleted='NA_NOT_APPLICABLE',operation='NA_NOT_APPLICABLE') {
  family<-if(category=='registered_hypothesis')e$multiplicity_family[match(endpoint,e$endpoint_id)]else 'NO_FORMAL_FAMILY'
  rows[[length(rows)+1L]]<<-data.frame(analysis_id=paste(endpoint,category,if(deleted!='NA_NOT_APPLICABLE')deleted else construct,sep='__'),endpoint_id=endpoint,population_reference_id=pop,category=category,construct_id=construct,members=paste(members,collapse=';'),parent_members=parent,deleted_member=deleted,multiplicity_family=family,operation_name=operation)
 }
 for(i in seq_len(nrow(e))) {
  r<-e[i,];pop<-s07_population(r$endpoint_id);m<-s07_members(pop,cfg)
  add(r$endpoint_id,pop,'registered_hypothesis',r$construct,if(grepl('_member_',r$endpoint_id))r$construct else m)
  if(r$contract=='C4') {add(r$endpoint_id,pop,'bridge_bootstrap_BCa_sensitivity',r$construct,m,operation=paste0(r$endpoint_id,'__bootstrap'));next}
  if(grepl('_member_',r$endpoint_id))next
  for(cat in c('standardized_member_mean','paired_Wilcoxon','paired_percentile_bootstrap','score_method_direction'))add(r$endpoint_id,pop,cat,r$construct,m,operation=if(cat=='paired_percentile_bootstrap')paste0(r$endpoint_id,'__paired_percentile_bootstrap')else 'NA_NOT_APPLICABLE')
  if(r$multiplicity_family=='PROTEIN_CORE') {
   for(deleted in m)add(r$endpoint_id,pop,'leave_one_gene_out',r$construct,m[m!=deleted],paste(m,collapse=';'),deleted)
   add(r$endpoint_id,pop,'member_direction_concordance',r$construct,cfg$programs[[paste0(if(grepl('_G$',r$endpoint_id))'G'else 'H','_core')]]$members)
   if(grepl('_G$',r$endpoint_id))for(submodule in c('substrate_assembly','oxidation_support'))add(paste0('PDC000614_G_',submodule),pop,'G_submodule',paste0('G_',submodule),cfg$coverage_contracts$G_submodules[[submodule]])
  }
 }
 out<-do.call(rbind,rows);s01_validate_identity_vector(out$analysis_id,'S07 scope analysis');out$scope_order<-seq_len(nrow(out))
 out$timing<-'PRESPECIFIED_RULE_IMPLEMENTATION_MATERIALIZED_AFTER_PRIOR_EFFECT_ACCESS'
 out
}

s07_capture <- function(expr) {
 warnings<-character();value<-tryCatch(withCallingHandlers(force(expr),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 list(value=value,warnings=warnings,error=if(inherits(value,'error'))conditionMessage(value)else 'NA_NOT_APPLICABLE')
}

s07_empty <- function(e) {
 for(k in c('effect','ci_lower','ci_upper','p_value','q_value','standard_error','statistic','model_degrees_of_freedom','numerator_df','denominator_df','g_av','g_av_variance','g_av_ci_lower','g_av_ci_upper','raw_member_paired_mean_difference'))e[[k]]<-NA_real_
 e$biological_n<-0L;e$member<-if(grepl('_member_',e$endpoint_id))e$construct else 'NA_NOT_APPLICABLE'
 e$estimability<-e$standard_error_status<-e$confidence_interval_status<-e$g_av_state<-'NA_NONESTIMABLE'
 e$missing_reason<-e$uncertainty_reason<-e$g_av_reason<-'NOT_EXECUTED'
 e$ql_df_status<-'NA_NOT_APPLICABLE';e$native_direction<-'NA_NONESTIMABLE';e$multiplicity_status<-'FROZEN_FAMILY_BH'
 e$timing<-'PRESPECIFIED_BEFORE_DISEASE_EFFECT_ACCESS'
 e$execution_timing<-'S07_AFTER_S05_AND_S06_DISEASE_EFFECT_ACCESS'
 e$dependence<-e$evidence_role<-if(e$contract=='C3')'paired_human_protein_orthogonal_evidence'else 'prespecified_exploratory_same_cohort_rna_protein_bridge'
 e$common_pair_set_id<-s07_population(e$endpoint_id);e$provenance_id<-s07_m('ENTRY_AUTHORIZATION.json');e
}

s07_score <- function(input,pop,layer,members=pop$members,method='primary_rank') {
 md<-pop$metadata[[layer]];x<-input$raw[[if(layer=='RNA')'rna'else layer]]$values[,md$sample_id,drop=FALSE]
 if(method=='standardized_member_mean')s01_standardized_member_mean(x,members)else s01_rank_program_score(x,members,pop$backgrounds[[layer]])
}

s07_estimate <- function(score,md) {
 patients<-unique(md$patient_id);s07_pair_index(data.frame(opaque_identifier=md$sample_id,original_matrix_column=md$source_record_id,authoritative_aliquot_identifier=md$aliquot_id,canonical_case_identifier=md$patient_id,canonical_tissue_label=md$tissue,excluded=FALSE),patients,md$layer[1])
 if(!isTRUE(score$estimable))return(list(estimable=FALSE,state='NA_NONESTIMABLE',reason=score$reason,biological_n=as.integer(length(patients)),effect=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,p_value=NA_real_,statistic=NA_real_,standard_error=NA_real_,model_degrees_of_freedom=NA_real_,g_av=NA_real_,g_av_variance=NA_real_,g_av_ci_lower=NA_real_,g_av_ci_upper=NA_real_,g_av_state='NA_NONESTIMABLE',g_av_reason=score$reason))
 s01_require(identical(score$sample_id,md$sample_id),'score/sample order mismatch')
 t<-score$score[match(paste(patients,'Tumor'),paste(md$patient_id,md$tissue))];n<-score$score[match(paste(patients,'NAT'),paste(md$patient_id,md$tissue))]
 a<-s01_paired_delta_t(t,n,patients);g<-s01_paired_g_av(t,n,patients)
 a$standard_error<-if(a$estimable)sd(t-n)/sqrt(length(patients))else NA_real_
 a$model_degrees_of_freedom<-if(a$estimable)as.double(length(patients)-1)else NA_real_
 a$g_av<-g$g_av;a$g_av_variance<-g$variance;a$g_av_ci_lower<-g$ci_lower;a$g_av_ci_upper<-g$ci_upper;a$g_av_state<-g$state;a$g_av_reason<-if(is.null(g$reason))'NA_NOT_APPLICABLE'else g$reason
 a$source_data<-data.frame(patient_id=patients,Tumor_score=t,NAT_score=n,delta=t-n);a
}

s07_apply <- function(row,a) {
 for(k in intersect(names(row),names(a)))if(length(a[[k]])==1L)row[[k]]<-a[[k]]
 row$estimability<-a$state;row$missing_reason<-if(a$estimable)'NA_NOT_APPLICABLE'else a$reason
 row$standard_error_status<-row$confidence_interval_status<-if(a$estimable)'AVAILABLE'else 'NA_NONESTIMABLE'
 row$uncertainty_reason<-if(a$estimable)'NA_NOT_APPLICABLE'else a$reason;row$native_direction<-s01_direction_label(row$effect);row
}

s07_pair_bootstrap <- function(d,seed,cfg=s07_cfg()$diagnostic_bootstrap$paired_program) {
 s01_validate_paired_inputs(d$Tumor_score,d$NAT_score,d$patient_id);s01_require(identical(d$delta,d$Tumor_score-d$NAT_score),'bootstrap source delta mismatch')
 n<-nrow(d);draws<-rep(NA_real_,cfg$replicates);indices<-matrix(integer(),nrow=0,ncol=0)
 if(n>=3L&&all(is.finite(d$delta))) {
  indices<-s01_with_production_rng(seed,replicate(cfg$replicates,sample.int(n,n,replace=TRUE)))
  draws<-colMeans(matrix(d$delta[indices],nrow=n))
 }
 valid<-sum(is.finite(draws));ok<-valid>=cfg$minimum_valid_replicates
 ci<-if(ok)quantile(draws[is.finite(draws)],c(.025,.975),type=cfg$quantile_type,names=FALSE)else c(NA_real_,NA_real_)
 list(estimable=ok,state=if(ok)'ESTIMABLE'else 'NA_NONESTIMABLE',reason=if(ok)'NA_NOT_APPLICABLE'else 'FEWER_THAN_REQUIRED_FINITE_REPLICATES_OR_INVALID_INPUT',biological_n=as.integer(n),effect=if(all(is.finite(d$delta)))mean(d$delta)else NA_real_,ci_lower=ci[1],ci_upper=ci[2],valid_replicates=as.integer(valid),seed=seed,draws=draws,indices=indices,resampling_unit='complete_patient_pair',quantile_type=cfg$quantile_type)
}

s07_wilcoxon <- function(d) {
 if(any(!is.finite(d$delta))||all(d$delta==0))return(list(estimable=FALSE,state='NA_NONESTIMABLE',reason='NONFINITE_OR_ALL_ZERO_PAIRED_DIFFERENCES',statistic=NA_real_,p_value=NA_real_,warnings=character()))
 z<-s07_capture(stats::wilcox.test(d$Tumor_score,d$NAT_score,paired=TRUE,exact=FALSE,correct=FALSE,alternative='two.sided'))
 s01_require(!inherits(z$value,'error')&&!length(z$warnings),paste('Unexpected Wilcoxon failure',z$error,paste(z$warnings,collapse=';')))
 list(estimable=TRUE,state='ESTIMABLE',reason='NA_NOT_APPLICABLE',biological_n=nrow(d),effect=mean(d$delta),statistic=unname(z$value$statistic),p_value=z$value$p.value,warnings=z$warnings)
}

s07_limma <- function(input,pop) {
 md<-pop$metadata$protein;x<-input$raw$protein$values[pop$backgrounds$protein,md$sample_id,drop=FALSE]
 p<-s01_prepare_paired_model(x,md,'patient_id','patient')
 fit<-s07_capture(limma::eBayes(limma::lmFit(p$values,p$design)))
 s01_require(!inherits(fit$value,'error')&&!length(fit$warnings),paste('Unexpected protein limma failure; owner review:',fit$error,paste(fit$warnings,collapse=';')))
 extracted<-s01_extract_limma_contract(fit$value,p,'expression~tissue+patient','publication_TMT_median_centered_log2_ratio_no_renormalization')
 j<-match(pop$members,extracted$feature_id);s01_require(!anyNA(j),'member absent from broad fitted background')
 list(fit=fit$value,prepared=p,extracted=extracted,member_index=j,members=pop$members,metadata=md,warnings=fit$warnings,background=pop$backgrounds$protein)
}

s07_seed <- function(operation) {
 m<-s07_read('config/implementation/s01_target_seed_map.tsv');i<-which(m$operation_name==operation)
 s01_require(length(i)==1L&&m$rng_kind[i]=="L'Ecuyer-CMRG",'unregistered stochastic operation');as.integer(m$seed[i])
}

s07_analyze_protein <- function(input,scope) {
 cfg<-s07_cfg();e<-s07_crosswalk();e<-e[e$contract=='C3',];scores<-estimates<-models<-diagnostics<-list();rows<-list()
 for(popid in names(input$raw$ordinary)) {
  pop<-input$populations[[popid]];scores[[popid]]<-s07_score(input,pop,'protein');estimates[[popid]]<-s07_estimate(scores[[popid]],pop$metadata$protein)
  if(startsWith(popid,'protein_'))models[[popid]]<-s07_limma(input,pop)
 }
 for(i in seq_len(nrow(e))) {
  r<-s07_empty(e[i,]);pop<-input$populations[[r$common_pair_set_id]];r$biological_n<-as.integer(length(pop$patients));r$members<-paste(pop$members,collapse=';');r$background_n<-length(pop$backgrounds$protein);r$background_order_sha256<-s07_lf_hash(pop$backgrounds$protein);r$patient_order_sha256<-s07_lf_hash(pop$patients)
  if(grepl('_member_',r$endpoint_id)) {
   m<-models[[r$common_pair_set_id]];a<-m$extracted;j<-match(r$member,a$feature_id)
   for(k in c('effect','standard_error','statistic','p_value','ci_lower','ci_upper'))r[[k]]<-a[[k]][j]
   r$model_degrees_of_freedom<-a$degrees_of_freedom[j];r$estimability<-a$state[j]
   r$missing_reason<-if(r$estimability=='ESTIMABLE')'NA_NOT_APPLICABLE'else 'INVALID_NATIVE_LIMMA_OUTPUT_OWNER_REVIEW'
   r$standard_error_status<-r$confidence_interval_status<-if(r$estimability=='ESTIMABLE')'AVAILABLE'else 'NA_NONESTIMABLE';r$uncertainty_reason<-r$missing_reason
   md<-pop$metadata$protein;v<-input$raw$protein$values[r$member,md$sample_id];r$raw_member_paired_mean_difference<-mean(v[md$tissue=='Tumor']-v[md$tissue=='NAT']);r$native_direction<-s01_direction_label(r$raw_member_paired_mean_difference)
   r$g_av_state<-r$g_av_reason<-'NA_NOT_APPLICABLE'
  } else r<-s07_apply(r,estimates[[r$common_pair_set_id]])
  rows[[i]]<-r
 }
 for(i in which(scope$multiplicity_family=='NO_FORMAL_FAMILY'&!startsWith(scope$population_reference_id,'bridge_'))) {
  r<-scope[i,];pop<-input$populations[[r$population_reference_id]];base<-estimates[[r$population_reference_id]];members<-strsplit(r$members,';',fixed=TRUE)[[1]];score<-NULL
  if(r$category%in%c('standardized_member_mean','leave_one_gene_out','G_submodule')) {
   score<-s07_score(input,pop,'protein',members,if(r$category=='standardized_member_mean')'standardized_member_mean'else 'primary_rank')
   a<-s07_estimate(score,pop$metadata$protein)
  } else if(r$category=='paired_Wilcoxon')a<-s07_wilcoxon(base$source_data)
  else if(r$category=='paired_percentile_bootstrap')a<-s07_pair_bootstrap(base$source_data,s07_seed(r$operation_name))
  else if(r$category=='member_direction_concordance') {
   md<-pop$metadata$protein;m<-input$raw$protein$values[pop$members,md$sample_id,drop=FALSE]
   effects<-rowMeans(m[,md$tissue=='Tumor',drop=FALSE]-m[,md$tissue=='NAT',drop=FALSE])
   denom<-cfg$diagnostics$member_direction_concordance[[if(pop$population_id=='protein_G_program')'denominator_G'else 'denominator_H']]
   a<-s01_member_direction_concordance(base$effect,effects,denom);a$estimable<-a$state=='ESTIMABLE';a$reason<-if(a$estimable)'NA_NOT_APPLICABLE'else if(length(effects)!=denom)'ALOX12_ABSENT_REGISTERED_G_DENOMINATOR_SEVEN_UNCOMPUTABLE'else 'FULL_EFFECT_NOT_DIRECTIONAL_OR_NONFINITE'
   a$denominator<-denom;a$observed_member_n<-length(effects);a$observed_member_effects<-effects;a$biological_n<-length(pop$patients)
  } else if(r$category=='score_method_direction') {
   alt<-Filter(function(d)d$registry$population_reference_id==r$population_reference_id&&d$registry$category=='standardized_member_mean',diagnostics)[[1]]$estimate
   a<-s01_score_method_sensitivity(base$effect,alt$effect,isTRUE(alt$estimable));a$estimable<-!a$score_method_nonestimable;a$state<-if(a$estimable)'ESTIMABLE'else 'NA_NONESTIMABLE';a$reason<-if(a$estimable)'NA_NOT_APPLICABLE'else 'STANDARDIZED_MEMBER_MEAN_NONESTIMABLE';a$biological_n<-length(pop$patients)
  } else stop('Unsupported S07 protein diagnostic')
  diagnostics[[r$analysis_id]]<-list(registry=r,score=score,estimate=a,primary_effect=base$effect)
 }
 list(rows=do.call(rbind,rows),scores=scores,estimates=estimates,models=models,diagnostics=diagnostics)
}

s07_bridge_inference <- function(d,seed,cfg=s07_cfg()$bridge_inference) {
 s01_validate_identity_vector(d$patient_id,'bridge patient');s01_require(all(c('RNA_delta','protein_delta')%in%names(d)),'bridge delta schema')
 n<-nrow(d);input_valid<-all(is.finite(c(d$RNA_delta,d$protein_delta)))
 a<-list(biological_n=as.integer(n),rho=NA_real_,effect=NA_real_,p_value=NA_real_,nominal_p_before_interval_gate=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,bca_ci_lower=NA_real_,bca_ci_upper=NA_real_,bca_state='NA_NONESTIMABLE',bca_reason='PRIMARY_NONESTIMABLE',valid_bootstrap_replicates=0L,seed=seed,required_input_valid=input_valid,estimable=FALSE,state='NA_NONESTIMABLE',reason='INSUFFICIENT_COMPLETE_CASES_OR_INVALID_REQUIRED_INPUT',warnings=character(),bca_warnings=character(),boot=NULL)
 if(n<cfg$minimum_complete_cases||!input_valid)return(a)
 tt<-s07_capture(stats::cor.test(d$RNA_delta,d$protein_delta,method='spearman',exact=FALSE,alternative='two.sided'))
 a$warnings<-tt$warnings
 s01_require(!inherits(tt$value,'error')&&all(tt$warnings%in%'the standard deviation is zero'),paste('Unexpected Spearman failure',tt$error,paste(tt$warnings,collapse=';')))
 a$rho<-a$effect<-as.numeric(tt$value$estimate);a$nominal_p_before_interval_gate<-tt$value$p.value
 if(!is.finite(a$rho)){a$reason<-'NONFINITE_SPEARMAN_RHO';return(a)}
 bt<-s07_capture(s01_with_production_rng(seed,boot::boot(d,statistic=s07_boot_spearman,R=cfg$replicates,sim='ordinary',stype='i')))
 s01_require(!inherits(bt$value,'error')&&all(bt$warnings%in%'the standard deviation is zero'),paste('Unexpected bridge bootstrap failure',bt$error,paste(unique(bt$warnings),collapse=';')))
 a$warnings<-c(a$warnings,bt$warnings);a$boot<-bt$value;v<-as.numeric(bt$value$t);a$valid_bootstrap_replicates<-as.integer(sum(is.finite(v)))
 if(a$valid_bootstrap_replicates<cfg$minimum_valid_replicates){a$reason<-'FEWER_THAN_9500_FINITE_BOOTSTRAP_REPLICATES';return(a)}
 ci<-quantile(v[is.finite(v)],c(.025,.975),type=7,names=FALSE);a$ci_lower<-ci[1];a$ci_upper<-ci[2]
 a$estimable<-TRUE;a$state<-'ESTIMABLE';a$reason<-'NA_NOT_APPLICABLE';a$p_value<-tt$value$p.value
 bca<-s07_capture(boot::boot.ci(bt$value,conf=.95,type='bca'));a$bca_warnings<-bca$warnings
 if(!inherits(bca$value,'error')&&!is.null(bca$value$bca)&&all(is.finite(bca$value$bca[4:5]))) {
  a$bca_ci_lower<-as.numeric(bca$value$bca[4]);a$bca_ci_upper<-as.numeric(bca$value$bca[5]);a$bca_state<-'ESTIMABLE';a$bca_reason<-'NA_NOT_APPLICABLE'
 } else a$bca_reason<-if(bca$error!='NA_NOT_APPLICABLE')bca$error else 'NATIVE_BOOT_BCA_INTERVAL_UNAVAILABLE_OR_NONFINITE'
 a
}

s07_bridge <- function(input,axis) {
 pop<-input$populations[[paste0('bridge_',axis)]];e<-s07_crosswalk();r<-s07_empty(e[e$endpoint_id==paste0('CPTAC_bridge_',axis),])
 scores<-lapply(c('RNA','protein'),function(layer)s07_score(input,pop,layer));names(scores)<-c('RNA','protein')
 s01_require(all(vapply(scores,function(s)isTRUE(s$estimable),logical(1))),'frozen bridge score inputs invalid')
 deltas<-lapply(names(scores),function(layer) {
  md<-pop$metadata[[layer]];s<-scores[[layer]];s01_require(identical(s$sample_id,md$sample_id),'bridge score identity mismatch')
  data.frame(patient_id=pop$patients,Tumor_score=s$score[md$tissue=='Tumor'],NAT_score=s$score[md$tissue=='NAT'])
 });names(deltas)<-names(scores)
 d<-data.frame(patient_id=pop$patients,RNA_Tumor_score=deltas$RNA$Tumor_score,RNA_NAT_score=deltas$RNA$NAT_score,protein_Tumor_score=deltas$protein$Tumor_score,protein_NAT_score=deltas$protein$NAT_score)
 d$RNA_delta<-d$RNA_Tumor_score-d$RNA_NAT_score;d$protein_delta<-d$protein_Tumor_score-d$protein_NAT_score
 a<-s07_bridge_inference(d,s07_seed(paste0(r$endpoint_id,'__bootstrap')))
 for(k in c('effect','biological_n','p_value','ci_lower','ci_upper'))r[[k]]<-a[[k]]
 for(k in c('rho','bca_ci_lower','bca_ci_upper','bca_state','bca_reason','valid_bootstrap_replicates','nominal_p_before_interval_gate','required_input_valid'))r[[k]]<-a[[k]]
 r$bootstrap_ci_lower<-r$ci_lower;r$bootstrap_ci_upper<-r$ci_upper;r$estimability<-a$state;r$missing_reason<-a$reason
 r$members<-paste(pop$members,collapse=';');r$seed<-a$seed;r$RNA_background_n<-length(pop$backgrounds$RNA);r$protein_background_n<-length(pop$backgrounds$protein)
 r$RNA_background_order_sha256<-s07_lf_hash(pop$backgrounds$RNA);r$protein_background_order_sha256<-s07_lf_hash(pop$backgrounds$protein);r$patient_order_sha256<-s07_lf_hash(pop$patients)
 r$standard_error_status<-'NA_NOT_APPLICABLE';r$confidence_interval_status<-if(a$estimable)'AVAILABLE'else 'NA_NONESTIMABLE';r$uncertainty_reason<-'PERCENTILE_BOOTSTRAP_NO_COEFFICIENT_STANDARD_ERROR';r$g_av_state<-r$g_av_reason<-'NA_NOT_APPLICABLE'
 r$native_direction<-s01_direction_label(a$rho)
 list(row=r,estimate=a,scores=scores,deltas=d,population_id=pop$population_id)
}

s07_adjust <- function(rows,cfg=s07_cfg()) {
 for(f in unique(rows$multiplicity_family)) {
  j<-which(rows$multiplicity_family==f);spec<-cfg$multiplicity$families[[f]]
  s01_require(!is.null(spec)&&length(j)==spec$frozen_size&&all(rows$frozen_family_size[j]==spec$frozen_size)&&identical(rows$family_hypothesis[j],spec$hypotheses),'fixed family membership/order mismatch')
  rows$q_value[j]<-unname(s01_frozen_bh(setNames(rows$p_value[j],rows$family_hypothesis[j]),f,cfg))
 };rows
}

s07_diagnostic_rows <- function(protein) {
 do.call(rbind,lapply(protein$diagnostics,function(d) {
  r<-d$registry;a<-d$estimate
  for(k in c('effect','ci_lower','ci_upper','p_value','statistic','standard_error','model_degrees_of_freedom','biological_n','g_av','g_av_variance','g_av_ci_lower','g_av_ci_upper','valid_replicates','seed','concordant_members','denominator','proportion','observed_member_n'))r[[k]]<-if(is.null(a[[k]]))NA_real_ else a[[k]]
  r$diagnostic_id<-r$analysis_id;r$diagnostic_state<-a$state;r$missing_reason<-if(is.null(a$reason))'NA_NOT_APPLICABLE'else a$reason;r$q_value<-NA_real_;r$multiplicity_status<-'NO_FORMAL_FAMILY';r$role<-r$dependence<-'same_data_robustness';r$provenance_id<-s07_m('ENTRY_AUTHORIZATION.json')
  r$primary_effect<-d$primary_effect;r$primary_direction<-s01_direction_label(d$primary_effect);r$alternative_direction<-s01_direction_label(r$effect)
  r$primary_exact_zero<-is.finite(r$primary_effect)&&r$primary_effect==0;r$alternative_exact_zero<-if(is.finite(r$effect))r$effect==0 else NA
  r$loo_direction_stable<-if(r$category=='leave_one_gene_out')s01_loo_direction_status(d$primary_effect,r$effect,isTRUE(a$estimable))else 'NA_NOT_APPLICABLE'
  r$g_av_state<-if(is.null(a$g_av_state))'NA_NOT_APPLICABLE'else a$g_av_state;r$g_av_reason<-if(is.null(a$g_av_reason))'NA_NOT_APPLICABLE'else a$g_av_reason
  for(k in c('primary_rank_unrounded_raw_paired_mean_difference','standardized_member_mean_unrounded_raw_paired_mean_difference'))r[[k]]<-if(is.null(a[[k]]))NA_real_ else a[[k]]
  for(k in c('primary_rank_direction','standardized_member_mean_direction','score_method_sensitive'))r[[k]]<-if(is.null(a[[k]]))'NA_NOT_APPLICABLE'else a[[k]]
  for(k in c('score_method_nonestimable','exact_zero_direction_flag'))r[[k]]<-if(is.null(a[[k]]))NA else a[[k]]
  r
 }))
}

s07_assemble <- function(protein,bridges) {
 p<-s07_adjust(protein$rows);b<-s07_adjust(do.call(rbind,lapply(bridges,`[[`,'row')))
 for(i in seq_len(nrow(b))) {
  z<-s01_bridge_state(b$biological_n[i],b$rho[i],b$p_value[i],b$ci_lower[i],b$ci_upper[i],b$q_value[i],b$valid_bootstrap_replicates[i],b$required_input_valid[i])
  b$bridge_state[i]<-z$state;b$bridge_q_ci_inconsistent[i]<-if(is.na(z$bridge_q_ci_inconsistent))'NA_NONESTIMABLE'else as.character(z$bridge_q_ci_inconsistent);b$bridge_nominal_p_ci_inconsistent[i]<-if(is.na(z$bridge_nominal_p_ci_inconsistent))'NA_NONESTIMABLE'else as.character(z$bridge_nominal_p_ci_inconsistent)
 }
 list(protein=p,bridge=b,diagnostics=s07_diagnostic_rows(protein))
}

s07_run_target <- function(kind,input_path,scope,axis=NULL) {
 s07_guard();resource<-s07_resource();start<-s07_utc();input<-read_rds(input_path)
 a<-if(kind=='protein')s07_analyze_protein(input,scope)else s07_bridge(input,axis)
 a$execution<-list(start_utc=start,end_utc=s07_utc(),resource=resource,workers=1L,threads=1L)
 s07_save(a,s07_o(if(kind=='protein')'protein_source.rds'else paste0('bridge_',axis,'_source.rds')))
}

s07_tables <- function(result,protein,bridges,input) {
 rows<-result$protein;member<-rows$multiplicity_family=='PROTEIN_CORE_MEMBERS';scores<-values<-deltas<-diagnostic_deltas<-bootstrap<-list()
 add_source<-function(id,pop,layer,sc,method,members=pop$members) {
  if(!isTRUE(sc$estimable))return(invisible(NULL))
  md<-pop$metadata[[layer]];bgh<-s07_lf_hash(pop$backgrounds[[layer]])
  scores[[length(scores)+1L]]<<-data.frame(analysis_id=id,population_reference_id=pop$population_id,score_method=method,members=paste(members,collapse=';'),background_n=length(pop$backgrounds[[layer]]),background_order_sha256=bgh,md,score=as.numeric(sc$score))
 }
 for(id in names(protein$scores)) {
  pop<-input$populations[[id]];endpoint<-rows$endpoint_id[rows$common_pair_set_id==id&rows$member=='NA_NOT_APPLICABLE']
  add_source(endpoint,pop,'protein',protein$scores[[id]],'primary_rank')
  deltas[[id]]<-cbind(endpoint_id=endpoint,population_reference_id=id,protein$estimates[[id]]$source_data)
 }
 for(d in protein$diagnostics) {
  pop<-input$populations[[d$registry$population_reference_id]]
  if(!is.null(d$score))add_source(d$registry$analysis_id,pop,'protein',d$score,d$registry$category,strsplit(d$registry$members,';',fixed=TRUE)[[1]])
  if(!is.null(d$estimate$source_data))diagnostic_deltas[[length(diagnostic_deltas)+1L]]<-cbind(analysis_id=d$registry$analysis_id,population_reference_id=pop$population_id,d$estimate$source_data)
  if(!is.null(d$estimate$draws))bootstrap[[length(bootstrap)+1L]]<-data.frame(operation_name=d$registry$operation_name,replicate=seq_along(d$estimate$draws),value=d$estimate$draws,finite=is.finite(d$estimate$draws),seed=d$estimate$seed)
 }
 for(b in bridges) {
  pop<-input$populations[[b$population_id]]
  for(layer in c('RNA','protein'))add_source(b$row$endpoint_id,pop,layer,b$scores[[layer]],'primary_rank')
  if(!is.null(b$estimate$boot)){v<-as.numeric(b$estimate$boot$t);bootstrap[[length(bootstrap)+1L]]<-data.frame(operation_name=paste0(b$row$endpoint_id,'__bootstrap'),replicate=seq_along(v),value=v,finite=is.finite(v),seed=b$estimate$seed)}
 }
 for(pop in input$populations)for(layer in names(pop$metadata)) {
  md<-pop$metadata[[layer]];x<-input$raw[[if(layer=='RNA')'rna'else layer]]$values
  for(g in pop$members)values[[length(values)+1L]]<-data.frame(population_reference_id=pop$population_id,member=g,background_order_sha256=s07_lf_hash(pop$backgrounds[[layer]]),md,value=as.numeric(x[g,md$sample_id]),input_scale=if(layer=='RNA')'UQ_normalized_log2_RSEM_plus_1'else 'TMT_median_centered_log2_ratio')
 }
 bridge_deltas<-do.call(rbind,lapply(bridges,function(b)cbind(endpoint_id=b$row$endpoint_id,population_reference_id=b$population_id,b$deltas)))
 for(layer in c('RNA','protein'))for(tissue in c('Tumor','NAT')) {
  bridge_deltas[[paste(layer,tissue,'record_id',sep='_')]]<-unlist(lapply(bridges,function(b){md<-input$populations[[b$population_id]]$metadata[[layer]];md$source_record_id[md$tissue==tissue]}),use.names=FALSE)
 }
 tables<-list(protein_estimates=rows,cptac_bridge_estimates=result$bridge,all_program_estimates=rows[!member,],all_member_estimates=rows[member,],nominal_protein_diagnostics=result$diagnostics,score_method_diagnostics=result$diagnostics[result$diagnostics$category=='score_method_direction',],patient_program_scores=do.call(rbind,scores),patient_member_values=do.call(rbind,values),protein_paired_deltas=do.call(rbind,deltas),diagnostic_paired_deltas=do.call(rbind,diagnostic_deltas),bridge_paired_deltas=bridge_deltas,bootstrap_replicates=do.call(rbind,bootstrap))
 for(k in names(tables))rownames(tables[[k]])<-NULL
 tables
}

s07_close <- function(a,b,label,tolerance=1e-12) {
 # Scientific identity is checked separately; names must not alter numeric NA masks.
 a<-as.numeric(a);b<-as.numeric(b)
 s01_require(length(a)==length(b)&&identical(is.na(a),is.na(b))&&all(is.finite(a[!is.na(a)]))&&all(is.finite(b[!is.na(b)]))&&all(abs(a[!is.na(a)]-b[!is.na(b)])<=tolerance*pmax(1,abs(b[!is.na(b)]))),paste(label,'independent numeric anchor mismatch'));invisible(TRUE)
}

s07_anchor_pair <- function(a,d,label) {
 s01_require(identical(d$delta,d$Tumor_score-d$NAT_score)&&a$biological_n==nrow(d),'paired delta or n changed')
 if(a$estimable) {
  z<-stats::t.test(d$delta,mu=0,alternative='two.sided');s07_close(c(a$effect,a$statistic,a$p_value,a$ci_lower,a$ci_upper),c(mean(d$delta),as.numeric(z$statistic),z$p.value,as.numeric(z$conf.int)),label)
 }
 if(a$g_av_state=='ESTIMABLE') {
  z<-metafor::escalc(measure='SMCRP',m1i=mean(d$Tumor_score),m2i=mean(d$NAT_score),sd1i=sd(d$Tumor_score),sd2i=sd(d$NAT_score),ni=nrow(d),ri=cor(d$Tumor_score,d$NAT_score),vtype='LS',correct=TRUE)
  s07_close(c(a$g_av,a$g_av_variance,a$g_av_ci_lower,a$g_av_ci_upper),c(z$yi,z$vi,z$yi-qnorm(.975)*sqrt(z$vi),z$yi+qnorm(.975)*sqrt(z$vi)),paste(label,'metafor'))
 };invisible(TRUE)
}

s07_direct_score <- function(input,pop,layer,members,method) {
 md<-pop$metadata[[layer]];x<-input$raw[[if(layer=='RNA')'rna'else layer]]$values
 if(!all(members%in%rownames(x)))return(NULL)
 if(method=='standardized_member_mean') {
  y<-x[members,md$sample_id,drop=FALSE];if(any(apply(y,1,sd)<=apply(y,1,function(v)sqrt(.Machine$double.eps)*max(1,abs(v)))))return(NULL)
  return(as.numeric(colMeans(t(scale(t(y),center=TRUE,scale=TRUE)))))
 }
 y<-x[pop$backgrounds[[layer]],md$sample_id,drop=FALSE]
 vapply(seq_len(ncol(y)),function(j)sum((rank(y[,j],ties.method='average')[match(members,rownames(y))]-1)/(nrow(y)-1))/length(members),0)
}

s07_validate_models <- function(result,protein,bridges,input,scope) {
 e<-s07_crosswalk();s07_compare(e[e$contract=='C3',],result$protein[,names(e)],'C3 registered order and metadata');s07_compare(e[e$contract=='C4',],result$bridge[,names(e)],'C4 registered order and metadata')
 s07_compare(result$protein,s07_adjust(protein$rows),'source to assembled C3')
 recon<-s07_assemble(protein,bridges);s07_compare(result$bridge,recon$bridge,'source to assembled C4');s07_compare(result$diagnostics,recon$diagnostics,'source to assembled diagnostics')
 s01_require(identical(result$diagnostics$analysis_id,scope$analysis_id[scope$multiplicity_family=='NO_FORMAL_FAMILY'&!startsWith(scope$population_reference_id,'bridge_')]),'complete diagnostic scope/order mismatch')
 for(id in names(protein$estimates)) {
  pop<-input$populations[[id]];direct<-s07_direct_score(input,pop,'protein',s07_members(id),'primary_rank')
  s07_close(protein$scores[[id]]$score,direct,paste(id,'literal rank'));a<-protein$estimates[[id]];md<-pop$metadata$protein
  s07_close(a$source_data$Tumor_score,direct[md$tissue=='Tumor'],paste(id,'tumor score source'));s07_close(a$source_data$NAT_score,direct[md$tissue=='NAT'],paste(id,'NAT score source'));s07_anchor_pair(a,a$source_data,id)
 }
 for(d in protein$diagnostics) {
  r<-d$registry;a<-d$estimate;pop<-input$populations[[r$population_reference_id]]
  if(r$category=='leave_one_gene_out') {
   s01_require(identical(r$parent_members,paste(pop$members,collapse=';'))&&identical(r$members,paste(pop$members[pop$members!=r$deleted_member],collapse=';')),'LOO exact scientific identity corrupted')
  }
  if(r$category=='G_submodule')s01_require(identical(r$members,paste(s07_cfg()$coverage_contracts$G_submodules[[sub('G_','',r$construct_id,fixed=TRUE)]],collapse=';')),'submodule identity corrupted')
  if(!is.null(d$score)) {
   direct<-s07_direct_score(input,pop,'protein',strsplit(r$members,';',fixed=TRUE)[[1]],r$category)
   if(is.null(direct))s01_require(!a$estimable,'missing-member or zero variance diagnostic substituted')else {
    s07_close(d$score$score,direct,paste(r$analysis_id,'independent score'));s07_anchor_pair(a,a$source_data,r$analysis_id)
   }
  }
  if(r$category=='paired_percentile_bootstrap') {
   dd<-protein$estimates[[pop$population_id]]$source_data;n<-nrow(dd)
   idx<-s01_with_production_rng(a$seed,replicate(s07_cfg()$diagnostic_bootstrap$paired_program$replicates,sample.int(n,n,replace=TRUE)))
   s01_require(identical(idx,a$indices),'complete-pair bootstrap indices changed')
   v<-apply(idx,2,function(i)mean(dd$Tumor_score[i]-dd$NAT_score[i]));s07_close(a$draws,v,'pair bootstrap draws');s01_require(a$valid_replicates==sum(is.finite(v)),'pair valid count')
   s07_close(c(a$ci_lower,a$ci_upper),quantile(v[is.finite(v)],c(.025,.975),type=7,names=FALSE),'pair bootstrap type7')
   anchor<-s01_pair_bootstrap(dd$Tumor_score,dd$NAT_score,dd$patient_id,a$seed);s07_close(c(a$ci_lower,a$ci_upper),c(anchor$ci_lower,anchor$ci_upper),'accepted pair helper')
  }
  if(r$category=='paired_Wilcoxon') {
   dd<-protein$estimates[[pop$population_id]]$source_data;z<-stats::wilcox.test(dd$Tumor_score,dd$NAT_score,paired=TRUE,exact=FALSE,correct=FALSE,alternative='two.sided');s07_close(c(a$statistic,a$p_value),c(as.numeric(z$statistic),z$p.value),'native Wilcoxon')
  }
 }
 for(id in names(protein$models)) {
  m<-protein$models[[id]];pop<-input$populations[[id]];x<-input$raw$protein$values[pop$backgrounds$protein,pop$metadata$protein$sample_id,drop=FALSE]
  s01_require(identical(x,m$prepared$values)&&identical(m$background,pop$backgrounds$protein),'broad limma background changed')
  md<-pop$metadata$protein;design<-model.matrix(~factor(tissue,levels=c('NAT','Tumor'))+factor(patient_id),md)
  s01_require(qr(design)$rank==length(pop$patients)+1L,'paired limma rank')
  fit<-limma::eBayes(limma::lmFit(x,design));tab<-limma::topTable(fit,coef=2,number=Inf,sort.by='none',confint=.95)
  for(g in pop$members) {
   j<-match(g,m$extracted$feature_id);r<-result$protein[result$protein$common_pair_set_id==id&result$protein$member==g,]
   s07_close(c(r$effect,r$statistic,r$p_value,r$ci_lower,r$ci_upper),as.numeric(tab[g,c('logFC','t','P.Value','CI.L','CI.R')]),paste(g,'native limma topTable'))
   z<-lm(as.numeric(x[g,])~factor(tissue,levels=c('NAT','Tumor'))+factor(patient_id),data=md)
   s07_close(r$effect,as.numeric(coef(z)[2]),paste(g,'independent paired lm coefficient'))
   s07_close(r$effect,mean(x[g,md$tissue=='Tumor']-x[g,md$tissue=='NAT']),paste(g,'raw paired member mean'))
   s07_close(r$standard_error,fit$stdev.unscaled[g,2]*sqrt(fit$s2.post[j]),paste(g,'moderated SE'))
   s07_close(r$model_degrees_of_freedom,fit$df.total[j],paste(g,'moderated df'))
  }
 }
 for(b in bridges) {
  pop<-input$populations[[b$population_id]];a<-b$estimate;d<-b$deltas
  s01_require(identical(d$patient_id,pop$patients),'quadruplet patient ordering')
  for(layer in c('RNA','protein')) {
   direct<-s07_direct_score(input,pop,layer,pop$members,'primary_rank');s07_close(b$scores[[layer]]$score,direct,'bridge independent layer rank')
   md<-pop$metadata[[layer]];s07_close(d[[paste0(layer,'_delta')]],direct[md$tissue=='Tumor']-direct[md$tissue=='NAT'],'bridge paired delta')
  }
  tt<-stats::cor.test(d$RNA_delta,d$protein_delta,method='spearman',exact=FALSE,alternative='two.sided')
  s07_close(c(a$rho,a$nominal_p_before_interval_gate),c(as.numeric(tt$estimate),tt$p.value),'native bridge Spearman')
  if(!is.null(a$boot)) {
   idx<-s01_with_production_rng(a$seed,boot::boot.array(a$boot,indices=TRUE))
   v<-apply(idx,1,function(i){x<-d$RNA_delta[i];y<-d$protein_delta[i];if(sd(x)==0||sd(y)==0)NA_real_ else cor(rank(x,ties.method='average'),rank(y,ties.method='average'))})
   s07_close(as.numeric(a$boot$t),v,'bridge complete quadruplet bootstrap ranks');s01_require(a$valid_bootstrap_replicates==sum(is.finite(v)),'bridge bootstrap count')
   if(a$estimable)s07_close(c(a$ci_lower,a$ci_upper),quantile(v[is.finite(v)],c(.025,.975),type=7,names=FALSE),'bridge percentile primary')
   if(a$bca_state=='ESTIMABLE') {z<-s07_capture(boot::boot.ci(a$boot,conf=.95,type='bca'));s07_close(c(a$bca_ci_lower,a$bca_ci_upper),as.numeric(z$value$bca[4:5]),'BCa native anchor')}
  }
 }
 for(r in list(result$protein,result$bridge))for(f in unique(r$multiplicity_family)) {
  j<-which(r$multiplicity_family==f);valid<-is.finite(r$p_value[j]);q<-rep(NA_real_,length(j));q[valid]<-p.adjust(r$p_value[j][valid],'BH',n=s07_cfg()$multiplicity$families[[f]]$frozen_size)
  s07_close(r$q_value[j],q,paste(f,'fixed BH'))
 }
 s01_validate_method_output_rows(result$protein, s07_read('config/hypothesis_estimand_crosswalk.tsv'))
 list(status='PASS',UTC=s07_utc(),C3_hypotheses=nrow(result$protein),C4_hypotheses=nrow(result$bridge),nominal_protein_diagnostics=nrow(result$diagnostics),exact_frozen_identities=TRUE,independent_rank_and_z_scores=TRUE,native_Student_t=TRUE,metafor_SMCRP_LS=TRUE,broad_limma_topTable_and_paired_lm=TRUE,complete_patient_bootstrap_reconstruction=TRUE,bridge_native_cor_test_and_boot_ci=TRUE,fixed_four_families=TRUE,scientific_metadata_reconciliation=TRUE)
}

s07_validate_table_keys <- function(tables) {
 keys<-list(protein_estimates='endpoint_id',cptac_bridge_estimates='endpoint_id',all_program_estimates='endpoint_id',all_member_estimates='endpoint_id',nominal_protein_diagnostics='analysis_id',score_method_diagnostics='analysis_id',patient_program_scores=c('analysis_id','layer','patient_id','tissue'),patient_member_values=c('population_reference_id','layer','member','patient_id','tissue'),protein_paired_deltas=c('endpoint_id','patient_id'),diagnostic_paired_deltas=c('analysis_id','patient_id'),bridge_paired_deltas=c('endpoint_id','patient_id'),bootstrap_replicates=c('operation_name','replicate'))
 for(k in names(keys))s01_assert_unique_key(tables[[k]],keys[[k]],paste('S07',k))
 for(k in c('protein_estimates','cptac_bridge_estimates','all_program_estimates','all_member_estimates')) {
  schema<-yaml::read_yaml(s07_path('config/implementation/s01_output_schemas.yml'))$outputs[[k]]
  s01_require(all(schema$required_columns%in%names(tables[[k]])),paste('required schema',k))
 }
 invisible(TRUE)
}
  run <- function() {
    context$guard();context$resource_check()
    input<-s07_save(s07_inputs(),s07_o('validated_inputs.rds'));scope<-s07_scope()
    protein_path<-s07_run_target('protein',input,scope)
    bridge_paths<-vapply(c('G','H'),function(axis)s07_run_target('bridge',input,scope,axis),'',USE.NAMES=FALSE)
    protein<-read_rds(protein_path);bridges<-lapply(bridge_paths,read_rds);inputs<-read_rds(input)
    result<-s07_assemble(protein,bridges);assembled<-s07_save(result,s07_o('assembled.rds'))
    validation<-s07_validate_models(result,protein,bridges,inputs,scope)
    tables<-s07_tables(result,protein,bridges,inputs);s07_validate_table_keys(tables)
    out<-vapply(names(tables),function(n){p<-context$output('results',paste0(n,'.tsv'));x<-tables[[n]];for(k in names(x))if(is.double(x[[k]]))x[[k]]<-ifelse(is.na(x[[k]]),NA_character_,sprintf('%.17g',x[[k]]));context$write_tsv(x,p);s07_compare(tables[[n]],s07_read(p),n);p},'')
    context$write_rds(validation,context$output('qc','numerical_validation.rds'))
    list(scope='original_protein',tables=out,objects=vapply(c(input,protein_path,bridge_paths,assembled),resolve_path,''),validation=validation)
  }
  list(run=run,input=s07_inputs,project=s07_tables,validate=s07_validate_models,
       kernels=list(score=s07_score,estimate=s07_estimate,pair_bootstrap=s07_pair_bootstrap,bridge=s07_bridge_inference),resolve=resolve_path)
}
