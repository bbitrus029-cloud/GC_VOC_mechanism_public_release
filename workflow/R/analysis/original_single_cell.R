# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Public adaptation: explicit inputs and outputs; original scientific rules retained.

original_single_cell_module <- function(context) {
  historical_namespace <- 's06_single_cell_effects'
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
  s06_path <- resolve_path
  s06_m <- function(p) file.path('data/manifests',historical_namespace,p)
  s06_o <- function(p) file.path('data/processed',historical_namespace,p)
  s06_read <- function(path) utils::read.delim(context$input(resolve_path(path)),quote='',comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
  s06_cfg <- function() context$config
  s06_utc <- utc
  s06_guard <- function() context$guard()
  s06_resource <- function() context$resource_check()
  s06_save <- write_rds
s06_native_members <- function(p,compartment,config=s06_cfg()) {
 cfg<-config$single_cell_member_model
 s01_require(as.character(packageVersion('edgeR'))==cfg$edgeR_version,'edgeR version mismatch')
 md<-p$sample_data;md$sample_id<-md$unit_id;md$donor_id<-md$donor
 md$source_tissue<-md$tissue;md$tissue<-ifelse(md$tissue=='Normal','NAT',md$tissue)
 prepared<-s01_prepare_paired_model(p$counts,md,'donor_id','donor',TRUE)
 s01_require(prepared$biological_n>=cfg$minimum_complete_pairs,'insufficient frozen member pairs')
 result<-s01_sc_member_rows(compartment,config)
 result$biological_n<-prepared$biological_n;result$estimability<-'NA_NONESTIMABLE';result$missing_reason<-'NOT_EVALUATED'
 for(k in c('effect','statistic','numerator_df','denominator_df','p_value','q_value','standard_error','ci_lower','ci_upper'))result[[k]]<-NA_real_
 result$standard_error_status<-result$confidence_interval_status<-'NA_UNAVAILABLE'
 result$uncertainty_reason<-cfg$uncertainty_reason;result$effect_direction<-'Tumor_minus_NAT';result$biological_unit<-'donor';result$ql_df_status<-'NA_NONESTIMABLE'
 warnings<-character()
 native<-tryCatch(withCallingHandlers({
  dge<-edgeR::normLibSizes(edgeR::DGEList(counts=p$counts),method=cfg$normalization)
  s01_require(all(is.finite(dge$samples$norm.factors)&dge$samples$norm.factors>0),'invalid TMM factors')
  fit<-edgeR::glmQLFit(dge,prepared$design,legacy=cfg$glmQLFit$legacy,robust=cfg$glmQLFit$robust,dispersion=cfg$glmQLFit$dispersion,top.proportion=cfg$glmQLFit$top.proportion,abundance.trend=cfg$glmQLFit$abundance.trend,covariate.trend=cfg$glmQLFit$covariate.trend,prior.count=cfg$glmQLFit$prior.count)
  test<-edgeR::glmQLFTest(fit,coef=prepared$coefficient)
  list(dge=dge,fit=fit,test=test)
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
 if(length(warnings))stop(paste('UNEXPECTED_EDGER_WARNING:',paste(warnings,collapse=' | ')))
 if(inherits(native,'error'))result$missing_reason<-paste('NATIVE_MODEL_OR_TEST_FAILURE:',conditionMessage(native)) else {
  extracted<-s06_extract_native(native$test,result$member)
  for(i in seq_len(nrow(result))) {
   d<-p$detection[match(result$member[i],p$detection$gene_id),]
   if(nrow(d)!=1L || is.na(d$available) || !d$available) {result$missing_reason[i]<-'MEMBER_DETECTION_FAILED';next}
   v<-unlist(extracted[i,c('effect','statistic','numerator_df','denominator_df','p_value')])
   if(all(is.finite(v))&&v[2]>=0&&v[3]>0&&v[4]>0&&v[5]>=0&&v[5]<=1) {
    result[i,names(v)]<-as.list(v);result$estimability[i]<-'ESTIMABLE';result$missing_reason[i]<-'NA_NOT_APPLICABLE';result$ql_df_status[i]<-'AVAILABLE'
   } else result$missing_reason[i]<-'INVALID_NATIVE_TEST_OUTPUT'
  }
 }
 s01_validate_sc_member_output(result,'unadjusted',compartment,config)
 list(results=result,native=native,design=prepared$design,coefficient=prepared$coefficient,sample_data=md,input_counts=p$counts,background=p$background,modeling_background=p$modeling_background,original_library_sizes=p$original_library_sizes,eligible_units=p$eligible_units,modeling_units=p$modeling_units,library_sizes=p$modeling_library_sizes,detection=p$detection[p$detection$gene_id%in%result$member,],warnings=warnings,scope='AUTHORIZED_S06_PRODUCTION_ADAPTER')
}

s06_inputs <- function() {
 b<-read_rds('data/processed/s03_single_cell_ingress/raw_pseudobulks.rds');p<-read_rds('data/processed/s03_single_cell_ingress/preparation_objects.rds')
 units<-s06_read('data/manifests/s03_single_cell_ingress/pseudobulk_units.tsv');mapping<-s06_read('data/manifests/s03_single_cell_ingress/gene_mapping.tsv')
 s01_assert_exact_order(units$unit_id,colnames(b$counts),'S06 pseudobulk sample')
 s01_assert_exact_order(mapping$feature_id,rownames(b$counts),'S06 original gene')
 s01_require(all(mapping$unique)&identical(mapping$feature_id,mapping$gene_id),'source gene mapping changed')
 s01_require(identical(dim(b$counts),c(25287L,336L))&&sum(units$cell_count)==92406,'fixed universe mismatch')
 s01_require(all(colSums(b$counts)==units$total_umi),'original UMI mismatch')
 s06_compare(b$units,units,'persisted unit metadata')
 ids<-s06_read('data/manifests/s03_single_cell_ingress/sample_identity.tsv')
 s01_assert_unique_key(ids,c('donor','tissue'),'sample donor tissue')
 i<-match(units$sample_id,ids$sample);s01_require(!anyNA(i)&&all(units$donor==ids$donor[i])&&all(units$tissue==ids$tissue[i])&&all(units$chemistry==ids$chemistry[i]),'source sample suffix bridge mismatch')
 s01_require(all(vapply(split(ids$chemistry,ids$donor),function(x)length(unique(x))==1L,TRUE)),'donor chemistry mismatch')
 backgrounds<-s06_read('data/manifests/s03_single_cell_ingress/ordered_state_member_backgrounds.tsv')
 detection<-s06_read('data/manifests/s03_single_cell_ingress/all_gene_detection.tsv');pairs<-s06_read('data/manifests/s03_single_cell_ingress/complete_pair_identities.tsv')
 for(comp in names(p)) {
  j<-which(units$compartment==comp);q<-s03_prepare_compartment(b$counts[,j,drop=FALSE],b$units[j,],s06_cfg()$single_cell_member_model)
  s01_require(identical(q,p[[comp]]),paste('S03 preparation mismatch',comp))
  z<-backgrounds[backgrounds$compartment==comp,];s01_assert_exact_order(z$gene_id,q$background,'state background')
  s01_assert_exact_order(z$gene_id[z$member_model_retained],q$modeling_background,'model background')
  s01_assert_exact_order(pairs$donor_id[pairs$compartment==comp],q$complete_pair_donors,'complete donors')
  d<-detection[detection$compartment==comp,setdiff(names(detection),'compartment')];rownames(d)<-NULL
  s06_compare(d,q$detection,'all-eligible-unit detection')
 }
 states<-s06_read('data/manifests/s03_single_cell_ingress/state_role_freeze.tsv')
 for(i in which(states$endpoint_eligibility=='ELIGIBLE')) {
  z<-states[i,];q<-p[[z$compartment]];members<-strsplit(z$available_members,',',fixed=TRUE)[[1]]
  expected<-s06_cfg()$programs[[paste0(z$program,'_core')]]$members;d<-q$detection
  s01_assert_exact_order(members,expected[expected%in%d$gene_id[d$available]],'frozen available member')
  s01_require(length(q$complete_pair_donors)==z$complete_paired_donors,'frozen state n')
 }
 markers<-s06_read('reports/qc/s03_single_cell_ingress/marker_decisions.tsv');s01_require(nrow(markers)==7L&&all(markers$state=='PASS'),'marker gate mismatch')
 list(raw=b,prepared=p,states=states,units=units,mapping=mapping)
}

s06_scope <- function() {
 e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv');states<-s06_read('data/manifests/s03_single_cell_ingress/state_role_freeze.tsv')
 out<-data.frame(analysis_id=e$endpoint_id,endpoint_id=e$endpoint_id,category='registered_C2',members='FROZEN_S03',method=e$primary_method)
 for(i in which(states$endpoint_eligibility=='ELIGIBLE')) {
  z<-states[i,];mem<-strsplit(z$available_members,',',fixed=TRUE)[[1]]
  cats<-c('standardized_member_mean','paired_Wilcoxon','paired_percentile_bootstrap',paste0('leave_one_gene_out__',mem),'member_direction_concordance')
  for(cat in cats)out<-rbind(out,data.frame(analysis_id=paste(z$endpoint_id,cat,sep='__'),endpoint_id=z$endpoint_id,category=cat,members=z$available_members,method='prespecified_nominal_diagnostic_NO_FORMAL_FAMILY'))
 }
 # Complete fixed context programs and G submodules remain nominal diagnostics.
 # Use only the four admitted compartments and their frozen available genes/pairs.
 for(comp in unique(states$compartment[states$endpoint_eligibility=='ELIGIBLE'])) {
  sets<-s06_cfg()$coverage_contracts$context_programs
  if(comp=='myeloid')sets<-c(sets,s06_cfg()$coverage_contracts$G_submodules[c('substrate_assembly','oxidation_support')])
  for(n in names(sets))out<-rbind(out,data.frame(analysis_id=paste('GSE206785',comp,n,sep='__'),endpoint_id=paste0('GSE206785_',comp,'_H'),category=paste0('fixed_diagnostic__',n),members=paste(sets[[n]],collapse=','),method='complete_fixed_members_paired_t_nominal_NO_FORMAL_FAMILY'))
 }
 s01_require(!anyDuplicated(out$analysis_id),'duplicate locked scope');out
}

s06_empty <- function(e) {
 # Only current result fields are projected; historical role freeze remains unchanged.
 e$estimates_status<-e$owner_acceptance<-NULL
 e$biological_n<-as.integer(e$structural_n);e$donor_count<-e$biological_n
 for(k in c('effect','ci_lower','ci_upper','p_value','q_value','g_av','g_av_ci_lower','g_av_ci_upper','standard_error','statistic','numerator_df','denominator_df','raw_paired_proportion_difference'))e[[k]]<-NA_real_
 e$estimability<-'NA_NONESTIMABLE';e$missing_reason<-e$reason
 e$multiplicity_status<-if(e$frozen_family_size>0)'FROZEN_FAMILY_BH'else 'NO_FORMAL_FAMILY'
 e$timing<-'PRESPECIFIED_BEFORE_DISEASE_EFFECT_ACCESS';e$dependence<-'donor_aware_cellular_context';e$provenance_id<-s06_m('ENTRY_AUTHORIZATION.json')
 e$member<-e$broad_label<-e$compartment<-e$composition_method<-'NA_NOT_APPLICABLE';e$marker_audit_state<-e$marker_state
 e$standard_error_status<-e$confidence_interval_status<-'NA_NONESTIMABLE';e$uncertainty_reason<-e$reason;e$ql_df_status<-'NA_NOT_APPLICABLE'
 e$native_direction<-'NA_NONESTIMABLE';e$g_av_state<-'NA_NOT_APPLICABLE';e$g_av_reason<-'NA_NOT_APPLICABLE';e
}

s06_apply <- function(e,a) {
 for(k in intersect(names(e),names(a)))if(length(a[[k]])==1L)e[[k]]<-a[[k]]
 e$estimability<-a$state;e$missing_reason<-if(isTRUE(a$estimable))'NA_NOT_APPLICABLE'else if(!is.null(a$reason))a$reason else 'REGISTERED_NONESTIMABILITY'
 e$native_direction<-s01_direction_label(e$effect)
 e$confidence_interval_status<-if(isTRUE(a$estimable))'AVAILABLE'else 'NA_NONESTIMABLE';e$uncertainty_reason<-e$missing_reason;e$standard_error_status<-if(is.finite(e$standard_error))"AVAILABLE"else "NA_NOT_APPLICABLE";e
}

s06_score_estimate <- function(score,p) {
 donors<-p$complete_pair_donors;md<-p$sample_data
 if(!isTRUE(score$estimable))return(list(estimable=FALSE,state='NA_NONESTIMABLE',reason=score$reason,biological_n=length(donors),effect=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,p_value=NA_real_,g_av=NA_real_,g_av_ci_lower=NA_real_,g_av_ci_upper=NA_real_))
 s01_assert_exact_order(score$sample_id,p$modeling_units,'score model units')
 ni<-match(paste(donors,'Normal'),paste(md$donor,md$tissue));ti<-match(paste(donors,'Tumor'),paste(md$donor,md$tissue));s01_require(!anyNA(c(ni,ti)),'program pairing mismatch')
 n<-score$score[ni];t<-score$score[ti];a<-s01_paired_delta_t(t,n,donors);g<-s01_paired_g_av(t,n,donors)
 a$g_av<-g$g_av;a$g_av_ci_lower<-g$ci_lower;a$g_av_ci_upper<-g$ci_upper;a$g_av_state<-g$state;a$g_av_reason<-if(is.null(g$reason))'NA_NOT_APPLICABLE'else g$reason
 a$standard_error<-if(isTRUE(a$estimable))sd(t-n)/sqrt(length(donors))else NA_real_
 a$source_data<-data.frame(donor=donors,NAT_unit=p$modeling_units[ni],Tumor_unit=p$modeling_units[ti],NAT_score=n,Tumor_score=t,delta=t-n);a
}

s06_state_target <- function(inputs,guard,scope) {
 s06_guard();s06_resource();x<-read_rds(inputs);e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv');e<-e[e$endpoint_id%in%x$states$endpoint_id,]
 estimates<-list();scores<-list();rows<-list();diagnostics<-list()
 for(i in seq_len(nrow(e))) {
  z<-x$states[match(e$endpoint_id[i],x$states$endpoint_id),];r<-s06_empty(e[i,]);r$broad_label<-r$compartment<-z$compartment
  if(z$endpoint_eligibility=='ELIGIBLE') {
   p<-x$prepared[[z$compartment]];mat<-x$raw$counts[p$background,p$modeling_units,drop=FALSE];mem<-strsplit(z$available_members,',',fixed=TRUE)[[1]]
   score<-s01_rank_program_score(mat,mem,p$background);a<-s06_score_estimate(score,p);scores[[z$endpoint_id]]<-score;estimates[[z$endpoint_id]]<-a;r<-s06_apply(r,a)
   sr<-scope[scope$endpoint_id==z$endpoint_id&scope$category!='registered_C2',]
   for(j in seq_len(nrow(sr))) {
    d<-sr[j,];cat<-d$category;members<-strsplit(d$members,',',fixed=TRUE)[[1]];alt<-NULL
    if(cat=='standardized_member_mean')alt<-s06_score_estimate(s01_standardized_member_mean(mat,mem),p)
    if(startsWith(cat,'leave_one_gene_out__'))alt<-s06_score_estimate(s01_rank_program_score(mat,setdiff(mem,sub('leave_one_gene_out__','',cat,fixed=TRUE)),p$background),p)
    if(startsWith(cat,'fixed_diagnostic__')) {
     available<-p$detection$gene_id[p$detection$available]
     sc<-if(all(members%in%available))s01_rank_program_score(mat,members,p$background)else s01_nonestimable(list(reason='COMPLETE_FIXED_DIAGNOSTIC_MEMBERS_NOT_AVAILABLE'))
     alt<-s06_score_estimate(sc,p)
    }
    if(cat=='paired_percentile_bootstrap') {
     seeds<-s06_read('config/implementation/s01_target_seed_map.tsv');seed<-seeds$seed[match(d$analysis_id,seeds$operation_name)]
     alt<-s01_pair_bootstrap(a$source_data$Tumor_score,a$source_data$NAT_score,a$source_data$donor,seed);alt$effect<-a$effect;alt$seed<-seed;alt$biological_n<-a$biological_n
    }
    if(cat=='paired_Wilcoxon') {
     # No blanket warning suppression in the S06 adapter.
     w<-stats::wilcox.test(a$source_data$Tumor_score,a$source_data$NAT_score,paired=TRUE,exact=FALSE,correct=FALSE,alternative='two.sided')
     alt<-list(estimable=is.finite(w$p.value),state=if(is.finite(w$p.value))'ESTIMABLE'else 'NA_NONESTIMABLE',effect=a$effect,p_value=w$p.value,statistic=unname(w$statistic),biological_n=a$biological_n)
    }
    if(cat=='member_direction_concordance') {
     full<-s06_cfg()$programs[[paste0(z$program,'_core')]]$members
     # The registered 7/4 denominator is not silently reduced for partial SC scores.
     if(!all(full%in%mem))alt<-list(estimable=FALSE,state='NA_NONESTIMABLE',reason='REGISTERED_FULL_7_OR_4_MEMBER_DIAGNOSTIC_UNAVAILABLE_IN_REDUCED_MEMBER_STATE')else {
      md<-p$sample_data;delta<-vapply(full,function(g)mean(mat[g,md$tissue=='Tumor']-mat[g,md$tissue=='Normal']),0)
      c<-s01_member_direction_concordance(a$effect,delta,length(full));alt<-c;alt$effect<-a$effect;alt$biological_n<-a$biological_n
     }
    }
    s01_require(!is.null(alt),paste('unhandled diagnostic',cat));diagnostics[[d$analysis_id]]<-list(registry=d,estimate=alt)
   }
  }
  rows[[r$endpoint_id]]<-r
 }
 s06_save(list(rows=do.call(rbind,rows),scores=scores,estimates=estimates,diagnostics=diagnostics),s06_o('states.rds'))
}

s06_member_target <- function(comp,inputs,guard) {
 s06_guard();s06_resource();x<-read_rds(inputs);a<-s06_native_members(x$prepared[[comp]],comp)
 s06_save(a,s06_o(paste0('members_',comp,'.rds')))
}

s06_composition_data <- function() {
 c<-s06_read('data/manifests/s03_single_cell_ingress/composition_counts.tsv');n<-s06_read('data/manifests/s03_single_cell_ingress/composition_denominators.tsv');sets<-s06_read('data/manifests/s03_single_cell_ingress/composition_donor_sets.tsv')
 s01_assert_unique_key(c,c('donor','tissue','compartment'),'composition counts');s01_assert_unique_key(n,c('donor','tissue'),'composition total')
 i<-match(paste(c$donor,c$tissue),paste(n$donor,n$tissue));s01_require(!anyNA(i),'composition denominator join');c$total_count<-n$retained_cells[i];c$compartment_count<-c$cell_count
 totals<-aggregate(cell_count~donor+tissue,c,sum);j<-match(paste(n$donor,n$tissue),paste(totals$donor,totals$tissue));s01_require(all(n$retained_cells==totals$cell_count[j]),'composition denominator sum mismatch')
 c$source_tissue<-c$tissue;c$tissue<-ifelse(c$tissue=='Normal','NAT','Tumor');c$proportion<-c$cell_count/c$total_count
 c$empirical_logit<-s01_empirical_logit(c$cell_count,c$total_count);c$CLR<-NA_real_
 for(key in unique(paste(c$donor,c$tissue))) {j<-which(paste(c$donor,c$tissue)==key);s01_require(setequal(c$compartment[j],s06_cfg()$composition_inference$compartments),'seven parts required');c$CLR[j]<-s01_clr(c$cell_count[j])}
 for(k in c('primary_100','strict_200','count_aware'))c[[k]]<-sets[[k]][match(c$donor,sets$donor)]
 mincount<-tapply(n$retained_cells,n$donor,min);s01_require(all(sets$primary_100==(mincount[sets$donor]>=100))&&all(sets$strict_200==(mincount[sets$donor]>=200))&&all(sets$count_aware),'composition donor sets mismatch')
 c[order(c$donor,match(c$tissue,c('NAT','Tumor')),match(c$compartment,s06_cfg()$composition_inference$compartments)),]
}

s06_composition_target <- function(inputs,guard) {
 s06_guard();s06_resource();source<-s06_composition_data();e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv');e<-e[grepl('^SC_COMPOSITION',e$multiplicity_family),]
 rows<-estimates<-deltas<-list()
 for(i in seq_len(nrow(e))) {
  r<-s06_empty(e[i,]);r$compartment<-r$family_hypothesis
  family<-r$multiplicity_family;flag<-switch(family,SC_COMPOSITION_PRIMARY='primary_100',SC_COMPOSITION_STRICT200='strict_200',SC_COMPOSITION_CLR='primary_100',SC_COMPOSITION_QUASIBINOMIAL='count_aware')
  d<-source[source$compartment==r$compartment & source[[flag]],];n<-d[d$tissue=='NAT',];t<-d[d$tissue=='Tumor',];s01_assert_exact_order(n$donor,t$donor,'composition pair')
  field<-if(family=='SC_COMPOSITION_CLR')'CLR'else 'empirical_logit';r$composition_method<-if(family=='SC_COMPOSITION_QUASIBINOMIAL')'quasibinomial_logit_tissue_plus_donor'else field
  delta<-data.frame(endpoint_id=r$endpoint_id,donor=n$donor,NAT_count=n$cell_count,Tumor_count=t$cell_count,NAT_total=n$total_count,Tumor_total=t$total_count,NAT_proportion=n$proportion,Tumor_proportion=t$proportion,proportion_delta=t$proportion-n$proportion,NAT_transformed=n[[field]],Tumor_transformed=t[[field]],delta=t[[field]]-n[[field]],transform=field)
  if(family=='SC_COMPOSITION_QUASIBINOMIAL')a<-s01_fit_quasibinomial(d)else a<-s01_paired_delta_t(t[[field]],n[[field]],n$donor)
  if(!isTRUE(a$estimable)&&is.null(a$reason))a$reason<-paste(names(a$flags)[a$flags],collapse=';')
  a$biological_n<-nrow(n);if(is.null(a$audit_trail))a$standard_error<-if(isTRUE(a$estimable))sd(delta$delta)/sqrt(nrow(n))else NA_real_;r<-s06_apply(r,a);r$donor_count<-nrow(n);r$raw_paired_proportion_difference<-mean(delta$proportion_delta)
  rows[[r$endpoint_id]]<-r;estimates[[r$endpoint_id]]<-a;deltas[[r$endpoint_id]]<-delta
 }
 s06_save(list(rows=do.call(rbind,rows),estimates=estimates,source=source,deltas=do.call(rbind,deltas)),s06_o('composition.rds'))
}

s06_compare <- function(x,y,label) {
 s01_require(is.data.frame(x)&&is.data.frame(y)&&identical(names(x),names(y))&&nrow(x)==nrow(y),paste(label,'schema mismatch'))
 for(k in names(x)) {
  if(is.numeric(x[[k]])) {
   u<-as.numeric(x[[k]]);v<-as.numeric(y[[k]])
   s01_require(identical(is.na(u),is.na(v))&&all(abs(u[!is.na(u)]-v[!is.na(v)])<=1e-13*pmax(1,abs(u[!is.na(u)]))),paste(label,k,'value mismatch'))
  }else s01_require(identical(as.character(x[[k]]),as.character(y[[k]])),paste(label,k,'identity mismatch'))
 };invisible(TRUE)
}

s06_adjust <- function(rows,config=s06_cfg()) {
 for(f in unique(rows$multiplicity_family[rows$frozen_family_size>0])) {
  spec<-config$multiplicity$families[[f]];j<-which(rows$multiplicity_family==f)
  s01_require(length(j)==spec$frozen_size&&all(rows$frozen_family_size[j]==spec$frozen_size),'frozen family size shrinkage')
  s01_require(setequal(rows$family_hypothesis[j],spec$hypotheses),'family identity mismatch')
  k<-j[match(spec$hypotheses,rows$family_hypothesis[j])]
  rows$q_value[k]<-unname(s01_frozen_bh_values(setNames(rows$p_value[k],spec$hypotheses),as.integer(spec$frozen_size),spec$hypotheses))
 };rows
}

s06_member_projection <- function(models) {
 e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv');e<-e[e$multiplicity_family%in%c('SC_G_MEMBERS','SC_H_MEMBERS'),]
 raw<-do.call(rbind,lapply(models,function(m)m$results));raw<-s01_adjust_sc_member_families(raw,s06_cfg())
 rows<-lapply(seq_len(nrow(e)),function(i) {
  r<-s06_empty(e[i,]);a<-raw[match(r$family_hypothesis,raw$family_hypothesis),]
  for(k in intersect(names(a),names(r)))r[[k]]<-a[[k]]
  r$broad_label<-a$compartment;r$native_direction<-s01_direction_label(r$effect);r
 });do.call(rbind,rows)
}

s06_diagnostic_table <- function(states) {
 out<-do.call(rbind,lapply(states$diagnostics,function(d) {
  a<-d$estimate;r<-d$registry
  for(k in c('effect','ci_lower','ci_upper','p_value','statistic','biological_n','g_av','g_av_ci_lower','g_av_ci_upper','valid_replicates','seed','concordant_members','denominator','proportion'))r[[k]]<-if(is.null(a[[k]]))NA_real_ else a[[k]]
  r$diagnostic_state<-a$state;r$missing_reason<-if(is.null(a$reason))'NA_NOT_APPLICABLE'else a$reason;r$q_value<-NA_real_;r$multiplicity_status<-r$multiplicity_family<-'NO_FORMAL_FAMILY'
  r$timing<-'PRESPECIFIED_BEFORE_DISEASE_EFFECT_ACCESS';r$dependence<-'same_data_robustness';r$loo_direction_stable<-'NA_NOT_APPLICABLE'
  primary<-states$estimates[[r$endpoint_id]]
  if(startsWith(r$category,'leave_one_gene_out__'))r$loo_direction_stable<-s01_loo_direction_status(primary$effect,a$effect,isTRUE(a$estimable))
  r
 }))
 s06_correct_metadata(out)
}

s06_assemble <- function(state_path,member_paths,comp_path,scope) {
 s<-read_rds(state_path);models<-lapply(member_paths,read_rds);c<-read_rds(comp_path)
 e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv')
 rows<-rbind(s$rows,s06_member_projection(models),c$rows);rows<-rows[match(e$endpoint_id,rows$endpoint_id),];rows<-s06_adjust(rows)
 s01_require(nrow(rows)==56L&&!anyDuplicated(rows$endpoint_id),'C2 population mismatch')
 diag<-s06_diagnostic_table(s)
 s01_require(setequal(diag$analysis_id,scope$analysis_id[scope$category!='registered_C2']),'locked diagnostic scope incomplete')
 source<-do.call(rbind,lapply(names(s$estimates),function(id)cbind(endpoint_id=id,s$estimates[[id]]$source_data)))
 alt_source<-do.call(rbind,lapply(names(s$diagnostics),function(id){a<-s$diagnostics[[id]]$estimate;if(is.null(a$source_data))return(NULL);cbind(analysis_id=id,a$source_data)}))
 methods<-do.call(rbind,lapply(names(s$estimates),function(id) {
  a<-s$estimates[[id]];z<-s$diagnostics[[paste0(id,'__standardized_member_mean')]]$estimate
  cbind(endpoint_id=id,as.data.frame(s01_score_method_sensitivity(a$effect,z$effect,isTRUE(z$estimable)),stringsAsFactors=FALSE))
 }))
 audit<-do.call(rbind,lapply(names(c$estimates),function(id){a<-c$estimates[[id]];if(is.null(a$audit_trail))return(NULL);z<-a$audit_trail;z$warning_messages<-jsonlite::toJSON(z$warning_messages,auto_unbox=FALSE);r<-as.data.frame(z,stringsAsFactors=FALSE);for(k in names(a$flags))r[[paste0('flag_',k)]]<-a$flags[[k]];cbind(endpoint_id=id,r)}))
 result<-list(rows=rows,diagnostics=diag,source=source,diagnostic_source=alt_source,score_methods=methods,quasibinomial_audit=audit,state_path=state_path,member_paths=member_paths,composition_path=comp_path)
 s06_save(result,s06_o('assembled.rds'))
}

s06_tables <- function(result) {
 f<-result$rows;state<-!grepl('^SC_COMPOSITION',f$multiplicity_family);member<-f$multiplicity_family%in%c('SC_G_MEMBERS','SC_H_MEMBERS');c<-read_rds(result$composition_path)
 models<-lapply(result$member_paths,read_rds)
 inputs<-do.call(rbind,lapply(models,function(m){d<-m$sample_data;d$compartment<-m$results$compartment[1];d$model_library_size<-as.numeric(m$library_sizes);d$TMM_factor<-m$native$dge$samples$norm.factors;d}))
 tables<-list(single_cell_state_estimates=f[state,],all_program_estimates=f[state&!member,],all_member_estimates=f[member,],captured_composition_estimates=f[!state,],donor_program_deltas=result$source,nominal_program_diagnostics=result$diagnostics,diagnostic_donor_deltas=result$diagnostic_source,score_method_diagnostics=result$score_methods,quasibinomial_audit=result$quasibinomial_audit,composition_source_data=c$source,composition_paired_deltas=c$deltas,member_model_input_provenance=inputs,hypothesis_coverage=f[,c('endpoint_id','biological_n','estimability','missing_reason','multiplicity_family','family_hypothesis','frozen_family_size')])
 schemas<-read_yaml('config/implementation/s01_output_schemas.yml')$outputs
 for(n in names(tables)) {
  x<-tables[[n]];if(n%in%names(schemas))s01_require(all(schemas[[n]]$required_columns%in%names(x)),paste('schema missing fields',n))
  keys<-intersect(c('endpoint_id','analysis_id','compartment','donor','tissue'),names(x))
  if(length(keys))x<-x[do.call(order,c(x[keys],list(method='radix'))),,drop=FALSE]
  rownames(x)<-NULL;tables[[n]]<-x
 };tables
}

s06_anchor_pair <- function(a,d,label,g=TRUE) {
 s01_require(identical(d$delta,d$Tumor_score-d$NAT_score),paste(label,'paired source delta'))
 delta<-d$delta;n<-nrow(d);s01_require(a$biological_n==n,paste(label,'biological n'))
 if(isTRUE(a$estimable)) {
  tt<-stats::t.test(delta,mu=0);u<-c(a$effect,a$statistic,a$p_value,a$ci_lower,a$ci_upper);v<-c(mean(delta),unname(tt$statistic),tt$p.value,unname(tt$conf.int))
  s01_require(isTRUE(all.equal(u,v,tolerance=1e-12)),paste(label,'native t anchor'))
 }
 if(g&&identical(a$g_av_state,'ESTIMABLE')) {
  z<-metafor::escalc(measure='SMCRP',m1i=mean(d$Tumor_score),m2i=mean(d$NAT_score),sd1i=sd(d$Tumor_score),sd2i=sd(d$NAT_score),ni=n,ri=cor(d$Tumor_score,d$NAT_score),vtype='LS',correct=TRUE)
  v<-c(as.numeric(z$yi),as.numeric(z$yi)-qnorm(.975)*sqrt(z$vi),as.numeric(z$yi)+qnorm(.975)*sqrt(z$vi))
  s01_require(max(abs(c(a$g_av,a$g_av_ci_lower,a$g_av_ci_upper)-v))<=1e-12,paste(label,'metafor anchor'))
 };invisible(TRUE)
}

s06_validate <- function(result) {
 s06_guard();x<-s06_inputs();s<-read_rds(result$state_path);c<-read_rds(result$composition_path);models<-lapply(result$member_paths,read_rds)
 recon<-rbind(s$rows,s06_member_projection(models),c$rows);e<-s06_read('data/manifests/s03_single_cell_ingress/endpoint_role_freeze.tsv');recon<-s06_adjust(recon[match(e$endpoint_id,recon$endpoint_id),])
 s06_compare(result$rows,recon,'model-to-assembled reconciliation')
 crosswalk<-s06_read('config/hypothesis_estimand_crosswalk.tsv');crosswalk<-crosswalk[crosswalk$contract=='C2',]
 s06_compare(crosswalk,result$rows[,names(crosswalk)],'crosswalk identity')
 for(id in names(s$estimates)) {
  z<-x$states[match(id,x$states$endpoint_id),];p<-x$prepared[[z$compartment]];mat<-x$raw$counts[p$background,p$modeling_units,drop=FALSE];members<-strsplit(z$available_members,',',fixed=TRUE)[[1]]
  # Direct rank formula, independently indexed from the saved scoring helper output.
  scores<-vapply(seq_len(ncol(mat)),function(j)mean((rank(mat[,j],ties.method='average')[match(members,rownames(mat))]-1)/(nrow(mat)-1)),0)
  s01_require(identical(unname(scores),s$scores[[id]]$score),paste(id,'direct rank reconstruction'))
  a<-s$estimates[[id]];s06_anchor_pair(a,a$source_data,id)
  d<-result$source[result$source$endpoint_id==id,setdiff(names(result$source),'endpoint_id')];s06_compare(a$source_data,d,paste(id,'donor projection'))
 }
 for(d in s$diagnostics) {
  a<-d$estimate;cat<-d$registry$category;id<-d$registry$endpoint_id;p<-x$prepared[[x$states$compartment[match(id,x$states$endpoint_id)]]]
  if(!is.null(a$source_data))s06_anchor_pair(a,a$source_data,d$registry$analysis_id)
  if(cat=='paired_percentile_bootstrap') {
   dd<-s$estimates[[id]]$source_data;cfg<-s06_cfg()$diagnostic_bootstrap$paired_program
   draws<-s01_with_production_rng(a$seed,vapply(seq_len(cfg$replicates),function(i)mean(dd$delta[sample.int(nrow(dd),nrow(dd),replace=TRUE)]),0))
   ci<-quantile(draws[is.finite(draws)],c(.025,.975),type=cfg$quantile_type,names=FALSE)
   s01_require(a$valid_replicates==sum(is.finite(draws))&&identical(c(a$ci_lower,a$ci_upper),ci),'bootstrap seed/percentile reconstruction')
  }
 }
 for(m in models) {
  p<-x$prepared[[m$results$compartment[1]]];s01_require(identical(m$input_counts,p$counts)&&identical(m$library_sizes,colSums(p$counts)),'member input matrix/library mismatch')
  s01_require(identical(m$eligible_units,p$eligible_units)&&identical(m$background,p$background),'narrowed member detection denominator')
  if(!inherits(m$native,'error')) {
   expected<-s06_extract_native(m$native$test,m$results$member)
   # Separate literal mapping ensures extraction cannot swap native F, df.test or df.total.
   native<-m$native$test;j<-match(m$results$member,rownames(native$table))
   s01_require(identical(expected$statistic,native$table$F[j])&&identical(expected$numerator_df,rep_len(native$df.test,nrow(native$table))[j])&&identical(expected$denominator_df,rep_len(native$df.total,nrow(native$table))[j]),'native QL field mapping')
   s06_compare(expected,m$results[,names(expected)],'native member output extraction')
   s01_require(identical(as.numeric(m$native$dge$samples$lib.size),as.numeric(colSums(p$counts))),'edgeR filtered library sizes')
  }
 }
 s06_compare(s06_composition_data(),c$source,'composition source reconstruction')
 for(id in names(c$estimates)) {
  a<-c$estimates[[id]];d<-c$deltas[c$deltas$endpoint_id==id,];r<-c$rows[c$rows$endpoint_id==id,]
  source<-c$source[c$source$donor%in%d$donor & c$source$compartment==r$compartment,];n<-source[source$tissue=='NAT',];t<-source[source$tissue=='Tumor',]
  s01_assert_exact_order(d$donor,n$donor,'composition saved donor order')
  s01_require(identical(d$NAT_count,n$cell_count)&&identical(d$Tumor_count,t$cell_count)&&identical(d$NAT_total,n$total_count)&&identical(d$Tumor_total,t$total_count),'composition source-to-pair counts')
  s01_require(identical(d$proportion_delta,t$proportion-n$proportion)&&r$raw_paired_proportion_difference==mean(d$proportion_delta),'descriptive raw proportion delta')
  s01_require(identical(d$NAT_transformed,n[[d$transform[1]]])&&identical(d$Tumor_transformed,t[[d$transform[1]]]),'composition transformed source-to-pair')
  if(is.null(a$audit_trail)) {
   dd<-data.frame(donor=d$donor,NAT_score=d$NAT_transformed,Tumor_score=d$Tumor_transformed,delta=d$delta);s06_anchor_pair(a,dd,id,FALSE)
  }else {
   # The saved audit is checked against direct native GLM facts and t calibration.
   fit<-withCallingHandlers(glm(cbind(compartment_count,total_count-compartment_count)~tissue+donor,data=transform(source,tissue=factor(tissue,c('NAT','Tumor')),donor=factor(donor)),family=quasibinomial('logit'),control=glm.control(epsilon=1e-8,maxit=25,trace=FALSE)),warning=function(w){s01_require(conditionMessage(w)%in%a$warnings,'new GLM warning during anchor');invokeRestart('muffleWarning')})
   sm<-summary(fit);beta<-sm$coefficients['tissueTumor','Estimate'];se<-sm$coefficients['tissueTumor','Std. Error']
   facts<-list(converged=fit$converged,boundary=fit$boundary,rank=fit$rank,columns=ncol(model.matrix(fit)),fitted=fitted(fit),coefficient=beta,se=se,dispersion=sm$dispersion,df_residual=fit$df.residual)
   adj<-s01_adjudicate_quasibinomial(facts,a$warnings);s01_require(identical(adj$flags,a$flags)&&identical(adj$audit_trail,a$audit_trail),'GLM audit reconstruction')
   if(a$estimable)s01_require(isTRUE(all.equal(c(a$effect,a$standard_error,a$p_value,a$ci_lower,a$ci_upper),c(beta,se,2*pt(-abs(beta/se),fit$df.residual),beta+qt(.025,fit$df.residual)*se,beta+qt(.975,fit$df.residual)*se),tolerance=1e-12)),'GLM t calibration')
  }
 }
 s06_reconcile_sources(result,x,s,c)
 # Reconcile every auxiliary projection, including diagnostics and audit facts.
 s06_compare(result$diagnostics,s06_diagnostic_table(s),'nominal diagnostic projection')
 s01_require(sum(result$rows$estimability=='NA_NONESTIMABLE')>=9L,'frozen non-estimability rows omitted')
 f<-result$rows;structural<-f$endpoint_eligibility=='NA_NONESTIMABLE';s01_require(sum(structural)==9L&&all(is.na(f$effect[structural]))&&all(is.na(f$p_value[structural])),'effects generated for frozen unavailable states')
 for(family in unique(f$multiplicity_family[f$frozen_family_size>0])) {
  z<-f[f$multiplicity_family==family,];ok<-is.finite(z$p_value);q<-rep(NA_real_,nrow(z));q[ok]<-p.adjust(z$p_value[ok],'BH',n=s06_cfg()$multiplicity$families[[family]]$frozen_size)
  s01_require(identical(q,z$q_value),paste('full-family BH anchor',family))
 }
 list(status='PASS',verified_utc=s06_utc(),C2_rows=56L,state_rows=14L,member_rows=14L,composition_rows=28L,nominal_rows=nrow(result$diagnostics),donor_delta_rows=nrow(result$source),direct_rank_and_t_reconstruction=TRUE,metafor_g_av_anchor=TRUE,bootstrap_reconstruction=TRUE,native_QL_extraction=TRUE,composition_source_and_GLM_anchors=TRUE,full_family_BH=TRUE)
}

s06_assert_estimate_row <- function(row,a,label) {
 fields<-intersect(c('effect','ci_lower','ci_upper','p_value','statistic','biological_n','standard_error','g_av','g_av_ci_lower','g_av_ci_upper'),names(a))
 expected<-as.data.frame(a[fields],stringsAsFactors=FALSE)
 s06_compare(expected,row[,fields,drop=FALSE],label)
 s01_require(identical(row$estimability,a$state),paste(label,'estimability mismatch'))
 if(identical(a$state,'ESTIMABLE'))s01_require(row$missing_reason=='NA_NOT_APPLICABLE',paste(label,'spurious failure reason'))
 invisible(TRUE)
}

s06_reconcile_sources <- function(result,x,s,c) {
 for(id in names(s$estimates)) {
  p<-x$prepared[[x$states$compartment[match(id,x$states$endpoint_id)]]]
  a<-s$estimates[[id]];direct<-s06_score_estimate(s$scores[[id]],p)
  s06_compare(a$source_data,direct$source_data,paste(id,'scores-to-paired-source'))
  s06_assert_estimate_row(s$rows[s$rows$endpoint_id==id,],a,paste(id,'state estimate-to-row'))
 }
 for(id in names(c$estimates))s06_assert_estimate_row(c$rows[c$rows$endpoint_id==id,],c$estimates[[id]],paste(id,'composition estimate-to-row'))
 for(d in s$diagnostics) {
  a<-d$estimate;id<-d$registry$endpoint_id;cat<-d$registry$category;z<-x$states[match(id,x$states$endpoint_id),];p<-x$prepared[[z$compartment]]
  mat<-x$raw$counts[p$background,p$modeling_units,drop=FALSE];mem<-strsplit(z$available_members,',',fixed=TRUE)[[1]]
  if(!is.null(a$source_data)) {
   members<-if(startsWith(cat,'leave_one_gene_out__'))setdiff(mem,sub('leave_one_gene_out__','',cat,fixed=TRUE))else if(startsWith(cat,'fixed_diagnostic__'))strsplit(d$registry$members,',',fixed=TRUE)[[1]]else mem
   if(cat=='standardized_member_mean') {
    values<-mat[members,,drop=FALSE];scaled<-sweep(sweep(values,1,rowMeans(values),'-'),1,apply(values,1,sd),'/');score<-colMeans(scaled)
   }else score<-vapply(seq_len(ncol(mat)),function(j)mean((rank(mat[,j],ties.method='average')[match(members,rownames(mat))]-1)/(nrow(mat)-1)),0)
   direct<-s06_score_estimate(list(estimable=TRUE,score=unname(score),sample_id=colnames(mat)),p)
   s06_compare(a$source_data,direct$source_data,paste(d$registry$analysis_id,'diagnostic formula source'))
  }
  if(cat=='paired_Wilcoxon') {
   dd<-s$estimates[[id]]$source_data;w<-wilcox.test(dd$Tumor_score,dd$NAT_score,paired=TRUE,exact=FALSE,correct=FALSE,alternative='two.sided')
   s01_require(identical(a$p_value,w$p.value)&&identical(a$statistic,unname(w$statistic)),'saved nominal Wilcoxon mismatch')
  }
 }
 alt_source<-do.call(rbind,lapply(names(s$diagnostics),function(id){a<-s$diagnostics[[id]]$estimate;if(is.null(a$source_data))return(NULL);cbind(analysis_id=id,a$source_data)}))
 s06_compare(result$diagnostic_source,alt_source,'complete diagnostic donor-source projection')
 methods<-do.call(rbind,lapply(names(s$estimates),function(id){a<-s$estimates[[id]];z<-s$diagnostics[[paste0(id,'__standardized_member_mean')]]$estimate;cbind(endpoint_id=id,as.data.frame(s01_score_method_sensitivity(a$effect,z$effect,isTRUE(z$estimable)),stringsAsFactors=FALSE))}))
 s06_compare(result$score_methods,methods,'score-method flags-to-source')
 audit<-do.call(rbind,lapply(names(c$estimates),function(id){a<-c$estimates[[id]];if(is.null(a$audit_trail))return(NULL);z<-a$audit_trail;z$warning_messages<-jsonlite::toJSON(z$warning_messages,auto_unbox=FALSE);r<-as.data.frame(z,stringsAsFactors=FALSE);for(k in names(a$flags))r[[paste0('flag_',k)]]<-a$flags[[k]];cbind(endpoint_id=id,r)}))
 s06_compare(result$quasibinomial_audit,audit,'quasibinomial audit-to-source')
 cc<-c$source;uu<-x$units;idx<-match(paste(cc$donor,cc$source_tissue,cc$compartment),paste(uu$donor,uu$tissue,uu$compartment))
 s01_require(!anyNA(idx)&&identical(cc$cell_count,uu$cell_count[idx]),'composition counts-to-original retained universe')
 invisible(TRUE)
}

s06_expected_metadata <- function() {
 states<-s06_read('data/manifests/s03_single_cell_ingress/state_role_freeze.tsv')
 cfg<-s06_cfg();out<-list()
 add<-function(id,category,parent,program,construct,population,members,deleted='NA_NOT_APPLICABLE',parent_members='NA_NOT_APPLICABLE') {
  out[[id]]<<-data.frame(analysis_id=id,category=category,endpoint_id=if(parent=='NA_NOT_APPLICABLE')id else parent,
   scientific_program=program,construct_id=construct,parent_state_id=parent,population_reference_state_id=population,
   deleted_member=deleted,members=paste(members,collapse=','),parent_members=parent_members,stringsAsFactors=FALSE)
 }
 for(i in which(states$endpoint_eligibility=='ELIGIBLE')) {
  z<-states[i,];mem<-strsplit(z$available_members,',',fixed=TRUE)[[1]];id<-z$endpoint_id
  for(cat in c('standardized_member_mean','paired_Wilcoxon','paired_percentile_bootstrap',paste0('leave_one_gene_out__',mem),'member_direction_concordance')) {
   deleted<-if(startsWith(cat,'leave_one_gene_out__'))substring(cat,nchar('leave_one_gene_out__')+1L)else 'NA_NOT_APPLICABLE'
   add(paste0(id,'__',cat),cat,id,paste0(z$program,'_core'),paste0(z$program,'_core'),id,mem[mem!=deleted],deleted,z$available_members)
  }
 }
 for(comp in c('epithelial','fibroblast','myeloid','B_plasma')) {
  sets<-cfg$coverage_contracts$context_programs
  if(comp=='myeloid')sets<-c(sets,cfg$coverage_contracts$G_submodules[c('substrate_assembly','oxidation_support')])
  for(n in names(sets)) {
   submodule<-n%in%names(cfg$coverage_contracts$G_submodules)
   add(paste('GSE206785',comp,n,sep='__'),paste0('fixed_diagnostic__',n),
    if(submodule)'GSE206785_myeloid_G'else 'NA_NOT_APPLICABLE',if(submodule)'G_core'else n,n,
    paste0('GSE206785_',comp,'_H'),sets[[n]])
  }
 }
 do.call(rbind,out)
}

s06_correct_metadata <- function(d) {
 e<-s06_expected_metadata();scope<-s06_scope()
 s01_require(nrow(d)==nrow(e)&&!anyDuplicated(d$analysis_id)&&setequal(d$analysis_id,e$analysis_id),'diagnostic population mismatch')
 j<-match(d$analysis_id,e$analysis_id);k<-match(d$analysis_id,scope$analysis_id)
 s01_require(!anyNA(k)&&all(d$category==e$category[j]),'diagnostic category mapping mismatch')
 for(n in setdiff(names(e),'analysis_id'))d[[n]]<-e[[n]][j]
 d$original_scope_row<-k;d$original_scope_analysis_id<-scope$analysis_id[k]
 d$original_scope_endpoint_id<-scope$endpoint_id[k];d$original_scope_members<-scope$members[k]
 d$metadata_provenance<-'RESULT_INFORMED_DELIVERY_CORRECTION_OF_HISTORICAL_SCOPE'
 d
}

s06_validate_metadata <- function(d,states=NULL,inputs=NULL) {
 e<-s06_expected_metadata();scope<-s06_scope()
 s01_require(nrow(d)==56L&&!anyDuplicated(d$analysis_id)&&setequal(d$analysis_id,e$analysis_id),'diagnostic population mismatch')
 j<-match(d$analysis_id,e$analysis_id);s06_compare(e[j,],d[,names(e)],'independent diagnostic semantics')
 k<-match(d$analysis_id,scope$analysis_id)
 s01_require(identical(as.integer(d$original_scope_row),as.integer(k))&&identical(d$original_scope_analysis_id,scope$analysis_id[k])&&identical(d$original_scope_endpoint_id,scope$endpoint_id[k])&&identical(d$original_scope_members,scope$members[k]),'original scope mapping mismatch')
 s01_require(all(is.na(d$q_value))&&all(d$multiplicity_family=='NO_FORMAL_FAMILY'),'diagnostic role mismatch')
 if(!is.null(states)) {
  s01_require(!anyDuplicated(names(states$diagnostics))&&setequal(names(states$diagnostics),e$analysis_id),'source diagnostic population mismatch')
  for(i in seq_len(nrow(d))) {
   r<-d[i,];a<-states$diagnostics[[r$analysis_id]]
   s06_compare(scope[k[i],],a$registry,'unchanged historical source registry')
   s01_require(identical(a$estimate$state,r$diagnostic_state),'source diagnostic state mismatch')
   p<-inputs$prepared[[inputs$states$compartment[match(r$population_reference_state_id,inputs$states$endpoint_id)]]]
   mem<-strsplit(r$members,',',fixed=TRUE)[[1]];mat<-inputs$raw$counts[p$background,p$modeling_units,drop=FALSE]
   if(startsWith(r$category,'fixed_diagnostic__')&&!all(mem%in%p$detection$gene_id[p$detection$available])) {
    s01_require(!isTRUE(a$estimate$estimable)&&r$missing_reason=='COMPLETE_FIXED_DIAGNOSTIC_MEMBERS_NOT_AVAILABLE','fixed diagnostic coverage state mismatch')
   }
   if(!is.null(a$estimate$source_data)) {
    if(r$category=='standardized_member_mean') {
     v<-mat[mem,,drop=FALSE];v<-sweep(sweep(v,1,rowMeans(v),'-'),1,apply(v,1,sd),'/');score<-colMeans(v)
    } else score<-vapply(seq_len(ncol(mat)),function(h)mean((rank(mat[,h],ties.method='average')[match(mem,rownames(mat))]-1)/(nrow(mat)-1)),0)
    direct<-s06_score_estimate(list(estimable=TRUE,score=unname(score),sample_id=colnames(mat)),p)
    s06_compare(direct$source_data,a$estimate$source_data,paste(r$analysis_id,'metadata-to-source computation'))
   }
  }
 }
 invisible(TRUE)
}
  run <- function() {
    context$guard();context$resource_check()
    input<-s06_save(s06_inputs(),s06_o('validated_inputs.rds'))
    scope<-s06_scope()
    states<-s06_state_target(input,TRUE,scope)
    members<-vapply(c('myeloid','fibroblast','B_plasma'),function(comp)s06_member_target(comp,input,TRUE),'',USE.NAMES=FALSE)
    composition<-s06_composition_target(input,TRUE)
    assembled<-s06_assemble(states,members,composition,scope);result<-read_rds(assembled)
    validation<-s06_validate(result);s06_validate_metadata(result$diagnostics,read_rds(states),read_rds(input))
    tables<-s06_tables(result)
    out<-vapply(names(tables),function(n){p<-context$output('results',paste0(n,'.tsv'));context$write_tsv(tables[[n]],p);s06_compare(tables[[n]],s06_read(p),n);p},'')
    context$write_rds(validation,context$output('qc','numerical_validation.rds'))
    list(scope='original_single_cell',tables=out,objects=vapply(c(input,states,members,composition,assembled),resolve_path,''),validation=validation)
  }
  list(run=run,input=s06_inputs,project=s06_tables,validate=s06_validate,
       kernels=list(native_members=s06_native_members,extract_native=s06_extract_native,score_estimate=s06_score_estimate),resolve=resolve_path)
}
