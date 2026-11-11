# Portable replay of registered PE2 computations from explicit released values.
# Accepted estimates are never inputs. Comparison is a separate subsequent command.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L || !args[1] %in% c('B','D','E','test_B','test_DE'))
  stop('Usage: second_post_effect.R B|D|E|test_B|test_DE INPUT_MANIFEST.tsv OUTPUT_DIRECTORY')
if(!file.exists('PUBLIC_MANIFEST.tsv'))stop('Run from the public release root.')
source('workflow/R/common/load_analysis.R');load_public_analysis()
source('workflow/R/methods/second_post_effect.R')
stage<-args[1];cx<-public_analysis_context('second_post_effect',args[2],args[3],create=FALSE)
cx$resource_check();runtime<-cx$runtime(if(stage %in% c('B','test_B'))'epithelial' else 'supplementary')
R<-'results/supplementary/second_post_effect/';D<-'results/diagnostics/supplementary/second_post_effect/'
tab<-function(name,diagnostic=FALSE)cx$read_tsv(paste0(if(diagnostic)D else R,name))
write_tsv<-function(x,n,role='results')cx$write_tsv(x,cx$output(role,n))
write_json<-function(x,n,role='qc')cx$write_json(x,cx$output(role,n))
same<-function(x,y)stopifnot(isTRUE(all.equal(unname(x),unname(y),tolerance=1e-12,check.attributes=FALSE)))
seeds<-cx$read_tsv('workflow/config/implementation/second_post_effect_seed_map.tsv')
seed<-function(k){i<-grep(paste0('^PE2_',k,'_'),seeds$operation_name);stopifnot(length(i)==1);seeds$seed[i]}
roster<-cx$read_tsv('workflow/config/second_post_effect_endpoints.tsv')
H<-c('ALDH2','ALDH3A1','ALDH3A2','ALDH7A1')
close_family<-function(x) {
 ids<-roster$hypothesis_id[roster$family==unique(x$family)]
 stopifnot(identical(x$hypothesis_id,ids),!anyDuplicated(ids))
 z<-pe_close_family(setNames(x$p_value,ids),x$state,ids);stopifnot(z$closed)
 x$q_value<-unname(z$q);x$family_closed<-TRUE;x
}
if(stage %in% c('test_B','test_DE')) {
 source('workflow/tests/test_second_post_effect.R')
} else if(stage=='B') {
 v<-tab('B_member_values.tsv');ref<-tab('B_reference.tsv',TRUE);norm<-tab('B_normalization.tsv',TRUE)
 ids<-sprintf('Pat%02d',c(1:20,22:24));units<-as.vector(rbind(paste0(ids,'|NAT'),paste0(ids,'|Tumor')))
 stopifnot(identical(ref$member,H),identical(norm$unit_id,units),nrow(v)==184L,
   all(ref$reference_units==46L),all(ref$SD_denominator==45L),all(norm$prior_count==2),
   all(norm$TMM_reference_unit=='Pat11|Tumor'),all(norm$effective_library_size>0))
 same(norm$effective_library_size,norm$library_size*norm$TMM_factor)
 x<-t(vapply(H,function(g){a<-v[v$member==g,];stopifnot(identical(a$unit_id,units));a$logCPM},numeric(46)))
 dimnames(x)<-list(H,units);same(rowMeans(x),ref$pooled_mean);same(apply(x,1,sd),ref$sample_SD)
 z<-sweep(sweep(x,1,ref$pooled_mean,'-'),1,ref$sample_SD,'/')
 for(i in 1:4)same(z[i,],v$z[v$member==H[i]])
 scores<-colMeans(z);n<-scores[seq(1,46,2)];t<-scores[seq(2,46,2)]
 est<-s01_paired_delta_t(t,n,ids);g<-s01_paired_g_av(t,n,ids);indices<-pe2_indices(ids,seed('B'))
 draws<-apply(indices,2,function(j)mean((t-n)[j]));ci<-pe2_interval(draws)
 write_tsv(data.frame(donor=ids,NAT_unit=units[seq(1,46,2)],Tumor_unit=units[seq(2,46,2)],NAT_score=n,Tumor_score=t,delta=t-n),'B_donor_values.tsv')
 write_tsv(data.frame(omitted_donor=ids,effect=vapply(seq_along(ids),function(i)mean((t-n)[-i]),0),refit_n=22L),'B_omissions.tsv','diagnostics')
 write_tsv(data.frame(draw=1:10000,effect=draws),'B_bootstrap_draws.tsv','diagnostics')
 result<-data.frame(hypothesis_id='PE2_B_H4_TMM_LOGCPM_Z',role='TECHNICAL_SCORING_SENSITIVITY',family='NO_FORMAL_FAMILY',biological_n=23L,effect=est$effect,ci_lower=est$ci_lower,ci_upper=est$ci_upper,p_value=est$p_value,state=est$state,bootstrap_lower=ci$lower,bootstrap_upper=ci$upper,valid_draws=ci$valid,bootstrap_state=ci$state,units='mean_member_z',timing='post_effect_exploratory')
 write_tsv(result,'B_diagnostic.tsv');write_json(g,'B_g_av.json','diagnostics')
} else if(stage=='D') {
 values<-tab('D_patient_deltas.tsv');moments<-tab('D_frozen_reference_moments.tsv',TRUE)
 stopifnot(nrow(values)==113L,identical(unique(values$source),c('GSE122401','TCGA_STAD')),nrow(moments)==14L)
 results<-loo<-corr<-checks<-list()
 for(src in c('GSE122401','TCGA_STAD')) {
  d<-values[values$source==src,];ids<-d$patient_id;h<-d$delta_H4_rank
  stopifnot(length(ids)==if(src=='GSE122401')80L else 33L,!anyDuplicated(ids),!anyDuplicated(d$adjacent_sample),!anyDuplicated(d$tumor_sample),!any(d$adjacent_sample %in% d$tumor_sample))
  delta<-t(as.matrix(d[c('D_pit4','D_gland2','PGC_z')]));colnames(delta)<-ids
  stopifnot(all(is.finite(delta)),all(is.finite(h)))
  c<-cor(t(delta),method='spearman');ij<-which(upper.tri(c),arr.ind=TRUE)
  corr[[src]]<-data.frame(source=src,component_1=rownames(c)[ij[,1]],component_2=colnames(c)[ij[,2]],rho=c[ij],biological_n=length(ids),family='NO_FORMAL_FAMILY')
  for(k in rownames(delta)) {
   id<-paste0('PE2_D_',src,'_',k,'_beta_component');fit<-pe2_slope(delta[k,],h,ids)
   results[[id]]<-data.frame(hypothesis_id=id,source=src,component=k,as.data.frame(fit[c('state','reason','biological_n','effect','standard_error','statistic','residual_df','ci_lower','ci_upper','p_value')]),family='PE2_D_COMPONENT_ASSOC',family_size=6L,units='H4_rank_change_per_fixed_component_unit',timing='post_effect_exploratory')
   loo[[id]]<-do.call(rbind,lapply(seq_along(ids),function(i){f<-pe2_slope(delta[k,-i],h[-i],ids[-i],minimum_n=1L);data.frame(hypothesis_id=id,omitted_patient=ids[i],effect=f$effect,state=f$state,reason=f$reason,refit_n=length(ids)-1L)}))
   checks[[id]]<-list(diagnostics=fit$diagnostics,excluded=fit$excluded)
  }
  cx$resource_check(6,src)
 }
 write_tsv(values,'D_patient_deltas.tsv');write_tsv(close_family(do.call(rbind,results)),'D_results.tsv')
 write_tsv(do.call(rbind,loo),'D_omissions.tsv','diagnostics');write_tsv(do.call(rbind,corr),'D_component_correlations.tsv','diagnostics')
 write_json(checks,'D_model_diagnostics.json','diagnostics');write_tsv(moments,'D_frozen_reference_moments.tsv','diagnostics')
 summary<-do.call(rbind,lapply(names(loo),function(id){a<-loo[[id]];v<-a$state=='ESTIMABLE'&is.finite(a$effect);main<-results[[id]]$effect;data.frame(hypothesis_id=id,omissions=nrow(a),valid_omissions=sum(v),invalid_omissions=sum(!v),minimum=if(any(v))min(a$effect[v]) else NA_real_,maximum=if(any(v))max(a$effect[v]) else NA_real_,direction_changes=sum(sign(a$effect[v])!=sign(main)))}))
 write_tsv(summary,'D_influence_summary.tsv','diagnostics')
} else if(stage=='E') {
 v<-tab('E_paired_values.tsv');join<-tab('E_join_population.tsv',TRUE);ratios<-tab('E_ratio_identity.tsv',TRUE);bridge<-tab('E_aliquot_RNA_bridge.tsv',TRUE)
 ids<-join$patient_id;stopifnot(length(ids)==80L,!anyDuplicated(ids),identical(join$analytical_order,1:80),all(join$join_state=='COMMON_COMPLETE_H4'),!any(join$excluded),sum(join$declared_technical_ratios)==130L, sum(join$declared_technical_ratios==2L)==50L,sum(join$declared_technical_ratios==1L)==30L,nrow(v)==320L,identical(unique(v$member),H),nrow(ratios)==130L)
 x<-y<-matrix(NA_real_,4,80,dimnames=list(H,ids))
 for(g in seq_along(H)) {
  a<-v[v$member==H[g],];stopifnot(identical(a$patient_id,ids),identical(a$technical_ratios,join$declared_technical_ratios),all(a$biological_unit=='patient'))
  same(a$RNA_delta,a$RNA_log2_UQ_plus1_Tumor-a$RNA_log2_UQ_plus1_NAT)
  x[g,]<-a$RNA_delta;y[g,]<-a$protein_mean_centered_log2_T_N
 }
 stopifnot(identical(as.integer(table(factor(ratios$patient_id,levels=ids))),join$declared_technical_ratios))
 for(i in seq_len(nrow(ratios))) {
  a<-ratios[i,];tt<-bridge[bridge$aliquot_submitter_id==a$tumor_aliquot,];nn<-bridge[bridge$aliquot_submitter_id==a$adjacent_aliquot,]
  stopifnot(nrow(tt)==1L,nrow(nn)==1L,tt$patient_id==a$patient_id,nn$patient_id==a$patient_id,tt$sample_type=='Tumor',nn$sample_type=='Solid Tissue Normal')
 }
 indices<-pe2_indices(ids,seed('E'));z<-pe2_bridge(x,y,ids,indices)
 z$results$hypothesis_id<-paste0('PE2_E_MUN_',z$results$member,'_rho');z$results$family<-'PE2_E_MUN_MEMBER_BRIDGE';z$results$family_size<-4L;z$results$units<-'Spearman_rho';z$results$timing<-'post_effect_exploratory'
 write_tsv(v,'E_paired_values.tsv');write_tsv(close_family(z$results),'E_results.tsv');write_tsv(z$omissions,'E_omissions.tsv','diagnostics');write_tsv(as.data.frame(z$draws),'E_bootstrap_draws.tsv','diagnostics')
 summary<-do.call(rbind,lapply(seq_len(4),function(g){a<-z$omissions[z$omissions$member==z$results$member[g],];good<-is.finite(a$rho);data.frame(member=z$results$member[g],omissions=nrow(a),valid_omissions=sum(good),invalid_omissions=sum(!good),minimum=if(any(good))min(a$rho[good]) else NA_real_,maximum=if(any(good))max(a$rho[good]) else NA_real_,direction_changes=sum(sign(a$rho[good])!=sign(z$results$effect[g])))}))
 write_tsv(summary,'E_influence_summary.tsv','diagnostics')
 write_tsv(join,'E_join_population.tsv','diagnostics');write_tsv(ratios,'E_ratio_identity.tsv','diagnostics');write_tsv(bridge,'E_aliquot_RNA_bridge.tsv','diagnostics')
}
write_json(list(status='COMPUTATION_COMPLETED_NOT_YET_COMPARED',stage=stage,runtime=runtime,session=capture.output(sessionInfo()),new_inference=FALSE,starting_point=if(stage=='B')'normalized_member_values' else if(stage=='D')'patient_component_and_H_deltas' else if(stage=='E')'paired_member_values' else 'synthetic'),paste0(stage,'_run.json'))
write_tsv(cx$verify_inputs(),'verified_input_identities.tsv','qc');cx$guard();cx$resource_check(6,'stage_end')
cat('Completed',stage,'; accepted reference comparison is a separate command.\n')
