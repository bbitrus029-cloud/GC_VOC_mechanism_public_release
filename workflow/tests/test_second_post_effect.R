# Focused independent anchors; only synthetic values, never outcome-tuned.
if(stage=='test_B') {
  counts<-outer(seq_len(30),seq_len(8),function(g,s)(g*17+s*13)%%97+1)
  counts<-counts*rep(c(1,7,2,4,9,1,5,3),each=30)
  dimnames(counts)<-list(paste0('g',1:30),paste0('s',1:8));members<-rownames(counts)[1:4]
  z<-pe2_normalized_h(counts,members)
  native<-edgeR::calcNormFactors(edgeR::DGEList(counts=counts));x<-edgeR::cpm(native,log=TRUE,prior.count=2)[members,]
  expected<-rowMeans(scale(t(x),center=TRUE,scale=TRUE))
  stopifnot(isTRUE(all.equal(z$score$score,unname(expected),tolerance=1e-13)),isTRUE(all.equal(z$normalization$TMM_factor,native$samples$norm.factors,tolerance=1e-13)))
  # Member-only normalization would change the answer and must be detectable.
  wrong<-edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=counts[members,])),log=TRUE,prior.count=2)
  stopifnot(max(abs(x-wrong))>1e-3)
  ids<-paste0('p',1:4);idx<-pe2_indices(ids,211203342);d<-z$score$score[c(2,4,6,8)]-z$score$score[c(1,3,5,7)]
  draws<-apply(idx,2,function(j)mean(d[j]));ci<-pe2_interval(draws)
  native_boot<-s01_pair_bootstrap(z$score$score[c(2,4,6,8)],z$score$score[c(1,3,5,7)],ids,211203342)
  stopifnot(identical(ci$lower,native_boot$ci_lower),identical(ci$upper,native_boot$ci_upper),ci$valid==10000L)
  tst<-s01_paired_delta_t(z$score$score[c(2,4,6,8)],z$score$score[c(1,3,5,7)],ids);anchor<-t.test(d)
  stopifnot(abs(tst$p_value-anchor$p.value)<1e-13,max(abs(c(tst$ci_lower,tst$ci_upper)-anchor$conf.int))<1e-13)
  write_json(list(status='PASS',checks=c('native_full_background_TMM_logCPM','sample_SD_equal_mean','reject_member_only_context','inherited_bootstrap_order_and_interval','paired_t_native_anchor')),'test_B_PASS.json','qc')
} else {
  x<-seq(-2,3,length.out=25);y<-.7*x+sin(1:25)*(1+x^2);ids<-paste0('p',1:25)
  z<-pe2_slope(x,y,ids);f<-lm(y~x);deleted<-vapply(1:25,function(i)coef(lm(y[-i]~x[-i])),numeric(2));dif<-sweep(deleted,1,coef(f),'-');v<-tcrossprod(dif)
  stopifnot(isTRUE(all.equal(unname(pe_hc3(f)$covariance),unname(v),tolerance=1e-12)),abs(z$p_value-2*pt(-abs(coef(f)[2]/sqrt(v[2,2])),23))<1e-13)
  stopifnot(pe2_slope(rep(1,25),y,ids)$state=='NA_NONESTIMABLE',pe2_slope(x[1:19],y[1:19],ids[1:19])$state=='NA_NONESTIMABLE')
  members<-c('GKN1','GKN2','MUC5AC','TFF1','MUC6','TFF2','PGC');m<-matrix(1:70,7,dimnames=list(members,paste0('s',1:10)))
  ref<-list(state='READY',members=members,metadata=data.frame(sample_id=colnames(m)),mean=setNames(1:7,members),sd=setNames(2:8,members));c<-pe2_components(m,ref)
  stopifnot(identical(rownames(m),members),identical(names(ref$mean),members),identical(names(ref$sd),members),identical(dim(c),c(3L,10L)),identical(rownames(c),c('D_pit4','D_gland2','PGC_z')),identical(colnames(c),ref$metadata$sample_id))
  stopifnot(max(abs(c['PGC_z',]-(m['PGC',]-7)/8))<1e-13,max(abs(c['D_pit4',]-colMeans((m[1:4,]-1:4)/(2:5))))<1e-13)
  # Deliberately shuffled donor labels: preserve the accepted vector order.
  ids<-paste0('p',c(21:25,1:20));rna<-rbind(round(x),x,x^2,sin(x));rownames(rna)<-paste0('g',1:4);colnames(rna)<-ids;protein<-rna+matrix(cos(1:100),4);dimnames(protein)<-dimnames(rna)
  idx<-pe2_indices(ids,134607028);bridge<-pe2_bridge(rna,protein,ids,idx)
  stopifnot(identical(dimnames(rna),list(paste0('g',1:4),ids)),identical(dimnames(protein),dimnames(rna)),identical(dim(idx),c(25L,10000L)),identical(dim(bridge$draws),c(10000L,4L)),identical(colnames(bridge$draws),rownames(rna)),identical(bridge$results$member,rownames(rna)))
  native_idx<-s01_with_production_rng(134607028,replicate(10000L,sample.int(length(ids),length(ids),replace=TRUE)))
  stopifnot(identical(idx,native_idx),identical(ids,paste0('p',c(21:25,1:20))))
  for(member in rownames(rna))stopifnot(identical(bridge$omissions$omitted_patient[bridge$omissions$member==member],ids))
  stopifnot(nrow(bridge$omissions)==100,all(bridge$results$valid_draws==10000))
  for(g in 1:4){j<-idx[,21];stopifnot(identical(unname(bridge$draws[21,g]),unname(cor(rank(rna[g,j]),rank(protein[g,j])))),abs(bridge$results$p_value[g]-cor.test(rna[g,],protein[g,],method='spearman',exact=FALSE)$p.value)<1e-13)}
  ci<-pe2_interval(c(rep(1,9499),rep(NA_real_,501)));stopifnot(ci$state=='NA_NONESTIMABLE',is.na(ci$lower))
  p<-setNames(c(.01,.04,NA,.5),letters[1:4]);z<-pe_close_family(p,c('ESTIMABLE','ESTIMABLE','NA_UNAVAILABLE','ESTIMABLE'),names(p));stopifnot(z$closed,identical(unname(z$q[c(1,2,4)]),p.adjust(p[c(1,2,4)],'BH',n=4)|>unname()),is.na(z$q[3]))
  stopifnot(!pe_close_family(p,c('ESTIMABLE','ESTIMABLE','PENDING','ESTIMABLE'),names(p))$closed)
  vals<-matrix(1:12,4,dimnames=list(paste0('g',1:4),c('r1','r2','r3')));agg<-pe_protein_patient_ratios(vals,data.frame(ratio_column=colnames(vals),patient_id=c('p2','p1','p2')))
  stopifnot(identical(agg$patient_ids,c('p1','p2')),identical(dim(agg$values),c(4L,2L)),identical(rownames(agg$values),rownames(vals)),identical(colnames(agg$values),agg$patient_ids))
  stopifnot(identical(unname(agg$values[,'p2']),unname(rowMeans(vals[,c(1,3)]))),identical(unname(agg$values[,'p1']),as.numeric(vals[,'r2'])))
  stopifnot(identical(agg$excluded,character()),identical(agg$rule,'all_declared_replicates_of_all_H4_required_equal_mean_of_log2_ratios'))
  write_json(list(status='PASS',checks=c('HC3_independent_delete_one_covariance_identity','HC3_t_n_minus_2','fixed_component_moments_and_labels','constant_and_small_n_failure','shared_donor_resampling_with_ties_and_nonlexical_order','exact_gene_patient_dimension_order','native_Spearman_P','invalid_CI','fixed_BH_and_pending_closure','equal_mean_all_technical_ratios_and_patient_assignment')),'test_DE_PASS.json','qc')
}
