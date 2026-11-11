# Registered PE2 methods. Pure functions; no real-data execution on source().
pe2_indices <- function(ids, seed) {
  s01_validate_identity_vector(ids, 'PE2 analytical donor order')
  s01_with_production_rng(seed, replicate(10000L, sample.int(length(ids), length(ids), replace=TRUE)))
}
pe2_interval <- function(draws) {
  ok <- is.finite(draws)
  ci <- if(sum(ok)>=9500L) quantile(draws[ok],c(.025,.975),type=7,names=FALSE) else c(NA_real_,NA_real_)
  list(valid=sum(ok),lower=ci[1],upper=ci[2],state=if(sum(ok)>=9500L)'ESTIMABLE' else 'NA_NONESTIMABLE')
}
pe2_normalized_h <- function(counts, members) {
  s01_validate_numeric_matrix(counts,'PE2 full modeling counts',TRUE)
  s01_require(length(members)==4L && all(members %in% rownames(counts)),'PE2 exact H4')
  y <- edgeR::DGEList(counts=counts)
  s01_require(all(y$samples$lib.size>0),'PE2 positive libraries')
  # Mirror edgeR's default reference selection only for audit, never override it.
  f75 <- edgeR:::.calcFactorQuantile(counts,lib.size=y$samples$lib.size,p=.75)
  ref <- if(median(f75)<1e-20) which.max(colSums(sqrt(counts))) else which.min(abs(f75-mean(f75)))
  y <- edgeR::calcNormFactors(y,method='TMM')
  eff <- y$samples$lib.size*y$samples$norm.factors
  s01_require(all(is.finite(eff)&eff>0),'PE2 effective libraries')
  logcpm <- edgeR::cpm(y,log=TRUE,prior.count=2,normalized.lib.sizes=TRUE)
  score <- s01_standardized_member_mean(logcpm,members)
  x <- logcpm[members,,drop=FALSE]; mu<-rowMeans(x); sig<-apply(x,1,sd)
  z <- if(score$estimable) sweep(sweep(x,1,mu,'-'),1,sig,'/') else x*NA_real_
  list(score=score,logcpm=x,z=z,mean=mu,sd=sig,normalization=data.frame(unit_id=colnames(counts),library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors,effective_library_size=eff,TMM_reference_unit=colnames(counts)[ref],prior_count=2),reference_index=ref)
}
pe2_components <- function(expression, reference) {
  s01_require(reference$state=='READY' && identical(colnames(expression),reference$metadata$sample_id),'PE2 frozen D7 population')
  x<-expression[reference$members,,drop=FALSE]
  s01_require(all(is.finite(x)) && all(is.finite(reference$sd)&reference$sd>0),'PE2 complete D7 measurements')
  z<-sweep(sweep(x,1,reference$mean,'-'),1,reference$sd,'/')
  rbind(D_pit4=colMeans(z[c('GKN1','GKN2','MUC5AC','TFF1'),,drop=FALSE]),D_gland2=colMeans(z[c('MUC6','TFF2'),,drop=FALSE]),PGC_z=z['PGC',])
}
pe2_slope <- function(x,y,ids,minimum_n=20L) {
  s01_validate_identity_vector(ids,'PE2 slope patients')
  s01_require(length(x)==length(ids)&&length(y)==length(ids),'PE2 slope dimensions')
  good<-is.finite(x)&is.finite(y); n<-sum(good)
  out<-list(state='NA_NONESTIMABLE',reason='INSUFFICIENT_COMPLETE_PAIRS',biological_n=n,effect=NA_real_,standard_error=NA_real_,statistic=NA_real_,residual_df=n-2L,ci_lower=NA_real_,ci_upper=NA_real_,p_value=NA_real_,excluded=ids[!good])
  if(n<minimum_n) return(out)
  x<-x[good];y<-y[good]
  if(s01_effective_zero(x)){out$reason<-'CONSTANT_COMPONENT';return(out)}
  fit<-lm(y~x); h<-pe_hc3(fit)
  if(h$state!='AVAILABLE'){out$reason<-h$reason;return(out)}
  out$effect<-unname(coef(fit)[2]);out$standard_error<-unname(h$se[2]);out$statistic<-out$effect/out$standard_error
  out$ci_lower<-unname(h$lower[2]);out$ci_upper<-unname(h$upper[2]);out$p_value<-2*pt(-abs(out$statistic),df=n-2)
  if(!all(is.finite(unlist(out[c('effect','standard_error','statistic','ci_lower','ci_upper','p_value')])))){out$reason<-'INVALID_HC3_INFERENCE';out$p_value<-NA_real_;return(out)}
  out$state<-'ESTIMABLE';out$reason<-'NONE';out$fit<-fit
  out$diagnostics<-list(condition_number=kappa(model.matrix(fit),exact=TRUE),standardized_condition_number=kappa(cbind(1,scale(x)),exact=TRUE),rank=fit$rank,max_leverage=max(hatvalues(fit)),zero_component_supported=min(x)<=0&&max(x)>=0)
  out
}
pe2_bridge <- function(rna,protein,ids,indices) {
  s01_validate_identity_vector(ids,'PE2 bridge donors')
  s01_require(is.matrix(rna)&&identical(dimnames(rna),dimnames(protein))&&identical(colnames(rna),ids)&&nrow(rna)==4L,'PE2 common bridge order')
  s01_require(all(is.finite(rna))&&all(is.finite(protein))&&ncol(rna)==length(ids),'PE2 complete common bridge')
  s01_require(identical(dim(indices),c(length(ids),10000L))&&all(indices>=1L&indices<=length(ids)),'PE2 shared index matrix')
  draws<-matrix(NA_real_,10000L,4L,dimnames=list(NULL,rownames(rna))); rows<-loo<-list()
  for(g in seq_len(4L)) {
    x<-rna[g,];y<-protein[g,];n<-length(ids)
    eligible<-n>=20L&&!s01_effective_zero(x)&&!s01_effective_zero(y)&&length(unique(x))>=5L&&length(unique(y))>=5L
    rho<-if(sd(x)>0&&sd(y)>0)cor(x,y,method='spearman') else NA_real_
    if(eligible) for(b in seq_len(10000L)) {j<-indices[,b];draws[b,g]<-suppressWarnings(cor(x[j],y[j],method='spearman'))}
    ci<-pe2_interval(draws[,g]); valid<-eligible&&ci$state=='ESTIMABLE'&&is.finite(rho)
    rows[[g]]<-data.frame(member=rownames(rna)[g],biological_n=n,effect=rho,ci_lower=ci$lower,ci_upper=ci$upper,p_value=if(valid)cor.test(x,y,method='spearman',exact=FALSE)$p.value else NA_real_,state=if(valid)'ESTIMABLE' else 'NA_NONESTIMABLE',reason=if(valid)'NONE' else if(!eligible)'INSUFFICIENT_VARIATION_OR_UNITS' else 'INVALID_BOOTSTRAP_CI_RHO_DESCRIPTIVE',valid_draws=ci$valid)
    loo[[g]]<-data.frame(member=rownames(rna)[g],omitted_patient=ids,rho=vapply(seq_along(ids),function(i)suppressWarnings(cor(x[-i],y[-i],method='spearman')),0),refit_n=n-1L)
    loo[[g]]$state<-ifelse(is.finite(loo[[g]]$rho),'ESTIMABLE','NA_NONESTIMABLE');loo[[g]]$reason<-ifelse(is.finite(loo[[g]]$rho),'NONE','CONSTANT_OMISSION_VARIABLE')
  }
  list(results=do.call(rbind,rows),draws=draws,omissions=do.call(rbind,loo))
}
