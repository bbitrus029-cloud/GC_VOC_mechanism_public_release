# Patient-score-to-estimate replay. Reference effect tables are never read here.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('Usage: workflow/environment/run-r workflow/pipelines/reproduce_core_bulk.R OUTPUT_DIRECTORY')
root <- normalizePath(getwd(), mustWork=TRUE)
if (!file.exists('PUBLIC_MANIFEST.tsv')) stop('Run from the public release root.')
for (p in c('contracts','statistics','public_paths')) source(paste0('workflow/R/common/',p,'.R'))
out <- resolve_output_directory(args[1],root)
read_input <- function(path) {
  if (!file.exists(path)) stop('Missing declared public input: ',path,'; no fallback or download.')
  m <- read.delim('PUBLIC_MANIFEST.tsv',quote='',colClasses='character')
  j <- match(path,m$path)
  if (is.na(j) || digest::digest(file=path,algo='sha256') != m$sha256[j]) stop('Public input integrity mismatch: ',path)
  read.delim(path,quote='',comment.char='',check.names=FALSE,colClasses='character')
}
scores <- read_input('results/original/bulk/patient_score_deltas.tsv')
cw <- read_input('workflow/config/original_hypothesis_estimand_crosswalk.tsv')
cfg <- yaml::read_yaml('workflow/config/original_analysis.yml')
spec <- data.frame(source=c('GSE29272','GSE29272','GSE27342_64','GSE27342_64'),
 endpoint_id=c('GSE29272_G','GSE29272_H','GSE27342_primary_G','GSE27342_primary_H'),
 n=c(112L,112L,49L,49L), family=rep(c('BULK_CORE_GSE29272','BULK_CORE_GSE27342_64'),each=2))
stopifnot(all(c('patient_id','NAT_score','Tumor_score','delta','source','analysis_id') %in% names(scores)))
rows <- list(); patients <- list(); checks <- list()
for (i in seq_len(nrow(spec))) {
  e <- spec[i,]; x <- scores[scores$source==e$source & scores$analysis_id==e$endpoint_id,,drop=FALSE]
  s01_validate_identity_vector(x$patient_id,'core patient')
  stopifnot(nrow(x)==e$n)
  for (k in c('NAT_score','Tumor_score','delta')) x[[k]] <- as.numeric(x[[k]])
  stopifnot(all(is.finite(as.matrix(x[c('NAT_score','Tumor_score','delta')]))))
  error <- abs(x$delta-(x$Tumor_score-x$NAT_score))
  stopifnot(all(error <= 1e-13*pmax(1,abs(x$Tumor_score),abs(x$NAT_score))))
  accession <- if(e$source=='GSE29272') 'GSE29272' else 'GSE27342'
  ledger <- read_input(paste0('data/metadata/',accession,'/qualification/',accession,'_pair_ledger.tsv'))
  eligible <- ledger[ledger$pair_excluded=='FALSE' & ledger$core_population=='TRUE',,drop=FALSE]
  stopifnot(!anyDuplicated(paste(eligible$patient_id,eligible$tissue)),
    all(table(eligible$patient_id)==2L),setequal(eligible$tissue,c('NAT','Tumor')),
    setequal(x$patient_id,eligible$patient_id))
  patients[[e$endpoint_id]] <- sort(x$patient_id,method='radix')
  r <- cw[cw$endpoint_id==e$endpoint_id,,drop=FALSE]
  stopifnot(nrow(r)==1L,r$source==e$source,r$multiplicity_family==e$family,
    r$frozen_family_size=='2',r$estimand=='paired_tumor_minus_NAT_program_score_difference')
  t <- s01_paired_delta_t(x$Tumor_score,x$NAT_score,x$patient_id)
  g <- s01_paired_g_av(x$Tumor_score,x$NAT_score,x$patient_id)
  rows[[i]] <- data.frame(endpoint_id=r$endpoint_id,hypothesis_id=r$hypothesis_id,source=e$source,
    analysis_id=e$endpoint_id,multiplicity_family=e$family,family_hypothesis=r$family_hypothesis,
    frozen_family_size=2L,biological_unit='patient',biological_n=t$biological_n,
    contrast='Tumor_minus_paired_NAT',effect=t$effect,standard_error=sd(x$Tumor_score-x$NAT_score)/sqrt(nrow(x)),
    statistic=t$statistic,degrees_of_freedom=t$degrees_of_freedom,ci_lower=t$ci_lower,ci_upper=t$ci_upper,
    p_value=t$p_value,q_value=NA_real_,estimability=t$state,g_av=g$g_av,g_av_variance=g$variance,
    g_av_ci_lower=g$ci_lower,g_av_ci_upper=g$ci_upper,g_av_state=g$state)
  checks[[i]] <- data.frame(endpoint_id=e$endpoint_id,unique_patients=nrow(x),pairing_ledger_match=TRUE,
    max_saved_delta_difference=max(error),contrast='Tumor_minus_paired_NAT',result='PASS')
}
stopifnot(identical(patients[[1]],patients[[2]]),identical(patients[[3]],patients[[4]]))
r <- do.call(rbind,rows)
for (f in unique(r$multiplicity_family)) {
  j <- which(r$multiplicity_family==f); frozen <- cfg$multiplicity$families[[f]]
  stopifnot(frozen$frozen_size==2L,identical(r$family_hypothesis[j],frozen$hypotheses))
  r$q_value[j] <- unname(s01_frozen_bh_values(setNames(r$p_value[j],r$family_hypothesis[j]),2L,frozen$hypotheses))
}
# The validated output directory is created only after all input and population checks.
if (!dir.exists(out) && !dir.create(out,recursive=TRUE,showWarnings=FALSE)) stop('Output path creation failed: ',out)
if (!dir.exists(out)) stop('Output path was not created: ',out)
write_precise <- function(x,path) {
  for (k in names(x)) if(is.double(x[[k]])) x[[k]] <- ifelse(is.na(x[[k]]),'NA',sprintf('%.17g',x[[k]]))
  write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
}
write_precise(r,file.path(out,'core_bulk_estimates.tsv'))
write_precise(do.call(rbind,checks),file.path(out,'input_checks.tsv'))
cat('Patient-score-to-estimate replay complete: four estimates; 112 and 49 patient pairs; two separate BH families.\n')
