# Render retained descriptive plots from bundled permitted tables. No fitting/scoring.
args <- commandArgs(TRUE)
if (length(args) != 2L || !args[1] %in% c('tissue_context','protein'))
  stop('Usage: Rscript --vanilla workflow/pipelines/render.R tissue_context|protein OUTPUT_DIRECTORY')
root <- normalizePath(getwd(), mustWork=TRUE)
if (!file.exists('PUBLIC_MANIFEST.tsv')) stop('Run from the release directory.')
source('workflow/R/methods/accepted_plotting.R')
accepted_render_validate <- function(scope, tables) {
  columns <- function(x, required) {
    if (!all(required %in% names(x))) stop('Missing accepted source columns: ',
      paste(setdiff(required, names(x)), collapse = ', '))
  }
  unique_key <- function(x) {
    if (anyNA(x) || anyDuplicated(x)) stop('Missing or duplicate accepted source identity')
  }
  finite <- function(x) if (any(!is.finite(as.matrix(x)))) stop('Nonfinite required plotting value')
  if (scope == 'bridge') {
    d <- tables$pairs; r <- tables$estimates
    columns(d, c('endpoint_id', 'patient_id', 'RNA_delta', 'protein_delta'))
    columns(r, c('endpoint_id', 'biological_n', 'rho', 'ci_lower', 'ci_upper', 'q_value'))
    unique_key(r$endpoint_id)
    stopifnot(setequal(r$endpoint_id, c('CPTAC_bridge_G', 'CPTAC_bridge_H')),
      setequal(d$endpoint_id, r$endpoint_id))
    for (id in r$endpoint_id) {
      x <- d[d$endpoint_id == id, ]; unique_key(x$patient_id)
      stopifnot(nrow(x) == r$biological_n[r$endpoint_id == id])
    }
    finite(d[, c('RNA_delta', 'protein_delta')])
  } else if (scope == 'epithelial') {
    columns(tables$program, c('source', 'endpoint', 'biological_n', 'state', 'effect', 'ci_lower', 'ci_upper'))
    columns(tables$member, c('endpoint', 'member', 'state', 'effect'))
    columns(tables$composition, c('endpoint', 'raw_proportion_bootstrap_state',
      'raw_paired_proportion_difference', 'raw_proportion_bootstrap_lower', 'raw_proportion_bootstrap_upper'))
    columns(tables$donors, c('donor', 'NAT_score', 'Tumor_score'))
    stopifnot(nrow(tables$program) == 3L, nrow(tables$member) == 8L,
      nrow(tables$composition) == 2L, nrow(tables$donors) == 23L)
    unique_key(tables$donors$donor)
    finite(tables$donors[, c('NAT_score', 'Tumor_score')])
  } else if (scope == 'tissue_context') {
    columns(tables$pairs, c('source', 'patient_id', 'delta_D', 'delta_H'))
    columns(tables$support, c('source', 'biological_n', 'zero_supported', 'intercept_label'))
    unique_key(tables$support$source)
    stopifnot(setequal(tables$support$source, c('GSE122401', 'TCGA_STAD')),
      setequal(tables$pairs$source, c('GSE29272', 'GSE27342_core', 'GSE122401', 'TCGA_STAD')),
      all(tables$support$source %in% tables$pairs$source))
    for (s in tables$support$source) {
      x <- tables$pairs[tables$pairs$source == s, ]; line <- tables[[s]]
      columns(line, c('source', 'delta_D', 'fit', 'ci_lower', 'ci_upper'))
      unique_key(x$patient_id)
      stopifnot(nrow(x) == tables$support$biological_n[tables$support$source == s],
        nrow(line) == 100L, all(line$source == s), all(diff(line$delta_D) > 0))
      finite(x[, c('delta_D', 'delta_H')]); finite(line[, c('delta_D', 'fit', 'ci_lower', 'ci_upper')])
    }
  } else if (scope == 'protein') {
    d <- tables$members
    columns(d, c('hypothesis_id', 'family_hypothesis', 'biological_n', 'effect', 'ci_lower', 'ci_upper', 'state'))
    unique_key(d$hypothesis_id)
    stopifnot(nrow(d) == 4L, all(d$biological_n == 80L), all(d$state == 'ESTIMABLE'))
    finite(d[, c('effect', 'ci_lower', 'ci_upper')])
  } else stop('Unknown retained plotting scope: ', scope)
  invisible(TRUE)
}


scope <- args[1]
paths <- switch(scope,
 tissue_context=c(pairs='results/supplementary/tissue_context/patient_paired_changes.tsv',
 support='results/diagnostics/supplementary/tissue_context/intercept_support.tsv',
 GSE122401='results/supplementary/tissue_context/GSE122401_fitted_line_source.tsv',
 TCGA_STAD='results/supplementary/tissue_context/TCGA_STAD_fitted_line_source.tsv'),
 protein=c(members='results/supplementary/protein/member_figure_source.tsv'))
manifest <- read.delim('PUBLIC_MANIFEST.tsv',check.names=FALSE,quote='',colClasses='character')
for (p in paths) {
 if (!file.exists(p)) stop('Missing bundled plotting table: ',p,'. See REPRODUCIBILITY.md; no fallback or download is available.')
 i <- match(p,manifest$path)
 if (is.na(i) || digest::digest(file=p,algo='sha256') != manifest$sha256[i]) stop('Plot input integrity mismatch: ',p)
}
d <- lapply(paths,function(p)read.delim(p,check.names=FALSE,comment.char='',quote='"',stringsAsFactors=FALSE))
accepted_render_validate(scope,d)
source('workflow/R/common/public_paths.R')
out <- resolve_output_directory(args[2], root)
if (!dir.exists(out) && !dir.create(out, recursive=TRUE, showWarnings=FALSE))
  stop('Output directory creation failed: ', out)
if (!dir.exists(out)) stop('Output directory was not created: ', out)
render <- function(stem,width,height,draw) {
 pdf <- file.path(out,paste0(stem,'.pdf')); png <- file.path(out,paste0(stem,'.png'))
 grDevices::cairo_pdf(pdf,width=width,height=height,pointsize=11)
 tryCatch(draw(),finally=grDevices::dev.off())
 grDevices::png(png,width=width*150,height=height*150,res=150,type='cairo')
 tryCatch(draw(),finally=grDevices::dev.off())
 stopifnot(all(file.info(c(pdf,png))$size > 1000))
}
if (scope=='tissue_context') for (s in d$support$source) {
 x <- d$pairs[d$pairs$source==s,]
 render(paste0(s,'_paired_change_scatter'),7.5,5.5,
   function()accepted_draw_tissue_context(x,d[[s]],d$support,s))
} else render('protein_centered_ratio_effects',7,4.5,function()accepted_draw_protein(d$members))
cat('Rendered accepted descriptive plots; no scientific recomputation.\n')
