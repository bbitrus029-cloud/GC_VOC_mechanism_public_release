#!/usr/bin/env Rscript
# Definitions, source readers and preparation tests only. No study effect fitting.
args <- commandArgs(TRUE)
if(length(args)!=1L) stop('Usage: check_reconstruction.R OUTSIDE_EMPTY_DIRECTORY')
source('workflow/R/common/load_analysis.R');load_public_analysis()
out <- resolve_output_directory(args[1L],getwd())
if(!dir.exists(out) && !dir.create(out,recursive=TRUE)) stop('Output path creation failed.')
manifest <- file.path(out,'synthetic_inputs.tsv')
writeLines('logical_input	local_path	sha256',manifest)
ctx <- public_preparation_context(public_analysis_context('synthetic_checks',manifest,file.path(out,'fixtures')))
g150 <- gcvoc_gse150290(ctx); pdc <- gcvoc_pdc000614(ctx)
# Resource state must never enter the protected input-hash registry.
ctx$resource_check(minimum=0);ctx$guard();stopifnot(nrow(ctx$verify_inputs())>0)
source('workflow/tests/test_readers_preparation.R')
for(name in c('original_bulk','original_single_cell','original_protein','original_integration','external_bulk','epithelial')) {
  mod <- get(paste0(name,'_module'))(ctx)
  stopifnot(is.list(mod),is.function(mod$run))
}
stopifnot(is.function(tissue_protein_module(ctx,out)$run))
stopifnot(is.function(bulk_preparation(ctx)$s02_array),is.function(protein_preparation(ctx)$prepare))
cat('PASS: seven module definitions, preparation factories and synthetic reader checks. No study analyses run.\n')

source('workflow/tests/test_input_bindings.R')
test_input_bindings(file.path(out,'input_binding_tests'))
