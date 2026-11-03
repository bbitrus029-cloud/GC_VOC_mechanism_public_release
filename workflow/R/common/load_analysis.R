# Import definitions only; no scientific reconstruction is run by this loader.
load_public_analysis <- function(root=getwd(), envir=parent.frame()) {
  paths <- c(paste0('common/',c('contracts','statistics','artifacts','public_paths','public_context'),'.R'),
    paste0('methods/',c('supplementary_statistics','protein_ratios','pseudobulk_preparation','native_single_cell'),'.R'),
    paste0('io/',c('gse122401','gse206785','gse150290','pdc000614'),'.R'),
    paste0('preparation/',c('bulk','single_cell','protein','supplementary'),'.R'),
    paste0('analysis/',c('original_bulk','original_single_cell','original_protein','original_integration','external_bulk','epithelial','tissue_protein'),'.R'))
  for(p in paths) sys.source(file.path(root,'workflow/R',p),envir=envir)
  invisible(paths)
}
public_preparation_context <- function(ctx) {
  ctx$supplementary_config <- yaml::read_yaml(ctx$input('config/post_effect_supplement.yml'))
  ctx$route_specification <- yaml::read_yaml(ctx$input('workflow/config/implementation/gse150290_route.yml'))
  ctx$broad_audit_configuration <- jsonlite::read_json(ctx$input('workflow/config/implementation/gse150290_audit_panels.json'),simplifyVector=TRUE)
  ctx$seed_map <- read.delim(ctx$input('workflow/config/implementation/gse150290_full_seed_map.tsv'),stringsAsFactors=FALSE)
  ctx
}
run_public_reconstruction <- function(module,manifest,out,A_directory=NULL) {
  ctx <- public_analysis_context(module,manifest,out,create=FALSE)
  required <- read.delim(ctx$input('workflow/config/reconstruction_inputs.tsv'),quote='',colClasses='character')
  keys <- required$logical_input[required$module==module]
  if(!length(keys)) stop('Unknown reconstruction module: ',module)
  # Fail before creating outputs when any declared prerequisite is unavailable.
  invisible(ctx$input(keys))
  if(module=='tissue_protein') {
    if(is.null(A_directory)) stop('Missing independently reconstructed A processed directory.')
    A_directory <- normalizePath(A_directory,mustWork=TRUE)
    for(n in c('rank_scores.rds','diagnostic_scores.tsv')) invisible(ctx$input(file.path(A_directory,n)))
  }
  ctx$runtime(if(module=='epithelial')'epithelial' else if(startsWith(module,'original_'))'original' else 'supplementary')
  factory <- get(paste0(module,'_module'),mode='function')
  mod <- if(module=='tissue_protein')factory(ctx,A_directory) else factory(ctx)
  result <- mod$run()
  ctx$write_tsv(ctx$verify_inputs(),ctx$output('qc','verified_input_identities.tsv'))
  cat('Reconstruction completed: ',module,'; outputs in ',ctx$output_root,'\n',sep='')
  invisible(result)
}
