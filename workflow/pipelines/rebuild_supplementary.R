args <- commandArgs(TRUE)
if(!length(args) || !args[1] %in% c('A','B','CD') || length(args)!=if(args[1]=='CD')4L else 3L)
 stop('Usage: workflow/environment/run-r workflow/pipelines/rebuild_supplementary.R A|B|CD INPUT_MANIFEST OUTPUT_DIRECTORY [A_PROCESSED_DIRECTORY required for CD]')
if(!file.exists('PUBLIC_MANIFEST.tsv'))stop('Run from the public release root.')
source('workflow/R/common/load_analysis.R');load_public_analysis()
module <- c(A='external_bulk',B='epithelial',CD='tissue_protein')[[args[1]]]
run_public_reconstruction(module,args[2],args[3],if(length(args)==4L)args[4]else NULL)
