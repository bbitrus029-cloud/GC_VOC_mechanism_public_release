# Parse/load the released function library; execute only the included synthetic test.
if (!file.exists('PUBLIC_MANIFEST.tsv')) stop('Run from the release directory.')
files <- c('workflow/R/common/contracts.R','workflow/R/common/statistics.R',
 'workflow/R/methods/supplementary_statistics.R','workflow/R/methods/protein_ratios.R',
 'workflow/R/methods/pseudobulk_preparation.R','workflow/R/methods/native_single_cell.R',
 'workflow/R/io/gse122401.R','workflow/R/methods/accepted_plotting.R')
for (p in files) { parse(p); source(p) }
source('workflow/tests/test_methods.R')
cat('Library parse/import and synthetic native-anchor test passed; no study data fitted.
')
