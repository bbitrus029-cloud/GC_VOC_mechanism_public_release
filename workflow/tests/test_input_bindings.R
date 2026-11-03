# Synthetic fixtures only. Run from check_reconstruction.R; all writes stay outside the release.
test_input_bindings <- function(fixtures) {
  stopifnot(!file.exists(fixtures));dir.create(fixtures,recursive=TRUE)
  hash <- function(p) digest::digest(file=p,algo='sha256')
  put_manifest <- function(name,keys=character(),paths=character()) {
    p <- file.path(fixtures,paste0(name,'.tsv'))
    write.table(data.frame(logical_input=keys,local_path=paths,sha256=vapply(paths,hash,'')),p,sep='\t',quote=FALSE,row.names=FALSE)
    p
  }
  empty <- put_manifest('empty')
  context <- function(name,m=empty) public_analysis_context(name,m,file.path(fixtures,name))
  bad <- file.path(fixtures,'synthetic_config.yml')
  txt <- readLines('workflow/config/original_analysis.yml')
  stopifnot(any(grepl('ACSL4',txt,fixed=TRUE)))
  writeLines(sub('ACSL4','SYNTHETIC_WRONG_MEMBER',txt,fixed=TRUE),bad)
  testthat::test_that('bundled controls and governing catalogs cannot be shadowed', {
    keys <- c('config/analysis_config.yml','./config/analysis_config.yml',
      'workflow/config/original_analysis.yml','./workflow//config/original_analysis.yml',
      file.path(getwd(),'workflow/config/original_analysis.yml'),
      'workflow/config/bundled_input_aliases.tsv','workflow/config/reconstruction_inputs.tsv',
      'workflow/config/external_input_identities.tsv','config/implementation/s01_target_seed_map.tsv')
    for(i in seq_along(keys)) {
      m <- put_manifest(paste0('shadow',i),keys[i],bad)
      testthat::expect_error(context(paste0('shadow_run',i),m),'cannot be shadowed')
      testthat::expect_false(dir.exists(file.path(fixtures,paste0('shadow_run',i))))
    }
    good <- context('bundled_positive')
    testthat::expect_identical(good$input('./config/analysis_config.yml'),normalizePath('workflow/config/original_analysis.yml'))
    testthat::expect_error(good$input('workflow/config/../config/original_analysis.yml'),'traversal')
  })
  testthat::test_that('a self-consistent artificial PDC report cannot enter C/D', {
    ids <- read.delim('provenance/eligibility/supplementary/supplementary_readiness/bindings/D_ratio_column_identity.tsv',check.names=FALSE)
    stopifnot(nrow(ids)==130L)
    fake <- data.frame(Gene=c('ALDH2','ALDH3A1','ALDH3A2','ALDH7A1'),check.names=FALSE)
    for(n in ids$ratio_column)fake[[n]]<-rep(1,4) # Artificial values; no study measurements.
    fake$NCBIGeneID <- 1:4;fake$Authority <- 'SYNTHETIC';fake$Organism <- 'SYNTHETIC'
    p <- file.path(fixtures,'artificial_report.tsv');write.table(fake,p,sep='\t',quote=FALSE,row.names=FALSE)
    key <- 'data/raw/post_effect_supplement/protein/Gastric_Cancer_Korea_Proteome.itraq.tsv'
    ctx <- context('bad_acquisition',put_manifest('bad_acquisition_manifest',key,p))
    testthat::expect_error(ctx$input(key),'Authoritative acquisition SHA-256 mismatch')
    mod <- tissue_protein_module(ctx,fixtures)
    actual_protein_input <- get('s13_protein_input',envir=environment(mod$run))
    testthat::expect_error(actual_protein_input(),'Authoritative acquisition SHA-256 mismatch')
    other <- read.delim('workflow/config/external_input_identities.tsv',stringsAsFactors=FALSE)
    other <- other[other$binding_rule=='EXACT_ACQUISITION' & other$object!='Gastric_Cancer_Korea_Proteome.itraq.tsv',]
    for(i in seq_len(nrow(other))) {
      k <- paste0('different/logical/prefix/',other$object[i])
      c2 <- context(paste0('other_acquisition_',i),put_manifest(paste0('other_manifest_',i),k,p))
      testthat::expect_error(c2$input(k),'Authoritative acquisition SHA-256 mismatch')
    }
  })
  testthat::test_that('declared reconstructed inputs and owned atomic outputs remain usable', {
    external <- file.path(fixtures,'rank_scores.rds');saveRDS(list(synthetic=TRUE,value=pi),external)
    key <- 'data/processed/post_effect_supplement/s10/GSE122401_adapter.rds'
    ctx <- context('owned',put_manifest('owned_manifest',c(key,external),c(external,external)))
    testthat::expect_identical(ctx$read_rds(key),list(synthetic=TRUE,value=pi))
    testthat::expect_identical(ctx$input(external),external) # Explicit external A-style dependency.
    sentinel <- file.path(fixtures,'unlisted.rds');saveRDS(list(sentinel='preserve'),sentinel);before<-hash(sentinel)
    testthat::expect_error(ctx$input(sentinel),'Missing declared input')
    testthat::expect_error(ctx$input(file.path(ctx$output_root,'..','unlisted.rds')),'traversal')
    link <- file.path(ctx$output_root,'escape.rds');stopifnot(file.symlink(sentinel,link))
    testthat::expect_error(ctx$input(link),'symbolic link')
    unregistered <- file.path(ctx$output_root,'unregistered.rds');saveRDS(list(synthetic=TRUE),unregistered)
    testthat::expect_error(ctx$read_rds(unregistered),'Unregistered')
    reserved <- ctx$allowed('reserved.rds');saveRDS(list(synthetic=TRUE),reserved)
    testthat::expect_error(ctx$read_rds(reserved),'Unregistered')
    p <- ctx$write_rds(list(synthetic=TRUE,value=2),ctx$output('processed','owned.rds'))
    testthat::expect_identical(ctx$read_rds(p),list(synthetic=TRUE,value=2))
    saveRDS(list(synthetic=TRUE,value=3),p)
    testthat::expect_error(ctx$read_rds(p),'Run-produced input integrity mismatch')
    ctx$write_rds(list(synthetic=TRUE,value=4),p)
    testthat::expect_identical(ctx$read_rds(p),list(synthetic=TRUE,value=4))
    for(name in c('original_bulk','original_single_cell','original_protein','original_integration')) {
      mod <- get(paste0(name,'_module'))(ctx)
      reader <- get('read_rds',envir=environment(mod$run))
      testthat::expect_error(reader(file.path(ctx$output_root,'..','unlisted.rds')),'traversal')
    }
    testthat::expect_identical(hash(sentinel),before)
    testthat::expect_identical(readRDS(sentinel),list(sentinel='preserve'))
  })
  testthat::test_that('protein preparation resolves its logical prefix before each low-level reader', {
    physical <- file.path(fixtures,'physical');dir.create(physical)
    meta <- file.path(physical,'meta.tsv');clinical <- file.path(physical,'clinical.tsv')
    write.table(data.frame(Sample_ID=c('C3L.00001.N','C3L.00001.T'),Aliquot_ID=c('CPT0000000001','CPT0000000002'),Case='C3L.00001',tissue_type=c('NAT','Tumor')),meta,sep='\t',quote=FALSE,row.names=FALSE)
    write.table(data.frame(Sample='synthetic',Case='C3L.00001'),clinical,sep='\t',quote=FALSE,row.names=FALSE)
    matrix <- data.frame(Gene=c('SYNTH_A','SYNTH_B'),check.names=FALSE)
    matrix[['C3L.00001.N']]<-c(1,2);matrix[['C3L.00001.T']]<-c(2,3)
    rna<-file.path(physical,'rna.tsv');protein<-file.path(physical,'protein.tsv')
    write.table(matrix,rna,sep='\t',quote=FALSE,row.names=FALSE);write.table(matrix,protein,sep='\t',quote=FALSE,row.names=FALSE)
    prefix <- 'synthetic/logical/source'
    names <- c('GC_metatable.tsi','GC_Jan2025_clinical_data.tsi','synthetic_RNA.cct','synthetic_protein.cct')
    paths <- c(meta,clinical,rna,protein);keys<-file.path(prefix,names)
    m <- put_manifest('protein_prep',keys,paths);ctx<-context('protein_prep_run',m)
    calls<-character();original<-ctx$input
    ctx$input<-function(p){q<-original(p);calls<<-c(calls,q);q}
    settings<-list(source_directory=prefix,rna_file=names[3],protein_file=names[4],rna_sha256=hash(rna),protein_sha256=hash(protein))
    # Actual readers accept all synthetic identities/matrices, then the fixed study-size guard rejects them. No model is fitted.
    testthat::expect_error(protein_preparation(ctx)$prepare(settings),'physical dimensions')
    testthat::expect_identical(calls,paths)
    ctx2<-context('protein_undeclared',put_manifest('protein_partial',keys[-2],paths[-2]))
    testthat::expect_error(protein_preparation(ctx2)$prepare(settings),'Missing declared input')
    ctx3<-context('protein_tampered',m);cat('\n',file=meta,append=TRUE)
    testthat::expect_error(protein_preparation(ctx3)$prepare(settings),'integrity mismatch|Input changed')
  })
  cat('PASS: input authorities, protein preparation paths and run-output dependencies. Synthetic fixtures only.\n')
}
