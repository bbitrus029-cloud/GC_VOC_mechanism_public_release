# Synthetic reader/preparation fixtures; ctx, g150 and pdc supplied by check_reconstruction.R.
testthat::test_that("continuous abundance parser excludes the numeric pseudoheader", {
  p <- ctx$allowed("pseudoheader_fixture.tsv.gz")
  con <- gzfile(p,"wt");writeLines(c("gene_id\t999.25","ENSG000001\t0.125","ENSG000002\t3.75"),con);close(con)
  schema <- gse122401_inspect_processed_schemas(p,expected_samples=1L,
    uq_normalization_documented=TRUE,rsem_documented=TRUE,
    documented_feature_count=2L,expected_feature_count=2L)
  actual <- gse122401_read_processed_gene_file(p,schema)
  testthat::expect_identical(actual$gene_id,c("ENSG000001","ENSG000002"))
  testthat::expect_identical(actual$uq_normalized_expression,c(.125,3.75))
  testthat::expect_identical(schema$skip_rows,1L)
  testthat::expect_identical(schema$processed_scale,"continuous_RSEM_upper_quartile_normalized_abundance_not_raw_counts")
  testthat::expect_false(gse122401_numeric_token("1junk"))
  con <- gzfile(p,"wt");writeLines(c("gene_id\t999.25","ENSG000001\t1junk","ENSG000002\t3.75"),con);close(con)
  testthat::expect_error(gse122401_read_processed_gene_file(p,schema),"non-numeric")
  unlink(p)
})

testthat::test_that("GSE206785 frozen universe and raw source identities remain exact", {
  expected <- gse206785_fixed_commit_contract()
  testthat::expect_identical(expected$commit,"dd0e37b1ab5dcf41caa2bce5eea08aa9446446da")
  testthat::expect_identical(expected$raw_fragment_names,sprintf("raw_counts.csv.gz%02d",0:8))
  testthat::expect_error(gse206785_validate_formal_matrix_source(expected$commit,"GSE206785_scgex.txt.gz"),"not used")
  manifest <- expected$files[match(expected$raw_fragment_names,expected$files$file_name),]
  testthat::expect_true(gse206785_validate_raw_fragment_manifest(manifest))
  bad <- manifest;bad$sha256[1] <- paste(rep("0",64),collapse="")
  testthat::expect_error(gse206785_validate_raw_fragment_manifest(bad),"SHA-256 mismatch")
  testthat::expect_error(gse206785_parse_count_gzip("synthetic"),"not formal raw UMI")
  testthat::expect_identical(gse206785_canonical_hash_contract()$expected_digest,
    "1500b3bf1e4d8ec87b1a8dc892ce6a119cfb304e4f80e521435b8b3a00454349")
})

testthat::test_that("protein readers keep source identity and complete paired units", {
  testthat::expect_identical(unname(vapply(c("C3L.01234.T","C3L-01234-N"),
    pdc$pdc_gate_canonical_sample_id,character(1))),
    c("C3L.01234.T","C3L.01234.N"))
  testthat::expect_error(pdc$pdc_gate_canonical_sample_id("unresolved"),"sample")
  x <- data.frame(canonical_case_identifier=rep(c("C3L.00001","C3L.00002","C3L.00003"),each=2),
    canonical_tissue_label=rep(c("Tumor","NAT"),3),
    original_matrix_column=paste0("sample",1:6),opaque_identifier=paste0("opaque",1:6),
    authoritative_aliquot_identifier=paste0("aliquot",1:6),excluded=FALSE,stringsAsFactors=FALSE)
  x$authoritative_sample_identifier <- paste(x$canonical_case_identifier,rep(c("T","N"),3),sep=".")
  # The identity validator is shared by pairing and masked readers.
  testthat::expect_true(pdc$pdc_gate_validate_identity_table(x))
  bad <- x;bad$authoritative_sample_identifier[2] <- bad$authoritative_sample_identifier[1]
  testthat::expect_error(pdc$pdc_gate_validate_identity_table(bad),"duplicate|unique|identity")
})

testthat::test_that("raw archive paths and effective feature identity fail closed", {
  entries <- c("raw_gene_bc_matrices","raw_gene_bc_matrices/hg19",
    paste0("raw_gene_bc_matrices/hg19/",c("matrix.mtx","genes.tsv","barcodes.tsv")))
  testthat::expect_identical(g150$trial_archive_layout(entries),"raw_gene_bc_matrices/hg19")
  testthat::expect_error(g150$trial_archive_layout(c(entries,"../outside")),"Unsafe")
  testthat::expect_error(g150$trial_archive_layout(c(entries,entries[3])),"duplicate")
  cw <- read.delim(ctx$input("workflow/config/implementation/gse150290_jchain_crosswalk.tsv"),stringsAsFactors=FALSE)
  f <- data.frame(feature_id=c("ENSG00000132465","ENSG_fixture"),source_symbol=c("IGJ","BG"),source_row=c(7669L,1L))
  out <- g150$ann_approved_features(f,cw)
  testthat::expect_identical(out[,names(f)],f)
  testthat::expect_identical(out$canonical_symbol,c("JCHAIN","BG"))
  testthat::expect_identical(out$historical_gene_version,c(6L,NA_integer_))
  f$source_symbol[2] <- "JCHAIN"
  testthat::expect_error(g150$ann_approved_features(f,cw),"colliding")
})

testthat::test_that("formal QC preserves original symbols and exact boundary rules", {
  f <- data.frame(feature_id=paste0("id",1:6),source_symbol=c("MT-CO1","HBB","HBB","ALDH2","GENE1","GENE2"),
    canonical_symbol=c("MT-CO1","HBB","HBB","ALDH2","GENE1","GENE2"),source_row=1:6)
  x <- Matrix::Matrix(matrix(c(10,20,30,40,50,50,0,0,0,0,0,0,0,0,0,100,100,100),6,
    dimnames=list(f$feature_id,c("a","b","c"))),sparse=TRUE)
  q <- g150$trial_metrics(x,f)
  testthat::expect_identical(q$original_UMI,c(200,0,300))
  testthat::expect_equal(q$combined_MT_Hb_fraction,c(.3,NA,0),tolerance=0)
  testthat::expect_identical(g150$trial_low_quality(q,c(TRUE,FALSE,TRUE))$keep,c(FALSE,FALSE,TRUE))
  x@x[1] <- .5
  testthat::expect_error(g150$trial_metrics(x,f),"integer UMI")
  q <- data.frame(original_UMI=c(100,200,200,200,20000,20001),detected_features=rep(30,6),combined_MT_Hb_fraction=c(0,0,.1,0,0,0))
  z <- g150$trial_low_quality(q,rep(TRUE,6))
  testthat::expect_false(z$flags$UMI_OUTSIDE_RANGE[1])
  testthat::expect_false(z$flags$UMI_OUTSIDE_RANGE[5])
  testthat::expect_true(z$flags$UMI_OUTSIDE_RANGE[6])
  testthat::expect_true(z$flags$MT_HB_AT_LEAST_10_PERCENT[3])
})

testthat::test_that("preparation preserves original rows and unmatched detection units", {
  cfg <- ctx$supplementary_config$B
  h <- g150$s10b_h4_members()
  units <- data.frame(unit_id=paste0("fixture_",1:17),donor=c(rep(sprintf("d%02d",1:8),each=2),"unmatched"),
    tissue=c(rep(c("Normal","Tumor"),8),"Normal"),cell_count=20L)
  genes <- c(h,sprintf("bg%04d",1:596))
  x <- matrix(30,600,17,dimnames=list(genes,units$unit_id))
  x["bg0009",] <- 0;x["bg0009",1:8] <- 1;x["bg0010",1:16] <- 0
  p <- g150$s10b_prepare(x,units,h,cfg)
  p16 <- g150$s10b_prepare(x[,1:16],units[1:16,],h,cfg)
  testthat::expect_false("bg0009" %in% p$background)
  testthat::expect_true("bg0009" %in% p16$background)
  testthat::expect_equal(p$detection$nonzero_fraction[p$detection$gene_id=="bg0009"],8/17,tolerance=0)
  testthat::expect_identical(p$original_library_sizes,colSums(x))
  testthat::expect_identical(p$modeling_library_sizes,colSums(p$counts))
  testthat::expect_identical(ncol(p$counts),16L)
  testthat::expect_identical(p$eligible_units,units$unit_id)
  testthat::expect_false("unmatched" %in% p$complete_pair_donors)
  bad <- units;bad$unit_id[2] <- bad$unit_id[1]
  testthat::expect_error(g150$s10b_prepare(x,bad,h,cfg),"unit|duplicat|order")
  ids <- sprintf("ENSG_fixture_%04d",1:600);rownames(x) <- ids
  features <- data.frame(feature_id=ids,source_symbol=genes,source_row=1:600,canonical_symbol=genes)
  mapped <- g150$prep_endpoint(x,units,features,cfg)
  testthat::expect_identical(mapped$member_qualification$feature_id,ids[1:4])
  features$canonical_symbol[4] <- "missing"
  partial <- g150$prep_endpoint(x,units,features,cfg)
  testthat::expect_identical(sum(partial$member_qualification$technical_eligible),3L)
  testthat::expect_false(partial$exact_H4_technical_eligible)
  testthat::expect_identical(partial$counts,mapped$counts)
})

testthat::test_that("canonical panel scoring retains ties and empty population semantics", {
  f <- data.frame(feature_id=c("id1","id2","id3"),source_symbol=c("OLD_A","B","C"),canonical_symbol=c("A","B","C"),source_row=1:3)
  x <- Matrix::Matrix(matrix(c(1,4,0,2,6,0),3,dimnames=list(f$feature_id,c("d1","d2"))),sparse=TRUE)
  z <- g150$trial_panel_scores(x,f,rep(TRUE,3),list(high=c("A","B"),low=c("A","C")))
  testthat::expect_equal(unname(z$scores),matrix(c(.75,.75,.25,.25),2),tolerance=0)
  testthat::expect_identical(g150$trial_panel_scores(x,f,rep(TRUE,3),list(missing=c("A","D")))$state,"NA_NONESTIMABLE")
  ties <- matrix(.5,8,2,dimnames=list(NULL,c("epithelial","other")))
  testthat::expect_identical(g150$trial_adjudicate(ties,minimum_donors=8L)$state,"FAIL")
  spec <- g150$audit_spec()
  genes <- unique(c(unlist(spec$annotation_qualification$assignment_panels),
    unlist(spec$annotation_qualification$broad_audit_panels),
    unlist(spec$annotation_qualification$lineage_assignment_panels),
    unlist(spec$annotation_qualification$lineage_audit_panels),paste0("BG",1:5)))
  f <- data.frame(feature_id=paste0("empty",seq_along(genes)),source_symbol=genes,canonical_symbol=genes,source_row=seq_along(genes))
  x <- Matrix::Matrix(matrix(numeric(),length(genes),0,dimnames=list(f$feature_id,character())),sparse=TRUE)
  z <- g150$ann_annotate(x,f,character(),character(),spec=spec)
  testthat::expect_identical(nrow(z$cells),0L)
  testthat::expect_identical(z$cells$author_malignancy_status,character())
  testthat::expect_identical(z$cells$proliferation_state,character())
})

testthat::test_that("full graph operation identities and current module loading are explicit", {
  text <- paste(deparse(g150$full_pooled_graph),collapse="\n")
  testthat::expect_match(text,"post_effect_B.GSE150290.full_v1.pooled_pca",fixed=TRUE)
  testthat::expect_match(text,"post_effect_B.GSE150290.full_v1.pooled_graph",fixed=TRUE)
  testthat::expect_false(grepl("trial_v1",text,fixed=TRUE))
  files <- c("workflow/R/io/gse150290.R",paste0("workflow/R/methods/gse150290_",c("qc","annotation","population"),".R"))
  current <- paste(unlist(lapply(files,readLines)),collapse="\n")
  testthat::expect_false(grepl("body\\([^)]*\\)\\s*<-|parse\\(text|environment\\([^)]*\\)\\s*<-",current))
  testthat::expect_false(grepl("src/R/",current,fixed=TRUE))
})

testthat::test_that("typed table output retains quoted strings, missing values and order", {
  x <- data.frame(identity=c("b","a","c","d"),value=c(pi,1e-22,0,NA_real_),
    flag=c(TRUE,FALSE,NA,TRUE),state=c("PASS","FAIL","UNKNOWN","NONESTIMABLE"),
    message=c('quoted "text"',"multiline\nmessage","tab\tmessage","literal \\ path"))
  p <- ctx$allowed("quoted_reader_fixture.tsv.gz")
  g150$trial_precise_table(x,p,gzip=TRUE)
  con <- gzfile(p,"rt")
  observed <- read.delim(con,colClasses=c("character","numeric","logical","character","character"))
  close(con)
  testthat::expect_identical(observed,x)
  testthat::expect_true(g150$ann_compare(observed,x))
  changed <- observed;changed$value[1] <- changed$value[1]+1
  testthat::expect_false(g150$ann_compare(changed,x))
  missing <- observed;missing$message[1] <- NA_character_
  testthat::expect_false(g150$ann_compare(missing,x))
  frozen_seed <- read.delim(ctx$input("workflow/config/implementation/gse150290_full_seed_map.tsv"),stringsAsFactors=FALSE)
  testthat::expect_identical(ctx$seed_map,frozen_seed)
  testthat::expect_identical(names(frozen_seed),c("derivation_version","master_seed","operation_name","seed","rng_kind","hash_prefix"))
  operation <- "post_effect_B.GSE150290.full_v1.pooled_pca"
  expected_seed <- frozen_seed$seed[match(operation,frozen_seed$operation_name)]
  rng <- g150$trial_seed(operation)
  testthat::expect_identical(rng$seed,expected_seed)
  testthat::expect_identical(RNGkind()[1],"L'Ecuyer-CMRG")
  testthat::expect_identical(BiocParallel::bpnworkers(rng$bp),1L)
})
