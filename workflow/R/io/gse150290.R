# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.
# Factory below provides explicit immutable configuration and output context.
# Helpers are sourced into a new module environment once; no function is replaced.
gcvoc_gse150290 <- function(context) {
  stopifnot(is.function(context$input), is.character(context$root),
    is.data.frame(context$seed_map), is.list(context$route_specification),
    is.list(context$supplementary_config), is.list(context$broad_audit_configuration))
  trial_spec <- function() context$route_specification
  trial_atomic <- function(path,writer,check) {
    stopifnot(is.function(context$atomic), is.function(context$guard))
    context$guard(); context$atomic(path,writer,check)
  }

trial_archive_layout <- function (entries) 
{
    normalized <- sub("/$", "", entries)
    trial_require(length(entries) > 0L && !anyNA(entries) && !anyDuplicated(normalized) && !any(grepl("^/|(^|/)\\.\\.?(/|$)|\\\\", 
        normalized)), "Unsafe or duplicate archive paths")
    matrix <- normalized[grepl("(^|/)raw_gene_bc_matrices/hg19/matrix[.]mtx$", normalized)]
    trial_require(length(matrix) == 1L, "Unresolved physical capture/archive identities")
    root <- sub("/hg19/matrix[.]mtx$", "", matrix)
    expected <- c(root, paste0(root, "/hg19"), paste0(root, "/hg19/", c("matrix.mtx", "genes.tsv", "barcodes.tsv")))
    trial_require(all(expected[3:5] %in% normalized) && all(normalized %in% expected), "Unresolved physical capture/archive identities")
    paste0(root, "/hg19")
}

trial_raw_capture <- function (archive, cache) 
{
    archive <- context$input(archive)
    cache <- context$allowed(cache)
    entries <- utils::untar(archive, list = TRUE)
    relative_prefix <- trial_archive_layout(entries)
    kinds <- system2("/usr/bin/tar", c("--list", "--verbose", "--gzip", "--file", shQuote(archive)), 
        stdout = TRUE)
    trial_require(is.null(attr(kinds, "status")) && length(kinds) == length(entries) && all(substr(kinds, 
        1, 1) %in% c("-", "d")), "Unsupported archive entry type")
    trial_require(!dir.exists(cache), "Unregistered existing capture cache")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    utils::untar(archive, exdir = cache)
    prefix <- file.path(cache, relative_prefix)
    features <- read.delim(file.path(prefix, "genes.tsv"), header = FALSE, colClasses = "character", 
        check.names = FALSE)
    trial_require(ncol(features) == 2L, "Source feature schema")
    names(features) <- c("feature_id", "source_symbol")
    barcodes <- readLines(file.path(prefix, "barcodes.tsv"))
    triplets <- Matrix::readMM(file.path(prefix, "matrix.mtx"))
    trial_require(all(is.finite(triplets@x)) && all(triplets@x >= 0) && all(triplets@x == round(triplets@x)), 
        "Noninteger original coordinate counts")
    trial_require(!anyDuplicated(as.double(triplets@i) + nrow(triplets) * as.double(triplets@j)), "Duplicate original matrix coordinates")
    trial_require(nrow(triplets) == nrow(features) && ncol(triplets) == length(barcodes) && !anyDuplicated(barcodes) && 
        all(nzchar(barcodes)) && !anyDuplicated(features$feature_id) && !anyNA(features) && all(nzchar(features$feature_id)) && 
        all(nzchar(features$source_symbol)), "Raw dimension or identity conflict")
    counts <- as(triplets, "CsparseMatrix")
    dimnames(counts) <- list(features$feature_id, barcodes)
    list(counts = counts, features = features, barcodes = barcodes, unpacked_files = trial_manifest(file.path(prefix, 
        c("matrix.mtx", "genes.tsv", "barcodes.tsv"))))
}

trial_complete_map <- function (metrics, identity, qc = NULL, ed = NULL, dbl = NULL) 
{
    n <- nrow(metrics)
    cellkey <- paste(paste0("GSE150290_", identity$gsm), identity$gsm, metrics$barcode, sep = "|")
    trial_require(!anyDuplicated(cellkey), "Duplicate full cell key")
    map <- data.frame(cell_key = cellkey, source = "GSE150290", sample_id = paste0("GSE150290_", identity$gsm), 
        capture_id = identity$gsm, donor = identity$donor, tissue_origin = identity$tissue, source_sample_label = identity$source_name, 
        original_source_label = "ORIGINAL_LABEL_NOT_PROVIDED", metrics, author_low_quality = NA, author_doublet = NA, 
        emptyDrops_PValue = NA_real_, emptyDrops_FDR = NA_real_, emptyDrops_Limited = NA, scDblFinder_score = NA_real_, 
        scDblFinder_class = "NOT_RUN", technical_state = "NOT_EVALUATED", provisional_cluster = NA_character_, 
        provisional_broad_label = "UNRESOLVED", provisional_lineage = "UNRESOLVED", epithelial_membership = "UNRESOLVED", 
        lineage = "UNKNOWN", author_malignancy_status = "UNKNOWN", author_malignancy_evidence = "UNKNOWN", 
        proliferation_state = "UNKNOWN", first_exclusion_reason = NA_character_, all_applicable_reasons = NA_character_)
    if (!is.null(ed)) {
        trial_require(nrow(ed$result) == n, "Complete emptyDrops map")
        map$emptyDrops_PValue <- ed$result$PValue
        map$emptyDrops_FDR <- ed$result$FDR
        map$emptyDrops_Limited <- ed$result$Limited
    }
    if (is.null(qc)) 
        return(map)
    map$author_low_quality <- qc$author_low
    map$author_doublet <- qc$author_doublet
    flags <- qc$flags
    map$technical_state <- ifelse(qc$keep, "PENDING_DOUBLET", "EXCLUDED_QC")
    if (!is.null(dbl)) {
        j <- which(qc$keep)
        trial_require(identical(rownames(dbl$result), map$barcode[j]), "Doublet mapping/order")
        map$scDblFinder_score[j] <- dbl$result$scDblFinder.score
        map$scDblFinder_class[j] <- as.character(dbl$result$scDblFinder.class)
        map$technical_state[j] <- ifelse(map$scDblFinder_class[j] == "singlet", "QC_SINGLET_PENDING_ANNOTATION", 
            "EXCLUDED_DOUBLET")
    }
    flags$COMPUTATIONAL_DOUBLET <- map$scDblFinder_class == "doublet"
    flags$UNRESOLVED_DOUBLET <- qc$keep & !map$scDblFinder_class %in% c("singlet", "doublet")
    precedence <- trial_spec()$cell_retention$reason_precedence
    flags <- flags[, intersect(precedence, names(flags)), drop = FALSE]
    map$first_exclusion_reason <- apply(flags, 1, function(z) {
        i <- which(z)
        if (length(i)) 
            names(flags)[i[1]]
        else NA_character_
    })
    map$all_applicable_reasons <- apply(flags, 1, function(z) paste(names(flags)[which(z)], collapse = ";"))
    cbind(map, setNames(flags, paste0("flag_", names(flags))))
}

trial_manifest <- function (paths) 
{
    paths <- sort(unique(paths), method = "radix")
    trial_require(all(file.exists(paths)), "Missing material")
    data.frame(path = paths, bytes = as.numeric(file.info(paths)$size), sha256 = vapply(paths, trial_hash, 
        ""), row.names = NULL)
}

trial_hash <- function (p) 
digest::digest(file = p, algo = "sha256", serialize = FALSE)

trial_require <- function (x, message) 
if (!isTRUE(x)) stop(message, call. = FALSE)

trial_precise_table <- function (x, path, gzip = FALSE) 
{
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(x) <- NULL
    for (k in names(x)) if (is.factor(x[[k]])) 
        x[[k]] <- as.character(x[[k]])
    classes <- vapply(x, function(z) if (is.integer(z)) 
        "integer"
    else if (is.double(z)) 
        "numeric"
    else if (is.logical(z)) 
        "logical"
    else "character", "")
    encoded <- x
    for (k in names(x)[classes == "numeric"]) {
        encoded[[k]] <- format(x[[k]], digits = 17, scientific = TRUE, trim = TRUE)
        encoded[[k]][is.na(x[[k]])] <- NA_character_
    }
    trial_atomic(path, function(p) {
        con <- if (gzip) 
            gzfile(p, "wt")
        else file(p, "wt")
        on.exit(close(con))
        write.table(encoded, con, sep = "\t", row.names = FALSE, quote = which(classes == "character"), 
            qmethod = "double", na = "NA")
    }, function(p) {
        con <- if (gzip) 
            gzfile(p, "rt")
        else file(p, "rt")
        on.exit(close(con))
        read <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, colClasses = classes)
        trial_require(identical(read, x), "Substantive table readback mismatch")
    })
}

trial_map_write <- function (map, path) 
{
    trial_require(!anyDuplicated(map$cell_key) && !anyNA(map$cell_key), "Exhaustive map key")
    trial_precise_table(map, path, gzip = TRUE)
}

trial_json <- function (x, p) 
{
    txt <- jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
    trial_atomic(p, function(f) writeLines(txt, f), function(f) trial_require(identical(jsonlite::fromJSON(f, 
        simplifyVector = FALSE), jsonlite::fromJSON(txt, simplifyVector = FALSE)), "JSON readback mismatch"))
}

trial_rds <- function (x, p) 
trial_atomic(p, function(f) saveRDS(x, f, version = 3, compress = FALSE), function(f) trial_require(identical(readRDS(f), 
    x), "RDS readback mismatch"))

  ann_need <- trial_require; mapping_need <- trial_require; s10b_require <- trial_require
  s10b_prepare_core <- s03_prepare_compartment
  for (name in c('gse150290_qc.R','gse150290_annotation.R','gse150290_population.R'))
    sys.source(file.path(context$root,'workflow/R/methods',name),envir=environment())
  exports <- ls(envir=environment(),pattern='^(trial_|ann_|audit_|full_|s10b_|prep_|mapping_)')
  mget(exports,envir=environment(),inherits=FALSE)
}
