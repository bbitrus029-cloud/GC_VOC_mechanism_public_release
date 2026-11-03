# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.

trial_sets_literal <- function (features, spec = trial_spec()) 
{
    a <- spec$annotation_qualification
    cfg <- context$supplementary_config
    protected <- unique(unlist(cfg$programs, use.names = FALSE))
    audit <- unique(unlist(c(a$broad_audit_panels, a$lineage_audit_panels), use.names = FALSE))
    symbol <- features$source_symbol
    trial_require(length(symbol) > 0L && !anyNA(symbol) && all(nzchar(symbol)) && !anyDuplicated(features$feature_id), 
        "Unresolved feature identity")
    ambiguous <- duplicated(symbol) | duplicated(symbol, fromLast = TRUE)
    mt <- startsWith(symbol, "MT-")
    hb <- symbol %in% spec$cell_retention$hemoglobin_symbols
    excluded <- symbol %in% c(protected, audit) | mt | hb
    list(protected = protected, audit = audit, ambiguous = ambiguous, mt = mt, hb = hb, doublet = !excluded, 
        annotation = !excluded & !ambiguous & !symbol %in% a$proliferation_markers_excluded_from_DR, 
        assignment = !excluded & !ambiguous & !symbol %in% a$proliferation_markers_excluded_from_DR, 
        audit_universe = !ambiguous & !symbol %in% protected, D7 = unlist(cfg$D7, use.names = FALSE))
}

trial_sets <- function (features, spec = trial_spec()) 
{
    s <- trial_sets_literal(features, spec)
    m <- ann_masks(features, spec = spec)
    s$ambiguous <- m$ambiguous
    s$mt <- m$mitochondrial
    s$hb <- m$hemoglobin
    s$doublet <- m$doublet_input
    s$annotation <- m$annotation_input
    s$assignment <- m$assignment_background
    s$audit_universe <- m$audit_background
    s
}

trial_metrics <- function (x, features, spec = trial_spec()) 
{
    trial_require(inherits(x, "sparseMatrix") && nrow(x) == nrow(features) && all(is.finite(x@x)) && 
        all(x@x >= 0) & all(x@x == round(x@x)), "Invalid original integer UMI")
    s <- trial_sets(features, spec)
    total <- Matrix::colSums(x)
    fraction <- function(mask) {
        v <- Matrix::colSums(x[mask, , drop = FALSE])/total
        v[total == 0] <- NA_real_
        v
    }
    data.frame(barcode = colnames(x), original_UMI = unname(total), detected_features = unname(Matrix::colSums(x > 
        0)), mitochondrial_fraction = unname(fraction(s$mt)), hemoglobin_fraction = unname(fraction(s$hb)), 
        combined_MT_Hb_fraction = unname(fraction(s$mt | s$hb)), row.names = NULL)
}

trial_emptydrops_check <- function (d, meta, totals, barcodes, spec = trial_spec()) 
{
    trial_require(nrow(d) == length(totals) && identical(rownames(d), barcodes) && all(c("Total", "LogProb", 
        "PValue", "Limited", "FDR") %in% names(d)), "emptyDrops barcode/schema integrity")
    trial_require(identical(as.numeric(d$Total), as.numeric(totals)), "emptyDrops original totals")
    lower <- spec$cell_retention$cell_calling$arguments$lower
    trial_require(length(meta$retain) == 1L && is.finite(meta$retain), "Missing/nonfinite emptyDrops automatic knee")
    above <- totals > lower
    trial_require(all(is.finite(d$PValue[above])) && all(is.finite(d$FDR[above])) && all(is.finite(d$LogProb[above])) && 
        !anyNA(d$Limited[above]), "Unexpected nonfinite emptyDrops output above lower")
    trial_require(all(d$PValue[above] >= 0 & d$PValue[above] <= 1) && all(d$FDR[above] >= 0 & d$FDR[above] <= 
        1), "Invalid technical probability")
    trial_require(!any(d$Limited[above] & d$FDR[above] > 0.001), "UNRESOLVED_MONTE_CARLO_RESOLUTION")
    above & is.finite(d$FDR) & d$FDR <= 0.001
}

trial_emptydrops <- function (x, operation, spec = trial_spec(), diagnostic_path = NULL) 
{
    rng <- trial_seed(operation)
    a <- spec$cell_retention$cell_calling$arguments
    a$niters <- as.integer(a$niters)
    a$alpha <- Inf
    result <- do.call(DropletUtils::emptyDrops, c(list(m = x), a, list(BPPARAM = rng$bp)))
    d <- as.data.frame(result)
    meta <- S4Vectors::metadata(result)
    output <- list(result = d, metadata = meta, arguments = a, rng = list(seed = rng$seed, before = rng$state, 
        after = .Random.seed, backend = "SerialParam", workers = BiocParallel::bpnworkers(rng$bp), RNGseed = BiocParallel::bpRNGseed(rng$bp)))
    if (!is.null(diagnostic_path)) 
        trial_rds(output, diagnostic_path)
    output$called <- trial_emptydrops_check(d, meta, Matrix::colSums(x), colnames(x), spec)
    output
}

trial_low_quality <- function (metrics, called, author_low = NULL, author_doublet = NULL, spec = trial_spec()) 
{
    n <- nrow(metrics)
    trial_require(length(called) == n && !anyNA(called), "Cell call shape")
    flags <- data.frame(ZERO_UMI = metrics$original_UMI == 0, NOT_CELL_CALLED = !called)
    if (is.null(author_low)) 
        author_low <- rep(NA, n)
    if (is.null(author_doublet)) 
        author_doublet <- rep(NA, n)
    trial_require(is.logical(author_low) && is.logical(author_doublet) && length(author_low) == n && 
        length(author_doublet) == n, "Author flag mapping")
    threshold <- function(v) {
        if (!any(called)) 
            return(NA_real_)
        z <- log2(v[called])
        trial_require(all(is.finite(z)), "NONFINITE_QC_METRIC in threshold population")
        median(z) - 3 * stats::mad(z, center = median(z), constant = 1.4825999999999999)
    }
    umi <- threshold(metrics$original_UMI)
    genes <- threshold(metrics$detected_features)
    flags$NONFINITE_QC_METRIC <- called & !(is.finite(metrics$original_UMI) & is.finite(metrics$detected_features) & 
        is.finite(metrics$combined_MT_Hb_fraction))
    flags$AUTHOR_LOW_QUALITY <- !is.na(author_low) & author_low
    flags$AUTHOR_DOUBLET <- !is.na(author_doublet) & author_doublet
    range <- spec$cell_retention$low_quality$UMI_range_inclusive
    flags$UMI_OUTSIDE_RANGE <- called & (metrics$original_UMI < range[1] | metrics$original_UMI > range[2])
    flags$MT_HB_AT_LEAST_10_PERCENT <- called & is.finite(metrics$combined_MT_Hb_fraction) & metrics$combined_MT_Hb_fraction >= 
        spec$cell_retention$low_quality$combined_MT_Hb_fraction_exclusive_upper
    flags$LOW_UMI_MAD <- if (any(called)) 
        called & log2(metrics$original_UMI) < umi
    else rep(FALSE, n)
    flags$LOW_DETECTED_FEATURES_MAD <- if (any(called)) 
        called & log2(metrics$detected_features) < genes
    else rep(FALSE, n)
    trial_require(!anyNA(flags), "Missing QC exclusion flag")
    list(flags = flags, keep = rowSums(flags) == 0, author_low = author_low, author_doublet = author_doublet, 
        thresholds = list(log2_UMI_lower = umi, log2_detected_lower = genes, population = sum(called), 
            empty_population = !any(called)))
}

trial_doublets <- function (x, features, operation, spec = trial_spec(), diagnostic_path = NULL) 
{
    sets <- trial_sets(features, spec)
    y <- x[sets$doublet, , drop = FALSE]
    means <- Matrix::rowMeans(y)
    square <- y
    square@x <- square@x^2
    variance <- Matrix::rowMeans(square) - means^2
    nfeatures <- min(2000L, sum(is.finite(variance) & variance > 0))
    n <- ncol(y)
    dbr <- 0.0080000000000000002 * n/1000
    trial_require(n >= 10 && nfeatures >= 3 && dbr > 0 && dbr < 1, "UNRESOLVED_DOUBLET_FEASIBILITY_BOUNDARY")
    dims <- min(20L, n - 1L, nfeatures - 1L)
    a <- spec$cell_retention$doublets$arguments
    a$artificialDoublets <- max(1500L, n)
    a$dbr <- dbr
    a$dbr.sd <- 0.14999999999999999 * dbr
    a$nfeatures <- as.integer(nfeatures)
    a$dims <- as.integer(dims)
    a$includePCs <- as.integer(min(19L, dims))
    a$k <- 5L
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = y))
    audit <- new.env(parent = emptyenv())
    audit$entered <- 0L
    audit$completed <- list()
    audit$predictions <- list()
    observer <- local({
        e <- audit
        function(stage, value = NULL, nthreads = NULL) {
            if (stage == "prediction") 
                e$predictions[[length(e$predictions) + 1L]] <- list(length = length(value), finite = is.numeric(value) && 
                  length(value) > 0 && all(is.finite(value)))
            else if (stage == "enter") 
                e$entered <- e$entered + 1L
            else if (inherits(value, "xgb.Booster")) {
                config <- xgboost::xgb.config(value)
                trial_require(is.list(config) && is.list(config$learner$generic_param), "Native xgboost configuration schema mismatch")
                e$completed[[length(e$completed) + 1L]] <- list(nthreads = nthreads, generic_parameters = config$learner$generic_param)
            }
        }
    })
    old <- getOption("s10_b_trial.xgb_observer")
    options(s10_b_trial.xgb_observer = observer)
    on.exit(options(s10_b_trial.xgb_observer = old), add = TRUE)
    trace(".xgbtrain", where = asNamespace("scDblFinder"), print = FALSE, tracer = quote(getOption("s10_b_trial.xgb_observer")("enter")), 
        exit = quote(getOption("s10_b_trial.xgb_observer")("exit", returnValue(), nthreads)))
    on.exit(untrace(".xgbtrain", where = asNamespace("scDblFinder")), add = TRUE)
    trace("predict.xgboost", where = asNamespace("xgboost"), print = FALSE, exit = quote(getOption("s10_b_trial.xgb_observer")("prediction", 
        returnValue())))
    on.exit(untrace("predict.xgboost", where = asNamespace("xgboost")), add = TRUE)
    rng <- trial_seed(operation)
    result <- do.call(scDblFinder::scDblFinder, c(list(sce = sce), a, list(BPPARAM = rng$bp)))
    d <- as.data.frame(SummarizedExperiment::colData(result))
    output <- list(result = d, arguments = a, features = rownames(y), native_metadata = S4Vectors::metadata(result), 
        xgboost_audit = as.list(audit), rng = list(seed = rng$seed, before = rng$state, after = .Random.seed, 
            backend = "SerialParam", workers = 1L, RNGseed = BiocParallel::bpRNGseed(rng$bp)))
    if (!is.null(diagnostic_path)) 
        trial_rds(output, diagnostic_path)
    trial_require(identical(rownames(d), colnames(y)) && all(c("scDblFinder.score", "scDblFinder.class") %in% 
        names(d)), "Missing doublet output or barcode identity")
    trial_require(all(is.finite(d$scDblFinder.score)) && !anyNA(d$scDblFinder.class) && all(d$scDblFinder.class %in% 
        c("singlet", "doublet")), "Unresolved native doublet score/class")
    trial_require(audit$entered == a$iter && length(audit$completed) == a$iter, "Native xgboost training failure swallowed by scDblFinder")
    trial_require(length(audit$predictions) == a$iter && all(vapply(audit$predictions, function(z) z$finite, 
        TRUE)), "Native xgboost prediction failure swallowed by scDblFinder")
    trial_require(all(vapply(audit$completed, function(z) z$nthreads == 1 && as.integer(z$generic_parameters$nthread) == 
        1, TRUE)), "Native xgboost numerical thread mismatch")
    output
}

trial_normalize <- function (x) 
{
    den <- Matrix::colSums(x)
    trial_require(all(is.finite(den) & den > 0), "Annotation normalization denominator")
    y <- x %*% Matrix::Diagonal(x = 10000/den)
    y@x <- log2(1 + y@x)
    dimnames(y) <- dimnames(x)
    y
}

full_pooled_graph <- function (x, features, donor, spec = trial_spec()) 
{
    s <- trial_sets(features, spec)
    y <- trial_normalize(x[s$annotation, , drop = FALSE])
    detection <- Matrix::rowSums(x[s$annotation, , drop = FALSE] > 0) >= 10
    y <- y[detection, , drop = FALSE]
    rng <- trial_seed("post_effect_B.GSE150290.full_v1.pooled_pca")
    v <- scran::modelGeneVar(y, block = donor, design = NULL, subset.row = NULL, subset.fit = NULL, equiweight = TRUE, 
        method = "fisher", BPPARAM = rng$bp)
    eligible <- which(is.finite(v$bio) & v$bio > 0)
    eligible <- eligible[order(-v$bio[eligible], rownames(y)[eligible], method = "radix")]
    selected <- head(eligible, 2000L)
    trial_require(length(selected) >= 3L && ncol(x) >= 3L, "Insufficient pooled PCA features/cells")
    rank <- min(30L, ncol(y) - 1L, length(selected) - 1L)
    rng <- trial_seed("post_effect_B.GSE150290.full_v1.pooled_pca")
    pca <- BiocSingular::runPCA(Matrix::t(y[selected, , drop = FALSE]), rank = as.integer(rank), center = TRUE, 
        scale = FALSE, get.rotation = TRUE, get.pcs = TRUE, BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE))
    ids <- rownames(y)[selected]
    for (j in seq_len(ncol(pca$rotation))) {
        a <- abs(pca$rotation[, j])
        largest <- which(a == max(a))
        largest <- largest[order(ids[largest], method = "radix")][1]
        if (pca$rotation[largest, j] < 0) {
            pca$rotation[, j] <- -pca$rotation[, j]
            pca$x[, j] <- -pca$x[, j]
        }
    }
    graph_rng <- trial_seed("post_effect_B.GSE150290.full_v1.pooled_graph")
    graph <- bluster::clusterRows(pca$x, full = TRUE, BLUSPARAM = bluster::SNNGraphParam(k = as.integer(min(20L, 
        ncol(x) - 1L)), type = "rank", BNPARAM = BiocNeighbors::KmknnParam(), num.threads = 1L, BPPARAM = NULL, 
        cluster.fun = "louvain", cluster.args = list(resolution = 1)))
    cluster <- as.character(graph$clusters)
    keys <- vapply(split(colnames(x), cluster), function(z) sort(z, method = "radix")[1], "")
    labels <- setNames(sprintf("cluster_%04d", seq_along(keys)), names(keys)[order(keys, method = "radix")])
    list(cluster = unname(labels[cluster]), selected_feature_ids = ids, pca = pca, center = Matrix::rowMeans(y[selected, 
        , drop = FALSE]), graph = graph, seed_pca = rng$seed, seed_graph = graph_rng$seed)
}

trial_units <- function (x, donor, label, permitted_rows) 
{
    trial_require(length(donor) == ncol(x) && length(label) == ncol(x) && !anyNA(c(donor, label)), "Pooled unit identities")
    key <- paste(donor, label, sep = "|")
    ids <- sort(unique(key), method = "radix")
    membership <- Matrix::sparseMatrix(i = seq_along(key), j = match(key, ids), x = 1, dims = c(length(key), 
        length(ids)))
    trial_require(length(permitted_rows) == nrow(x) && is.logical(permitted_rows) && !anyNA(permitted_rows), 
        "Explicit annotation aggregation universe required")
    counts <- x[permitted_rows, , drop = FALSE] %*% membership
    colnames(counts) <- ids
    original_UMI <- as.numeric(Matrix::colSums(x) %*% membership)
    units <- data.frame(unit_id = ids, donor = vapply(strsplit(ids, "|", fixed = TRUE), `[`, "", 1L), 
        label = vapply(strsplit(ids, "|", fixed = TRUE), `[`, "", 2L), cells = as.integer(table(factor(key, 
            levels = ids))), UMI = original_UMI)
    units$eligible <- units$cells >= 20L & units$UMI >= 10000
    list(counts = counts, units = units)
}

trial_adjudicate <- function (scores, expected = NULL, minimum_donors) 
{
    if (!is.null(scores)) 
        trial_require(all(is.finite(scores)), "Nonfinite marker evidence cannot be omitted")
    if (is.null(scores) || nrow(scores) < minimum_donors) 
        return(list(state = "NA_NONESTIMABLE", reason = "INSUFFICIENT_ELIGIBLE_DONORS_OR_MARKER_MAPPING", 
            assignment = "UNRESOLVED", donors = if (is.null(scores)) 0L else nrow(scores)))
    maxima <- apply(scores, 1, function(x) {
        j <- which(x == max(x))
        if (length(j) == 1L) 
            colnames(scores)[j]
        else "AMBIGUOUS"
    })
    medians <- apply(scores, 2, median)
    j <- which(medians == max(medians))
    winner <- if (length(j) == 1L) 
        names(medians)[j]
    else "AMBIGUOUS"
    candidate <- if (is.null(expected)) 
        winner
    else expected
    pass <- candidate %in% colnames(scores) && winner == candidate && mean(maxima == candidate) >= 0.59999999999999998
    list(state = if (pass) "PASS" else "FAIL", reason = if (pass) "NONE" else "UNIQUE_SELECTION_OR_STRICT_MEDIAN_RULE_FAILED", 
        assignment = if (pass) candidate else "UNRESOLVED", donors = nrow(scores), unique_maxima = maxima, 
        median_scores = medians)
}

trial_seed <- function (operation) 
{
    x <- context$seed_map
    j <- match(operation, x$operation_name)
    ann_need(!is.na(j), "Unfrozen full-cohort operation seed")
    seed <- x$seed[j]
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    list(seed = seed, bp = BiocParallel::SerialParam(RNGseed = seed, stop.on.error = TRUE, progressbar = FALSE), 
        state = .Random.seed)
}
