# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.
full_gate_policy <- "GSE150290_EFFECTIVE_ANNOTATION_GATE_08_V1"

ann_compare <- function(x,y) {
  if(is.data.frame(x)&&is.data.frame(y)) {
    if(!identical(names(x),names(y))||nrow(x)!=nrow(y))return(FALSE)
    return(all(vapply(names(x),function(k) {
      a<-x[[k]];b<-y[[k]]
      if(!identical(is.na(a),is.na(b)))return(FALSE)
      if(all(is.na(a))||!length(a))return(TRUE)
      if(is.numeric(a)&&is.numeric(b))return(identical(as.double(a),as.double(b)))
      identical(a,b)
    },TRUE)))
  }
  identical(x,y)
}

mapping_gtf_attribute <- function (x, key, required = TRUE) 
{
    pattern <- paste0("(^|;)[[:space:]]*", key, "[[:space:]]+\"[^\"]*\"")
    hits <- regmatches(x, gregexpr(pattern, x, perl = TRUE))
    mapping_need(all(lengths(hits) <= 1L), paste("Duplicate GTF attribute", key))
    if (required) 
        mapping_need(all(lengths(hits) == 1L), paste("Missing GTF attribute", key))
    vapply(hits, function(z) if (!length(z)) 
        NA_character_
    else sub(".*\"([^\"]*)\"$", "\\1", z), "")
}

mapping_gtf_genes <- function (lines, line_numbers = seq_along(lines)) 
{
    keep_comments <- !startsWith(lines, "#")
    lines <- lines[keep_comments]
    line_numbers <- line_numbers[keep_comments]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    mapping_need(all(lengths(parts) == 9L), "GTF must have nine fields")
    keep <- vapply(parts, function(z) z[3] == "gene", TRUE)
    parts <- parts[keep]
    line_numbers <- line_numbers[keep]
    if (!length(parts)) 
        return(NULL)
    z <- do.call(rbind, parts)
    attrs <- z[, 9]
    data.frame(feature_id = mapping_gtf_attribute(attrs, "gene_id"), source_symbol = mapping_gtf_attribute(attrs, 
        "gene_name", required = FALSE), version = mapping_gtf_attribute(attrs, "gene_version", required = FALSE), 
        biotype = mapping_gtf_attribute(attrs, "gene_biotype"), seq_region = z[, 1], start = as.numeric(z[, 
            4]), end = as.numeric(z[, 5]), strand = z[, 7], gtf_line = line_numbers, stringsAsFactors = FALSE)
}

mapping_feature_masks <- function (features, identity, spec, cfg) 
{
    mapping_need(length(identity) == nrow(features) && !anyNA(identity), "Feature identity vector")
    a <- spec$annotation_qualification
    protected <- unique(unlist(cfg$programs, use.names = FALSE))
    audit <- unique(unlist(c(a$broad_audit_panels, a$lineage_audit_panels), use.names = FALSE))
    ambiguous <- duplicated(identity) | duplicated(identity, fromLast = TRUE)
    mt <- startsWith(features$source_symbol, "MT-")
    hb <- features$source_symbol %in% spec$cell_retention$hemoglobin_symbols
    excluded <- identity %in% c(protected, audit) | mt | hb
    assignment <- !excluded & !ambiguous & !identity %in% a$proliferation_markers_excluded_from_DR
    data.frame(source_row = features$source_row, feature_id = features$feature_id, protected = identity %in% 
        protected, held_out = identity %in% audit, mitochondrial = mt, hemoglobin = hb, proliferation = identity %in% 
        a$proliferation_markers_excluded_from_DR, D7 = identity %in% unlist(cfg$D7), ambiguous = ambiguous, 
        doublet_input = !excluded, annotation_input = assignment, assignment_background = assignment, 
        audit_background = !ambiguous & !identity %in% protected, stringsAsFactors = FALSE)
}

ann_approved_features <- function (features, crosswalk, assembly = "GRCh37", reference = "10x_hg19_1.2.0_GTF82_supported") 
{
    required <- c("feature_id", "source_symbol", "source_row")
    ann_need(all(required %in% names(features)) && !anyNA(features[, required]) && !anyDuplicated(features$source_row), 
        "Invalid source tuple keys")
    tuple <- list(original_feature_id = "ENSG00000132465", original_symbol = "IGJ", original_source_row = 7669L, 
        historical_gene_version = 6L, canonical_symbol = "JCHAIN", canonical_hgnc = "HGNC:5713", assembly = "GRCh37", 
        reference_scope = "10x_hg19_1.2.0_GTF82_supported")
    ann_need(nrow(crosswalk) == 1L && all(names(tuple) %in% names(crosswalk)), "Only one approved tuple is allowed")
    for (k in names(tuple)) ann_need(identical(as.character(crosswalk[[k]]), as.character(tuple[[k]])), 
        paste("Unapproved tuple", k))
    ann_need(isTRUE(crosswalk$active) && crosswalk$action == "APPROVED_SAME_GENE_RECONCILIATION" && assembly == 
        tuple$assembly && reference == tuple$reference_scope, "Unapproved activation or context")
    stable <- sub("[.][0-9]+$", "", features$feature_id)
    version <- ifelse(grepl("[.][0-9]+$", features$feature_id), sub(".*[.]", "", features$feature_id), 
        "")
    j <- which(stable == tuple$original_feature_id)
    ann_need(length(j) == 1L && !anyDuplicated(stable) && features$source_symbol[j] == tuple$original_symbol && 
        features$source_row[j] == tuple$original_source_row && version[j] %in% c("", "6") && sum(features$source_symbol == 
        "IGJ") == 1L && !any(features$source_symbol == "JCHAIN"), "Source tuple missing, ambiguous, versioned incorrectly or canonically colliding")
    out <- features
    out$canonical_symbol <- out$source_symbol
    out$canonical_hgnc <- NA_character_
    out$historical_gene_version <- NA_integer_
    out$reconciliation_state <- "UNCHANGED_LITERAL_SOURCE"
    out$canonical_symbol[j] <- "JCHAIN"
    out$canonical_hgnc[j] <- "HGNC:5713"
    out$historical_gene_version[j] <- 6L
    out$reconciliation_state[j] <- "APPROVED_UNIQUE_SAME_GENE_IDENTITY"
    ann_need(identical(out[, names(features)], features), "Source metadata changed")
    out
}

ann_masks <- function (features, identity = features$canonical_symbol, spec = trial_spec()) 
mapping_feature_masks(features, identity, spec, context$supplementary_config)

ann_sets <- function (features, spec = trial_spec()) 
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

ann_check_reuse <- function (features, spec = trial_spec()) 
{
    literal <- ann_masks(features, features$source_symbol, spec)
    canonical <- ann_masks(features, features$canonical_symbol, spec)
    flags <- setdiff(names(literal), c("source_row", "feature_id"))
    evidence <- data.frame(mask = flags, original_rows = vapply(literal[flags], sum, 0L), canonical_rows = vapply(canonical[flags], 
        sum, 0L), changed_rows = vapply(flags, function(k) sum(literal[[k]] != canonical[[k]]), 0L))
    ann_need(all(evidence$changed_rows == 0L), "Changed reuse mask: upstream rerun not authorized")
    evidence
}

ann_panel_scores <- function (counts, features, universe, panels) 
{
    ann_need(nrow(counts) == nrow(features) && length(universe) == nrow(features) && is.logical(universe) && 
        !anyNA(universe), "Panel universe identity")
    symbols <- features$canonical_symbol
    ann_need(length(symbols) == nrow(features) && !anyNA(symbols), "Missing canonical scoring interface")
    genes <- unlist(panels, use.names = FALSE)
    mapped <- vapply(genes, function(g) sum(symbols == g) == 1L && sum(symbols[universe] == g) == 1L, 
        TRUE)
    if (!all(mapped)) 
        return(list(state = "NA_NONESTIMABLE", reason = "REQUIRED_MARKER_MAPPING_FAILED", missing = unique(genes[!mapped]), 
            scores = NULL))
    if (!ncol(counts)) 
        return(list(state = "NA_NONESTIMABLE", reason = "NO_ELIGIBLE_DONOR_UNITS", missing = character(), 
            scores = matrix(numeric(), 0, length(panels), dimnames = list(NULL, names(panels)))))
    ann_need(sum(universe) >= 2L, "Insufficient marker rank universe")
    z <- as.matrix(counts[universe, , drop = FALSE])
    lib <- colSums(z)
    ann_need(all(lib > 0 & is.finite(lib)), "Invalid pooled marker normalization")
    z <- log2(1 + sweep(z, 2, lib, "/") * 1000000)
    ranks <- apply(z, 2, function(v) (rank(v, ties.method = "average") - 1)/(length(v) - 1))
    if (is.null(dim(ranks))) 
        ranks <- matrix(ranks, ncol = ncol(z))
    rownames(ranks) <- symbols[universe]
    scores <- vapply(panels, function(panel) colMeans(ranks[panel, , drop = FALSE]), numeric(ncol(z)))
    if (is.null(dim(scores))) 
        scores <- matrix(scores, nrow = ncol(z), dimnames = list(NULL, names(panels)))
    rownames(scores) <- colnames(counts)
    list(state = "SCORED", reason = "NONE", missing = character(), scores = scores)
}

trial_panel_scores <- function (counts, features, universe, panels) 
{
    ann_need(nrow(counts) == nrow(features) && length(universe) == nrow(features) && is.logical(universe) && 
        !anyNA(universe), "Panel universe identity")
    symbols <- features$canonical_symbol
    ann_need(length(symbols) == nrow(features) && !anyNA(symbols), "Missing canonical scoring interface")
    genes <- unlist(panels, use.names = FALSE)
    mapped <- vapply(genes, function(g) sum(symbols == g) == 1L && sum(symbols[universe] == g) == 1L, 
        TRUE)
    if (!all(mapped)) 
        return(list(state = "NA_NONESTIMABLE", reason = "REQUIRED_MARKER_MAPPING_FAILED", missing = unique(genes[!mapped]), 
            scores = NULL))
    if (!ncol(counts)) 
        return(list(state = "NA_NONESTIMABLE", reason = "NO_ELIGIBLE_DONOR_UNITS", missing = character(), 
            scores = matrix(numeric(), 0, length(panels), dimnames = list(NULL, names(panels)))))
    ann_need(sum(universe) >= 2L, "Insufficient marker rank universe")
    z <- as.matrix(counts[universe, , drop = FALSE])
    lib <- colSums(z)
    ann_need(all(lib > 0 & is.finite(lib)), "Invalid pooled marker normalization")
    z <- log2(1 + sweep(z, 2, lib, "/") * 1000000)
    ranks <- apply(z, 2, function(v) (rank(v, ties.method = "average") - 1)/(length(v) - 1))
    if (is.null(dim(ranks))) 
        ranks <- matrix(ranks, ncol = ncol(z))
    rownames(ranks) <- symbols[universe]
    scores <- vapply(panels, function(panel) colMeans(ranks[panel, , drop = FALSE]), numeric(ncol(z)))
    if (is.null(dim(scores))) 
        scores <- matrix(scores, nrow = ncol(z), dimnames = list(NULL, names(panels)))
    rownames(scores) <- colnames(counts)
    list(state = "SCORED", reason = "NONE", missing = character(), scores = scores)
}

trial_score_groups <- function (x, features, donor, label, universe, panels, groups = sort(unique(label)), expected = NULL, 
    minimum_donors = 2L) 
{
    u <- trial_units(x, donor, label, universe)
    f <- features[universe, , drop = FALSE]
    calls <- scores <- list()
    used <- character()
    for (g in groups) {
        keep <- which(u$units$eligible & u$units$label == g)
        s <- trial_panel_scores(u$counts[, keep, drop = FALSE], f, rep(TRUE, nrow(f)), panels)
        target <- if (is.null(expected)) 
            NULL
        else unname(expected[g])
        if (!is.null(target) && (!length(target) || is.na(target) || !target %in% names(panels))) {
            decision <- list(state = "NOT_APPLICABLE", reason = "NO_EXPECTED_PANEL_FOR_THIS_LABEL", assignment = "UNRESOLVED", 
                donors = length(keep))
        }
        else decision <- trial_adjudicate(s$scores, expected = target, minimum_donors = minimum_donors)
        calls[[length(calls) + 1L]] <- data.frame(group = g, expected = if (is.null(target)) 
            NA_character_
        else target, state = decision$state, reason = if (s$state == "NA_NONESTIMABLE") 
            paste(decision$reason, s$reason, sep = ";")
        else decision$reason, assignment = decision$assignment, eligible_donors = length(keep), scored_donors = if (is.null(s$scores)) 
            0L
        else nrow(s$scores), minimum_donors = minimum_donors, marker_mapping = s$state, missing_markers = paste(s$missing, 
            collapse = ";"))
        if (!is.null(s$scores) && nrow(s$scores)) {
            used <- union(used, unlist(panels, use.names = FALSE))
            for (j in seq_len(ncol(s$scores))) scores[[length(scores) + 1L]] <- data.frame(group = g, 
                unit_id = u$units$unit_id[keep], donor = u$units$donor[keep], panel = colnames(s$scores)[j], 
                score = as.numeric(s$scores[, j]))
        }
    }
    emptycalls <- data.frame(group = character(), expected = character(), state = character(), reason = character(), 
        assignment = character(), eligible_donors = integer(), scored_donors = integer(), minimum_donors = integer(), 
        marker_mapping = character(), missing_markers = character())
    emptyscores <- data.frame(group = character(), unit_id = character(), donor = character(), panel = character(), 
        score = numeric())
    list(calls = if (length(calls)) do.call(rbind, calls) else emptycalls, scores = if (length(scores)) do.call(rbind, 
        scores) else emptyscores, units = u$units, markers_actually_scored = used)
}

ann_annotate <- function (x, features, donor, cluster, spec = trial_spec()) 
{
    trial_require(ncol(x) == length(cluster) && length(donor) == length(cluster) && !anyNA(c(donor, cluster)) && 
        !anyDuplicated(colnames(x)), "Pooled annotation identities")
    sets <- trial_sets(features, spec)
    a <- spec$annotation_qualification
    broad <- trial_score_groups(x, features, donor, cluster, sets$assignment, a$assignment_panels)
    assigned <- setNames(broad$calls$assignment, broad$calls$group)
    broad_cell <- unname(assigned[cluster])
    trial_require(!anyNA(broad_cell), "Complete broad assignment mapping")
    cluster_audit <- trial_score_groups(x, features, donor, cluster, sets$audit_universe, a$broad_audit_panels, 
        groups = names(assigned), expected = assigned, minimum_donors = 8L)
    supported <- broad_cell != "UNRESOLVED"
    overall_groups <- sort(unique(c("epithelial", broad_cell[supported])))
    broad_audit <- trial_score_groups(x[, supported, drop = FALSE], features, donor[supported], broad_cell[supported], 
        sets$audit_universe, a$broad_audit_panels, groups = overall_groups, expected = setNames(overall_groups, 
            overall_groups), minimum_donors = 8L)
    epithelial <- broad_cell == "epithelial"
    fine <- trial_score_groups(x[, epithelial, drop = FALSE], features, donor[epithelial], cluster[epithelial], 
        sets$assignment, a$lineage_assignment_panels)
    fine_map <- setNames(fine$calls$assignment, fine$calls$group)
    proposed <- rep("UNKNOWN", ncol(x))
    proposed[supported & !epithelial] <- "NOT_EPITHELIAL"
    proposed[epithelial] <- unname(fine_map[cluster[epithelial]])
    proposed[proposed == "UNRESOLVED"] <- "other_unassigned"
    trial_require(!anyNA(proposed), "Complete finer assignment mapping")
    named_lineage <- epithelial & proposed != "other_unassigned"
    fine_groups <- sort(unique(proposed[named_lineage]))
    fine_audit <- trial_score_groups(x[, named_lineage, drop = FALSE], features, donor[named_lineage], 
        proposed[named_lineage], sets$audit_universe, a$lineage_audit_panels, groups = fine_groups, expected = setNames(fine_groups, 
            fine_groups), minimum_donors = 8L)
    final <- proposed
    passed <- fine_audit$calls$group[fine_audit$calls$state == "PASS"]
    final[epithelial & !proposed %in% passed] <- "other_unassigned"
    cells <- data.frame(cell_key = colnames(x), donor = donor, provisional_cluster = cluster, provisional_broad_label = broad_cell, 
        provisional_lineage = proposed, lineage_after_audit = final, epithelial_membership = ifelse(epithelial, 
            "PROVISIONAL_EPITHELIAL_PENDING_FORMAL_AUDIT", ifelse(supported, "PROVISIONAL_NON_EPITHELIAL", 
                "UNRESOLVED")), fine_audit_reason = ifelse(epithelial & !proposed %in% passed, "FINER_AUDIT_NOT_PASSED_OR_NONESTIMABLE", 
            "NOT_APPLICABLE_OR_PASSED"), author_malignancy_status = rep("UNKNOWN", ncol(x)), proliferation_state = rep("UNKNOWN", 
            ncol(x)))
    actual <- unique(c(broad$markers_actually_scored, cluster_audit$markers_actually_scored, broad_audit$markers_actually_scored, 
        fine$markers_actually_scored, fine_audit$markers_actually_scored))
    list(cells = cells, broad_assignment = broad, cluster_heldout_audit = cluster_audit, broad_heldout_audit = broad_audit, 
        lineage_assignment = fine, lineage_heldout_audit = fine_audit, D7 = data.frame(gene = sets$D7, 
            actually_scored = sets$D7 %in% actual, use = ifelse(sets$D7 %in% unlist(a$lineage_assignment_panels), 
                "LINEAGE_ASSIGNMENT", "HELD_OUT_LINEAGE_AUDIT"), interpretation = "ANNOTATION_LINKED_CONSISTENCY_NOT_INDEPENDENT_D7_VALIDATION"))
}

ann_reconcile <- function (map, annotation) 
{
    keep <- map$technical_state == "QC_SINGLET_PENDING_ANNOTATION"
    j <- match(map$cell_key[keep], annotation$cells$cell_key)
    trial_require(!anyNA(j) && !anyDuplicated(j), "Annotation/barcode bridge mismatch")
    a <- annotation$cells[j, , drop = FALSE]
    for (k in c("provisional_cluster", "provisional_broad_label", "provisional_lineage", "epithelial_membership")) map[keep, 
        k] <- a[[k]]
    map$lineage[keep] <- a$lineage_after_audit
    map$fine_audit_reason <- "NOT_EVALUATED"
    map$fine_audit_reason[keep] <- a$fine_audit_reason
    reasons <- ifelse(a$provisional_broad_label == "UNRESOLVED", "UNRESOLVED_EPITHELIAL_IDENTITY", ifelse(a$provisional_broad_label == 
        "epithelial", "RETAINED_EPITHELIAL", "RETAINED_NON_EPITHELIAL"))
    map$technical_state[keep] <- ifelse(reasons == "UNRESOLVED_EPITHELIAL_IDENTITY", "EXCLUDED_UNRESOLVED_BROAD_IDENTITY", 
        paste0(reasons, "_PROVISIONAL"))
    map$first_exclusion_reason[keep] <- reasons
    map$all_applicable_reasons[keep] <- reasons
    for (reason in c("UNRESOLVED_EPITHELIAL_IDENTITY", "RETAINED_EPITHELIAL", "RETAINED_NON_EPITHELIAL")) {
        map[[paste0("flag_", reason)]] <- FALSE
        map[keep, paste0("flag_", reason)] <- reasons == reason
    }
    map
}

ann_map_preservation <- function (original, map) 
{
    keep <- original$technical_state == "QC_SINGLET_PENDING_ANNOTATION"
    changed <- c("provisional_cluster", "provisional_broad_label", "provisional_lineage", "epithelial_membership", 
        "lineage", "technical_state", "first_exclusion_reason", "all_applicable_reasons")
    ann_need(ann_compare(map[, setdiff(names(original), changed)], original[, setdiff(names(original), 
        changed)]) && ann_compare(map[!keep, names(original)], original[!keep, ]), "Original upstream values changed")
    ann_need(identical(map$cell_key, original$cell_key) && !anyDuplicated(map$cell_key), "Barcode identity/order changed")
    invisible(TRUE)
}

audit_spec <- function (config = context$broad_audit_configuration) 
{
    p <- list(epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"), fibroblast = c("COL1A1", "COL1A2", "COL3A1", 
        "DCN"), endothelial = c("VWF", "EMCN", "KDR"), myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS"), 
        T_NK = c("CD3D", "CD3E", "NKG7"), B_plasma = c("CD79A", "CD74", "MS4A1", "MZB1"), mast = c("KIT", 
            "TPSAB1", "CPA3"))
    ann_need(identical(config$broad_audit_panels, p) && isTRUE(config$active) && identical(as.integer(config$complete_roster_sizes), 
        c(4L, 4L, 3L, 4L, 3L, 4L, 3L)) && config$minimum_cells == 20 && config$minimum_original_UMI == 
        10000 && config$minimum_assignment_donors == 2 && config$minimum_audit_donors == 8 && config$unique_maximum_fraction == 
        0.59999999999999998 && isTRUE(config$strict_median) && isTRUE(config$common_ranking_background) && 
        isTRUE(config$equal_donor_weight), "Unapproved audit roster or rule")
    s <- trial_spec()
    s$annotation_qualification$broad_audit_panels <- p
    s
}

audit_masks <- function (f, new = audit_spec(), old = trial_spec()) 
{
    x <- ann_masks(f, spec = old)
    y <- ann_masks(f, spec = new)
    ann_need(identical(x$source_row, y$source_row) && identical(x$feature_id, y$feature_id), "Changed source order")
    flags <- setdiff(names(x), c("source_row", "feature_id"))
    d <- data.frame(mask = flags, old_rows = vapply(x[flags], sum, 0L), new_rows = vapply(y[flags], sum, 
        0L), changed_rows = vapply(flags, function(k) sum(x[[k]] != y[[k]]), 0L))
    ann_need(all(d$changed_rows == 0L), "Changed upstream/ranking mask; reuse blocked")
    d
}

audit_reason <- function (calls) 
{
    vapply(seq_len(nrow(calls)), function(i) {
        if (calls$state[i] == "NOT_APPLICABLE") 
            return("NO_EXPECTED_PANEL_FOR_THIS_LABEL")
        z <- character()
        if (calls$marker_mapping[i] != "SCORED" && nzchar(calls$missing_markers[i])) 
            z <- c(z, "REQUIRED_MARKER_MAPPING_FAILED")
        if (calls$eligible_donors[i] == 0L) 
            z <- c(z, "NO_ELIGIBLE_DONOR_UNITS")
        if (calls$eligible_donors[i] < calls$minimum_donors[i]) 
            z <- c(z, "INSUFFICIENT_ELIGIBLE_DONORS_BELOW_8")
        if (!length(z)) 
            if (calls$state[i] == "PASS") 
                "NONE"
            else calls$reason[i]
        else paste(z, collapse = ";")
    }, "")
}

audit_diagnostics <- function (evidence) 
{
    out <- list()
    for (kind in c("cluster_heldout_audit", "broad_heldout_audit", "lineage_heldout_audit")) {
        e <- evidence[[kind]]
        for (g in e$calls$group) {
            z <- e$scores[e$scores$group == g, , drop = FALSE]
            call <- e$calls[e$calls$group == g, , drop = FALSE]
            if (!nrow(z)) 
                next
            donors <- sort(unique(z$donor), method = "radix")
            panels <- unique(z$panel)
            m <- matrix(NA_real_, length(donors), length(panels), dimnames = list(donors, panels))
            m[cbind(match(z$donor, donors), match(z$panel, panels))] <- z$score
            ann_need(!anyNA(m), "Incomplete competitor coverage")
            winners <- apply(m, 1, function(v) {
                j <- which(v == max(v))
                if (length(j) == 1L) 
                  colnames(m)[j]
                else "AMBIGUOUS"
            })
            target <- call$expected
            med <- apply(m, 2, median)
            valid <- !is.na(target) && target %in% panels
            frac <- if (valid) 
                mean(winners == target)
            else NA_real_
            strict <- if (valid) 
                all(med[target] > med[panels != target])
            else NA
            out[[length(out) + 1L]] <- data.frame(stage = kind, group = g, donor = donors, unique_maximum = winners, 
                expected = target, expected_unique_fraction = frac, strict_expected_median = strict, 
                formal_state = call$state, interpretation = "COMPUTABLE_SCORE_DIAGNOSTIC_ONLY; formal eight-donor gate unchanged")
        }
    }
    if (length(out)) 
        do.call(rbind, out)
    else data.frame(stage = character(), group = character(), donor = character(), unique_maximum = character(), 
        expected = character(), expected_unique_fraction = numeric(), strict_expected_median = logical(), 
        formal_state = character(), interpretation = character())
}

audit_evaluate <- function (x, features, donor, cluster, spec = audit_spec()) 
{
    a <- ann_annotate(x, features, donor, cluster, spec)
    for (k in c("cluster_heldout_audit", "broad_heldout_audit", "lineage_heldout_audit")) a[[k]]$calls$nonestimability_reason <- audit_reason(a[[k]]$calls)
    cluster_calls <- a$cluster_heldout_audit$calls
    label_calls <- a$broad_heldout_audit$calls
    j <- match(a$cells$provisional_cluster, cluster_calls$group)
    b <- match(a$cells$provisional_broad_label, label_calls$group)
    a$cells$cluster_formal_audit_state <- cluster_calls$state[j]
    a$cells$cluster_formal_audit_reason <- cluster_calls$nonestimability_reason[j]
    a$cells$broad_formal_audit_state <- ifelse(is.na(b), "NOT_APPLICABLE", label_calls$state[b])
    a$cells$broad_formal_audit_reason <- ifelse(is.na(b), "UNRESOLVED_BROAD_IDENTITY", label_calls$nonestimability_reason[b])
    a$cells$overall_epithelial_formal_audit_state <- ifelse(a$cells$provisional_broad_label == "epithelial", 
        label_calls$state[match("epithelial", label_calls$group)], "NOT_APPLICABLE")
    a$cells$audit_qualified_epithelial_membership <- ifelse(a$cells$provisional_broad_label == "epithelial", 
        ifelse(a$cells$overall_epithelial_formal_audit_state == "PASS", "INCLUDED_FORMAL_AUDIT_PASS", 
            "UNRESOLVED_PENDING_FORMAL_AUDIT"), ifelse(a$cells$provisional_broad_label == "UNRESOLVED", 
            "UNRESOLVED_BROAD_IDENTITY", "PROVISIONAL_NON_EPITHELIAL"))
    a$diagnostics <- audit_diagnostics(a)
    a
}

audit_reconcile <- function (original, a) 
{
    m <- ann_reconcile(original, a)
    ann_map_preservation(original, m)
    keep <- original$technical_state == "QC_SINGLET_PENDING_ANNOTATION"
    j <- match(original$cell_key[keep], a$cells$cell_key)
    fields <- c("cluster_formal_audit_state", "cluster_formal_audit_reason", "broad_formal_audit_state", 
        "broad_formal_audit_reason", "overall_epithelial_formal_audit_state", "audit_qualified_epithelial_membership")
    for (k in fields) {
        m[[k]] <- rep("NOT_EVALUATED_UPSTREAM_EXCLUDED", nrow(m))
        m[[k]][keep] <- a$cells[[k]][j]
    }
    m
}

full_graph_serial <- function (g, keys) 
list(cell_keys = keys, cluster = g$cluster, selected_feature_ids = g$selected_feature_ids, rotation = g$pca$rotation, 
    PC_scores = g$pca$x, centering = g$center, edges = igraph::as_edgelist(g$graph$objects$graph, names = FALSE), 
    edge_weights = igraph::E(g$graph$objects$graph)$weight, seed_pca = g$seed_pca, seed_graph = g$seed_graph)

full_gate_need <- function (ok, message) 
if (!isTRUE(ok)) stop(message, call. = FALSE)

full_gate_digest <- function (x) 
digest::digest(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null", digits = NA), algo = "sha256", 
    serialize = FALSE)

full_gate_context <- function (overall_calls, cluster_calls, source_bindings) 
list(call_fingerprint = full_gate_digest(list(overall = overall_calls, clusters = cluster_calls)), source_bindings = source_bindings)

full_gate_validate_calls <- function (x, minimum_audit_donors) 
{
    fields <- c("group", "expected", "state", "eligible_donors", "scored_donors", "minimum_donors", "marker_mapping")
    full_gate_need(is.data.frame(x) && all(fields %in% names(x)), "Gate call schema incomplete")
    full_gate_need(!anyNA(x[, fields]) && !anyDuplicated(x$group) && all(nzchar(x$group)), "Gate calls have missing or duplicate identities")
    full_gate_need(all(x$state %in% c("PASS", "FAIL", "NA_NONESTIMABLE", "NOT_APPLICABLE")), "Unknown formal audit state")
    for (field in c("eligible_donors", "scored_donors", "minimum_donors")) full_gate_need(is.numeric(x[[field]]) && 
        all(is.finite(x[[field]])) && all(x[[field]] >= 0 & x[[field]] == floor(x[[field]])), "Invalid donor accounting")
    full_gate_need(all(x$minimum_donors == minimum_audit_donors) && all(x$scored_donors <= x$eligible_donors), 
        "Contradictory donor gates")
    estimable <- x$state %in% c("PASS", "FAIL")
    full_gate_need(all(x$eligible_donors[estimable] >= minimum_audit_donors) && all(x$scored_donors[estimable] == 
        x$eligible_donors[estimable]) && all(x$marker_mapping[estimable] == "SCORED"), "Contradictory estimable audit")
    invisible(TRUE)
}

full_effective_admission_gate <- function (overall_calls, cluster_calls, context, dispositions = list(), disposition_registry = list(), 
    minimum_audit_donors = 8L) 
{
    full_gate_validate_calls(overall_calls, minimum_audit_donors)
    full_gate_validate_calls(cluster_calls, minimum_audit_donors)
    full_gate_need(identical(context$call_fingerprint, full_gate_context(overall_calls, cluster_calls, 
        context$source_bindings)$call_fingerprint), "Stale or mismatched audit context")
    full_gate_need(length(context$source_bindings) > 0L, "Missing saved-evidence binding")
    e <- overall_calls[overall_calls$group == "epithelial", , drop = FALSE]
    full_gate_need(nrow(e) <= 1L && (!nrow(e) || e$expected == "epithelial"), "Contradictory overall epithelial identity")
    overall_state <- if (nrow(e)) 
        e$state
    else "NA_NONESTIMABLE"
    groups <- vapply(dispositions, function(d) d$group, "")
    full_gate_need(!anyDuplicated(groups), "Duplicate disposition identity")
    failed <- cluster_calls$group[cluster_calls$state == "FAIL"]
    full_gate_need(all(groups %in% failed), "Stale disposition does not match a failed cluster")
    full_gate_need(identical(sort(as.character(names(disposition_registry))), sort(groups)), "Disposition registry incomplete or mismatched")
    findings <- cluster_calls[, c("group", "expected", "state", "eligible_donors", "scored_donors"), 
        drop = FALSE]
    findings$disposition_status <- rep("NOT_REQUIRED", nrow(findings))
    findings$epithelial_boundary_relevance <- rep("NOT_ADJUDICATED_BY_GATE", nrow(findings))
    findings$blocking <- rep(FALSE, nrow(findings))
    findings$reason <- rep("NONE", nrow(findings))
    for (i in seq_len(nrow(findings))) {
        if (findings$state[i] == "NA_NONESTIMABLE") {
            findings$reason[i] <- "CLUSTER_AUDIT_NONESTIMABLE_RETAINED_NOT_A_SEPARATE_EIGHT_DONOR_VETO"
            next
        }
        if (findings$state[i] != "FAIL") 
            next
        j <- match(findings$group[i], groups)
        if (is.na(j)) {
            findings$disposition_status[i] <- "MISSING_REQUIRES_REVIEW"
            findings$epithelial_boundary_relevance[i] <- "UNKNOWN"
            findings$blocking[i] <- TRUE
            findings$reason[i] <- "ESTIMABLE_CLUSTER_DISAGREEMENT_WITHOUT_BOUND_DISPOSITION"
            next
        }
        d <- dispositions[[j]]
        full_gate_need(identical(disposition_registry[[d$group]], full_gate_digest(d)), "Incorrectly authenticated disposition")
        full_gate_need(identical(d$audit_fingerprint, context$call_fingerprint) && identical(d$source_bindings, 
            context$source_bindings) && identical(d$expected, findings$expected[i]), "Stale or mismatched disposition evidence")
        full_gate_need(is.character(d$reason) && length(d$reason) == 1L && nzchar(d$reason) && is.character(d$authority) && 
            length(d$authority) == 1L && nzchar(d$authority), "Incomplete disposition rationale/authority")
        full_gate_need(d$relevance %in% c("EPITHELIAL_BOUNDARY", "UNRELATED", "UNKNOWN") && d$status %in% 
            c("PENDING", "DOCUMENTED_UNRELATED"), "Unsupported disposition state")
        unrelated <- d$relevance == "UNRELATED" && d$status == "DOCUMENTED_UNRELATED"
        full_gate_need((unrelated && d$expected != "epithelial" && isFALSE(d$pending_adjudication)) || 
            (!unrelated && d$status == "PENDING" && d$relevance != "UNRELATED" && isTRUE(d$pending_adjudication)), 
            "Contradictory disposition clearance")
        findings$disposition_status[i] <- d$status
        findings$epithelial_boundary_relevance[i] <- d$relevance
        findings$blocking[i] <- !unrelated
        findings$reason[i] <- d$reason
    }
    blocks <- findings[findings$blocking, , drop = FALSE]
    permit <- identical(overall_state, "PASS") && nrow(blocks) == 0L
    reason <- if (!identical(overall_state, "PASS")) 
        "OVERALL_EPITHELIAL_AUDIT_NOT_PASS"
    else if (nrow(blocks)) 
        blocks$reason[1]
    else "NONE"
    list(schema_version = 1L, policy_id = full_gate_policy, overall_epithelial_state = overall_state, 
        subsidiary_findings = findings, epithelial_boundary_adjudication = if (nrow(blocks)) "PENDING" else "NO_PENDING_BOUNDARY_CONFLICT_IN_BOUND_EVIDENCE", 
        annotation_freeze_allowed = permit, effective_annotation_freeze_permission = permit, population_freeze_allowed = permit, 
        reason = reason, blocking_findings = blocks, input_context = context, disposition_bindings = disposition_registry, 
        scope = "Gate permission only; no owner acceptance, population freeze or effect release")
}

full_gate_assert_permission <- function (gate) 
{
    full_gate_need(identical(gate$policy_id, full_gate_policy), "Legacy component gate is not current authority")
    flags <- gate[c("annotation_freeze_allowed", "effective_annotation_freeze_permission", "population_freeze_allowed")]
    full_gate_need(length(flags) == 3L && all(vapply(flags, isTRUE, TRUE)), paste("Effective population freeze blocked:", 
        gate$reason))
    full_gate_need(identical(gate$overall_epithelial_state, "PASS") && identical(gate$reason, "NONE") && 
        nrow(gate$blocking_findings) == 0L && !identical(gate$epithelial_boundary_adjudication, "PENDING"), 
        "Contradictory effective permission")
    invisible(TRUE)
}

full_gate_require_saved <- function (gate_path, gate_identity, recompute, verify_identity) 
{
    full_gate_need(identical(gate_identity$path, gate_path), "Wrong effective gate path")
    verify_identity(gate_identity)
    saved <- jsonlite::read_json(gate_path, simplifyVector = FALSE)
    full_gate_need(identical(saved$policy_id, full_gate_policy), "Legacy component gate cannot authorize preparation")
    full_gate_need(length(saved$input_bindings) > 0L && length(saved$implementation_bindings) > 0L, "Incomplete saved gate input/implementation bindings")
    for (binding in c(saved$input_bindings, saved$implementation_bindings)) verify_identity(binding)
    current <- recompute()
    canonical <- jsonlite::fromJSON(jsonlite::toJSON(current, auto_unbox = TRUE, null = "null", na = "null", 
        digits = NA), simplifyVector = FALSE)
    full_gate_need(identical(saved$result, canonical), "Stale or contradictory current effective gate")
    full_gate_assert_permission(current)
}
