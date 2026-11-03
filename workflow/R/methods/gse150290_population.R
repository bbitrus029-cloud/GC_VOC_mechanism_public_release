# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.

s10b_cell_keys <- function (sample_id, library_id, barcode) 
{
    s10b_require(length(sample_id) == length(barcode) && length(library_id) == length(barcode), "Cell identity lengths")
    s10b_require(!anyNA(c(sample_id, library_id, barcode)) && all(nzchar(c(sample_id, library_id, barcode))) && 
        !any(grepl("[|\t\r\n]", c(sample_id, library_id, barcode))), "Cell key components")
    key <- paste(sample_id, library_id, barcode, sep = "|")
    s10b_require(!anyDuplicated(key), "Duplicated physical cell identity or technical duplicate")
    key
}

s10b_validate_cell_map <- function (map, keys) 
{
    required <- c("cell_key", "sample_id", "technical_library", "source_barcode", "donor", "tissue_origin", 
        "original_source_label", "qc_pass", "epithelial_membership", "inclusion_reason", "lineage", "author_malignancy_status", 
        "author_malignancy_evidence", "proliferation_state")
    s10b_require(all(required %in% names(map)), "Incomplete annotation map fields")
    s10b_require(identical(as.character(map$cell_key), as.character(keys)), "Cell map is not exhaustive and ordered")
    s10b_require(identical(map$cell_key, s10b_cell_keys(map$sample_id, map$technical_library, map$source_barcode)), 
        "Physical cell key mismatch")
    s10b_require(all(map$tissue_origin %in% c("NAT", "Tumor")) && !anyNA(map$qc_pass) && is.logical(map$qc_pass), 
        "Unresolved tissue or QC")
    for (k in c("donor", "epithelial_membership", "inclusion_reason", "lineage", "author_malignancy_status", 
        "author_malignancy_evidence", "proliferation_state")) s10b_require(!anyNA(map[[k]]) && all(nzchar(map[[k]])), 
        paste("Missing annotation field", k))
    s10b_require(all(map$epithelial_membership %in% c("INCLUDED", "EXCLUDED", "UNRESOLVED")), "Invalid epithelial membership")
    s10b_require(all(map$lineage[map$epithelial_membership == "INCLUDED"] != "UNKNOWN"), "Finer uncertainty must be other_unassigned epithelial")
    libs <- unique(map[, c("technical_library", "sample_id", "donor", "tissue_origin")])
    s10b_require(!anyDuplicated(libs$technical_library), "Technical library assigned to conflicting identities")
    s10b_require(!any(map$author_malignancy_status != "UNKNOWN" & map$author_malignancy_evidence %in% 
        c("UNKNOWN", "TISSUE_ORIGIN_ONLY")), "Malignancy requires author evidence")
    invisible(TRUE)
}

s10b_aggregate <- function (counts, map, qualification, endpoint = "overall_epithelial") 
{
    s10b_require(identical(qualification$QC, "VALIDATED") && identical(qualification$annotation, "VALIDATED") && 
        nzchar(qualification$specification_sha256), "QC and annotation qualification required before populations or H coverage")
    s10b_validate_cell_map(map, colnames(counts))
    values <- if (inherits(counts, "sparseMatrix")) 
        counts@x
    else counts
    s10b_require(!is.null(rownames(counts)) && !anyNA(rownames(counts)) && all(nzchar(rownames(counts))) && 
        !anyDuplicated(rownames(counts)) && all(is.finite(values)) && all(values >= 0) && all(values == 
        round(values)), "Invalid aggregation UMI or feature identity")
    keep <- map$qc_pass & map$epithelial_membership == "INCLUDED"
    if (endpoint != "overall_epithelial") 
        keep <- keep & map$lineage == endpoint
    if (!any(keep)) 
        return(list(counts = matrix(numeric(), nrow(counts), 0L, dimnames = list(rownames(counts), character())), 
            units = data.frame(unit_id = character(), donor = character(), tissue = character(), cell_count = integer(), 
                technical_libraries = character()), biological_unit = "donor", original_counts = TRUE, 
            endpoint = endpoint, population_reason = "NO_RETAINED_EPITHELIAL_CELLS"))
    z <- map[keep, , drop = FALSE]
    unit <- paste(z$donor, z$tissue_origin, sep = "|")
    ids <- sort(unique(unit), method = "radix")
    membership <- Matrix::sparseMatrix(i = seq_along(unit), j = match(unit, ids), x = 1, dims = c(length(unit), 
        length(ids)))
    summed <- as.matrix(counts[, keep, drop = FALSE] %*% membership)
    colnames(summed) <- ids
    units <- data.frame(unit_id = ids, donor = vapply(strsplit(ids, "|", fixed = TRUE), `[`, "", 1L), 
        tissue = vapply(strsplit(ids, "|", fixed = TRUE), `[`, "", 2L), cell_count = as.integer(table(factor(unit, 
            levels = ids))), technical_libraries = vapply(ids, function(id) paste(sort(unique(z$technical_library[unit == 
            id]), method = "radix"), collapse = "|"), ""))
    units$tissue[units$tissue == "NAT"] <- "Normal"
    list(counts = summed, units = units, biological_unit = "donor", original_counts = TRUE)
}

s10b_h4_members <- function () 
{
    h <- unlist(context$supplementary_config$programs$H_core4, use.names = FALSE)
    s10b_require(is.character(h) && length(h) == 4L && !anyNA(h) && !anyDuplicated(h) && all(nzchar(h)), 
        "Invalid authoritative H4 configuration")
    h
}

s10b_prepare <- function (counts, units, members, cfg) 
{
    h <- s10b_h4_members()
    s10b_require(!missing(members) && is.character(members) && length(members) == length(h) && !anyNA(members) && 
        !anyDuplicated(members) && setequal(members, h), "MALFORMED_H4_ROSTER: request the complete unique configured four-member roster")
    requested_order <- members
    members <- h
    s10b_require(is.matrix(counts) && is.numeric(counts) && nrow(counts) > 0L, "Preparation counts must be a numeric features-by-units matrix")
    s01_validate_identity_vector(rownames(counts), "preparation gene row names")
    s10b_require(is.data.frame(units) && all(c("unit_id", "donor", "tissue", "cell_count") %in% names(units)) && 
        nrow(units) == ncol(counts), "Invalid preparation unit schema or dimensions")
    s10b_require(is.numeric(units$cell_count) && all(is.finite(units$cell_count)) && all(units$cell_count >= 
        0) && all(units$cell_count == round(units$cell_count)), "Preparation cell_count must contain finite nonnegative integers")
    if (ncol(counts) > 0L) {
        s01_validate_numeric_matrix(counts, "B preparation", require_counts = TRUE)
        s01_assert_exact_order(units$unit_id, colnames(counts), "preparation unit")
        s01_assert_unique_key(units, c("donor", "tissue"), "preparation donor tissue")
        s10b_require(all(units$tissue %in% c("Normal", "Tumor")), "Invalid preparation tissue labels")
    }
    original <- colSums(counts)
    eligible <- units$cell_count >= cfg$minimum_cells & original >= cfg$minimum_original_UMI
    if (!any(eligible)) {
        p <- list(original_library_sizes = original, eligible = eligible, eligible_units = character(), 
            complete_pair_donors = character(), modeling_units = character(), background = character(), 
            modeling_background = character(), counts = counts[integer(), integer(), drop = FALSE], modeling_library_sizes = numeric(), 
            sample_data = units[FALSE, , drop = FALSE], design = NULL, detection = data.frame(gene_id = rownames(counts), 
                eligible_units = 0L, nonzero_units = 0L, nonzero_fraction = NA_real_, paired_both_nonzero = 0L, 
                available = FALSE), model_fitted = FALSE, normalization_run = FALSE)
        structural <- if (ncol(counts) == 0L) 
            "NO_RETAINED_EPITHELIAL_CELLS"
        else "NO_ELIGIBLE_DONOR_TISSUE_UNITS"
    }
    else {
        p <- s10b_prepare_core(counts, units, cfg)
        structural <- if (length(p$complete_pair_donors) < cfg$minimum_pairs) 
            "INSUFFICIENT_COMPLETE_PAIRS"
        else "NONE"
    }
    p$input_units <- units
    p$detection$fraction_state <- if (any(eligible)) 
        "DEFINED"
    else "UNDEFINED"
    p$detection$fraction_reason <- if (any(eligible)) 
        "NONE"
    else "ZERO_ELIGIBLE_UNIT_DENOMINATOR"
    index <- match(members, p$detection$gene_id)
    available <- !is.na(index) & p$detection$available[index] & members %in% p$modeling_background
    available[is.na(available)] <- FALSE
    member_state <- data.frame(member = members, technical_eligible = structural == "NONE" & available, 
        reason = if (structural != "NONE") 
            structural
        else ifelse(available, "RUNTIME_MODEL_PENDING", "MEMBER_MAPPING_OR_DETECTION_FAILED"))
    p$member_qualification <- member_state
    p$requested_member_order <- requested_order
    p$member_roster_status <- "VALID_COMPLETE_H4_CANONICAL_ORDER"
    p$exact_H4_technical_eligible <- structural == "NONE" && all(available)
    p$structural_reason <- structural
    p$effect_runtime_state <- "PENDING_EFFECT_RELEASE"
    p
}

s10b_composition_population <- function (map, cfg) 
{
    s10b_validate_cell_map(map, map$cell_key)
    z <- map[map$qc_pass & map$epithelial_membership == "INCLUDED", ]
    count <- as.data.frame(table(factor(z$donor, levels = sort(unique(z$donor), method = "radix")), factor(z$tissue_origin, 
        levels = c("NAT", "Tumor"))))
    names(count) <- c("donor", "tissue", "retained_epithelial_cells")
    count$eligible <- count$retained_epithelial_cells >= cfg$composition_minimum_epithelial_cells_per_tissue
    pairs <- sort(names(which(table(as.character(count$donor[count$eligible])) == 2L)), method = "radix")
    list(counts = count, complete_pair_donors = pairs, structural_eligible = length(pairs) >= cfg$minimum_pairs, 
        effect_runtime_state = "PENDING_EFFECT_RELEASE")
}

s10b_lineage_containment <- function (lineage_donors, overall_donors) 
{
    s10b_require(!anyDuplicated(lineage_donors) && !anyDuplicated(overall_donors), "Duplicated donor identity")
    list(available = all(lineage_donors %in% overall_donors), unmatched_donors = setdiff(lineage_donors, 
        overall_donors), scores_computed = FALSE)
}

s10b_roster <- function (source, lineages, categories, lineage_eligibility) 
{
    s10b_require(!is.null(lineages) && !is.null(categories), "Unresolved lineage or composition taxonomy cannot become L=0 or K=1")
    s10b_require(!anyDuplicated(lineages) && !anyDuplicated(categories) && length(categories) > 0L, "Invalid frozen taxonomy")
    s10b_require(length(lineage_eligibility) == length(lineages) && !anyNA(lineage_eligibility), "Lineage eligibility roster mismatch")
    h <- s10b_h4_members()
    make <- function(family, endpoint, member = "H_core4", technical = "PENDING_TECHNICAL_CHECK") data.frame(hypothesis_id = paste("HYP_PE_B", 
        source, endpoint, member, sep = "_"), family = family, endpoint = endpoint, member = member, 
        technical_state = technical, inference_state = "PENDING_EFFECT_RELEASE", timing = "post_effect_exploratory")
    rows <- list(make("PEB_EPI_H4", "overall_epithelial"))
    for (i in seq_along(lineages)) rows[[length(rows) + 1L]] <- make("PEB_LINEAGE_H4", lineages[i], technical = if (lineage_eligibility[i]) 
        "TECHNICALLY_ELIGIBLE"
    else "TECHNICALLY_UNAVAILABLE_PLANNED_SLOT_RETAINED")
    for (e in c("overall_epithelial", lineages)) for (g in h) rows[[length(rows) + 1L]] <- make("PEB_H_MEMBERS", 
        e, g)
    for (k in categories) rows[[length(rows) + 1L]] <- make("PEB_CAPTURED_COMPOSITION", paste0("composition_", 
        k), "captured_category")
    out <- do.call(rbind, rows)
    sizes <- c(PEB_EPI_H4 = 1L, PEB_LINEAGE_H4 = length(lineages), PEB_H_MEMBERS = 4L * (length(lineages) + 
        1L), PEB_CAPTURED_COMPOSITION = length(categories))
    out$family_size <- unname(sizes[out$family])
    s10b_require(!anyDuplicated(out$hypothesis_id), "Hypothesis identity collision")
    list(hypotheses = out, L = length(lineages), K = length(categories), L_H4_eligible = sum(lineage_eligibility), 
        family_sizes = sizes)
}

s10b_apply_endpoint_qualifications <- function (roster, qualifications) 
{
    out <- roster
    rows <- out$hypotheses
    expression <- rows$family != "PEB_CAPTURED_COMPOSITION"
    endpoints <- unique(rows$endpoint[expression])
    s10b_require(!is.null(names(qualifications)) && !anyDuplicated(names(qualifications)) && setequal(names(qualifications), 
        endpoints), "Incomplete endpoint qualification join")
    rows$technical_reason <- "PENDING_TECHNICAL_CHECK"
    for (e in endpoints) {
        p <- qualifications[[e]]
        m <- p$member_qualification
        s10b_require(identical(m$member, s10b_h4_members()) && nrow(m) == 4L && !anyNA(m$technical_eligible) && 
            !anyNA(m$reason), "Invalid endpoint member qualification")
        j <- which(expression & rows$endpoint == e)
        for (i in j) {
            program <- rows$member[i] == "H_core4"
            eligible <- if (program) 
                p$exact_H4_technical_eligible
            else m$technical_eligible[match(rows$member[i], m$member)]
            reason <- if (program) {
                if (p$structural_reason != "NONE") 
                  p$structural_reason
                else if (eligible) 
                  "RUNTIME_PROGRAM_PENDING"
                else "INCOMPLETE_H4_MEASUREMENT"
            }
            else m$reason[match(rows$member[i], m$member)]
            s10b_require(length(eligible) == 1L && !is.na(eligible) && length(reason) == 1L && !is.na(reason), 
                "Invalid output mapping")
            rows$technical_state[i] <- if (eligible) 
                "TECHNICALLY_ELIGIBLE"
            else "TECHNICALLY_UNAVAILABLE_PLANNED_SLOT_RETAINED"
            rows$technical_reason[i] <- reason
        }
    }
    out$hypotheses <- rows
    out$L_H4_eligible <- sum(vapply(qualifications[setdiff(endpoints, "overall_epithelial")], function(p) isTRUE(p$exact_H4_technical_eligible), 
        TRUE))
    out
}

prep_feature_bridge <- function (features, members = s10b_h4_members()) 
{
    trial_require(all(c("feature_id", "source_symbol", "source_row", "canonical_symbol") %in% names(features)) && 
        !anyDuplicated(features$feature_id) && !anyNA(features$feature_id), "Invalid original feature identities")
    do.call(rbind, lapply(members, function(g) {
        i <- which(!is.na(features$canonical_symbol) & features$canonical_symbol == g)
        data.frame(member = g, matches = length(i), mapping_state = if (!length(i)) 
            "MISSING"
        else if (length(i) > 1L) 
            "AMBIGUOUS"
        else "UNIQUE", feature_id = if (length(i) == 1L) 
            features$feature_id[i]
        else NA_character_, source_symbol = if (length(i) == 1L) 
            features$source_symbol[i]
        else NA_character_, source_row = if (length(i) == 1L) 
            features$source_row[i]
        else NA_integer_)
    }))
}

prep_endpoint <- function (counts, units, features, cfg) 
{
    trial_require(identical(rownames(counts), features$feature_id), "Count feature order changed")
    original <- colSums(counts)
    p <- s10b_prepare(counts, units, s10b_h4_members(), cfg)
    bridge <- prep_feature_bridge(features)
    index <- match(bridge$feature_id, p$detection$gene_id)
    detection_ok <- !is.na(index) & p$detection$available[index] & bridge$feature_id %in% p$modeling_background
    detection_ok[is.na(detection_ok)] <- FALSE
    available <- bridge$mapping_state == "UNIQUE" & detection_ok
    reason <- if (p$structural_reason != "NONE") 
        rep(p$structural_reason, 4L)
    else ifelse(bridge$mapping_state != "UNIQUE", paste0("MEMBER_IDENTITY_", bridge$mapping_state), ifelse(available, 
        "RUNTIME_MODEL_PENDING", "MEMBER_DETECTION_OR_PAIRED_NONZERO_FAILED"))
    p$member_qualification <- cbind(bridge, technical_eligible = p$structural_reason == "NONE" & available, 
        reason = reason, eligible_units = ifelse(is.na(index), sum(p$eligible), p$detection$eligible_units[index]), 
        nonzero_units = ifelse(is.na(index), NA_integer_, p$detection$nonzero_units[index]), nonzero_fraction = ifelse(is.na(index), 
            NA_real_, p$detection$nonzero_fraction[index]), paired_both_nonzero = ifelse(is.na(index), 
            NA_integer_, p$detection$paired_both_nonzero[index]))
    p$exact_H4_technical_eligible <- p$structural_reason == "NONE" && all(available)
    p$identifier_adapter <- "Exact canonical member to one original feature row; all original rows and library totals retained"
    trial_require(identical(original, p$original_library_sizes), "Identifier adaptation changed original totals")
    p
}

prep_composition <- function (map, roster, categories, cfg) 
{
    s10b_validate_cell_map(map, map$cell_key)
    trial_require(!anyDuplicated(categories) && length(categories) > 0L, "Invalid frozen K")
    donors <- sort(unique(roster$donor), method = "radix")
    units <- expand.grid(tissue_origin = c("NAT", "Tumor"), donor = donors, stringsAsFactors = FALSE)
    units <- units[, c("donor", "tissue_origin")]
    keys <- paste(units$donor, units$tissue_origin, sep = "|")
    sample_keys <- paste(roster$donor, roster$tissue, sep = "|")
    trial_require(!anyDuplicated(sample_keys), "Conflicting source donor/tissue identities")
    units$sample_observed <- keys %in% sample_keys
    z <- map[map$qc_pass & map$epithelial_membership == "INCLUDED", , drop = FALSE]
    trial_require(all(z$lineage %in% categories), "K is not exhaustive")
    counts <- as.integer(table(factor(paste(z$donor, z$tissue_origin, sep = "|"), levels = keys)))
    units$retained_epithelial_cells <- ifelse(units$sample_observed, counts, NA_integer_)
    units$eligible <- units$sample_observed & !is.na(units$retained_epithelial_cells) & units$retained_epithelial_cells >= 
        cfg$composition_minimum_epithelial_cells_per_tissue
    units$status <- ifelse(!units$sample_observed, "MISSING_TISSUE_SAMPLE", ifelse(units$eligible, "ELIGIBLE", 
        "BELOW_COMPOSITION_CELL_GATE"))
    pairs <- donors[vapply(donors, function(d) all(units$eligible[units$donor == d]), TRUE)]
    catcounts <- do.call(rbind, lapply(categories, function(k) {
        n <- as.integer(table(factor(paste(z$donor[z$lineage == k], z$tissue_origin[z$lineage == k], 
            sep = "|"), levels = keys)))
        data.frame(units, category = k, category_cells = ifelse(units$sample_observed, n, NA_integer_))
    }))
    sums <- vapply(keys, function(key) sum(catcounts$category_cells[paste(catcounts$donor, catcounts$tissue_origin, 
        sep = "|") == key]), 0)
    trial_require(identical(unname(is.na(sums)), !units$sample_observed) && all(sums[units$sample_observed] == 
        units$retained_epithelial_cells[units$sample_observed]), "K denominator accounting")
    list(units = units, category_counts = catcounts, complete_pair_donors = pairs, structural_eligible = length(pairs) >= 
        cfg$minimum_pairs, structural_reason = if (length(pairs) >= cfg$minimum_pairs) "NONE" else "INSUFFICIENT_COMPOSITION_PAIRS", 
        effect_runtime_state = "PENDING_EFFECT_RELEASE", proportions_computed = FALSE, degeneracy_assessed = FALSE)
}

prep_map <- function (saved, roster) 
{
    split <- strsplit(saved$cell_key, "|", fixed = TRUE)
    trial_require(all(lengths(split) == 3L), "Malformed original cell key")
    capture <- vapply(split, `[`, "", 2L)
    i <- match(capture, roster$gsm)
    trial_require(!anyNA(i) && all(saved$donor == roster$donor[i]), "Original cell donor/capture mismatch")
    included <- saved$audit_qualified_epithelial_membership == "INCLUDED_FORMAL_AUDIT_PASS"
    target <- saved$provisional_cluster == "cluster_0004"
    map <- data.frame(cell_key = saved$cell_key, sample_id = vapply(split, `[`, "", 1L), technical_library = capture, 
        source_barcode = vapply(split, `[`, "", 3L), donor = saved$donor, tissue_origin = roster$tissue[i], 
        original_source_label = "ORIGINAL_LABEL_NOT_PROVIDED", source_sample_label = roster$source_name[i], 
        qc_pass = TRUE, epithelial_membership = ifelse(included, "INCLUDED", "EXCLUDED"), inclusion_reason = ifelse(included, 
            "SUPPORTED_EPITHELIAL_FORMAL_BROAD_AND_CLUSTER_AUDIT_PASS", ifelse(target, "OWNER_ADJUDICATED_EXCLUSION_BIOLOGICAL_IDENTITY_UNRESOLVED", 
                ifelse(saved$provisional_broad_label == "UNRESOLVED", "UNRESOLVED_BROAD_IDENTITY_EXCLUDED", 
                  "PROVISIONAL_NON_EPITHELIAL_EXCLUDED"))), lineage = ifelse(included, saved$lineage_after_audit, 
            "NOT_IN_EPITHELIAL_UNIVERSE"), author_malignancy_status = saved$author_malignancy_status, 
        author_malignancy_evidence = "UNKNOWN", proliferation_state = saved$proliferation_state, saved[, 
            setdiff(names(saved), c("cell_key", "donor", "epithelial_membership", "author_malignancy_status", 
                "proliferation_state")), drop = FALSE], check.names = FALSE)
    s10b_validate_cell_map(map, map$cell_key)
    trial_require(nrow(map) == 191821L && sum(included) == 33013L && sum(target) == 7281L && all(!included[target]) && 
        setequal(unique(map$provisional_cluster[included]), c("cluster_0006", "cluster_0008", "cluster_0011")) && 
        sum(map$lineage == "other_unassigned") == 30756L && sum(map$lineage == "enteroendocrine") == 
        2257L, "Frozen epithelial identity anchors changed")
    map
}

prep_units_grid <- function (map, roster, endpoint) 
{
    z <- map[map$epithelial_membership == "INCLUDED" & (endpoint == "overall_epithelial" | map$lineage == 
        endpoint), , drop = FALSE]
    u <- roster[, c("gsm", "donor", "tissue")]
    u$unit_id <- paste(u$donor, u$tissue, sep = "|")
    u$tissue[u$tissue == "NAT"] <- "Normal"
    u$cell_count <- as.integer(table(factor(z$technical_library, levels = u$gsm)))
    u$technical_libraries <- u$gsm
    u[order(u$unit_id, method = "radix"), ]
}

prep_capture_aggregate <- function (x, map, features, identity, expected_cells) 
{
    trial_require(inherits(x$counts, "sparseMatrix") && identical(rownames(x$counts), features$feature_id) && 
        identical(x$features, features[, c("feature_id", "source_symbol")]) && all(is.finite(x$counts@x)) && 
        all(x$counts@x >= 0) && all(x$counts@x == round(x$counts@x)), "Original count scale/feature order mismatch")
    trial_require(identical(x$identity$gsm, identity$gsm) && identical(x$identity$donor, identity$donor) && 
        ncol(x$counts) == expected_cells && setequal(colnames(x$counts), map$cell_key), "Count/capture/cell identity mismatch")
    map <- map[match(colnames(x$counts), map$cell_key), , drop = FALSE]
    s10b_validate_cell_map(map, colnames(x$counts))
    original <- Matrix::colSums(x$counts)
    result <- setNames(lapply(c("overall_epithelial", "enteroendocrine"), function(e) s10b_aggregate(x$counts, 
        map, list(QC = "VALIDATED", annotation = "VALIDATED", specification_sha256 = "BOUND_PRE_H4_FREEZE"), 
        e)), c("overall_epithelial", "enteroendocrine"))
    for (e in names(result)) {
        keep <- map$epithelial_membership == "INCLUDED" & (e == "overall_epithelial" | map$lineage == 
            e)
        trial_require(sum(result[[e]]$counts) == sum(original[keep]), "Original cell-to-pseudobulk UMI conservation failed")
    }
    list(aggregates = result, original_UMI = original, cell_totals = data.frame(cell_key = map$cell_key, 
        original_UMI = as.numeric(original)))
}

prep_nom_candidates <- function (route) 
{
    a <- names(route$annotation_qualification$lineage_assignment_panels)
    b <- names(route$annotation_qualification$lineage_audit_panels)
    stopifnot(length(a) > 0L, !anyDuplicated(a), identical(a, b), !"other_unassigned" %in% a)
    a
}

prep_nom_evidence <- function (path) 
{
    columns <- names(read.delim(path, nrows = 0L, check.names = FALSE))
    stopifnot(sum(columns == "panel") == 1L)
    classes <- rep("NULL", length(columns))
    classes[columns == "panel"] <- "character"
    unique(read.delim(path, colClasses = classes, check.names = FALSE)$panel)
}

prep_nom_validate <- function (table, route, assignment = NULL, audit = NULL) 
{
    expected <- prep_nom_candidates(route)
    stopifnot(identical(table$candidate, c(expected, "other_unassigned")), !anyDuplicated(table$candidate), 
        is.logical(table$planned_L_endpoint), !anyNA(table$planned_L_endpoint))
    if (!is.null(assignment)) 
        stopifnot(setequal(assignment, expected))
    if (!is.null(audit)) 
        stopifnot(setequal(audit, expected))
    invisible(TRUE)
}

prep_nom_successor <- function (original, route, assignment, audit) 
{
    stopifnot(sum(original$candidate == "intestinal_absorptive") == 1L)
    successor <- original
    successor$candidate[successor$candidate == "intestinal_absorptive"] <- "intestinal_metaplastic"
    prep_nom_validate(successor, route, assignment, audit)
    stopifnot(identical(successor[, setdiff(names(successor), "candidate"), drop = FALSE], original[, 
        setdiff(names(original), "candidate"), drop = FALSE]))
    successor
}
