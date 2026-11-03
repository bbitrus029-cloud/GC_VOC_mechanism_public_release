# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.
# Initialize explicitly with gcvoc_pdc000614(context); no historical gate runs.
gcvoc_pdc000614 <- function(context) {
  stopifnot(is.function(context$input), is.character(context$root))
  separate_identity <- function() {
    x <- context$read_tsv('inputs/pdc_separate_exclusion.tsv')
    stopifnot(nrow(x)==1L, all(c('sample_id','case_id') %in% names(x)))
    x
  }
  pdc_gate_expected_root <- normalizePath(context$root,mustWork=TRUE)

pdc_gate_assess_authoritative_protein_provenance <- function (root, protein, reconciliation) 
{
    evidence <- pdc_gate_collect_authoritative_protein_evidence(root, protein, reconciliation)
    pdc_gate_validate_authoritative_protein_evidence(evidence)
    evidence$authoritative_post_publication_qc_matrix_verified <- TRUE
    evidence
}

pdc_gate_assess_bridge_endpoint <- function (rna, protein, members, endpoint_id, minimum_quadruplets = 15L) 
{
    if (!all(members %in% rownames(rna$values)) || !all(members %in% rownames(protein$values))) {
        pdc_gate_fail(paste0("Fixed bridge member coverage failed for ", endpoint_id))
    }
    cases <- sort(unique(c(rna$identity$canonical_case_identifier, protein$identity$canonical_case_identifier)), 
        method = "radix")
    rows <- vector("list", length(cases))
    retained <- list()
    retained_index <- 0L
    for (index in seq_along(cases)) {
        case <- cases[[index]]
        r <- rna$identity[rna$identity$canonical_case_identifier == case, , drop = FALSE]
        p <- protein$identity[protein$identity$canonical_case_identifier == case, , drop = FALSE]
        rt <- r[r$canonical_tissue_label == "Tumor", , drop = FALSE]
        rn <- r[r$canonical_tissue_label == "NAT", , drop = FALSE]
        pt <- p[p$canonical_tissue_label == "Tumor", , drop = FALSE]
        pn <- p[p$canonical_tissue_label == "NAT", , drop = FALSE]
        unique_structure <- all(c(nrow(rt), nrow(rn), nrow(pt), nrow(pn)) == 1L)
        excluded <- unique_structure && any(c(rt$excluded, rn$excluded, pt$excluded, pn$excluded))
        finite <- unique_structure && !excluded && all(is.finite(rna$values[members, c(rt$opaque_identifier, 
            rn$opaque_identifier), drop = FALSE])) && all(is.finite(protein$values[members, c(pt$opaque_identifier, 
            pn$opaque_identifier), drop = FALSE]))
        status <- if (!unique_structure) {
            "not_applicable_incomplete_quadruplet"
        }
        else if (excluded) {
            "excluded_whole_quadruplet"
        }
        else if (!finite) {
            "excluded_nonfinite_complete_case"
        }
        else {
            "retained_complete_quadruplet"
        }
        if (identical(status, "retained_complete_quadruplet")) {
            retained_index <- retained_index + 1L
            retained[[retained_index]] <- list(case = case, rna = c(rt$opaque_identifier, rn$opaque_identifier), 
                protein = c(pt$opaque_identifier, pn$opaque_identifier))
        }
        rows[[index]] <- data.frame(patient_id = case, rna_tumor_record_id = if (nrow(rt) == 1L) {
            rt$original_matrix_column
        }
        else "not_available", rna_nat_record_id = if (nrow(rn) == 1L) {
            rn$original_matrix_column
        }
        else "not_available", protein_tumor_record_id = if (nrow(pt) == 1L) {
            pt$original_matrix_column
        }
        else "not_available", protein_nat_record_id = if (nrow(pn) == 1L) {
            pn$original_matrix_column
        }
        else "not_available", endpoint_id = endpoint_id, exclusion_propagation_status = status, complete_case_status = if (identical(status, 
            "retained_complete_quadruplet")) 
            "complete"
        else "incomplete", reason = status, evidence_ids = "EVIDENCE_003;EVIDENCE_026;EVIDENCE_027;EVIDENCE_028", 
            stringsAsFactors = FALSE)
    }
    ledger <- do.call(rbind, rows)
    rownames(ledger) <- NULL
    if (length(retained) < minimum_quadruplets) {
        return(list(endpoint_id = endpoint_id, members = members, ledger = ledger, retained = retained, 
            complete_n = as.integer(length(retained)), coverage_and_complete_case_pass = FALSE, rna_background = character(), 
            protein_background = character(), protein_background_payload = NULL, rna_background_payload = NULL))
    }
    rna_opaque <- unlist(lapply(retained, `[[`, "rna"), use.names = FALSE)
    protein_opaque <- unlist(lapply(retained, `[[`, "protein"), use.names = FALSE)
    rna_background <- rownames(rna$values)[rowSums(is.finite(rna$values[, rna_opaque, drop = FALSE])) == 
        length(rna_opaque)]
    protein_background <- rownames(protein$values)[rowSums(is.finite(protein$values[, protein_opaque, 
        drop = FALSE])) == length(protein_opaque)]
    if (length(rna_background) < 2L || length(protein_background) < 2L || !all(members %in% rna_background) || 
        !all(members %in% protein_background)) {
        pdc_gate_fail("Frozen bridge background is incomplete")
    }
    list(endpoint_id = endpoint_id, members = members, ledger = ledger, retained = retained, complete_n = as.integer(length(retained)), 
        coverage_and_complete_case_pass = TRUE, rna_background = sort(rna_background, method = "radix"), 
        protein_background = sort(protein_background, method = "radix"), protein_background_payload = pdc_gate_background_payload(protein_background), 
        rna_background_payload = pdc_gate_background_payload(rna_background))
}

pdc_gate_assess_pair_endpoint <- function (matrix_input, members, endpoint_id, minimum_pairs = 3L) 
{
    if (!all(members %in% rownames(matrix_input$values))) {
        pdc_gate_fail(paste0("Fixed member coverage failed for ", endpoint_id))
    }
    identity <- matrix_input$identity
    cases <- sort(unique(identity$canonical_case_identifier), method = "radix")
    rows <- vector("list", length(cases))
    retained_cases <- character()
    retained_opaque <- character()
    for (index in seq_along(cases)) {
        case <- cases[[index]]
        z <- identity[identity$canonical_case_identifier == case, , drop = FALSE]
        tumor <- z[z$canonical_tissue_label == "Tumor", , drop = FALSE]
        nat <- z[z$canonical_tissue_label == "NAT", , drop = FALSE]
        unique_structure <- nrow(tumor) == 1L && nrow(nat) == 1L
        excluded <- unique_structure && any(c(tumor$excluded, nat$excluded))
        finite <- unique_structure && !excluded && all(is.finite(matrix_input$values[members, c(tumor$opaque_identifier, 
            nat$opaque_identifier), drop = FALSE]))
        status <- if (!unique_structure) {
            "not_applicable_incomplete_structure"
        }
        else if (excluded) {
            "excluded_whole_pair"
        }
        else if (!finite) {
            "excluded_nonfinite_complete_case"
        }
        else {
            "retained_complete_pair"
        }
        if (identical(status, "retained_complete_pair")) {
            retained_cases <- c(retained_cases, case)
            retained_opaque <- c(retained_opaque, tumor$opaque_identifier, nat$opaque_identifier)
        }
        rows[[index]] <- data.frame(patient_id = case, rna_tumor_record_id = "not_applicable", rna_nat_record_id = "not_applicable", 
            protein_tumor_record_id = if (nrow(tumor) == 1L) {
                tumor$original_matrix_column
            }
            else "not_available", protein_nat_record_id = if (nrow(nat) == 1L) {
                nat$original_matrix_column
            }
            else "not_available", endpoint_id = endpoint_id, exclusion_propagation_status = status, complete_case_status = if (identical(status, 
                "retained_complete_pair")) {
                "complete"
            }
            else "incomplete", reason = status, evidence_ids = "EVIDENCE_003;EVIDENCE_026;EVIDENCE_028", 
            stringsAsFactors = FALSE)
    }
    ledger <- do.call(rbind, rows)
    rownames(ledger) <- NULL
    if (length(retained_cases) < minimum_pairs) {
        pdc_gate_fail(paste0("Fewer than ", minimum_pairs, " complete protein pairs"))
    }
    pdc_gate_validate_paired_design(retained_cases)
    background <- rownames(matrix_input$values)[rowSums(is.finite(matrix_input$values[, retained_opaque, 
        drop = FALSE])) == length(retained_opaque)]
    if (length(background) < 2L || !all(members %in% background)) {
        pdc_gate_fail("Frozen protein rank background is incomplete")
    }
    list(endpoint_id = endpoint_id, members = members, ledger = ledger, retained_cases = retained_cases, 
        retained_opaque = retained_opaque, complete_n = as.integer(length(retained_cases)), background = sort(background, 
            method = "radix"), background_payload = pdc_gate_background_payload(background))
}

pdc_gate_attach_exclusions <- function (matrix_input, exclusion_samples, exclusion_cases) 
{
    matrix_input$identity$excluded <- matrix_input$identity$original_matrix_column %in% exclusion_samples | 
        matrix_input$identity$canonical_case_identifier %in% exclusion_cases
    matrix_input
}

pdc_gate_background_payload <- function (features) 
{
    if (!is.character(features) || !length(features) || anyNA(features) || any(!nzchar(features)) || 
        anyDuplicated(features)) {
        pdc_gate_fail("Background features must be unique nonempty symbols; duplicate detected")
    }
    features <- sort(features, method = "radix")
    preimage <- paste0(paste(features, collapse = "\n"), "\n")
    list(features_json = as.character(jsonlite::toJSON(features, auto_unbox = FALSE, pretty = FALSE)), 
        sha256 = digest::digest(charToRaw(enc2utf8(preimage)), algo = "sha256", serialize = FALSE))
}

pdc_gate_read_clinical_identity <- function (path) 
{
    x <- pdc_gate_read_allowlisted_tabular_identity(path, c("Sample", "Case"))
    category <- x$Sample == "CAT" | x$Case == "CAT"
    if (any(category)) {
        if (sum(category) != 1L || !all(x[category, ] == "CAT")) {
            pdc_gate_fail("LinkedOmics clinical category row is ambiguous")
        }
        x <- x[!category, , drop = FALSE]
    }
    cases <- vapply(x$Case, pdc_gate_canonical_case_id, character(1))
    if (anyDuplicated(cases) || anyNA(cases)) {
        pdc_gate_fail("LinkedOmics clinical case identities are non-unique")
    }
    sort(cases, method = "radix")
}

pdc_gate_read_linkedomics_identity <- function (path) 
{
    required <- c("Sample_ID", "Aliquot_ID", "Case", "tissue_type")
    x <- pdc_gate_read_allowlisted_tabular_identity(path, required)
    category <- x$Sample_ID == "CAT" | x$Aliquot_ID == "CAT" | x$Case == "CAT"
    if (any(category)) {
        if (sum(category) != 1L || !all(x[category, c("Sample_ID", "Aliquot_ID", "Case")] == "CAT")) {
            pdc_gate_fail("LinkedOmics metatable category row is ambiguous")
        }
        x <- x[!category, , drop = FALSE]
    }
    canonical_sample <- vapply(x$Sample_ID, pdc_gate_canonical_sample_id, character(1))
    canonical_case <- vapply(x$Case, pdc_gate_canonical_case_id, character(1))
    canonical_aliquot <- vapply(x$Aliquot_ID, pdc_gate_canonical_aliquot_id, character(1))
    tissue_token <- tolower(trimws(x$tissue_type))
    canonical_tissue <- ifelse(tissue_token %in% c("tumor", "tumour"), "Tumor", ifelse(tissue_token %in% 
        c("nat", "normal"), "NAT", NA_character_))
    if (anyNA(canonical_tissue) || anyDuplicated(canonical_sample) || anyDuplicated(canonical_aliquot) || 
        any(canonical_case != sub("[.][TN]$", "", canonical_sample)) || any(canonical_tissue != ifelse(endsWith(canonical_sample, 
        ".T"), "Tumor", "NAT"))) {
        pdc_gate_fail("LinkedOmics identity mapping is incomplete or conflicting")
    }
    out <- data.frame(original_matrix_column = canonical_sample, authoritative_aliquot_identifier = canonical_aliquot, 
        canonical_case_identifier = canonical_case, canonical_tissue_label = canonical_tissue, stringsAsFactors = FALSE)
    out <- out[order(out$original_matrix_column, method = "radix"), , drop = FALSE]
    rownames(out) <- NULL
    out
}

pdc_gate_read_masked_matrix <- function (path, identity, layer) 
{
    if (!layer %in% c("RNA", "protein")) {
        pdc_gate_fail("Masked matrix layer must be RNA or protein")
    }
    first <- readLines(path, n = 1L, warn = FALSE)
    if (length(first) != 1L) 
        pdc_gate_fail("Matrix header is absent")
    header <- strsplit(first, "\t", fixed = TRUE)[[1L]]
    if (length(header) < 2L || anyNA(header) || any(!nzchar(header)) || anyDuplicated(header)) {
        pdc_gate_fail("Matrix header is invalid")
    }
    sample_columns <- header[-1L]
    if (!setequal(sample_columns, identity$original_matrix_column) && !all(sample_columns %in% identity$original_matrix_column)) {
        pdc_gate_fail("Matrix sample identities are unmatched to allowlisted metadata")
    }
    metadata <- identity[match(sample_columns, identity$original_matrix_column), , drop = FALSE]
    if (anyNA(metadata) || anyDuplicated(paste(metadata$canonical_case_identifier, metadata$canonical_tissue_label, 
        sep = "\r"))) {
        pdc_gate_fail("Matrix identity mapping is missing, duplicate, or many-to-many")
    }
    opaque <- sprintf("opaque_%s_%04d", tolower(layer), seq_along(sample_columns))
    classes <- c("character", rep("numeric", length(sample_columns)))
    x <- utils::read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", 
        "NaN", ""), quote = "", comment.char = "", colClasses = classes)
    if (!identical(names(x), header)) 
        pdc_gate_fail("Matrix header readback differs")
    features <- trimws(as.character(x[[1L]]))
    if (anyNA(features) || any(!nzchar(features)) || anyDuplicated(features)) {
        pdc_gate_fail("Matrix canonical feature identities are non-unique")
    }
    values <- as.matrix(x[, -1L, drop = FALSE])
    if (!is.numeric(values)) 
        pdc_gate_fail("Matrix assay fields are not numeric")
    rownames(values) <- features
    colnames(values) <- opaque
    metadata$opaque_identifier <- opaque
    metadata <- metadata[, c("opaque_identifier", "original_matrix_column", "authoritative_aliquot_identifier", 
        "canonical_case_identifier", "canonical_tissue_label"), drop = FALSE]
    list(layer = layer, values = values, identity = metadata, feature_n = as.integer(nrow(values)), sample_n = as.integer(ncol(values)), 
        sha256 = pdc_gate_sha256(path), bytes = as.integer(file.info(path)$size), file_name = basename(path))
}

pdc_gate_reconcile_official_exclusions <- function (root) 
{
    source_root <- file.path(root, "cache", "pre_amendment_audit", "PDC000614_TMT13_gate_01", "source")
    paths <- list(design = file.path(source_root, "PDC_GraphQL_PDC000614_experimental_design_response_01.json"), 
        biospecimen = file.path(source_root, "PDC_GraphQL_PDC000614_biospecimen_response_01.json"), case_identity = file.path(source_root, 
            "PDC_GraphQL_PDC000614_plex13_case_identity_response_01.json"), table_s1 = file.path(source_root, 
            "Elsevier_S2666379126001734_mmc2_response_01.xlsx"))
    paths <- lapply(paths, context$input)
    channels <- pdc_gate_read_pdc_design_channels(paths$design)
    biological_identity <- pdc_gate_flatten_pdc_case_identity(paths$case_identity)
    reference_identity <- pdc_gate_read_pdc_reference_identity(paths$biospecimen)
    selected_biological <- biological_identity[biological_identity$source_identifier %in% channels$source_identifier, 
        , drop = FALSE]
    biospecimens <- rbind(selected_biological, reference_identity)
    classified <- pdc_gate_classify_pdc_channels(channels, biospecimens)
    pdc_gate_validate_pdc_tmt13_design(classified)
    pdc_identity <- classified[classified$channel_class == "biological_cpt_aliquot", c("canonical_aliquot_id", 
        "canonical_case_identifier", "canonical_tissue_label"), drop = FALSE]
    table_s1_all <- pdc_gate_read_table_s1_identity(paths$table_s1)
    table_identity <- pdc_gate_validate_table_s1_tmt13(table_s1_all)
    pdc_gate_validate_correction03_identity_contract(pdc_identity, table_identity)
    pdc_gate_validate_publication_prose_discrepancy(pdc_gate_contract()$publication_prose_reported_tmt13_exclusions)
    tmt_rows <- table_s1_all[as.character(table_s1_all$TMT18_set) == "13", , drop = FALSE]
    tmt_rows$canonical_sample <- vapply(tmt_rows$Sample_ID, pdc_gate_canonical_sample_id, character(1))
    tmt_rows$canonical_aliquot <- vapply(tmt_rows$Aliquot_ID, pdc_gate_canonical_aliquot_id, character(1))
    pdc_index <- match(tmt_rows$canonical_aliquot, classified$canonical_aliquot_id)
    if (anyNA(pdc_index)) 
        pdc_gate_fail("Table S1 aliquot is unmatched in PDC")
    tmt_rows$canonical_case <- vapply(tmt_rows$Case, pdc_gate_canonical_case_id, character(1))
    tmt_rows$canonical_tissue <- ifelse(endsWith(tmt_rows$canonical_sample, ".T"), "Tumor", "NAT")
    order_index <- order(tmt_rows$canonical_sample, method = "radix")
    tmt_rows <- tmt_rows[order_index, , drop = FALSE]
    pdc_index <- pdc_index[order_index]
    separate_index <- which(trimws(table_s1_all$Sample_ID) == separate_identity()$sample_id)
    if (length(separate_index) != 1L || as.character(table_s1_all$TMT18_set[[separate_index]]) == "13") {
        pdc_gate_fail("Separate unresolved-tissue exclusion is not uniquely resolved")
    }
    separate <- table_s1_all[separate_index, , drop = FALSE]
    exclusion_input <- data.frame(authoritative_sample_identifier = c(tmt_rows$canonical_sample, separate_identity()$sample_id), 
        canonical_case_identifier = c(tmt_rows$canonical_case, separate_identity()$case_id), authoritative_aliquot_identifier = c(tmt_rows$canonical_aliquot, 
            pdc_gate_canonical_aliquot_id(separate$Aliquot_ID[[1L]])), canonical_tissue_label = c(tmt_rows$canonical_tissue, 
            "NAT"), plex_identifier = c(rep("13CPTAC_GC_Proteome_JHU_Ascend", 16L), paste0("TMT18_set_", 
            trimws(separate$TMT18_set[[1L]]))), channel_identifier = c(classified$channel_identifier[pdc_index], 
            trimws(separate$TMT18_channel[[1L]])), exclusion_reason = c(rep("structured_official_TMT13_exclusion", 
            16L), "unresolved_tissue_assignment"), evidence_ids = c(rep("EVIDENCE_020;EVIDENCE_021;EVIDENCE_023;EVIDENCE_026", 
            16L), "EVIDENCE_001;EVIDENCE_026"), stringsAsFactors = FALSE)
    pdc_gate_validate_publication_qc_exclusions(exclusion_input)
    ledger <- data.frame(record_id = c(sprintf("TMT13_%02d", seq_len(16L)), "SEPARATE_01"), authoritative_sample_identifier = exclusion_input$authoritative_sample_identifier, 
        canonical_case_identifier = exclusion_input$canonical_case_identifier, canonical_tissue_label = exclusion_input$canonical_tissue_label, 
        plex_identifier = exclusion_input$plex_identifier, channel_identifier = exclusion_input$channel_identifier, 
        exclusion_status = "excluded", exclusion_reason = exclusion_input$exclusion_reason, evidence_ids = exclusion_input$evidence_ids, 
        identity_status = c(rep("resolved_exact_PDC_TableS1_aliquot_and_case_tissue_match", 16L), "resolved_separate_unresolved_tissue_assignment_exclusion"), 
        stringsAsFactors = FALSE)
    list(ledger = ledger, exclusions = exclusion_input, hashes = pdc_gate_identity_set_hashes(table_identity), 
        pdc_design = classified, table_s1_identity = table_identity)
}

pdc_gate_structural_pair_count <- function (identity) 
{
    cases <- unique(identity$canonical_case_identifier)
    sum(vapply(cases, function(case) {
        z <- identity[identity$canonical_case_identifier == case, , drop = FALSE]
        identical(sort(z$canonical_tissue_label, method = "radix"), c("NAT", "Tumor")) && !any(z$excluded)
    }, logical(1)))
}

pdc_gate_structural_quadruplet_count <- function (rna_identity, protein_identity) 
{
    cases <- unique(c(rna_identity$canonical_case_identifier, protein_identity$canonical_case_identifier))
    sum(vapply(cases, function(case) {
        r <- rna_identity[rna_identity$canonical_case_identifier == case, , drop = FALSE]
        p <- protein_identity[protein_identity$canonical_case_identifier == case, , drop = FALSE]
        identical(sort(r$canonical_tissue_label, method = "radix"), c("NAT", "Tumor")) && identical(sort(p$canonical_tissue_label, 
            method = "radix"), c("NAT", "Tumor")) && !any(c(r$excluded, p$excluded))
    }, logical(1)))
}

pdc_gate_build_complete_pairs <- function (records, minimum_pairs = 3L) 
{
    required <- c("patient_id", "sample_id", "tissue", "excluded", "finite")
    if (!is.data.frame(records) || !all(required %in% names(records)) || anyNA(records[required])) {
        pdc_gate_fail("Protein pair records are incomplete")
    }
    if (!is.logical(records$excluded) || !is.logical(records$finite)) {
        pdc_gate_fail("Protein pair exclusion and finite flags must be logical")
    }
    patients <- sort(unique(records$patient_id), method = "radix")
    rows <- lapply(patients, function(patient) {
        z <- records[records$patient_id == patient, , drop = FALSE]
        tissue_ok <- identical(sort(z$tissue, method = "radix"), c("NAT", "Tumor"))
        complete <- nrow(z) == 2L && tissue_ok && !any(z$excluded) && all(z$finite) && !anyDuplicated(z$sample_id)
        if (!complete) 
            return(NULL)
        data.frame(patient_id = patient, tumor_sample_id = z$sample_id[z$tissue == "Tumor"], nat_sample_id = z$sample_id[z$tissue == 
            "NAT"], stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    out <- if (length(rows)) 
        do.call(rbind, rows)
    else data.frame(patient_id = character(), tumor_sample_id = character(), nat_sample_id = character(), 
        stringsAsFactors = FALSE)
    rownames(out) <- NULL
    if (nrow(out) < as.integer(minimum_pairs)) {
        pdc_gate_fail(paste0("Fewer than ", as.integer(minimum_pairs), " complete protein pairs"))
    }
    out
}

pdc_gate_complete_quadruplets <- function (x, minimum_n = 15L) 
{
    required <- c("patient_id", "rna_tumor", "rna_nat", "protein_tumor", "protein_nat", "excluded", "ambiguous", 
        "finite")
    if (!is.data.frame(x) || !all(required %in% names(x))) {
        pdc_gate_fail("Bridge quadruplet input is incomplete")
    }
    record_fields <- c("rna_tumor", "rna_nat", "protein_tumor", "protein_nat")
    present <- stats::complete.cases(x[, record_fields, drop = FALSE]) & apply(x[, record_fields, drop = FALSE], 
        1L, function(z) all(nzchar(z)))
    keep <- present & !x$excluded & !x$ambiguous & x$finite
    out <- x[keep, , drop = FALSE]
    out <- out[order(out$patient_id, method = "radix"), , drop = FALSE]
    rownames(out) <- NULL
    if (nrow(out) < as.integer(minimum_n)) {
        pdc_gate_fail("Fewer than 15 bridge quadruplets remain complete")
    }
    out
}

pdc_gate_collect_authoritative_protein_evidence <- function (root, protein, reconciliation) 
{
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    bridge_source <- file.path(root, "cache", "pre_amendment_audit", "CPTAC_same_patient_bridge_gate_01", 
        "source")
    gate_source <- file.path(root, "cache", "pre_amendment_audit", "PDC000614_TMT13_gate_01", "source")
    publication_path <- file.path(bridge_source, "PMC13198309.html")
    catalog_path <- file.path(bridge_source, "CPTAC-GC_index.html")
    table_s1_path <- file.path(gate_source, "Elsevier_S2666379126001734_mmc2_response_01.xlsx")
    protein_path <- file.path(bridge_source, protein$file_name)
    publication_path <- context$input(publication_path)
    catalog_path <- context$input(catalog_path)
    table_s1_path <- context$input(table_s1_path)
    protein_path <- context$input(protein_path)
    paths <- c(publication_path, catalog_path, table_s1_path, protein_path)
    if (!all(file.exists(paths)) || any(dir.exists(paths))) {
        pdc_gate_fail("Authoritative protein-provenance source is absent")
    }
    publication <- readLines(publication_path, warn = FALSE, encoding = "UTF-8")
    publication_processed_matrices_linked <- any(grepl(paste("Processed data matrices utilized for this publication can be", 
        "accessed via LinkedOmics:"), publication, fixed = TRUE))
    publication_tmt_reprocessed_without_tmt13 <- any(grepl(paste("For TMT-labeled datasets (proteomics, phosphoproteomics,", 
        "glycoproteomics, and acetylomics), data were reprocessed without TMT 13."), publication, fixed = TRUE))
    publication_analyses_use_linkedomics_processed_matrices <- any(grepl(paste("All analyses were performed using the processed data matrices described", 
        "above and deposited on LinkedOmics, unless otherwise specified"), publication, fixed = TRUE))
    catalog <- readLines(catalog_path, warn = FALSE, encoding = "UTF-8")
    protein_file_name <- "GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor_NAT.cct"
    link_token <- paste0("href=\"./", protein_file_name, "\"")
    link_index <- grep(link_token, catalog, fixed = TRUE)
    if (length(link_index) != 1L) {
        pdc_gate_fail("Official catalog protein-file link is absent or ambiguous")
    }
    catalog_row_index <- seq.int(max(1L, link_index - 1L), min(length(catalog), link_index + 7L))
    catalog_row <- paste(catalog[catalog_row_index], collapse = " ")
    catalog_display_feature_n <- if (grepl("<td> 11293</td>", catalog_row, fixed = TRUE)) 
        11293L
    else NA_integer_
    catalog_sample_n <- if (grepl("<td>189</td>", catalog_row, fixed = TRUE)) 
        189L
    else NA_integer_
    catalog_scale_verified <- grepl("Proteome Gene TMT (median centered log2 ratio, tumor and normal)", 
        catalog_row, fixed = TRUE) && grepl("Expression (median centerted log2 ratio)", catalog_row, 
        fixed = TRUE)
    table_s1 <- pdc_gate_read_table_s1_identity(table_s1_path)
    passed <- toupper(trimws(table_s1$passed_QC)) == "TRUE"
    tissue <- tolower(trimws(table_s1$tissue_type))
    passed_biological <- passed & tissue %in% c("tumor", "nat")
    required_identity <- c("Sample_ID", "Aliquot_ID", "Case", "tissue_type")
    if (anyNA(table_s1[passed_biological, required_identity, drop = FALSE]) || any(!nzchar(as.matrix(table_s1[passed_biological, 
        required_identity, drop = FALSE])))) {
        pdc_gate_fail("Official Table S1 passed-QC identity is incomplete")
    }
    table_canonical_sample_ids <- vapply(table_s1$Sample_ID[passed_biological], pdc_gate_canonical_sample_id, 
        character(1))
    active_canonical_sample_ids <- vapply(protein$identity$original_matrix_column, pdc_gate_canonical_sample_id, 
        character(1))
    failed_ids <- reconciliation$exclusions$authoritative_sample_identifier
    if (length(failed_ids) != 17L || anyDuplicated(failed_ids)) {
        pdc_gate_fail("Authoritative failed-QC identity set is invalid")
    }
    list(publication_sha256 = pdc_gate_sha256(publication_path), catalog_sha256 = pdc_gate_sha256(catalog_path), 
        table_s1_sha256 = pdc_gate_sha256(table_s1_path), publication_processed_matrices_linked = publication_processed_matrices_linked, 
        publication_tmt_reprocessed_without_tmt13 = publication_tmt_reprocessed_without_tmt13, publication_analyses_use_linkedomics_processed_matrices = publication_analyses_use_linkedomics_processed_matrices, 
        catalog_protein_href = paste0("./", protein_file_name), catalog_display_feature_n = catalog_display_feature_n, 
        catalog_sample_n = catalog_sample_n, catalog_scale_verified = catalog_scale_verified, protein_file_name = protein$file_name, 
        protein_official_url = paste0("https://linkedomics.org/data_download/CPTAC-GC/", protein$file_name), 
        protein_bytes = as.integer(file.info(protein_path)$size), protein_sha256 = pdc_gate_sha256(protein_path), 
        table_canonical_sample_ids = table_canonical_sample_ids, active_canonical_sample_ids = active_canonical_sample_ids, 
        table_sample_n = as.integer(length(table_canonical_sample_ids)), table_tumor_n = as.integer(sum(tissue[passed_biological] == 
            "tumor")), table_nat_n = as.integer(sum(tissue[passed_biological] == "nat")), active_sample_n = as.integer(length(active_canonical_sample_ids)), 
        active_tumor_n = as.integer(sum(protein$identity$canonical_tissue_label == "Tumor")), active_nat_n = as.integer(sum(protein$identity$canonical_tissue_label == 
            "NAT")), table_sample_set_sha256 = pdc_gate_lf_set_sha256(table_canonical_sample_ids), active_sample_set_sha256 = pdc_gate_lf_set_sha256(active_canonical_sample_ids), 
        sample_set_sha256 = pdc_gate_lf_set_sha256(table_canonical_sample_ids), sample_sets_equal = identical(sort(table_canonical_sample_ids, 
            method = "radix"), sort(active_canonical_sample_ids, method = "radix")), failed_qc_record_n = as.integer(length(failed_ids)), 
        failed_qc_records_present_n = as.integer(sum(failed_ids %in% active_canonical_sample_ids)), failed_qc_records_absent_n = as.integer(sum(!failed_ids %in% 
            active_canonical_sample_ids)), physical_feature_n = as.integer(protein$feature_n), physical_unique_feature_n = as.integer(length(unique(rownames(protein$values)))), 
        physical_features_with_any_finite_n = as.integer(sum(rowSums(is.finite(protein$values)) > 0L)), 
        semantic_version_supplied = FALSE, second_reference_matrix_required = FALSE)
}

pdc_gate_validate_authoritative_protein_evidence <- function (evidence) 
{
    required <- c("publication_sha256", "catalog_sha256", "table_s1_sha256", "publication_processed_matrices_linked", 
        "publication_tmt_reprocessed_without_tmt13", "publication_analyses_use_linkedomics_processed_matrices", 
        "catalog_protein_href", "catalog_display_feature_n", "catalog_sample_n", "catalog_scale_verified", 
        "protein_file_name", "protein_official_url", "protein_bytes", "protein_sha256", "table_canonical_sample_ids", 
        "active_canonical_sample_ids", "table_sample_n", "table_tumor_n", "table_nat_n", "active_sample_n", 
        "active_tumor_n", "active_nat_n", "table_sample_set_sha256", "active_sample_set_sha256", "sample_set_sha256", 
        "sample_sets_equal", "failed_qc_record_n", "failed_qc_records_present_n", "failed_qc_records_absent_n", 
        "physical_feature_n", "physical_unique_feature_n", "physical_features_with_any_finite_n", "semantic_version_supplied", 
        "second_reference_matrix_required")
    if (!is.list(evidence) || !identical(names(evidence), required)) {
        pdc_gate_fail("Authoritative protein-provenance evidence schema is invalid")
    }
    if (!identical(evidence$publication_sha256, "3ba5e42f0bb3d02e780db4e4ee6f74e41a39771f53339cf19e43af187b60665d") || 
        !isTRUE(evidence$publication_processed_matrices_linked) || !isTRUE(evidence$publication_tmt_reprocessed_without_tmt13) || 
        !isTRUE(evidence$publication_analyses_use_linkedomics_processed_matrices)) {
        pdc_gate_fail("Primary-publication processed-matrix linkage is invalid")
    }
    expected_file <- "GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor_NAT.cct"
    if (!identical(evidence$catalog_sha256, "9deea7e50af4a486430ed514356261b2e3893986881f4e128b0fd759a69960ed") || 
        !identical(evidence$catalog_protein_href, paste0("./", expected_file)) || !identical(evidence$catalog_display_feature_n, 
        11293L) || !identical(evidence$catalog_sample_n, 189L) || !isTRUE(evidence$catalog_scale_verified)) {
        pdc_gate_fail("Official catalog linkage or exact 11,293 display metadata is invalid")
    }
    if (!identical(evidence$table_s1_sha256, "4ae1d97ea02d8691698f6782781b508f4b988bb5e016ece51991f733afacf901")) {
        pdc_gate_fail("Official Table S1 identity snapshot is invalid")
    }
    if (!identical(evidence$protein_file_name, expected_file) || !identical(evidence$protein_official_url, 
        paste0("https://linkedomics.org/data_download/CPTAC-GC/", expected_file)) || !identical(evidence$protein_bytes, 
        35726655L) || !identical(evidence$protein_sha256, "f826ecbc37646e0a012b0a0ea9855b0b791fc456f0c919274a04098ef1456700")) {
        pdc_gate_fail("Exact official protein file identity is invalid")
    }
    table_ids <- evidence$table_canonical_sample_ids
    active_ids <- evidence$active_canonical_sample_ids
    if (!is.character(table_ids) || !is.character(active_ids) || anyNA(table_ids) || anyNA(active_ids) || 
        any(!nzchar(table_ids)) || any(!nzchar(active_ids)) || anyDuplicated(table_ids) || anyDuplicated(active_ids)) {
        pdc_gate_fail("Official or active sample identity is duplicate or invalid")
    }
    if (!identical(evidence$table_sample_n, 189L) || !identical(evidence$active_sample_n, 189L) || length(table_ids) != 
        189L || length(active_ids) != 189L || !identical(evidence$table_tumor_n, 159L) || !identical(evidence$table_nat_n, 
        30L) || !identical(evidence$active_tumor_n, 159L) || !identical(evidence$active_nat_n, 30L)) {
        pdc_gate_fail("Official and active 189-sample 159/30 sample universe is invalid")
    }
    if (!isTRUE(evidence$sample_sets_equal) || !identical(sort(table_ids, method = "radix"), sort(active_ids, 
        method = "radix"))) {
        pdc_gate_fail("Official passed-QC and active protein sample sets are unequal")
    }
    expected_sample_hash <- "d1101c667e517734908465cf6fb6e3f8690b9283c1dabcc46e0f0688fc2fbea7"
    observed_table_hash <- pdc_gate_lf_set_sha256(table_ids)
    observed_active_hash <- pdc_gate_lf_set_sha256(active_ids)
    if (!identical(evidence$table_sample_set_sha256, expected_sample_hash) || !identical(evidence$active_sample_set_sha256, 
        expected_sample_hash) || !identical(evidence$sample_set_sha256, expected_sample_hash) || !identical(observed_table_hash, 
        expected_sample_hash) || !identical(observed_active_hash, expected_sample_hash)) {
        pdc_gate_fail("Official and active sample-set SHA-256 is invalid")
    }
    if (!identical(evidence$failed_qc_record_n, 17L) || !identical(evidence$failed_qc_records_present_n, 
        0L) || !identical(evidence$failed_qc_records_absent_n, 17L)) {
        pdc_gate_fail("Authoritative failed-QC exclusion absence is invalid")
    }
    if (!identical(evidence$physical_feature_n, 11294L) || !identical(evidence$physical_unique_feature_n, 
        11294L) || !identical(evidence$physical_features_with_any_finite_n, 11294L)) {
        pdc_gate_fail("Exact 11,294 physical feature identity is invalid")
    }
    if (!identical(evidence$semantic_version_supplied, FALSE) || !identical(evidence$second_reference_matrix_required, 
        FALSE)) {
        pdc_gate_fail("Authoritative route must not require a semantic version or second matrix")
    }
    TRUE
}

pdc_gate_fail <- function (message) 
{
    stop(message, call. = FALSE)
}

pdc_gate_validate_paired_design <- function (patient_ids) 
{
    if (!is.character(patient_ids) || anyNA(patient_ids) || any(!nzchar(patient_ids)) || anyDuplicated(patient_ids)) {
        pdc_gate_fail("Paired design is not full-rank because patient IDs duplicate")
    }
    if (length(patient_ids) < 3L) {
        pdc_gate_fail("At least three complete protein pairs are required")
    }
    patient <- factor(rep(patient_ids, each = 2L), levels = patient_ids)
    tissue <- factor(rep(c("NAT", "Tumor"), length(patient_ids)), levels = c("NAT", "Tumor"))
    design <- stats::model.matrix(~0 + patient + tissue)
    if (qr(design)$rank != ncol(design)) {
        pdc_gate_fail("Paired design is not full-rank")
    }
    TRUE
}

pdc_gate_canonical_case_id <- function (value) 
{
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        pdc_gate_fail("Case identifier must be one non-missing string")
    }
    value <- trimws(value)
    match <- regexec("^(C3[LN])[._-]([0-9]{5})$", value, perl = TRUE)
    pieces <- regmatches(value, match)[[1L]]
    if (length(pieces) != 3L) {
        pdc_gate_fail(paste0("Invalid case identifier: ", value))
    }
    paste(pieces[2L], pieces[3L], sep = ".")
}

pdc_gate_read_allowlisted_tabular_identity <- function (path, required) 
{
    if (!file.exists(path) || dir.exists(path)) {
        pdc_gate_fail("Allowlisted identity table is absent")
    }
    first <- readLines(path, n = 1L, warn = FALSE)
    if (length(first) != 1L) 
        pdc_gate_fail("Identity-table header is absent")
    header <- strsplit(first, "\t", fixed = TRUE)[[1L]]
    if (anyNA(header) || any(!nzchar(header)) || !all(required %in% header) || any(vapply(required, function(value) sum(header == 
        value), integer(1)) != 1L)) {
        pdc_gate_fail("Identity-table allowlist fields are incomplete or ambiguous")
    }
    classes <- ifelse(header %in% required, "character", "NULL")
    x <- utils::read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character(), 
        quote = "", comment.char = "", colClasses = classes)
    x <- x[, required, drop = FALSE]
    if (!identical(names(x), required)) {
        pdc_gate_fail("Identity-table field-selective readback is invalid")
    }
    x
}

pdc_gate_canonical_aliquot_id <- function (value) 
{
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        pdc_gate_fail("PDC aliquot identifier must be one non-missing string")
    }
    value <- trimws(value)
    compact <- "^CPT[0-9]{10}$"
    spaced <- "^(CPT[0-9]{6}) ([0-9]{4})$"
    if (grepl(compact, value, perl = TRUE)) {
        return(value)
    }
    match <- regexec(spaced, value, perl = TRUE)
    pieces <- regmatches(value, match)[[1L]]
    if (length(pieces) != 3L) {
        pdc_gate_fail(paste0("Invalid PDC aliquot identifier: ", value))
    }
    paste0(pieces[2L], pieces[3L])
}

pdc_gate_canonical_sample_id <- function (value) 
{
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        pdc_gate_fail("Sample identifier must be one non-missing string")
    }
    value <- trimws(value)
    pattern <- "^(C3[LN])[._-]([0-9]{5})[._-]([TN])$"
    match <- regexec(pattern, value, perl = TRUE)
    pieces <- regmatches(value, match)[[1L]]
    if (length(pieces) != 4L) {
        pdc_gate_fail(paste0("Invalid sample identifier: ", value))
    }
    paste(pieces[2L], pieces[3L], pieces[4L], sep = ".")
}

pdc_gate_sha256 <- function (path) 
{
    if (!file.exists(path) || dir.exists(path)) {
        pdc_gate_fail(paste0("Required regular file is absent: ", path))
    }
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

pdc_gate_classify_pdc_channels <- function (channels, biospecimens) 
{
    channel_fields <- c("channel_identifier", "source_identifier")
    biospecimen_fields <- c("source_identifier", "aliquot_is_ref", "canonical_case_identifier", "canonical_tissue_label")
    if (!is.data.frame(channels) || !all(channel_fields %in% names(channels)) || !is.data.frame(biospecimens) || 
        !all(biospecimen_fields %in% names(biospecimens)) || !nrow(channels) || anyNA(channels[channel_fields]) || 
        any(!nzchar(trimws(as.matrix(channels[channel_fields]))))) {
        pdc_gate_fail("PDC channel or biospecimen identity input is incomplete")
    }
    if (anyDuplicated(channels$channel_identifier) || anyDuplicated(channels$source_identifier) || anyDuplicated(biospecimens$source_identifier)) {
        pdc_gate_fail("PDC channel identity contains a duplicate or many-to-many match")
    }
    match_index <- match(channels$source_identifier, biospecimens$source_identifier)
    matched <- !is.na(match_index)
    reference_token <- rep(NA_character_, nrow(channels))
    reference_token[matched] <- as.character(biospecimens$aliquot_is_ref[match_index[matched]])
    valid_reference_token <- is.na(reference_token) | reference_token %in% c("Yes", "No")
    if (!all(valid_reference_token)) {
        pdc_gate_fail("PDC aliquot_is_ref contains an invalid Boolean token")
    }
    is_cpt <- grepl("^CPT[0-9]{10}$", trimws(channels$source_identifier), perl = TRUE)
    channel_class <- rep("invalid_or_unresolved", nrow(channels))
    channel_class[matched & reference_token == "Yes"] <- "official_reference"
    channel_class[matched & reference_token == "No" & is_cpt] <- "biological_cpt_aliquot"
    canonical_aliquot <- rep("not_available", nrow(channels))
    biological <- channel_class == "biological_cpt_aliquot"
    canonical_aliquot[biological] <- vapply(channels$source_identifier[biological], pdc_gate_canonical_aliquot_id, 
        character(1))
    canonical_aliquot[channel_class == "official_reference"] <- "not_applicable"
    canonical_case <- rep("not_available", nrow(channels))
    canonical_tissue <- rep("not_available", nrow(channels))
    canonical_case[biological] <- as.character(biospecimens$canonical_case_identifier[match_index[biological]])
    canonical_tissue[biological] <- as.character(biospecimens$canonical_tissue_label[match_index[biological]])
    canonical_case[channel_class == "official_reference"] <- "not_applicable"
    canonical_tissue[channel_class == "official_reference"] <- "not_applicable"
    data.frame(channel_identifier = as.character(channels$channel_identifier), source_identifier = as.character(channels$source_identifier), 
        channel_class = channel_class, canonical_aliquot_id = canonical_aliquot, canonical_case_identifier = canonical_case, 
        canonical_tissue_label = canonical_tissue, stringsAsFactors = FALSE)
}

pdc_gate_contract <- function (root = pdc_gate_expected_root) 
{
    pdc_gate_project_root(root)
    yaml::read_yaml(context$input("config/pdc000614_tmt13_gate_contract_v1.0.yml"))
}

pdc_gate_flatten_pdc_case_identity <- function (path) 
{
    parsed <- jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
    cases <- parsed$data
    if (!is.list(cases) || !length(cases)) {
        pdc_gate_fail("PDC case identity response is invalid")
    }
    rows <- list()
    index <- 0L
    for (case in cases) {
        if (!is.data.frame(case) || nrow(case) != 1L || !all(c("case_submitter_id", "samples") %in% names(case))) {
            pdc_gate_fail("PDC case identity record is malformed")
        }
        sample_rows <- case$samples[[1L]]
        if (!is.data.frame(sample_rows) || !all(c("sample_submitter_id", "tissue_type", "aliquots") %in% 
            names(sample_rows))) {
            pdc_gate_fail("PDC sample identity record is malformed")
        }
        for (sample_index in seq_len(nrow(sample_rows))) {
            aliquots <- sample_rows$aliquots[[sample_index]]
            if (!is.data.frame(aliquots) || !all(c("aliquot_submitter_id", "aliquot_is_ref") %in% names(aliquots))) {
                pdc_gate_fail("PDC aliquot identity record is malformed")
            }
            for (aliquot_index in seq_len(nrow(aliquots))) {
                index <- index + 1L
                rows[[index]] <- data.frame(source_identifier = trimws(aliquots$aliquot_submitter_id[[aliquot_index]]), 
                  aliquot_is_ref = as.character(aliquots$aliquot_is_ref[[aliquot_index]]), canonical_case_identifier = pdc_gate_canonical_case_id(case$case_submitter_id[[1L]]), 
                  canonical_tissue_label = switch(trimws(sample_rows$tissue_type[[sample_index]]), Tumor = "Tumor", 
                    Normal = "NAT", pdc_gate_fail("PDC tissue identity is unsupported")), stringsAsFactors = FALSE)
            }
        }
    }
    out <- do.call(rbind, rows)
    if (anyNA(out) || any(!nzchar(as.matrix(out))) || anyDuplicated(out$source_identifier)) {
        pdc_gate_fail("PDC biological aliquot identity is duplicate or ambiguous")
    }
    out
}

pdc_gate_identity_set_hashes <- function (x) 
{
    required <- c("canonical_aliquot_id", "canonical_case_identifier", "canonical_tissue_label")
    if (!is.data.frame(x) || !all(required %in% names(x)) || nrow(x) != 16L || anyNA(x[required]) || 
        any(!nzchar(trimws(as.matrix(x[required])))) || anyDuplicated(x$canonical_aliquot_id) || anyDuplicated(paste(x$canonical_case_identifier, 
        x$canonical_tissue_label, sep = "\r"))) {
        pdc_gate_fail("Identity-set hash input is incomplete, duplicate, or ambiguous")
    }
    list(canonical_aliquot_set_sha256 = pdc_gate_lf_set_sha256(x$canonical_aliquot_id), canonical_case_plus_tissue_set_sha256 = pdc_gate_lf_set_sha256(paste(x$canonical_case_identifier, 
        x$canonical_tissue_label, sep = "|")))
}

pdc_gate_read_pdc_design_channels <- function (path) 
{
    parsed <- jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
    design <- parsed$data$uiStudyExperimentalDesign
    if (!is.data.frame(design)) 
        pdc_gate_fail("PDC design response is invalid")
    selected <- design[design$plex_dataset_name == "13CPTAC_GC_Proteome_JHU_Ascend", , drop = FALSE]
    if (nrow(selected) != 1L || selected$experiment_number[[1L]] != 2L) {
        pdc_gate_fail("Exact PDC TMT13 plex or experiment identity is unresolved")
    }
    channel_fields <- grep("^tmt_", names(selected), value = TRUE)
    if (length(channel_fields) != 18L) {
        pdc_gate_fail("PDC TMT13 design must contain exactly 18 channels")
    }
    rows <- lapply(channel_fields, function(field) {
        cell <- selected[[field]][[1L]]
        if (!is.data.frame(cell) || nrow(cell) != 1L || !"aliquot_submitter_id" %in% names(cell) || is.na(cell$aliquot_submitter_id[[1L]]) || 
            !nzchar(trimws(cell$aliquot_submitter_id[[1L]]))) {
            pdc_gate_fail("PDC channel identity is absent or ambiguous")
        }
        data.frame(channel_identifier = field, source_identifier = trimws(cell$aliquot_submitter_id[[1L]]), 
            stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    if (anyDuplicated(out$channel_identifier) || anyDuplicated(out$source_identifier)) {
        pdc_gate_fail("PDC TMT13 channel mapping is duplicate or many-to-many")
    }
    out
}

pdc_gate_read_pdc_reference_identity <- function (path) 
{
    parsed <- jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
    x <- parsed$data$uiBiospecimenPerStudy
    required <- c("aliquot_submitter_id", "aliquot_is_ref")
    if (!is.data.frame(x) || !all(required %in% names(x))) {
        pdc_gate_fail("PDC biospecimen response is invalid")
    }
    x <- x[x$aliquot_submitter_id %in% c("NCI7-9", "Ref-13"), required, drop = FALSE]
    if (nrow(x) != 2L || anyDuplicated(x$aliquot_submitter_id) || any(x$aliquot_is_ref != "Yes")) {
        pdc_gate_fail("PDC reference-channel status is unresolved")
    }
    data.frame(source_identifier = x$aliquot_submitter_id, aliquot_is_ref = x$aliquot_is_ref, canonical_case_identifier = "not_applicable", 
        canonical_tissue_label = "not_applicable", stringsAsFactors = FALSE)
}

pdc_gate_read_table_s1_identity <- function (path) 
{
    required <- c("Sample_ID", "Aliquot_ID", "Case", "tissue_type", "passed_QC", "TMT18_set", "TMT18_set_position", 
        "TMT18_channel")
    if (!file.exists(path) || dir.exists(path)) {
        pdc_gate_fail("Official Table S1 snapshot is absent")
    }
    header <- names(readxl::read_excel(path, sheet = "Derived_Features", n_max = 0L))
    if (anyDuplicated(header) || !all(required %in% header)) {
        pdc_gate_fail("Official Table S1 identity/QC schema is incomplete")
    }
    column_types <- ifelse(header %in% required, "text", "skip")
    x <- readxl::read_excel(path, sheet = "Derived_Features", col_types = column_types)
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(x), required)) {
        pdc_gate_fail("Table S1 field-selective parser returned an unsafe schema")
    }
    x
}

pdc_gate_validate_correction03_identity_contract <- function (pdc, table_s1, expected_aliquot_set_sha256 = pdc_gate_contract()$structured_official_tmt13_exclusions$canonical_aliquot_set_sha256, 
    expected_case_tissue_set_sha256 = pdc_gate_contract()$structured_official_tmt13_exclusions$canonical_case_plus_tissue_set_sha256, 
    identity_sources = unname(unlist(pdc_gate_contract()$structured_official_tmt13_exclusions$identity_sources, 
        use.names = FALSE))) 
{
    required_sources <- c("official_primary_supplement_table_s1", "official_pdc_plex_and_biospecimen_metadata")
    if (!is.character(identity_sources) || anyNA(identity_sources) || !identical(identity_sources, required_sources)) {
        pdc_gate_fail("Correction 03 requires exactly two structured official identity sources")
    }
    pdc_gate_assert_identity_sets_equal(pdc, table_s1)
    for (x in list(pdc, table_s1)) {
        pdc_gate_assert_structured_count_contract(sum(x$canonical_tissue_label == "Tumor"), sum(x$canonical_tissue_label == 
            "NAT"))
    }
    pdc_hashes <- pdc_gate_identity_set_hashes(pdc)
    table_hashes <- pdc_gate_identity_set_hashes(table_s1)
    if (!identical(pdc_hashes, table_hashes)) {
        pdc_gate_fail("Structured official identity-set hashes disagree")
    }
    if (!identical(pdc_hashes$canonical_aliquot_set_sha256, expected_aliquot_set_sha256)) {
        pdc_gate_fail("Correction 03 canonical aliquot-set hash mismatch")
    }
    if (!identical(pdc_hashes$canonical_case_plus_tissue_set_sha256, expected_case_tissue_set_sha256)) {
        pdc_gate_fail("Correction 03 canonical case-plus-tissue-set hash mismatch")
    }
    TRUE
}

pdc_gate_validate_pdc_tmt13_design <- function (classified) 
{
    required <- c("channel_identifier", "source_identifier", "channel_class", "canonical_aliquot_id", 
        "canonical_case_identifier", "canonical_tissue_label")
    if (!is.data.frame(classified) || !all(required %in% names(classified)) || anyNA(classified[required])) {
        pdc_gate_fail("Classified PDC TMT13 design is incomplete")
    }
    biological <- classified$channel_class == "biological_cpt_aliquot"
    reference <- classified$channel_class == "official_reference"
    if (any(classified$channel_class == "invalid_or_unresolved") || sum(reference) != 2L || !identical(sort(classified$source_identifier[reference], 
        method = "radix"), c("NCI7-9", "Ref-13"))) {
        pdc_gate_fail("PDC TMT13 contains an invalid, unresolved, or misclassified reference channel")
    }
    if (sum(biological) != 16L) {
        pdc_gate_fail("PDC TMT13 must contain exactly 16 biological CPT aliquots")
    }
    aliquots <- classified$canonical_aliquot_id[biological]
    cases <- classified$canonical_case_identifier[biological]
    tissues <- classified$canonical_tissue_label[biological]
    if (anyDuplicated(aliquots) || anyDuplicated(paste(cases, tissues, sep = "\r"))) {
        pdc_gate_fail("PDC biological identities contain a duplicate or many-to-many match")
    }
    if (any(!tissues %in% c("Tumor", "NAT")) || sum(tissues == "Tumor") != 6L || sum(tissues == "NAT") != 
        10L) {
        pdc_gate_fail("PDC structured identity set must contain 6 tumor and 10 NAT records")
    }
    TRUE
}

pdc_gate_validate_publication_prose_discrepancy <- function (x) 
{
    expected <- list(tumor_samples = 10L, nat_samples = 6L, total_samples = 16L, identifiers_resolved = FALSE, 
        operative_exclusion_identity_source = FALSE, discrepancy_retained = TRUE)
    if (!is.list(x) || !identical(names(x), names(expected)) || !identical(x, expected)) {
        pdc_gate_fail("Publication prose 10/6 discrepancy must be retained as non-operative evidence")
    }
    TRUE
}

pdc_gate_validate_publication_qc_exclusions <- function (x) 
{
    required <- c("authoritative_sample_identifier", "canonical_case_identifier", "authoritative_aliquot_identifier", 
        "canonical_tissue_label", "plex_identifier", "channel_identifier", "exclusion_reason", "evidence_ids")
    if (!is.data.frame(x) || !identical(names(x), required) || nrow(x) != 17L || anyNA(x) || any(!nzchar(trimws(as.matrix(x))))) {
        pdc_gate_fail("Publication-QC exclusion ledger must contain exactly 16 TMT13 rows and one separate row")
    }
    pdc_gate_validate_identity_table(x)
    aliquots <- vapply(x$authoritative_aliquot_identifier, pdc_gate_canonical_aliquot_id, character(1))
    if (anyDuplicated(aliquots)) {
        pdc_gate_fail("Publication-QC exclusions contain a duplicate or many-to-many aliquot")
    }
    separate <- x[x$authoritative_sample_identifier == separate_identity()$sample_id, , drop = FALSE]
    if (nrow(separate) != 1L || separate$canonical_tissue_label[[1L]] != "NAT" || separate$exclusion_reason[[1L]] != 
        "unresolved_tissue_assignment") {
        pdc_gate_fail("The separately identified record must retain its separate unresolved tissue assignment reason")
    }
    is_tmt13 <- x$authoritative_sample_identifier != separate_identity()$sample_id & x$exclusion_reason == "structured_official_TMT13_exclusion"
    if (sum(is_tmt13) != 16L || sum(x$canonical_tissue_label[is_tmt13] == "Tumor") != 6L || sum(x$canonical_tissue_label[is_tmt13] == 
        "NAT") != 10L) {
        pdc_gate_fail("Structured official evidence must identify exactly 6 tumor and 10 NAT TMT13 exclusions")
    }
    if (length(unique(x$plex_identifier[is_tmt13])) != 1L || anyDuplicated(x$channel_identifier[is_tmt13])) {
        pdc_gate_fail("TMT13 plex or channel identity is incomplete or many-to-many")
    }
    if (sum(!is_tmt13) != 1L) {
        pdc_gate_fail("Publication exclusions contain an unexpected non-TMT13 row")
    }
    TRUE
}

pdc_gate_validate_table_s1_tmt13 <- function (x) 
{
    required <- c("Sample_ID", "Aliquot_ID", "Case", "tissue_type", "passed_QC", "TMT18_set", "TMT18_set_position", 
        "TMT18_channel")
    if (!is.data.frame(x) || !all(required %in% names(x))) {
        pdc_gate_fail("Table S1 identity/QC input is missing required fields")
    }
    tmt13 <- x[as.character(x$TMT18_set) == "13", required, drop = FALSE]
    if (nrow(tmt13) != 16L || anyNA(tmt13) || any(!nzchar(trimws(as.matrix(tmt13))))) {
        pdc_gate_fail("Table S1 TMT18_set 13 must contain exactly 16 complete identity rows")
    }
    passed_qc <- if (is.logical(tmt13$passed_QC)) {
        tmt13$passed_QC
    }
    else {
        token <- tolower(trimws(as.character(tmt13$passed_QC)))
        if (any(!token %in% c("true", "false"))) {
            pdc_gate_fail("Table S1 passed_QC contains an invalid Boolean token")
        }
        token == "true"
    }
    if (any(passed_qc)) {
        pdc_gate_fail("All Table S1 TMT18_set 13 records must have passed_QC false")
    }
    canonical_sample <- vapply(tmt13$Sample_ID, pdc_gate_canonical_sample_id, character(1))
    canonical_aliquot <- vapply(tmt13$Aliquot_ID, pdc_gate_canonical_aliquot_id, character(1))
    canonical_case <- vapply(tmt13$Case, pdc_gate_canonical_case_id, character(1))
    canonical_tissue <- as.character(tmt13$tissue_type)
    if (anyDuplicated(canonical_sample) || anyDuplicated(canonical_aliquot) || anyDuplicated(tmt13$TMT18_set_position) || 
        anyDuplicated(tmt13$TMT18_channel)) {
        pdc_gate_fail("Table S1 TMT13 identities and channels must be unique; duplicate found")
    }
    expected_case <- sub("[.][TN]$", "", canonical_sample)
    expected_tissue <- ifelse(endsWith(canonical_sample, ".T"), "Tumor", "NAT")
    if (any(canonical_case != expected_case) || any(canonical_tissue != expected_tissue)) {
        pdc_gate_fail("Table S1 contains a sample/case/tissue identity conflict")
    }
    if (sum(canonical_tissue == "Tumor") != 6L || sum(canonical_tissue == "NAT") != 10L) {
        pdc_gate_fail("Table S1 structured identity set must contain 6 tumor and 10 NAT records")
    }
    data.frame(canonical_aliquot_id = canonical_aliquot, canonical_case_identifier = canonical_case, 
        canonical_tissue_label = canonical_tissue, stringsAsFactors = FALSE)
}

pdc_gate_lf_set_sha256 <- function (values) 
{
    if (!is.character(values) || !length(values) || anyNA(values) || any(!nzchar(values)) || anyDuplicated(values)) {
        pdc_gate_fail("Identity-set hash input must be unique nonempty UTF-8 strings")
    }
    preimage <- paste0(paste(sort(enc2utf8(values), method = "radix"), collapse = "\n"), "\n")
    digest::digest(charToRaw(enc2utf8(preimage)), algo = "sha256", serialize = FALSE)
}

pdc_gate_project_root <- function (root = pdc_gate_expected_root) 
{
    resolved <- normalizePath(root, winslash = "/", mustWork = TRUE)
    if (!identical(resolved, pdc_gate_expected_root)) {
        pdc_gate_fail("PDC000614 TMT13 Gate requires the exact project root")
    }
    resolved
}

pdc_gate_assert_identity_sets_equal <- function (pdc, table_s1) 
{
    required <- c("canonical_aliquot_id", "canonical_case_identifier", "canonical_tissue_label")
    inputs <- list(PDC = pdc, Table_S1 = table_s1)
    for (name in names(inputs)) {
        x <- inputs[[name]]
        if (!is.data.frame(x) || !all(required %in% names(x)) || nrow(x) != 16L || anyNA(x[required]) || 
            anyDuplicated(x$canonical_aliquot_id) || anyDuplicated(paste(x$canonical_case_identifier, 
            x$canonical_tissue_label, sep = "\r"))) {
            pdc_gate_fail(paste0(name, " identity set is incomplete or duplicate"))
        }
    }
    if (!identical(sort(pdc$canonical_aliquot_id, method = "radix"), sort(table_s1$canonical_aliquot_id, 
        method = "radix"))) {
        pdc_gate_fail("PDC and Table S1 fail exact canonical aliquot-set equality")
    }
    pdc_case_tissue <- paste(pdc$canonical_case_identifier, pdc$canonical_tissue_label, sep = "\r")
    table_case_tissue <- paste(table_s1$canonical_case_identifier, table_s1$canonical_tissue_label, sep = "\r")
    if (!identical(sort(pdc_case_tissue, method = "radix"), sort(table_case_tissue, method = "radix"))) {
        pdc_gate_fail("PDC and Table S1 fail exact case-plus-tissue-set equality")
    }
    TRUE
}

pdc_gate_assert_structured_count_contract <- function (structured_tumor_n, structured_nat_n) 
{
    counts <- c(structured_tumor_n, structured_nat_n)
    if (!is.numeric(counts) || length(counts) != 2L || anyNA(counts) || any(counts < 0) || any(counts != 
        as.integer(counts))) {
        pdc_gate_fail("Structured identity counts must be nonnegative integers")
    }
    if (identical(as.integer(counts), c(6L, 10L))) {
        return(TRUE)
    }
    pdc_gate_fail("Structured official identity set must contain exactly 6 tumor and 10 NAT records")
}

pdc_gate_validate_identity_table <- function (x) 
{
    required <- c("authoritative_sample_identifier", "canonical_case_identifier", "canonical_tissue_label")
    if (!is.data.frame(x) || !all(required %in% names(x))) {
        pdc_gate_fail("Identity table is missing required fields")
    }
    if (!nrow(x) || anyNA(x[required]) || any(!nzchar(trimws(as.matrix(x[required]))))) {
        pdc_gate_fail("Identity table contains a missing identity")
    }
    canonical_sample <- vapply(x$authoritative_sample_identifier, pdc_gate_canonical_sample_id, character(1))
    if (anyDuplicated(canonical_sample)) {
        pdc_gate_fail("Identity table contains a duplicate or many-to-many sample")
    }
    tissues <- as.character(x$canonical_tissue_label)
    if (any(!tissues %in% c("Tumor", "NAT"))) {
        pdc_gate_fail("Identity table contains an ambiguous tissue label")
    }
    parsed_case <- sub("[.][TN]$", "", canonical_sample)
    parsed_tissue <- ifelse(endsWith(canonical_sample, ".T"), "Tumor", "NAT")
    if (any(parsed_case != x$canonical_case_identifier) || any(parsed_tissue != tissues)) {
        pdc_gate_fail("Identity table contains a sample/case/tissue conflict")
    }
    key <- paste(x$canonical_case_identifier, tissues, sep = "\r")
    if (anyDuplicated(key)) {
        pdc_gate_fail("Identity table contains a duplicate or many-to-many identity")
    }
    invisible(TRUE)
}

  list(pdc_gate_assess_authoritative_protein_provenance=pdc_gate_assess_authoritative_protein_provenance,
    pdc_gate_assess_bridge_endpoint=pdc_gate_assess_bridge_endpoint,
    pdc_gate_assess_pair_endpoint=pdc_gate_assess_pair_endpoint,
    pdc_gate_attach_exclusions=pdc_gate_attach_exclusions,
    pdc_gate_background_payload=pdc_gate_background_payload,
    pdc_gate_read_clinical_identity=pdc_gate_read_clinical_identity,
    pdc_gate_read_linkedomics_identity=pdc_gate_read_linkedomics_identity,
    pdc_gate_read_masked_matrix=pdc_gate_read_masked_matrix,
    pdc_gate_reconcile_official_exclusions=pdc_gate_reconcile_official_exclusions,
    pdc_gate_structural_pair_count=pdc_gate_structural_pair_count,
    pdc_gate_structural_quadruplet_count=pdc_gate_structural_quadruplet_count,
    pdc_gate_build_complete_pairs=pdc_gate_build_complete_pairs,
    pdc_gate_complete_quadruplets=pdc_gate_complete_quadruplets,
    pdc_gate_collect_authoritative_protein_evidence=pdc_gate_collect_authoritative_protein_evidence,
    pdc_gate_validate_authoritative_protein_evidence=pdc_gate_validate_authoritative_protein_evidence,
    pdc_gate_fail=pdc_gate_fail,
    pdc_gate_validate_paired_design=pdc_gate_validate_paired_design,
    pdc_gate_canonical_case_id=pdc_gate_canonical_case_id,
    pdc_gate_read_allowlisted_tabular_identity=pdc_gate_read_allowlisted_tabular_identity,
    pdc_gate_canonical_aliquot_id=pdc_gate_canonical_aliquot_id,
    pdc_gate_canonical_sample_id=pdc_gate_canonical_sample_id,
    pdc_gate_sha256=pdc_gate_sha256,
    pdc_gate_classify_pdc_channels=pdc_gate_classify_pdc_channels,
    pdc_gate_contract=pdc_gate_contract,
    pdc_gate_flatten_pdc_case_identity=pdc_gate_flatten_pdc_case_identity,
    pdc_gate_identity_set_hashes=pdc_gate_identity_set_hashes,
    pdc_gate_read_pdc_design_channels=pdc_gate_read_pdc_design_channels,
    pdc_gate_read_pdc_reference_identity=pdc_gate_read_pdc_reference_identity,
    pdc_gate_read_table_s1_identity=pdc_gate_read_table_s1_identity,
    pdc_gate_validate_correction03_identity_contract=pdc_gate_validate_correction03_identity_contract,
    pdc_gate_validate_pdc_tmt13_design=pdc_gate_validate_pdc_tmt13_design,
    pdc_gate_validate_publication_prose_discrepancy=pdc_gate_validate_publication_prose_discrepancy,
    pdc_gate_validate_publication_qc_exclusions=pdc_gate_validate_publication_qc_exclusions,
    pdc_gate_validate_table_s1_tmt13=pdc_gate_validate_table_s1_tmt13,
    pdc_gate_lf_set_sha256=pdc_gate_lf_set_sha256,
    pdc_gate_project_root=pdc_gate_project_root,
    pdc_gate_assert_identity_sets_equal=pdc_gate_assert_identity_sets_equal,
    pdc_gate_assert_structured_count_contract=pdc_gate_assert_structured_count_contract,
    pdc_gate_validate_identity_table=pdc_gate_validate_identity_table)
}
