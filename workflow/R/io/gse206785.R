# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.
# All raw/source paths are explicit arguments. The native parser compiles
# only when requested; TMPDIR must be inside the current isolated work area.
.gse206785_raw_umi_parser_environment <- new.env(parent=globalenv())

gse206785_spec <- function () 
{
    list(accession = "GSE206785", chemistry_authority = paste0("fixed_commit_cell_metadata.csv:Platform_with_exact_GEO_metadata_", 
        "Platform_corroboration"), metadata_schema = c("Sample", "Patient", "Tissue", "Platform", "Subtype", 
        "Type", "Annotation"), tissues = c("Normal", "Tumor"), platforms = c("SC3P", "SC5P"), compartments = c("epithelial", 
        "fibroblast", "endothelial", "myeloid", "T_NK", "B_plasma", "mast_other_immune"), programs = list(G = c("ACSL4", 
        "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B", "POR"), H = c("ALDH2", "ALDH3A1", "ALDH3A2", 
        "ALDH7A1")), qc_rules = gse206785_frozen_qc_rules(), type_to_compartment = gse206785_frozen_type_mapping(), 
        unit_rules = list(cells_min_inclusive = 20L, umi_min_inclusive = 10000, paired_donors_min_inclusive = 8L, 
            member_nonzero_fraction_min_inclusive = 0.5, member_nonzero_both_tissues_paired_donors_min_inclusive = 3L), 
        expected = list(metadata_rows = 111140L, samples = 48L, patients = 24L, chemistry_donors = c(SC3P = 14L, 
            SC5P = 10L), genes = 25287L, cells = 111140L, historical_fragments = 9L, fibroblast_raw_pseudobulk_sha256 = "1500b3bf1e4d8ec87b1a8dc892ce6a119cfb304e4f80e521435b8b3a00454349", 
            fixed_universe = c(numeric_qc = 107175L, subtype_only = 97914L, intersection = 95668L, mural_glial = 3262L, 
                seven_compartment = 92406L), formal_paired_donors = c(myeloid_G = 10L, myeloid_H = 10L, 
                fibroblast_H = 11L, B_plasma_H = 16L), composition = c(primary_100 = 20L, strict_200 = 18L, 
                count_aware = 24L)), disclosures = list(published_source_context_displayed = TRUE, published_source_context_classification = paste0("DOCUMENTED_PUBLISHED_SOURCE_CONTEXT_DISPLAY_WITHOUT_DEFINED_", 
            "PROJECT_DISEASE_EFFECT_ACCESS"), source_role_correction_classification = paste0("DOCUMENTED_STAGE_C2_FIELD_SPECIFIC_SOURCE_ROLE_ERROR_WITHOUT_", 
            "PROJECT_EFFECT_ACCESS"), raw_umi_input_correction_classification = paste0("DOCUMENTED_STAGE_C2_WRONG_PROCESSED_MATRIX_INPUT_WITHOUT_", 
            "PROJECT_DISEASE_EFFECT_ACCESS"), project_derived_effect_accessed = FALSE, g_or_h_effect_accessed = FALSE, 
            composition_effect_calculated_or_viewed = FALSE, result_ledger_action = "remain_header_only", 
            osf_prior_access_disclosure_required = TRUE, marker_audit_executed = FALSE, disease_effect_accessed = FALSE, 
            result_ledger_appended = FALSE))
}

gse206785_reconcile_metadata_sources <- function (author, geo, expected_rows = 111140L) 
{
    author_schema <- c("Sample", "Patient", "Tissue", "Platform", "Type", "Subtype")
    geo_schema <- c("Sample", "Patient", "Tissue", "Platform", "Subtype", "Type", "Annotation")
    if (!identical(names(author), author_schema)) {
        gse206785_abort("fixed-commit cell_metadata.csv schema mismatch")
    }
    if (!identical(names(geo), geo_schema)) {
        gse206785_abort("GEO processed metadata schema mismatch")
    }
    if (nrow(author) != expected_rows || nrow(geo) != expected_rows) {
        gse206785_abort("author and GEO metadata row-count mismatch")
    }
    common <- author_schema
    mismatch <- vapply(common, function(field) {
        !identical(as.character(author[[field]]), as.character(geo[[field]]))
    }, logical(1))
    if (any(mismatch)) {
        gse206785_abort(paste0("row-for-row metadata mismatch on: ", paste(common[mismatch], collapse = ",")))
    }
    combined <- author
    combined$Annotation <- geo$Annotation
    combined <- combined[c("Sample", "Patient", "Tissue", "Platform", "Subtype", "Type", "Annotation")]
    list(metadata = combined, common_fields = common, row_count = nrow(combined), rowwise_mismatch_count = 0L, 
        chemistry_field_match = TRUE)
}

gse206785_validate_metadata <- function (metadata, expected) 
{
    schema <- gse206785_spec()$metadata_schema
    if (!identical(names(metadata), schema)) {
        missing <- setdiff(schema, names(metadata))
        if (length(missing)) {
            gse206785_abort(paste("metadata is missing", paste(missing, collapse = ", ")))
        }
        gse206785_abort("metadata schema or column order mismatch")
    }
    identity_fields <- c("Sample", "Patient", "Tissue", "Platform")
    missing_identity <- vapply(metadata[identity_fields], function(x) any(is.na(x) | !nzchar(trimws(as.character(x)))), 
        logical(1))
    if (any(missing_identity)) {
        gse206785_abort("metadata contains missing identity values")
    }
    if (nrow(metadata) != expected$rows) {
        gse206785_abort("metadata row count mismatch")
    }
    unsupported <- setdiff(unique(metadata$Platform), c("SC3P", "SC5P"))
    if (length(unsupported)) {
        gse206785_abort("Platform vocabulary must be restricted to SC3P and SC5P")
    }
    for (sample in unique(metadata$Sample)) {
        z <- metadata[metadata$Sample == sample, identity_fields, drop = FALSE]
        for (field in c("Patient", "Tissue", "Platform")) {
            if (length(unique(z[[field]])) != 1L) {
                gse206785_abort(paste("sample has multiple", field, "assignments"))
            }
        }
    }
    sample_ledger <- unique(metadata[identity_fields])
    if (nrow(sample_ledger) != expected$samples || length(unique(sample_ledger$Sample)) != expected$samples) {
        gse206785_abort("metadata sample count mismatch")
    }
    if (length(unique(sample_ledger$Patient)) != expected$patients) {
        gse206785_abort("metadata patient count mismatch")
    }
    if (!setequal(unique(sample_ledger$Tissue), expected$tissues)) {
        gse206785_abort("metadata tissue vocabulary mismatch")
    }
    pair_counts <- table(sample_ledger$Patient, sample_ledger$Tissue)
    if (!identical(sort(colnames(pair_counts)), sort(expected$tissues)) || any(pair_counts != 1L)) {
        gse206785_abort("each patient-tissue identity must have exactly one sample")
    }
    suffix_expected <- ifelse(sample_ledger$Tissue == "Normal", "N", "T")
    if (any(!endsWith(sample_ledger$Sample, suffix_expected))) {
        gse206785_abort("sample suffix and tissue identity conflict")
    }
    sample_base <- sub("[NT]$", "", sample_ledger$Sample)
    base_bridge <- unique(data.frame(sample_base = sample_base, Patient = sample_ledger$Patient, stringsAsFactors = FALSE))
    if (anyDuplicated(base_bridge$sample_base) || anyDuplicated(base_bridge$Patient) || nrow(base_bridge) != 
        expected$patients) {
        gse206785_abort("sample-base to patient bridge is not one-to-one")
    }
    donor_platform <- unique(sample_ledger[c("Patient", "Platform")])
    mismatch <- names(which(table(donor_platform$Patient) != 1L))
    if (length(mismatch)) {
        gse206785_abort("authoritative within-donor chemistry mismatch")
    }
    chemistry <- table(donor_platform$Platform)
    chemistry_counts <- setNames(as.integer(chemistry[c("SC3P", "SC5P")]), c("SC3P", "SC5P"))
    chemistry_counts[is.na(chemistry_counts)] <- 0L
    if (!identical(chemistry_counts, expected$platform_donors)) {
        gse206785_abort("authoritative donor chemistry allocation mismatch")
    }
    tissue_order <- match(sample_ledger$Tissue, expected$tissues)
    sample_ledger <- sample_ledger[order(sample_ledger$Patient, tissue_order, sample_ledger$Sample), 
        , drop = FALSE]
    names(sample_ledger) <- c("sample", "donor", "tissue", "chemistry")
    rownames(sample_ledger) <- NULL
    list(sample_ledger = sample_ledger, chemistry_counts = chemistry_counts, within_donor_mismatches = 0L, 
        chemistry_eligible_donors = expected$patients)
}

gse206785_parse_soft_samples <- function (path) 
{
    if (!file.exists(path)) 
        gse206785_abort("family SOFT is missing")
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    samples <- list()
    supplementary <- character()
    current <- NULL
    flush <- function() {
        if (!is.null(current)) 
            samples[[length(samples) + 1L]] <<- current
    }
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (!length(line)) 
            break
        if (startsWith(line, "^SAMPLE = ")) {
            flush()
            current <- list(gsm = sub("^\\^SAMPLE = ", "", line), title = NA_character_, library_strategy = NA_character_)
            next
        }
        if (startsWith(line, "^")) {
            flush()
            current <- NULL
            next
        }
        if (startsWith(line, "!Series_supplementary_file")) {
            supplementary <- c(supplementary, basename(sub("^[^=]*= ?", "", line)))
        }
        if (!is.null(current) && startsWith(line, "!Sample_geo_accession")) {
            current$gsm <- sub("^[^=]*= ?", "", line)
        }
        if (!is.null(current) && startsWith(line, "!Sample_title")) {
            current$title <- sub("^[^=]*= ?", "", line)
        }
        if (!is.null(current) && startsWith(line, "!Sample_library_strategy")) {
            current$library_strategy <- sub("^[^=]*= ?", "", line)
        }
    }
    flush()
    sample_table <- do.call(rbind, lapply(samples, as.data.frame, stringsAsFactors = FALSE))
    rownames(sample_table) <- NULL
    list(samples = sample_table, series_supplementary_files = unique(supplementary))
}

gse206785_reconcile_soft_metadata <- function (soft_samples, chemistry_ledger, expected_samples) 
{
    if (nrow(soft_samples) != expected_samples || anyDuplicated(soft_samples$gsm) || anyDuplicated(soft_samples$title) || 
        any(is.na(soft_samples$gsm) | is.na(soft_samples$title))) {
        gse206785_abort("SOFT sample records are not unique and complete")
    }
    normalized_sample <- gse206785_soft_title_to_sample(soft_samples$title)
    if (!identical(sort(normalized_sample), sort(chemistry_ledger$sample))) {
        gse206785_abort("SOFT GSM-title population differs from processed metadata")
    }
    bridge <- data.frame(gsm = soft_samples$gsm, source_title = soft_samples$title, sample = normalized_sample, 
        stringsAsFactors = FALSE)
    bridge <- bridge[order(bridge$sample), , drop = FALSE]
    rownames(bridge) <- NULL
    invisible(bridge)
}

gse206785_parse_raw_umi_gzip <- function (path, expected_genes, expected_cells, core_members = unlist(gse206785_spec()$programs, use.names = FALSE)) 
{
    if (!file.exists(path)) 
        gse206785_abort("raw-UMI gzip is missing")
    if (length(expected_genes) != 1L || length(expected_cells) != 1L || expected_genes < 1L || expected_cells < 
        1L) {
        gse206785_abort("raw-UMI expected dimensions are invalid")
    }
    gse206785_compile_raw_umi_parser()
    parser <- get("gse206785_parse_raw_umi_cpp", envir = .gse206785_raw_umi_parser_environment, inherits = FALSE)
    parsed <- parser(normalizePath(path, mustWork = TRUE), as.integer(expected_genes), as.integer(expected_cells), 
        core_members)
    cell_ids <- as.character(parsed$cell_ids)
    suffix <- sub("^.*-", "", cell_ids)
    format_pass <- grepl(".+-.+", cell_ids) & nzchar(suffix)
    structure <- list(cell_count = as.integer(parsed$cell_count), gene_names = as.character(parsed$gene_names), 
        gene_count = as.integer(parsed$gene_count), cell_ids = cell_ids, cell_ids_unique = !anyDuplicated(cell_ids), 
        sample_suffix = suffix, cell_id_format_pass = all(format_pass), suffix_population = sort(unique(suffix)), 
        header_count = as.integer(parsed$header_count), value_encoding = "strict_nonnegative_integer_raw_UMI", 
        orientation = "gene_rows_by_cell_columns", strict_integer_tokens = TRUE, barcode_key_available = FALSE)
    metrics <- data.frame(detected_genes = as.integer(parsed$detected_genes), total_umi = as.numeric(parsed$total_umi), 
        mitochondrial_fraction = as.numeric(parsed$mitochondrial_fraction), as.data.frame(parsed$core_counts, 
            check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
    list(structure = structure, metrics = metrics)
}

gse206785_reconcile_release <- function (count_structure, metadata, expected_genes, expected_cells, expected_payloads) 
{
    required <- c("cell_count", "gene_names", "gene_count", "cell_ids_unique", "sample_suffix", "cell_id_format_pass", 
        "suffix_population", "header_count", "physical_payload_count", "documented_chain", "barcode_key_available", 
        "orientation", "strict_integer_tokens", "combined_bytes", "combined_sha256")
    if (!all(required %in% names(count_structure))) {
        gse206785_abort("count structure is incomplete")
    }
    if (count_structure$gene_count != expected_genes || length(count_structure$gene_names) != expected_genes) {
        gse206785_abort("count gene count mismatch")
    }
    if (anyDuplicated(count_structure$gene_names)) {
        gse206785_abort("count file contains duplicate gene rows")
    }
    if (count_structure$cell_count != expected_cells || length(count_structure$sample_suffix) != expected_cells) {
        gse206785_abort("count cell count mismatch")
    }
    if (!isTRUE(count_structure$cell_ids_unique)) {
        gse206785_abort("count file contains duplicate cell identifiers")
    }
    if (!isTRUE(count_structure$cell_id_format_pass)) {
        gse206785_abort("count sample-suffix reconstruction failed")
    }
    if (!identical(count_structure$orientation, "gene_rows_by_cell_columns") || !isTRUE(count_structure$strict_integer_tokens)) {
        gse206785_abort("formal raw-UMI orientation or integer-token contract failed")
    }
    expected_population <- sort(unique(metadata$Sample))
    if (!identical(sort(count_structure$suffix_population), expected_population)) {
        gse206785_abort("count-column sample-suffix population mismatch")
    }
    if (!identical(count_structure$sample_suffix, metadata$Sample)) {
        gse206785_abort("count-column and metadata order mismatch")
    }
    if (count_structure$physical_payload_count != expected_payloads || count_structure$header_count != 
        1L || count_structure$documented_chain != "nine_fixed_commit_byte_fragments_one_logical_gzip_stream") {
        gse206785_abort("deposited fragment chain mismatch")
    }
    contract <- gse206785_fixed_commit_contract()
    if (!identical(as.numeric(count_structure$combined_bytes), contract$combined_bytes) || !identical(as.character(count_structure$combined_sha256), 
        contract$combined_sha256)) {
        gse206785_abort("combined raw-UMI stream identity mismatch")
    }
    if (isTRUE(count_structure$barcode_key_available)) {
        gse206785_abort("a barcode-keyed join claim is prohibited")
    }
    list(gene_count = as.integer(count_structure$gene_count), cell_count = as.integer(count_structure$cell_count), 
        physical_payload_count = as.integer(count_structure$physical_payload_count), logical_gzip_stream_count = 1L, 
        header_count = as.integer(count_structure$header_count), documented_chain = as.character(count_structure$documented_chain), 
        sample_suffix_population = length(count_structure$suffix_population), order_suffix_reconciliation_pass = TRUE, 
        unique_cell_identifiers_pass = TRUE, unique_gene_identifiers_pass = TRUE, barcode_key_available = FALSE, 
        identity_chain = paste("fixed_author_row_order", "unique_count_identifiers", "per_row_sample_suffix_agreement", 
            sep = ";"), historical_fragment_expectation = "9_authoritative_fixed_commit_contract", historical_fragment_expectation_reconciled = identical(as.integer(count_structure$physical_payload_count), 
            9L))
}

gse206785_apply_cell_universe <- function (metrics, metadata, rules = gse206785_frozen_qc_rules(), type_mapping = gse206785_frozen_type_mapping(), 
    expected_counts = NULL) 
{
    gse206785_validate_qc_rules(rules)
    gse206785_validate_type_mapping(type_mapping)
    required <- c("detected_genes", "total_umi", "mitochondrial_fraction")
    if (!all(required %in% names(metrics)) || nrow(metrics) != nrow(metadata)) {
        gse206785_abort("cell-QC metric schema or row count mismatch")
    }
    values <- unlist(metrics[required], use.names = FALSE)
    if (any(!is.finite(values)) || any(values < 0)) {
        gse206785_abort("cell-QC metrics must be finite and nonnegative")
    }
    unknown_types <- setdiff(unique(metadata$Type), names(type_mapping))
    if (length(unknown_types)) {
        gse206785_abort("unexpected author broad Type")
    }
    numeric_qc_pass <- metrics$detected_genes >= rules$detected_genes_min_inclusive & metrics$detected_genes <= 
        rules$detected_genes_max_inclusive & metrics$total_umi >= rules$umi_min_inclusive & metrics$total_umi <= 
        rules$umi_max_inclusive & metrics$mitochondrial_fraction < rules$mitochondrial_fraction_max_exclusive
    subtype_pass <- !metadata$Subtype %in% rules$low_quality_subtypes
    intersection_pass <- numeric_qc_pass & subtype_pass
    descriptive_exclusion <- metadata$Type %in% rules$descriptive_exclusions
    compartment <- unname(type_mapping[metadata$Type])
    seven_compartment_pass <- intersection_pass & !descriptive_exclusion & !is.na(compartment)
    counts <- c(numeric_qc = sum(numeric_qc_pass), subtype_only = sum(subtype_pass), intersection = sum(intersection_pass), 
        mural_glial = sum(intersection_pass & descriptive_exclusion), seven_compartment = sum(seven_compartment_pass))
    storage.mode(counts) <- "integer"
    if (!is.null(expected_counts) && !identical(counts, expected_counts)) {
        gse206785_abort("unexplained fixed-QC count mismatch")
    }
    list(flags = data.frame(numeric_qc_pass = numeric_qc_pass, subtype_pass = subtype_pass, intersection_pass = intersection_pass, 
        descriptive_exclusion = descriptive_exclusion, seven_compartment_pass = seven_compartment_pass, 
        compartment = compartment, stringsAsFactors = FALSE), counts = counts)
}

gse206785_aggregate_raw_umi <- function (path, group_index, group_names, expected_genes, expected_cells, scheme_names = NULL) 
{
    if (is.null(dim(group_index))) {
        group_index <- matrix(group_index, ncol = 1L)
    }
    if (!is.matrix(group_index) || nrow(group_index) != expected_cells || !length(group_names) || anyDuplicated(group_names) || 
        any(is.na(group_index)) || any(!group_index %in% c(-1L, seq_along(group_names) - 1L))) {
        gse206785_abort("raw-UMI pseudobulk grouping contract mismatch")
    }
    if (is.null(scheme_names)) {
        scheme_names <- paste0("scheme_", seq_len(ncol(group_index)))
    }
    if (length(scheme_names) != ncol(group_index) || anyDuplicated(scheme_names)) {
        gse206785_abort("raw-UMI pseudobulk scheme-name contract mismatch")
    }
    gse206785_compile_raw_umi_parser()
    aggregator <- get("gse206785_aggregate_raw_umi_cpp", envir = .gse206785_raw_umi_parser_environment, 
        inherits = FALSE)
    result <- aggregator(normalizePath(path, mustWork = TRUE), as.integer(expected_genes), as.integer(expected_cells), 
        matrix(as.integer(group_index), nrow = nrow(group_index)), as.integer(length(group_names)))
    counts <- result$counts
    split_columns <- split(seq_len(ncol(counts)), rep(scheme_names, each = length(group_names)))
    out <- lapply(split_columns, function(index) {
        z <- counts[, index, drop = FALSE]
        dimnames(z) <- list(as.character(result$genes), group_names)
        z
    })
    if (length(out) == 1L) 
        out[[1L]]
    else out[scheme_names]
}

gse206785_apply_unit_eligibility <- function (units, rules = gse206785_frozen_unit_rules()) 
{
    gse206785_validate_unit_rules(rules)
    if (!all(c("cell_count", "total_umi") %in% names(units))) {
        gse206785_abort("structural unit schema is incomplete")
    }
    units$cell_count >= rules$cells_min_inclusive & units$total_umi >= rules$umi_min_inclusive
}

gse206785_aggregate_units <- function (metadata, metrics, flags) 
{
    members <- unlist(gse206785_spec()$programs, use.names = FALSE)
    if (!all(c("total_umi", members) %in% names(metrics)) || nrow(metadata) != nrow(metrics) || nrow(flags) != 
        nrow(metadata)) {
        gse206785_abort("cell-to-unit aggregation inputs are inconsistent")
    }
    keep <- flags$seven_compartment_pass
    if (!any(keep)) 
        gse206785_abort("seven-compartment universe is empty")
    key_data <- data.frame(donor = metadata$Patient[keep], tissue = metadata$Tissue[keep], compartment = flags$compartment[keep], 
        stringsAsFactors = FALSE)
    key <- interaction(key_data$donor, key_data$tissue, key_data$compartment, drop = TRUE, lex.order = TRUE)
    numeric_data <- cbind(total_umi = metrics$total_umi[keep], as.matrix(metrics[keep, members, drop = FALSE]))
    sums <- rowsum(numeric_data, key, reorder = FALSE)
    first <- !duplicated(key)
    keys <- key_data[first, , drop = FALSE]
    counts <- as.integer(table(key)[rownames(sums)])
    out <- cbind(keys, cell_count = counts, as.data.frame(sums), stringsAsFactors = FALSE)
    out$eligible <- gse206785_apply_unit_eligibility(out)
    out <- out[order(out$compartment, out$donor, match(out$tissue, c("Normal", "Tumor"))), , drop = FALSE]
    rownames(out) <- NULL
    out
}

gse206785_canonical_fibroblast_pseudobulk_hash <- function (path, metadata, final_flag, expected_genes, expected_cells, expected_columns, expected_digest, 
    contract = gse206785_canonical_hash_contract()) 
{
    gse206785_validate_canonical_hash_contract(contract)
    required <- c("Patient", "Tissue", "Type")
    if (!is.data.frame(metadata) || !all(required %in% names(metadata)) || nrow(metadata) != expected_cells || 
        length(final_flag) != expected_cells || !is.logical(final_flag) || anyNA(final_flag) || anyNA(metadata[, 
        required, drop = FALSE])) {
        gse206785_abort("canonical pseudobulk metadata contract mismatch")
    }
    donors <- sort(unique(as.character(metadata$Patient)), method = "radix")
    key <- do.call(rbind, lapply(donors, function(donor) {
        data.frame(Patient = rep(donor, 2L), Tissue = c("Normal", "Tumor"), stringsAsFactors = FALSE)
    }))
    if (nrow(key) != expected_columns || expected_genes < 1L || any(grepl("\r", c(metadata$Patient, metadata$Tissue), 
        fixed = TRUE))) {
        gse206785_abort("canonical pseudobulk dimension or identity contract mismatch")
    }
    observed_key <- unique(metadata[, c("Patient", "Tissue"), drop = FALSE])
    observed_key$Patient <- as.character(observed_key$Patient)
    observed_key$Tissue <- as.character(observed_key$Tissue)
    key_token <- paste(key$Patient, key$Tissue, sep = "\r")
    observed_token <- paste(observed_key$Patient, observed_key$Tissue, sep = "\r")
    if (anyDuplicated(key_token) || !setequal(observed_token, key_token)) {
        gse206785_abort("canonical pseudobulk donor/tissue key contract mismatch")
    }
    selected <- final_flag & as.character(metadata$Type) == "Fibroblast"
    cell_token <- paste(metadata$Patient, metadata$Tissue, sep = "\r")
    group_index <- rep.int(-1L, expected_cells)
    group_index[selected] <- match(cell_token[selected], key_token) - 1L
    if (anyNA(group_index[selected]) || !any(selected)) {
        gse206785_abort("canonical fibroblast cell selection is empty or unresolved")
    }
    if (length(expected_digest) != 1L || !grepl("^[0-9a-f]{64}$", expected_digest)) {
        gse206785_abort("canonical pseudobulk expected digest is invalid")
    }
    gse206785_compile_raw_umi_parser()
    writer <- get("gse206785_write_canonical_pseudobulk_preimage_cpp", envir = .gse206785_raw_umi_parser_environment, 
        inherits = FALSE)
    preimage_path <- tempfile("gse206785_canonical_preimage_", fileext = ".bin")
    if (file.exists(preimage_path)) {
        gse206785_abort("unexpected canonical preimage path exists")
    }
    on.exit(if (file.exists(preimage_path)) unlink(preimage_path), add = TRUE)
    result <- writer(normalizePath(path, mustWork = TRUE), preimage_path, as.integer(expected_genes), 
        as.integer(expected_cells), as.integer(group_index), as.integer(expected_columns))
    preimage_bytes <- as.numeric(file.info(preimage_path)$size)
    if (!is.finite(preimage_bytes) || !identical(preimage_bytes, as.numeric(result$preimage_bytes))) {
        gse206785_abort("canonical pseudobulk preimage byte-count mismatch")
    }
    observed_digest <- gse206785_sha256(preimage_path)
    if (!identical(observed_digest, expected_digest)) {
        gse206785_abort("canonical fibroblast pseudobulk digest mismatch")
    }
    list(reference_implementation = contract$reference_implementation, reference_implementation_sha256 = contract$reference_implementation_sha256, 
        reference_state = contract$reference_state, compartment = contract$compartment, cell_selection_rule = contract$cell_selection_rule, 
        donor_ordering_rule = contract$donor_ordering_rule, tissue_ordering_rule = contract$tissue_ordering_rule, 
        gene_ordering_rule = contract$gene_ordering_rule, gene_symbol_encoding = contract$gene_symbol_encoding, 
        integer_encoding = contract$integer_encoding, matrix_rows = as.integer(result$matrix_rows), matrix_columns = as.integer(result$matrix_columns), 
        integer_nonnegative = isTRUE(result$integer_nonnegative), matrix_persisted = FALSE, preimage_bytes = preimage_bytes, 
        observed_digest = observed_digest, expected_digest = expected_digest, checkpoint_status = "pass")
}

gse206785_detect_members <- function (units) 
{
    spec <- gse206785_spec()
    all_members <- unlist(spec$programs, use.names = FALSE)
    required <- c("donor", "tissue", "compartment", "eligible", all_members)
    if (!all(required %in% names(units))) {
        gse206785_abort("member-detection unit schema is incomplete")
    }
    if (anyDuplicated(units[c("donor", "tissue", "compartment")])) {
        gse206785_abort("duplicate structural pseudobulk unit")
    }
    out <- list()
    k <- 0L
    for (compartment in unique(units$compartment)) {
        eligible <- units[units$compartment == compartment & units$eligible, , drop = FALSE]
        for (program in names(spec$programs)) {
            for (member in spec$programs[[program]]) {
                k <- k + 1L
                nonzero <- eligible[[member]] > 0
                fraction <- if (nrow(eligible)) 
                  sum(nonzero)/nrow(eligible)
                else 0
                donor_pairs <- split(eligible, eligible$donor)
                both <- sum(vapply(donor_pairs, function(z) {
                  nrow(z) == 2L && setequal(z$tissue, spec$tissues) && all(z[[member]] > 0)
                }, logical(1)))
                fraction_pass <- fraction >= spec$unit_rules$member_nonzero_fraction_min_inclusive
                pair_pass <- both >= spec$unit_rules$member_nonzero_both_tissues_paired_donors_min_inclusive
                out[[k]] <- data.frame(compartment = compartment, program = program, member = member, 
                  eligible_pseudobulk_units = nrow(eligible), nonzero_pseudobulk_units = sum(nonzero), 
                  nonzero_unit_fraction = fraction, nonzero_both_tissues_paired_donors = as.integer(both), 
                  fraction_rule_pass = fraction_pass, paired_both_tissues_rule_pass = pair_pass, available = fraction_pass && 
                    pair_pass, stringsAsFactors = FALSE)
            }
        }
    }
    do.call(rbind, out)
}

gse206785_build_eligibility_matrix <- function (unit_counts, detection, chemistry_eligible_donors, published_context = NULL) 
{
    spec <- gse206785_spec()
    required_units <- c("compartment", "eligible_tumor_units", "eligible_normal_units", "complete_paired_donors")
    required_detection <- c("compartment", "program", "member", "available")
    if (!all(required_units %in% names(unit_counts)) || !all(required_detection %in% names(detection))) {
        gse206785_abort("eligibility-matrix inputs are incomplete")
    }
    if (!identical(unit_counts$compartment, spec$compartments)) {
        gse206785_abort("unit-count compartment order mismatch")
    }
    all_members <- unlist(spec$programs, use.names = FALSE)
    rows <- list()
    k <- 0L
    for (compartment in spec$compartments) {
        unit <- unit_counts[unit_counts$compartment == compartment, , drop = FALSE]
        for (program in c("G", "H")) {
            k <- k + 1L
            members <- spec$programs[[program]]
            z <- detection[detection$compartment == compartment & detection$program == program & detection$member %in% 
                members, , drop = FALSE]
            available <- members[members %in% z$member[z$available]]
            missing <- setdiff(members, available)
            required_count <- if (program == "G") 
                5L
            else 3L
            member_count_pass <- length(available) >= required_count
            coverage <- gse206785_program_coverage(program, available)
            anchor_pass <- coverage
            paired_pass <- unit$complete_paired_donors >= spec$unit_rules$paired_donors_min_inclusive
            estimable <- coverage && paired_pass
            plan_role <- gse206785_plan_role(compartment, program)
            gate_role <- if (!estimable) {
                "non_estimable_under_frozen_rules"
            }
            else if (plan_role == "formal_candidate") {
                "structurally_eligible_pending_marker_audit"
            }
            else if (plan_role == "prespecified_exploratory_candidate") {
                "exploratory_structurally_eligible_pending_marker_audit"
            }
            else {
                "descriptive_structurally_estimable"
            }
            reason <- if (!paired_pass) {
                "fewer_than_eight_complete_paired_donors"
            }
            else if (!coverage) {
                "frozen_program_coverage_failed"
            }
            else if (plan_role == "formal_candidate") {
                "formal_candidate_structurally_eligible_marker_audit_not_run"
            }
            else if (plan_role == "prespecified_exploratory_candidate") {
                "exploratory_candidate_structurally_eligible_marker_audit_not_run"
            }
            else {
                "descriptive_only_under_prespecified_plan_role"
            }
            row <- data.frame(compartment = compartment, program = program, chemistry_eligible_donors = as.integer(chemistry_eligible_donors), 
                eligible_tumor_units = as.integer(unit$eligible_tumor_units), eligible_normal_units = as.integer(unit$eligible_normal_units), 
                paired_pseudobulk_eligible_donors = as.integer(unit$complete_paired_donors), complete_paired_donors = as.integer(unit$complete_paired_donors), 
                eligible_unit_rule = "cells>=20_and_total_UMI>=10000", available_member_count = as.integer(length(available)), 
                required_member_count = required_count, available_members = paste(available, collapse = ","), 
                missing_members = paste(missing, collapse = ","), member_detection_rule_pass = member_count_pass, 
                program_coverage_pass = coverage, structural_requirement_pass = anchor_pass, minimum_eight_paired_donors_pass = paired_pass, 
                state_estimability = if (estimable) 
                  "structurally_estimable"
                else "non_estimable", prespecified_plan_role = plan_role, gate_role = gate_role, gate_role_reason = reason, 
                marker_audit_status = "NOT_RUN_IN_STAGE_C2", stringsAsFactors = FALSE)
            for (member in all_members) {
                value <- if (member %in% members) 
                  member %in% available
                else NA
                row[[paste0("member_", member, "_available")]] <- value
            }
            rows[[k]] <- row
        }
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

gse206785_unit_counts <- function (units) 
{
    spec <- gse206785_spec()
    rows <- lapply(spec$compartments, function(compartment) {
        z <- units[units$compartment == compartment & units$eligible, , drop = FALSE]
        paired <- split(z, z$donor)
        complete <- sum(vapply(paired, function(x) {
            nrow(x) == 2L && setequal(x$tissue, spec$tissues)
        }, logical(1)))
        data.frame(compartment = compartment, eligible_tumor_units = sum(z$tissue == "Tumor"), eligible_normal_units = sum(z$tissue == 
            "Normal"), complete_paired_donors = as.integer(complete), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

gse206785_donor_compartment_ledger <- function (units, donors) 
{
    grid <- expand.grid(donor = sort(donors), compartment = gse206785_spec()$compartments, KEEP.OUT.ATTRS = FALSE, 
        stringsAsFactors = FALSE)
    rows <- lapply(seq_len(nrow(grid)), function(i) {
        z <- units[units$donor == grid$donor[[i]] & units$compartment == grid$compartment[[i]], , drop = FALSE]
        one <- function(tissue, field, default) {
            value <- z[z$tissue == tissue, field]
            if (length(value)) 
                value[[1L]]
            else default
        }
        normal_eligible <- isTRUE(one("Normal", "eligible", FALSE))
        tumor_eligible <- isTRUE(one("Tumor", "eligible", FALSE))
        data.frame(donor = grid$donor[[i]], compartment = grid$compartment[[i]], normal_retained_cells = as.integer(one("Normal", 
            "cell_count", 0L)), tumor_retained_cells = as.integer(one("Tumor", "cell_count", 0L)), normal_total_umi = as.numeric(one("Normal", 
            "total_umi", 0)), tumor_total_umi = as.numeric(one("Tumor", "total_umi", 0)), normal_unit_eligible = normal_eligible, 
            tumor_unit_eligible = tumor_eligible, complete_pair_eligible = normal_eligible && tumor_eligible, 
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

gse206785_composition_input <- function (metadata, flags, donors) 
{
    retained <- metadata[flags$seven_compartment_pass, c("Patient", "Tissue"), drop = FALSE]
    observed <- as.data.frame(table(factor(retained$Patient, levels = sort(donors)), factor(retained$Tissue, 
        levels = c("Normal", "Tumor"))), stringsAsFactors = FALSE)
    names(observed) <- c("donor", "tissue", "retained_cells")
    observed$donor <- as.character(observed$donor)
    observed$tissue <- as.character(observed$tissue)
    observed$retained_cells <- as.integer(observed$retained_cells)
    observed$chemistry_eligible <- TRUE
    observed
}

gse206785_build_composition_eligibility <- function (counts) 
{
    required <- c("donor", "tissue", "retained_cells", "chemistry_eligible")
    if (!all(required %in% names(counts))) {
        gse206785_abort("composition input schema is incomplete")
    }
    if (anyDuplicated(counts[c("donor", "tissue")])) {
        gse206785_abort("duplicate donor-tissue composition denominator")
    }
    donors <- sort(unique(counts$donor))
    out <- lapply(donors, function(donor) {
        z <- counts[counts$donor == donor, , drop = FALSE]
        pair_valid <- nrow(z) == 2L && setequal(z$tissue, c("Normal", "Tumor")) && all(is.finite(z$retained_cells)) && 
            all(z$retained_cells >= 0)
        chemistry <- pair_valid && all(z$chemistry_eligible)
        minimum <- if (pair_valid) 
            min(z$retained_cells)
        else NA_integer_
        primary <- chemistry && minimum >= 100L
        strict <- chemistry && minimum >= 200L
        count_aware <- chemistry
        data.frame(donor = donor, minimum_retained_cells_across_tissues = minimum, primary_100 = primary, 
            strict_200 = strict, count_aware = count_aware, primary_reason = if (primary) 
                "included"
            else "fewer_than_100_retained_cells_in_at_least_one_tissue", strict_reason = if (strict) 
                "included"
            else "fewer_than_200_retained_cells_in_at_least_one_tissue", count_aware_reason = if (count_aware) 
                "included"
            else "invalid_or_chemistry_ineligible_pair", stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

gse206785_write_tsv <- function (x, path) 
{
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA", eol = "\n", 
        fileEncoding = "UTF-8")
    invisible(path)
}

gse206785_cell_universe_table <- function (observed) 
{
    expected <- gse206785_spec()$expected$fixed_universe
    data.frame(stage = c("numeric_qc_pass", "subtype_only_pass", "numeric_subtype_intersection", "mural_glial_exclusions", 
        "seven_compartment_universe"), expected = as.integer(expected), observed = as.integer(observed), 
        status = ifelse(as.integer(expected) == as.integer(observed), "pass", "fail"), stringsAsFactors = FALSE)
}

gse206785_validate_raw_fragment_manifest <- function (manifest, contract = gse206785_fixed_commit_contract()) 
{
    required <- c("file_name", "bytes", "sha256")
    if (!is.data.frame(manifest) || !all(required %in% names(manifest))) {
        gse206785_abort("raw-UMI fragment manifest schema is incomplete")
    }
    expected <- contract$files[match(contract$raw_fragment_names, contract$files$file_name), required, 
        drop = FALSE]
    observed <- manifest[, required, drop = FALSE]
    if (nrow(observed) != 9L || anyDuplicated(observed$file_name) || !identical(observed$file_name, expected$file_name)) {
        gse206785_abort("formal input requires exactly nine ordered raw-UMI fragments")
    }
    if (!identical(as.numeric(observed$bytes), as.numeric(expected$bytes)) || !identical(as.character(observed$sha256), 
        expected$sha256)) {
        gse206785_abort("raw-UMI fragment byte count or SHA-256 mismatch")
    }
    invisible(TRUE)
}

gse206785_concatenate_fragments <- function (paths, output_path, contract = gse206785_fixed_commit_contract()) 
{
    manifest <- data.frame(file_name = basename(paths), bytes = as.numeric(file.info(paths)$size), sha256 = vapply(paths, 
        gse206785_sha256, character(1)), stringsAsFactors = FALSE)
    gse206785_validate_raw_fragment_manifest(manifest, contract)
    validate_combined <- function(path) {
        file.exists(path) && identical(as.numeric(file.info(path)$size), contract$combined_bytes) && 
            identical(gse206785_sha256(path), contract$combined_sha256)
    }
    if (file.exists(output_path)) {
        if (!validate_combined(output_path)) {
            gse206785_abort("existing combined raw-UMI derivative is stale or invalid")
        }
        return(normalizePath(output_path, mustWork = TRUE))
    }
    partial <- paste0(output_path, ".part.", Sys.getpid())
    if (file.exists(partial)) {
        gse206785_abort("unexpected combined raw-UMI partial derivative exists")
    }
    to <- file(partial, open = "wb")
    promoted <- FALSE
    on.exit({
        try(close(to), silent = TRUE)
        if (!promoted && file.exists(partial)) unlink(partial)
    }, add = TRUE)
    for (path in paths) {
        from <- file(path, open = "rb")
        repeat {
            chunk <- readBin(from, what = "raw", n = 8L * 1024L * 1024L)
            if (!length(chunk)) 
                break
            writeBin(chunk, to, useBytes = TRUE)
        }
        close(from)
    }
    close(to)
    if (!validate_combined(partial)) {
        gse206785_abort("combined raw-UMI byte-chain identity mismatch")
    }
    if (!file.rename(partial, output_path)) {
        gse206785_abort("combined raw-UMI derivative promotion failed")
    }
    promoted <- TRUE
    normalizePath(output_path, mustWork = TRUE)
}

gse206785_validate_formal_matrix_source <- function (repository_commit, source_name) 
{
    contract <- gse206785_fixed_commit_contract()
    valid <- identical(repository_commit, contract$commit) && source_name %in% contract$raw_fragment_names
    if (!valid) {
        gse206785_abort(paste("GSE206785_scgex.txt.gz is supplementary processed expression and is", 
            "not used for formal raw UMI QC, pseudobulk, or member detection"))
    }
    invisible(TRUE)
}

gse206785_parse_count_gzip <- function (path) 
{
    gse206785_abort("processed or non-integer matrix is not formal raw UMI")
}

gse206785_validate_eligibility_matrix <- function (matrix, unit_counts, detection, chemistry_eligible_donors) 
{
    expected_keys <- data.frame(compartment = rep(gse206785_spec()$compartments, each = 2L), program = rep(c("G", 
        "H"), times = 7L), stringsAsFactors = FALSE)
    if (nrow(matrix) != 14L || !identical(matrix[c("compartment", "program")], expected_keys)) {
        gse206785_abort("eligibility matrix must contain exactly 14 ordered rows")
    }
    if (any(matrix$marker_audit_status != "NOT_RUN_IN_STAGE_C2")) {
        gse206785_abort("marker audit must remain NOT_RUN_IN_STAGE_C2")
    }
    recomputed <- gse206785_build_eligibility_matrix(unit_counts, detection, chemistry_eligible_donors)
    if (!identical(matrix, recomputed)) {
        gse206785_abort("eligibility matrix differs from recomputed roles or facts")
    }
    invisible(TRUE)
}

gse206785_validate_composition_eligibility <- function (ledger, counts) 
{
    expected <- gse206785_build_composition_eligibility(counts)
    if (!identical(ledger, expected)) {
        gse206785_abort("composition-denominator drift from recomputation")
    }
    invisible(TRUE)
}

gse206785_frozen_qc_rules <- function () 
{
    list(detected_genes_min_inclusive = 200L, detected_genes_max_inclusive = 7500L, umi_min_inclusive = 500, 
        umi_max_inclusive = 75000, mitochondrial_fraction_max_exclusive = 0.10000000000000001, low_quality_subtypes = c("F04", 
            "F15", "B03", "B08", "T03", "T04", "T05", "T16", "T21", "EN08", "EN09", "EN11", "EN13", "EP12", 
            "M15"), descriptive_exclusions = c("Mural", "Glial"))
}

gse206785_frozen_type_mapping <- function () 
{
    c(Epithelial = "epithelial", Fibroblast = "fibroblast", Endothelial = "endothelial", Myeloid = "myeloid", 
        `CD8+ T` = "T_NK", `CD4+ T` = "T_NK", `Innate lymphoid` = "T_NK", B = "B_plasma", Plasma = "B_plasma", 
        Mast = "mast_other_immune", Mural = NA_character_, Glial = NA_character_)
}

gse206785_abort <- function (message) 
{
    stop(message, call. = FALSE)
}

gse206785_soft_title_to_sample <- function (title) 
{
    pattern <- "^[A-Za-z][0-9]+_[A-Za-z]([0-9]{6})_([NT])$"
    matched <- grepl(pattern, title, perl = TRUE)
    out <- title
    out[matched] <- sub(pattern, "\\1\\2", title[matched], perl = TRUE)
    out
}

gse206785_compile_raw_umi_parser <- function () 
{
    if (exists("gse206785_parse_raw_umi_cpp", envir = .gse206785_raw_umi_parser_environment, inherits = FALSE)) {
        return(invisible(TRUE))
    }
    code <- "\n    #include <Rcpp.h>\n    #include <zlib.h>\n    #include <unordered_map>\n    #include <unordered_set>\n    #include <vector>\n    #include <string>\n    #include <limits>\n    #include <cctype>\n    #include <cstdint>\n    #include <fstream>\n    using namespace Rcpp;\n\n    static std::vector<std::string> raw_split_csv(const std::string& line) {\n      std::vector<std::string> fields;\n      size_t start = 0;\n      for (size_t i = 0; i <= line.size(); ++i) {\n        if (i == line.size() || line[i] == \",\"[0]) {\n          fields.push_back(line.substr(start, i - start));\n          start = i + 1;\n        }\n      }\n      return fields;\n    }\n\n    static unsigned long long raw_strict_uint(const std::string& token) {\n      if (token.empty()) stop(\"raw UMI token must be a strict nonnegative integer\");\n      unsigned long long value = 0;\n      for (size_t i = 0; i < token.size(); ++i) {\n        unsigned char ch = static_cast<unsigned char>(token[i]);\n        if (!std::isdigit(ch)) {\n          stop(\"raw UMI token must be a strict nonnegative integer\");\n        }\n        unsigned int digit = token[i] - \"0\"[0];\n        if (value > (std::numeric_limits<unsigned long long>::max() - digit) / 10) {\n          stop(\"raw UMI integer token overflow\");\n        }\n        value = value * 10 + digit;\n      }\n      if (value > 9007199254740991ULL) stop(\"raw UMI integer token overflow\");\n      return value;\n    }\n\n    static bool raw_mt_prefix(const std::string& gene) {\n      return gene.size() >= 3 &&\n        std::toupper(static_cast<unsigned char>(gene[0])) == \"M\"[0] &&\n        std::toupper(static_cast<unsigned char>(gene[1])) == \"T\"[0] &&\n        gene[2] == \"-\"[0];\n    }\n\n    // [[Rcpp::export]]\n    Rcpp::List gse206785_parse_raw_umi_cpp(\n        std::string path, int expected_genes, int expected_cells,\n        Rcpp::CharacterVector core_members) {\n      gzFile handle = gzopen(path.c_str(), \"rb\");\n      if (handle == NULL) stop(\"unable to open raw-UMI gzip\");\n      const int buffer_size = 16 * 1024 * 1024;\n      std::vector<char> buffer(buffer_size);\n      char* got = gzgets(handle, buffer.data(), buffer_size);\n      if (got == NULL) {\n        gzclose(handle);\n        stop(\"raw-UMI gzip is empty\");\n      }\n      std::string header(buffer.data());\n      if (header.empty() || header.back() != \"\\n\"[0]) {\n        gzclose(handle);\n        stop(\"raw-UMI header exceeds bounded parser buffer or is truncated\");\n      }\n      header.pop_back();\n      if (!header.empty() && header.back() == \"\\r\"[0]) header.pop_back();\n      std::vector<std::string> fields = raw_split_csv(header);\n      if (fields.size() != static_cast<size_t>(expected_cells + 1) ||\n          fields[0] != \"Gene\") {\n        gzclose(handle);\n        stop(\"raw-UMI orientation must be gene rows by cell columns\");\n      }\n      std::vector<std::string> cell_ids(fields.begin() + 1, fields.end());\n      std::unordered_set<std::string> cell_seen;\n      for (size_t i = 0; i < cell_ids.size(); ++i) {\n        if (cell_ids[i].empty() || !cell_seen.insert(cell_ids[i]).second) {\n          gzclose(handle);\n          stop(\"raw-UMI header contains duplicate or missing cell identifiers\");\n        }\n      }\n      std::unordered_map<std::string, int> core_lookup;\n      for (int i = 0; i < core_members.size(); ++i) {\n        core_lookup[as<std::string>(core_members[i])] = i;\n      }\n      std::vector<bool> core_found(core_members.size(), false);\n      std::vector<int> detected(expected_cells, 0);\n      std::vector<double> totals(expected_cells, 0);\n      std::vector<double> mitochondrial(expected_cells, 0);\n      NumericMatrix core_matrix(expected_cells, core_members.size());\n      std::vector<std::string> gene_names;\n      gene_names.reserve(expected_genes);\n      std::unordered_set<std::string> gene_seen;\n      int header_count = 1;\n      while ((got = gzgets(handle, buffer.data(), buffer_size)) != NULL) {\n        std::string line(buffer.data());\n        bool ended = !line.empty() && line.back() == \"\\n\"[0];\n        if (!ended && !gzeof(handle)) {\n          gzclose(handle);\n          stop(\"raw-UMI row exceeds bounded parser buffer\");\n        }\n        if (ended) line.pop_back();\n        if (!line.empty() && line.back() == \"\\r\"[0]) line.pop_back();\n        if (line.empty()) {\n          gzclose(handle);\n          stop(\"raw-UMI file contains an empty row\");\n        }\n        std::vector<std::string> row = raw_split_csv(line);\n        if (row.size() != static_cast<size_t>(expected_cells + 1)) {\n          gzclose(handle);\n          stop(\"raw-UMI row field count mismatch\");\n        }\n        std::string gene = row[0];\n        if (gene == \"Gene\") {\n          header_count += 1;\n          gzclose(handle);\n          stop(\"raw-UMI fragment chain contains a repeated logical header\");\n        }\n        if (gene.empty() || !gene_seen.insert(gene).second) {\n          gzclose(handle);\n          stop(\"raw-UMI file contains duplicate or missing gene identifiers\");\n        }\n        gene_names.push_back(gene);\n        bool mitochondrial_gene = raw_mt_prefix(gene);\n        int core_index = -1;\n        std::unordered_map<std::string, int>::const_iterator it =\n          core_lookup.find(gene);\n        if (it != core_lookup.end()) {\n          core_index = it->second;\n          core_found[core_index] = true;\n        }\n        for (int cell = 0; cell < expected_cells; ++cell) {\n          unsigned long long value = raw_strict_uint(row[cell + 1]);\n          double numeric_value = static_cast<double>(value);\n          if (value > 0) detected[cell] += 1;\n          totals[cell] += numeric_value;\n          if (!R_finite(totals[cell])) {\n            gzclose(handle);\n            stop(\"raw-UMI cell total overflow\");\n          }\n          if (mitochondrial_gene) mitochondrial[cell] += numeric_value;\n          if (core_index >= 0) core_matrix(cell, core_index) = numeric_value;\n        }\n      }\n      int close_status = gzclose(handle);\n      if (close_status != Z_OK) stop(\"raw-UMI gzip integrity failure\");\n      if (static_cast<int>(gene_names.size()) != expected_genes) {\n        stop(\"raw-UMI gene-row count mismatch\");\n      }\n      for (size_t i = 0; i < core_found.size(); ++i) {\n        if (!core_found[i]) stop(\"raw-UMI file is missing a frozen core member\");\n      }\n      NumericVector mt_fraction(expected_cells);\n      for (int i = 0; i < expected_cells; ++i) {\n        mt_fraction[i] = totals[i] > 0 ? mitochondrial[i] / totals[i] : 0;\n      }\n      colnames(core_matrix) = core_members;\n      return List::create(\n        _[\"gene_names\"] = wrap(gene_names),\n        _[\"cell_ids\"] = wrap(cell_ids),\n        _[\"gene_count\"] = static_cast<int>(gene_names.size()),\n        _[\"cell_count\"] = expected_cells,\n        _[\"header_count\"] = header_count,\n        _[\"detected_genes\"] = wrap(detected),\n        _[\"total_umi\"] = wrap(totals),\n        _[\"mitochondrial_fraction\"] = mt_fraction,\n        _[\"core_counts\"] = core_matrix\n      );\n    }\n\n    // [[Rcpp::export]]\n    Rcpp::List gse206785_aggregate_raw_umi_cpp(\n        std::string path, int expected_genes, int expected_cells,\n        Rcpp::IntegerMatrix group_index, int group_count) {\n      int scheme_count = group_index.ncol();\n      if (group_index.nrow() != expected_cells || scheme_count < 1 ||\n          group_count < 1) {\n        stop(\"raw-UMI aggregation group contract mismatch\");\n      }\n      for (int i = 0; i < expected_cells; ++i) {\n        for (int scheme = 0; scheme < scheme_count; ++scheme) {\n          if (group_index(i, scheme) < -1 ||\n              group_index(i, scheme) >= group_count) {\n            stop(\"raw-UMI aggregation group index is invalid\");\n          }\n        }\n      }\n      gzFile handle = gzopen(path.c_str(), \"rb\");\n      if (handle == NULL) stop(\"unable to open raw-UMI gzip for aggregation\");\n      const int buffer_size = 16 * 1024 * 1024;\n      std::vector<char> buffer(buffer_size);\n      char* got = gzgets(handle, buffer.data(), buffer_size);\n      if (got == NULL) {\n        gzclose(handle);\n        stop(\"raw-UMI gzip is empty during aggregation\");\n      }\n      std::string header(buffer.data());\n      if (header.empty() || header.back() != \"\\n\"[0]) {\n        gzclose(handle);\n        stop(\"raw-UMI aggregation header is truncated\");\n      }\n      header.pop_back();\n      if (!header.empty() && header.back() == \"\\r\"[0]) header.pop_back();\n      std::vector<std::string> fields = raw_split_csv(header);\n      if (fields.size() != static_cast<size_t>(expected_cells + 1) ||\n          fields[0] != \"Gene\") {\n        gzclose(handle);\n        stop(\"raw-UMI aggregation orientation mismatch\");\n      }\n      NumericMatrix aggregated(expected_genes, group_count * scheme_count);\n      CharacterVector genes(expected_genes);\n      std::unordered_set<std::string> gene_seen;\n      int gene_index = 0;\n      while ((got = gzgets(handle, buffer.data(), buffer_size)) != NULL) {\n        std::string line(buffer.data());\n        bool ended = !line.empty() && line.back() == \"\\n\"[0];\n        if (!ended && !gzeof(handle)) {\n          gzclose(handle);\n          stop(\"raw-UMI aggregation row exceeds bounded parser buffer\");\n        }\n        if (ended) line.pop_back();\n        if (!line.empty() && line.back() == \"\\r\"[0]) line.pop_back();\n        std::vector<std::string> row = raw_split_csv(line);\n        if (row.size() != static_cast<size_t>(expected_cells + 1) ||\n            gene_index >= expected_genes) {\n          gzclose(handle);\n          stop(\"raw-UMI aggregation row contract mismatch\");\n        }\n        if (row[0].empty() || !gene_seen.insert(row[0]).second) {\n          gzclose(handle);\n          stop(\"raw-UMI aggregation gene identity mismatch\");\n        }\n        genes[gene_index] = row[0];\n        for (int cell = 0; cell < expected_cells; ++cell) {\n          unsigned long long value = raw_strict_uint(row[cell + 1]);\n          for (int scheme = 0; scheme < scheme_count; ++scheme) {\n            int group = group_index(cell, scheme);\n            if (group >= 0) {\n              int column = scheme * group_count + group;\n              aggregated(gene_index, column) += static_cast<double>(value);\n              if (!R_finite(aggregated(gene_index, column))) {\n                gzclose(handle);\n                stop(\"raw-UMI pseudobulk aggregation overflow\");\n              }\n            }\n          }\n        }\n        gene_index += 1;\n      }\n      int close_status = gzclose(handle);\n      if (close_status != Z_OK) stop(\"raw-UMI aggregation gzip integrity failure\");\n      if (gene_index != expected_genes) {\n        stop(\"raw-UMI aggregation gene-row count mismatch\");\n      }\n      return List::create(\n        _[\"genes\"] = genes,\n        _[\"counts\"] = aggregated\n      );\n    }\n\n    // [[Rcpp::export]]\n    Rcpp::List gse206785_write_canonical_pseudobulk_preimage_cpp(\n        std::string raw_path, std::string preimage_path,\n        int expected_genes, int expected_cells,\n        Rcpp::IntegerVector group_index, int group_count) {\n      if (group_index.size() != expected_cells || group_count < 1) {\n        stop(\"canonical pseudobulk grouping contract mismatch\");\n      }\n      for (int cell = 0; cell < expected_cells; ++cell) {\n        if (group_index[cell] < -1 || group_index[cell] >= group_count) {\n          stop(\"canonical pseudobulk group index is invalid\");\n        }\n      }\n      gzFile handle = gzopen(raw_path.c_str(), \"rb\");\n      if (handle == NULL) stop(\"unable to open raw-UMI gzip for canonical hash\");\n      std::ofstream output(\n        preimage_path.c_str(), std::ios::binary | std::ios::trunc\n      );\n      if (!output.is_open()) {\n        gzclose(handle);\n        stop(\"unable to create canonical hash preimage\");\n      }\n      const int buffer_size = 16 * 1024 * 1024;\n      std::vector<char> buffer(buffer_size);\n      char* got = gzgets(handle, buffer.data(), buffer_size);\n      if (got == NULL) {\n        gzclose(handle);\n        stop(\"raw-UMI gzip is empty during canonical hash\");\n      }\n      std::string header(buffer.data());\n      if (header.empty() || header.back() != \"\\n\"[0]) {\n        gzclose(handle);\n        stop(\"canonical hash header is truncated\");\n      }\n      header.pop_back();\n      if (!header.empty() && header.back() == \"\\r\"[0]) header.pop_back();\n      std::vector<std::string> fields = raw_split_csv(header);\n      if (fields.size() != static_cast<size_t>(expected_cells + 1) ||\n          fields[0] != \"Gene\") {\n        gzclose(handle);\n        stop(\"canonical hash raw-UMI orientation mismatch\");\n      }\n      std::unordered_set<std::string> gene_seen;\n      int gene_count = 0;\n      unsigned long long bytes_written = 0;\n      const unsigned long long signed_max =\n        static_cast<unsigned long long>(std::numeric_limits<int64_t>::max());\n      while ((got = gzgets(handle, buffer.data(), buffer_size)) != NULL) {\n        std::string line(buffer.data());\n        bool ended = !line.empty() && line.back() == \"\\n\"[0];\n        if (!ended && !gzeof(handle)) {\n          gzclose(handle);\n          stop(\"canonical hash raw-UMI row exceeds bounded buffer\");\n        }\n        if (ended) line.pop_back();\n        if (!line.empty() && line.back() == \"\\r\"[0]) line.pop_back();\n        std::vector<std::string> row = raw_split_csv(line);\n        if (row.size() != static_cast<size_t>(expected_cells + 1) ||\n            gene_count >= expected_genes) {\n          gzclose(handle);\n          stop(\"canonical hash raw-UMI row contract mismatch\");\n        }\n        const std::string& gene = row[0];\n        if (gene.empty() || gene == \"Gene\" ||\n            !gene_seen.insert(gene).second) {\n          gzclose(handle);\n          stop(\"canonical hash gene identity mismatch\");\n        }\n        std::vector<unsigned long long> sums(group_count, 0);\n        for (int cell = 0; cell < expected_cells; ++cell) {\n          unsigned long long value = raw_strict_uint(row[cell + 1]);\n          int group = group_index[cell];\n          if (group >= 0) {\n            if (value > signed_max - sums[group]) {\n              gzclose(handle);\n              stop(\"canonical pseudobulk signed int64 overflow\");\n            }\n            sums[group] += value;\n          }\n        }\n        output.write(gene.data(), static_cast<std::streamsize>(gene.size()));\n        output.put(static_cast<char>(0));\n        for (int group = 0; group < group_count; ++group) {\n          uint64_t value = static_cast<uint64_t>(sums[group]);\n          for (int byte = 0; byte < 8; ++byte) {\n            output.put(static_cast<char>((value >> (8 * byte)) & 0xff));\n          }\n        }\n        if (!output.good()) {\n          gzclose(handle);\n          stop(\"canonical hash preimage write failed\");\n        }\n        bytes_written += static_cast<unsigned long long>(gene.size()) + 1ULL +\n          static_cast<unsigned long long>(group_count) * 8ULL;\n        gene_count += 1;\n      }\n      int close_status = gzclose(handle);\n      output.close();\n      if (close_status != Z_OK) {\n        stop(\"canonical hash raw-UMI gzip integrity failure\");\n      }\n      if (output.fail()) stop(\"canonical hash preimage close failed\");\n      if (gene_count != expected_genes) {\n        stop(\"canonical hash gene-row count mismatch\");\n      }\n      return List::create(\n        _[\"matrix_rows\"] = gene_count,\n        _[\"matrix_columns\"] = group_count,\n        _[\"integer_nonnegative\"] = true,\n        _[\"preimage_bytes\"] = static_cast<double>(bytes_written)\n      );\n    }\n  "
    old_pkg_libs <- Sys.getenv("PKG_LIBS", unset = NA_character_)
    on.exit({
        if (is.na(old_pkg_libs)) Sys.unsetenv("PKG_LIBS") else Sys.setenv(PKG_LIBS = old_pkg_libs)
    }, add = TRUE)
    Sys.setenv(PKG_LIBS = paste(Sys.getenv("PKG_LIBS"), "-lz"))
    Rcpp::sourceCpp(code = code, env = .gse206785_raw_umi_parser_environment, rebuild = FALSE, showOutput = FALSE, 
        verbose = FALSE)
    invisible(TRUE)
}

gse206785_fixed_commit_contract <- function () 
{
    file_name <- c("cell_metadata.csv", "README.md", sprintf("raw_counts.csv.gz%02d", 0:8))
    bytes <- c(4375642, 1150, rep(20000000, 8L), 14415491)
    sha256 <- c("d8936a22460b27348c7eb966b04f6a492ed9a05163ec7f663b2a02f3e91644be", "2d3b47fb920db7e068e6d19141a73ff237abc9f307b6d92d07d83c89b18ed0bd", 
        "832b05fc400af6194d330ef14dd1f545bab260be5ac38987f6b30d7a7ec4bebc", "adddaf9d2aa0da9b75bb7506bf5f20ae61e0421e9f2ec72d9646dcbe70652441", 
        "d4742a66a9c384880e508f88bcbf31c0862d7b8cb83bdcf8a39fc7603cce3fa7", "649464f7ac59f68da0ac4c1c2892ef06c0117a673795a4c88a50a26e3bb5a09a", 
        "902cb9723684deffe18d1b89470742371fafeeb49fe1dc2248b3c511d2087445", "1cec8d41d51cd24744330349980a4f2255ac3736dc5ef5ab5dc3d1dc398a466d", 
        "8f45232cc443ac667f06d62f7911891744cbf94296cae221589158a7f20708b6", "1b223bd271c671bc5910c6e904e278ad2a8c8da2787ce9e76a9a676701c9175b", 
        "5b1650d052c6f1802ed33ccd5538fea7092dcfb367cbc0b5be3da80874953117")
    list(repository = "chriscainx/gastric-cancer", commit = "dd0e37b1ab5dcf41caa2bce5eea08aa9446446da", 
        directory = "00.Processed Data", base_url = paste0("https://raw.githubusercontent.com/chriscainx/gastric-cancer/", 
            "dd0e37b1ab5dcf41caa2bce5eea08aa9446446da/", "00.Processed%20Data"), files = data.frame(file_name = file_name, 
            bytes = as.numeric(bytes), sha256 = sha256, stringsAsFactors = FALSE), raw_fragment_names = sprintf("raw_counts.csv.gz%02d", 
            0:8), combined_bytes = 174415491, combined_sha256 = "016d71dc669b44d82db49575686e5898db4f124af77179292ed85702e4a010a9", 
        combined_derivative_name = "raw_counts.csv.gz.combined_derivative")
}

gse206785_validate_qc_rules <- function (rules) 
{
    if (!identical(rules, gse206785_frozen_qc_rules())) {
        if (!identical(rules$low_quality_subtypes, gse206785_frozen_qc_rules()$low_quality_subtypes)) {
            gse206785_abort("low-quality subtype-list drift")
        }
        gse206785_abort("numeric QC-rule drift")
    }
    invisible(TRUE)
}

gse206785_validate_type_mapping <- function (mapping) 
{
    expected <- gse206785_frozen_type_mapping()
    if (!is.na(mapping[["Mural"]]) || !is.na(mapping[["Glial"]])) {
        bad <- if (!is.na(mapping[["Mural"]])) 
            "Mural"
        else "Glial"
        gse206785_abort(paste(bad, "relabeling is prohibited"))
    }
    if (!identical(mapping, expected)) {
        gse206785_abort("author broad-compartment mapping drift")
    }
    invisible(TRUE)
}

gse206785_frozen_unit_rules <- function () 
{
    gse206785_spec()$unit_rules
}

gse206785_validate_unit_rules <- function (rules) 
{
    if (!identical(rules, gse206785_frozen_unit_rules())) {
        gse206785_abort("pseudobulk eligibility-rule drift")
    }
    invisible(TRUE)
}

gse206785_canonical_hash_contract <- function () 
{
    list(reference_implementation = "run_sc_result_blind_stream.py", reference_implementation_sha256 = "7b08c149294e81d30020e37fe74b3a5c78a62ca913979a093655758019d17382", 
        reference_state = "result-blind Stage 8 technical audit", compartment = "fibroblast", cell_selection_rule = paste0("final_equals_numeric_qc_pass_AND_subtype_pass_AND_mapped_to_one_of_", 
            "seven_compartments_AND_Type_equals_Fibroblast"), donor_ordering_rule = "ascending_bytewise_ordinary_character", 
        tissue_ordering_rule = "Normal_then_Tumor", gene_ordering_rule = "original_raw_umi_matrix_row_order", 
        gene_symbol_encoding = "UTF-8_then_single_NUL_byte", integer_encoding = "signed_int64_little_endian", 
        row_delimiter = "none", column_metadata_in_preimage = "none", final_terminator = "none", expected_digest = "1500b3bf1e4d8ec87b1a8dc892ce6a119cfb304e4f80e521435b8b3a00454349", 
        matrix_rows = 25287L, matrix_columns = 48L, matrix_persisted = FALSE)
}

gse206785_sha256 <- function (path) 
{
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

gse206785_validate_canonical_hash_contract <- function (contract) 
{
    if (!identical(contract, gse206785_canonical_hash_contract())) {
        gse206785_abort("canonical pseudobulk hash contract drift")
    }
    invisible(TRUE)
}

gse206785_plan_role <- function (compartment, program) 
{
    key <- paste(compartment, program, sep = "_")
    if (key %in% c("myeloid_G", "myeloid_H", "fibroblast_H", "B_plasma_H")) {
        return("formal_candidate")
    }
    if (identical(key, "epithelial_H")) {
        return("prespecified_exploratory_candidate")
    }
    "descriptive_if_estimable_otherwise_non_estimable"
}

gse206785_program_coverage <- function (program, available_members) 
{
    available_members <- unique(available_members)
    if (identical(program, "G")) {
        members <- gse206785_spec()$programs$G
        available_members <- intersect(available_members, members)
        return(length(available_members) >= 5L && any(c("ACSL4", "LPCAT3") %in% available_members) && 
            sum(c("ALOX5", "ALOX12", "ALOX15", "ALOX15B", "POR") %in% available_members) >= 2L)
    }
    if (identical(program, "H")) {
        members <- gse206785_spec()$programs$H
        available_members <- intersect(available_members, members)
        return(length(available_members) >= 3L && any(c("ALDH3A1", "ALDH3A2") %in% available_members) && 
            any(c("ALDH2", "ALDH7A1") %in% available_members))
    }
    gse206785_abort("unsupported program")
}
