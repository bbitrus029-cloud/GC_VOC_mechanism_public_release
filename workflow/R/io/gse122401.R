# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.

gse122401_frozen_members <- function () 
{
    data.frame(program = c(rep("G-core", 7L), rep("H-core", 4L)), gene_symbol = c("ACSL4", "LPCAT3", 
        "ALOX5", "ALOX12", "ALOX15", "ALOX15B", "POR", "ALDH2", "ALDH3A1", "ALDH3A2", "ALDH7A1"), ensembl_gene_id = c("ENSG00000068366", 
        "ENSG00000111684", "ENSG00000012779", "ENSG00000108839", "ENSG00000161905", "ENSG00000179593", 
        "ENSG00000127948", "ENSG00000111275", "ENSG00000108602", "ENSG00000072210", "ENSG00000164904"), 
        stringsAsFactors = FALSE)
}

gse122401_numeric_token <- function (x) 
{
    grepl("^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][+-]?[0-9]+)?$", x) | x %in% c("NA", "NaN", "Inf", 
        "+Inf", "-Inf")
}

gse122401_inspect_processed_schemas <- function (paths, expected_samples = 160L, uq_normalization_documented, rsem_documented, exact_normalization_command_documented = TRUE, 
    documented_feature_count = 34908L, expected_feature_count = 34908L) 
{
    if (length(paths) != expected_samples || any(!file.exists(paths))) {
        stop("All formal processed files are required for schema inspection", call. = FALSE)
    }
    if (!isTRUE(uq_normalization_documented) || !isTRUE(rsem_documented) || !isTRUE(exact_normalization_command_documented) || 
        !identical(as.integer(documented_feature_count), as.integer(expected_feature_count)) || length(expected_feature_count) != 
        1L || is.na(expected_feature_count) || expected_feature_count < 1L) {
        stop("The exact RSEM UQ processed-expression documentation is unresolved", call. = FALSE)
    }
    signatures <- vapply(paths, function(path) {
        con <- gzfile(path, open = "rt")
        lines <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
        if (length(lines) != expected_feature_count + 1L) {
            stop("A processed gene table must contain exactly 34,908 valid data rows after one pseudoheader", 
                call. = FALSE)
        }
        fields <- strsplit(lines, "\t", fixed = TRUE)
        if (any(lengths(fields) != 2L)) {
            stop("Every processed row must contain exactly two positional fields", call. = FALSE)
        }
        first <- fields[[1L]]
        data_fields <- fields[-1L]
        ids <- vapply(data_fields, `[[`, character(1), 1L)
        tokens <- vapply(data_fields, `[[`, character(1), 2L)
        if (length(first) != 2L || !identical(first[[1L]], "gene_id") || !gse122401_numeric_token(first[[2L]]) || 
            any(!grepl("^ENSG[0-9]+([.][0-9]+)?$", ids))) {
            stop("The processed gene-table schema is not the frozen two-column RSEM UQ form", call. = FALSE)
        }
        if (any(!gse122401_numeric_token(tokens))) {
            stop("A processed gene table contains a nonnumeric processed token", call. = FALSE)
        }
        paste0("two_tab_fields|one_pseudoheader_gene_id_numeric_sentinel|", "ensembl_id|uq_rsem_abundance|feature_count=", 
            expected_feature_count)
    }, character(1))
    if (length(unique(signatures)) != 1L) {
        stop("All formal processed files must have an identical schema", call. = FALSE)
    }
    descriptor <- paste0("gene_id<TAB>uq_normalized_expression;skip_pseudoheader=1;", "rsem_command=exact_rsem_1.2.15;", 
        "normalization_command=quartile_norm.pl_-c_5_-q_75_-t_1000;", "feature_count=", expected_feature_count, 
        ";raw_counts=false")
    list(fields = c("gene_id", "uq_normalized_expression"), gene_column = "gene_id", expression_column = "uq_normalized_expression", 
        skip_rows = 1L, canonical_feature_count = as.integer(expected_feature_count), processed_input_documentation_status = paste0("exact_rsem_expression_calculation_quartile_norm_c5_q75_t1000_", 
            "and_34908_gene_statement_pass"), schema_hash = digest::digest(descriptor, algo = "sha256", 
            serialize = FALSE), processed_scale = "continuous_RSEM_upper_quartile_normalized_abundance_not_raw_counts")
}

gse122401_read_processed_gene_file <- function (path, schema) 
{
    required <- c("gene_column", "expression_column", "skip_rows")
    if (!all(required %in% names(schema)) || !file.exists(path)) {
        stop("Processed reader lacks its frozen schema or physical file", call. = FALSE)
    }
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    table <- utils::read.delim(con, header = FALSE, skip = schema$skip_rows, sep = "\t", quote = "", 
        comment.char = "", colClasses = "character", col.names = c(schema$gene_column, schema$expression_column), 
        check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(table) != 2L || !nrow(table)) {
        stop("Processed gene table is empty or does not have two columns", call. = FALSE)
    }
    if (!"canonical_feature_count" %in% names(schema) || nrow(table) != schema$canonical_feature_count) {
        stop("Processed gene table does not contain the frozen canonical feature count", call. = FALSE)
    }
    tokens <- table[[schema$expression_column]]
    if (any(!gse122401_numeric_token(tokens))) {
        stop("Processed expression column contains a non-numeric token", call. = FALSE)
    }
    table[[schema$expression_column]] <- suppressWarnings(as.numeric(tokens))
    table
}
