# Extracted maintained implementation; stable machine names and seeds retained.
pe_protein_ratio_identities <- function(columns, biospecimens) {
    s01_require(all(grepl("^Unshared Log T[0-9]+([.][12])?/N[0-9]+([.][12])?$", columns)), "invalid paired protein ratio label")
    s01_validate_identity_vector(columns, "protein ratio column")
    s01_validate_identity_vector(biospecimens$aliquot_submitter_id, "protein aliquot")
    labels <- sub("^Unshared Log ", "", columns)
    parts <- strsplit(labels, "/", fixed = TRUE)
    t <- vapply(parts, `[`, "", 1L)
    n <- vapply(parts, `[`, "", 2L)
    it <- match(sub("[.][12]$", "", t), biospecimens$aliquot_submitter_id)
    inat <- match(sub("[.][12]$", "", n), biospecimens$aliquot_submitter_id)
    s01_require(!anyNA(c(it, inat)), "unresolved protein aliquot")
    s01_require(all(biospecimens$case_submitter_id[it] == biospecimens$case_submitter_id[inat]) && all(biospecimens$sample_type[it] == 
        "Tumor") && all(biospecimens$sample_type[inat] == "Solid Tissue Normal"), "protein numerator denominator pairing")
    rt <- sub("^[^.]+", "", t)
    rn <- sub("^[^.]+", "", n)
    s01_require(identical(rt, rn), "protein replicate mismatch")
    data.frame(ratio_column = columns, patient_id = paste0("Pat", biospecimens$case_submitter_id[it]), tumor_aliquot = biospecimens$aliquot_submitter_id[it], 
        adjacent_aliquot = biospecimens$aliquot_submitter_id[inat], replicate = ifelse(rt == "", "single", rt))
}

pe_protein_patient_ratios <- function(x, identities) {
    s01_require(identical(colnames(x), identities$ratio_column), "protein ratio column ordering")
    s01_validate_identity_vector(rownames(x), "protein member")
    patients <- sort(unique(identities$patient_id), method = "radix")
    finite <- vapply(patients, function(p) all(is.finite(x[, identities$patient_id == p, drop = FALSE])), TRUE)
    complete <- patients[finite]
    out <- vapply(complete, function(p) rowMeans(x[, identities$patient_id == p, drop = FALSE]), numeric(nrow(x)))
    dimnames(out) <- list(rownames(x), complete)
    list(values = out, patient_ids = complete, excluded = patients[!finite], rule = "all_declared_replicates_of_all_H4_required_equal_mean_of_log2_ratios")
}

