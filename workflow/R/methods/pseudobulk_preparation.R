# Maintained reader/preparation functions extracted from accepted implementation.
# Stable function and randomization identities are deliberately retained.

s03_prepare_compartment <- function (counts, units, cfg) 
{
    s01_validate_numeric_matrix(counts, "preparation counts", TRUE)
    s01_assert_exact_order(units$unit_id, colnames(counts), "preparation unit")
    s01_assert_unique_key(units, c("donor", "tissue"), "preparation donor tissue")
    s01_require(all(units$tissue %in% c("Normal", "Tumor")), "preparation tissue identity")
    original <- colSums(counts)
    eligible <- units$cell_count >= cfg$minimum_cells & original >= cfg$minimum_original_UMI
    donors <- sort(names(which(table(units$donor[eligible]) == 2L)), method = "radix")
    j <- which(eligible & units$donor %in% donors)
    j <- j[order(units$donor[j], match(units$tissue[j], c("Normal", "Tumor")), method = "radix")]
    fraction <- if (any(eligible)) 
        rowMeans(counts[, eligible, drop = FALSE] > 0)
    else rep(0, nrow(counts))
    bg <- which(fraction >= cfg$detection_fraction)
    bg <- bg[order(rownames(counts)[bg], method = "radix")]
    modeling <- bg[rowSums(counts[bg, j, drop = FALSE]) > 0]
    mat <- counts[modeling, j, drop = FALSE]
    both <- vapply(seq_len(nrow(counts)), function(i) if (length(j)) 
        sum(vapply(split(counts[i, j] > 0, units$donor[j]), all, TRUE))
    else 0L, integer(1))
    design <- NULL
    if (length(donors) >= 2L) {
        design <- model.matrix(~tissue + donor, data = data.frame(tissue = factor(units$tissue[j], c("Normal", 
            "Tumor")), donor = factor(units$donor[j], levels = donors)))
        rownames(design) <- units$unit_id[j]
        s01_require(qr(design)$rank == ncol(design), "preparation design rank failure")
    }
    list(original_library_sizes = original, eligible = eligible, eligible_units = units$unit_id[eligible], 
        complete_pair_donors = donors, modeling_units = units$unit_id[j], background = rownames(counts)[bg], 
        modeling_background = rownames(counts)[modeling], counts = mat, modeling_library_sizes = colSums(mat), 
        sample_data = units[j, , drop = FALSE], design = design, detection = data.frame(gene_id = rownames(counts), 
            eligible_units = sum(eligible), nonzero_units = rowSums(counts[, eligible, drop = FALSE] > 
                0), nonzero_fraction = fraction, paired_both_nonzero = both, available = fraction >= 
                cfg$detection_fraction & both >= cfg$minimum_paired_donors_nonzero_both_tissues), model_fitted = FALSE, 
        normalization_run = FALSE)
}
