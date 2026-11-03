# Extracted maintained implementation; stable machine names and seeds retained.
epithelial_module <- function(context) {
  s12_root <- context$root
  s12_sha <- function(p) digest::digest(file=context$input(p),algo="sha256")
  s12_now <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")
  s12_rds <- function(x,p) context$write_rds(x,p)
  s12_json <- function(x,p) context$write_json(x,p)
  s12_text <- function(x,p) context$write_text(x,p)
  s12_guard <- function(...) context$guard()
  s12_effect_guard <- function(...) context$guard()
  s12_memory <- function(stage,minimum=8,prestart=FALSE) context$resource_check(minimum,stage)
  s12_hashes <- function(paths) data.frame(path=paths,bytes=file.info(paths)$size,sha256=vapply(paths,s12_sha,""))
  s12_atomic <- function(p,writer,check) context$atomic(p,writer,check)
  s12_m <- "data/manifests/post_effect_supplement/s12"
  s12_p <- context$output("processed"); s12_r <- context$output("results"); s12_q <- context$output("qc")
  s12_b <- "data/processed/post_effect_supplement/s10_b_preparation/prepare_01"
s12_read <- function(p) read.delim(context$input(p), stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")

s12_runtime <- function(...) context$runtime("epithelial")

s12_roster <- function() s12_read(file.path(s12_b, "complete_hypothesis_roster.tsv"))

s12_cfg <- function() yaml::read_yaml(context$input("config/post_effect_supplement.yml"))

s12_seed <- function(operation) {
    s <- s12_read("config/implementation/s12_seed_map.tsv")
    i <- match(operation, s$operation_name)
    s01_require(!is.na(i), "Unfrozen S12 seed")
    as.integer(s$seed[i])
}

s12_tsv <- function(x, p) {
    s12_atomic(p, function(f) {
        old <- options(digits = 17)
        on.exit(options(old))
        write.table(x, f, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
    }, function(f) {
        classes <- vapply(x, function(z) if (is.numeric(z)) 
            "numeric"
        else if (is.logical(z)) 
            "logical"
        else "character", "")
        z <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = classes)
        s01_require(identical(names(x), names(z)) && nrow(x) == nrow(z), "TSV schema readback")
        for (k in names(x)) s01_require(identical(is.na(x[[k]]), is.na(z[[k]])) && if (is.numeric(x[[k]])) 
            isTRUE(all.equal(as.numeric(x[[k]]), z[[k]], tolerance = 1e-13))
        else identical(as.character(x[[k]]), as.character(z[[k]])), paste("TSV readback", k))
    })
}

s12_members <- c("ALDH2", "ALDH3A1", "ALDH3A2", "ALDH7A1")
s12_bind <- function(p, raw) {
    s01_validate_numeric_matrix(raw$counts, "S12 original pseudobulk counts", TRUE)
    s01_validate_numeric_matrix(p$counts, "S12 paired modeling counts", TRUE)
    s01_validate_identity_vector(p$background, "S12 detection background")
    s01_validate_identity_vector(p$modeling_background, "S12 modeling background")
    s01_validate_identity_vector(p$modeling_units, "S12 model units")
    s01_validate_identity_vector(p$eligible_units, "S12 detection units")
    s01_validate_identity_vector(p$complete_pair_donors, "S12 donors")
    s01_require(setequal(rownames(p$counts), p$modeling_background) && setequal(colnames(p$counts), p$modeling_units), 
        "S12 frozen model dimensions/identities")
    s01_require(setequal(raw$units$unit_id, colnames(raw$counts)) && !anyDuplicated(raw$units$unit_id), "S12 original unit bridge")
    u <- raw$units[match(colnames(raw$counts), raw$units$unit_id), , drop = FALSE]
    s01_assert_unique_key(u, c("donor", "tissue"), "S12 raw donor/tissue")
    s01_require(all(u$tissue %in% c("Normal", "Tumor")), "S12 source tissue identity")
    cfg <- s12_cfg()$B
    totals <- colSums(raw$counts)
    s01_require(identical(totals[p$input_units$unit_id], p$original_library_sizes), "S12 original library totals")
    eligible <- u$cell_count >= cfg$minimum_cells & totals >= cfg$minimum_original_UMI
    s01_require(setequal(u$unit_id[eligible], p$eligible_units), "S12 detection-unit population")
    donors <- sort(names(which(table(u$donor[eligible]) == 2L)), method = "radix")
    s01_require(identical(donors, p$complete_pair_donors), "S12 complete-pair population")
    s01_require(all(p$background %in% rownames(raw$counts)), "S12 missing detection background")
    fraction <- rowMeans(raw$counts[, p$eligible_units, drop = FALSE] > 0)
    bg <- sort(rownames(raw$counts)[fraction >= cfg$detection_fraction], method = "radix")
    s01_require(identical(bg, p$background), "S12 all-eligible detection denominator/background")
    sc <- raw$counts[p$background, p$modeling_units, drop = FALSE]
    modeling <- p$background[rowSums(sc) > 0]
    s01_require(identical(modeling, p$modeling_background), "S12 all-zero modeling row rule")
    p$counts <- p$counts[p$modeling_background, p$modeling_units, drop = FALSE]
    s01_require(identical(p$counts, sc[p$modeling_background, , drop = FALSE]) && identical(colSums(p$counts), 
        p$modeling_library_sizes), "S12 full-background matrix/library totals")
    md <- p$sample_data
    s01_require(setequal(md$unit_id, p$modeling_units) && !anyDuplicated(md$unit_id), "S12 model metadata key")
    md <- md[match(p$modeling_units, md$unit_id), , drop = FALSE]
    j <- match(md$unit_id, u$unit_id)
    s01_require(identical(md$donor, u$donor[j]) && identical(md$tissue, u$tissue[j]), "S12 donor/tissue alignment")
    s01_require(identical(md$donor, rep(donors, each = 2L)) && identical(md$tissue, rep(c("Normal", "Tumor"), length(donors))), 
        "S12 frozen paired ordering")
    p$sample_data <- md
    b <- p$member_qualification
    s01_require(identical(b$member, s12_members) && all(b$mapping_state == "UNIQUE") && all(b$matches == 1L) && 
        !anyDuplicated(b$feature_id), "S12 unique frozen member bridge")
    f <- raw$feature_identity
    s01_require(setequal(f$feature_id, rownames(raw$counts)) && !anyDuplicated(f$feature_id), "S12 original feature identity")
    j <- match(b$feature_id, f$feature_id)
    s01_require(!anyNA(j) && identical(b$source_symbol, f$source_symbol[j]) && identical(as.integer(b$source_row), 
        as.integer(f$source_row[j])) && identical(b$member, f$canonical_symbol[j]), "S12 original feature bridge identity")
    di <- match(b$feature_id, p$detection$gene_id)
    s01_require(!anyNA(di), "S12 member detection identity")
    both <- vapply(b$feature_id, function(g) sum(raw$counts[g, p$modeling_units[seq(1L, ncol(sc), 2L)]] > 0 & raw$counts[g, 
        p$modeling_units[seq(2L, ncol(sc), 2L)]] > 0), 0L)
    ok <- fraction[b$feature_id] >= cfg$detection_fraction & both >= cfg$minimum_paired_donors_nonzero_both_tissues & 
        b$feature_id %in% p$modeling_background
    s01_require(identical(unname(b$technical_eligible), unname(ok & p$structural_reason == "NONE")) && all(p$detection$available[di] == 
        ok), "S12 frozen member eligibility")
    p$score_counts <- sc
    md$sample_id <- md$unit_id
    md$donor_id <- md$donor
    md$tissue <- ifelse(md$tissue == "Normal", "NAT", md$tissue)
    p$model_prepared <- s01_prepare_paired_model(p$counts, md, "donor_id", "donor", TRUE)
    p
}

s12_pairs <- function(p, values) {
    s01_require(identical(values$sample_id, p$modeling_units), "S12 score-unit alignment")
    list(ids = p$complete_pair_donors, n = values$score[seq(1L, length(values$score), 2L)], t = values$score[seq(2L, 
        length(values$score), 2L)])
}

s12_program <- function(p, method = "rank") {
    if (p$structural_reason != "NONE" || !p$exact_H4_technical_eligible) 
        return(list(state = "NA_UNAVAILABLE", reason = if (p$structural_reason != "NONE") p$structural_reason else "INCOMPLETE_H4_MEASUREMENT", 
            biological_n = length(p$complete_pair_donors)))
    ids <- p$member_qualification$feature_id
    sc <- if (method == "rank") 
        s01_rank_program_score(p$score_counts, ids, p$background)
    else s01_standardized_member_mean(p$score_counts, ids)
    if (!isTRUE(sc$estimable)) 
        return(list(state = "NA_NONESTIMABLE", reason = if (is.null(sc$reason)) "INVALID_PROGRAM_SCORE" else sc$reason, 
            biological_n = length(p$complete_pair_donors), score = sc))
    v <- s12_pairs(p, sc)
    a <- s01_paired_delta_t(v$t, v$n, v$ids)
    g <- s01_paired_g_av(v$t, v$n, v$ids)
    if (is.null(a$reason)) 
        a$reason <- if (a$state == "ESTIMABLE") 
            "NONE"
        else paste(names(a$flags)[a$flags], collapse = ";")
    a$standard_error <- if (isTRUE(a$estimable)) 
        sd(v$t - v$n)/sqrt(length(v$ids))
    else NA_real_
    list(state = a$state, reason = a$reason, biological_n = length(v$ids), estimate = a, g_av = g, score = sc, 
        source = data.frame(donor = v$ids, NAT_unit = p$modeling_units[seq(1L, ncol(p$score_counts), 2L)], Tumor_unit = p$modeling_units[seq(2L, 
            ncol(p$score_counts), 2L)], NAT_score = v$n, Tumor_score = v$t, delta = v$t - v$n), method = method)
}

s12_native <- function(p, fit_function = pe_sc_members) {
    b <- p$member_qualification
    x <- fit_function(p, b$feature_id, if (p$structural_reason == "NONE") 
        NULL
    else p$structural_reason)
    r <- x$results
    s01_require(identical(r$member, b$feature_id), "S12 native feature roster")
    r$feature_id <- r$member
    r$member <- b$member
    r$source_symbol <- b$source_symbol
    r$source_row <- b$source_row
    unavailable <- !b$technical_eligible
    r$state[unavailable] <- "NA_UNAVAILABLE"
    r$reason[unavailable] <- b$reason[unavailable]
    x$results <- r
    x$design <- p$model_prepared$design
    x$coefficient <- p$model_prepared$coefficient
    x$sample_data <- p$sample_data
    x$library_sizes <- p$modeling_library_sizes
    x$background <- p$modeling_background
    x$member_bridge <- b
    x$biological_n <- length(p$complete_pair_donors)
    x
}

s12_composition <- function(c, category, seed) {
    ids <- c$complete_pair_donors
    x <- c$category_counts
    s01_assert_unique_key(x, c("donor", "tissue_origin", "category"), "S12 composition key")
    s01_require(category %in% x$category, "S12 unknown composition category")
    z <- x[x$category == category & x$donor %in% ids, , drop = FALSE]
    n <- z[match(paste(ids, "NAT"), paste(z$donor, z$tissue_origin)), , drop = FALSE]
    t <- z[match(paste(ids, "Tumor"), paste(z$donor, z$tissue_origin)), , drop = FALSE]
    s01_require(!anyNA(n$donor) && !anyNA(t$donor) && all(n$eligible & t$eligible) && all(c(n$retained_epithelial_cells, 
        t$retained_epithelial_cells) >= 100), "S12 composition population/denominator")
    guard <- pe_composition_guard(t$category_cells, t$retained_epithelial_cells, n$category_cells, n$retained_epithelial_cells)
    src <- data.frame(category = category, donor = ids, NAT_count = n$category_cells, Tumor_count = t$category_cells, 
        NAT_total = n$retained_epithelial_cells, Tumor_total = t$retained_epithelial_cells, NAT_proportion = n$category_cells/n$retained_epithelial_cells, 
        Tumor_proportion = t$category_cells/t$retained_epithelial_cells, raw_delta = guard$raw_delta, NAT_logit = NA_real_, 
        Tumor_logit = NA_real_, logit_delta = NA_real_)
    if (guard$degenerate) {
        a <- list(state = "NA_NONESTIMABLE", reason = guard$reason, effect = NA_real_, p_value = NA_real_, ci_lower = NA_real_, 
            ci_upper = NA_real_)
        boot <- list(state = "NA_UNAVAILABLE", reason = "DEGENERATE_RAW_CATEGORY_NO_BOOTSTRAP_INTERVAL", ci_lower = NA_real_, 
            ci_upper = NA_real_, valid_replicates = 0L)
    }
    else {
        src$NAT_logit <- s01_empirical_logit(src$NAT_count, src$NAT_total)
        src$Tumor_logit <- s01_empirical_logit(src$Tumor_count, src$Tumor_total)
        src$logit_delta <- src$Tumor_logit - src$NAT_logit
        a <- s01_paired_delta_t(src$Tumor_logit, src$NAT_logit, ids)
        a$reason <- if (a$state == "ESTIMABLE") 
            "NONE"
        else paste(names(a$flags)[a$flags], collapse = ";")
        a$standard_error <- if (isTRUE(a$estimable)) 
            sd(src$logit_delta)/sqrt(length(ids))
        else NA_real_
        boot <- s01_pair_bootstrap(src$Tumor_proportion, src$NAT_proportion, ids, seed)
    }
    list(state = a$state, reason = a$reason, biological_n = length(ids), estimate = a, guard = guard, source = src, 
        bootstrap = boot, seed = seed, raw_mean_delta = mean(src$raw_delta))
}

s12_join_close <- function(roster, rows) {
    s01_require(!anyDuplicated(rows$hypothesis_id) && nrow(rows) == nrow(roster) && setequal(rows$hypothesis_id, 
        roster$hypothesis_id), "S12 complete-roster join")
    rows <- rows[match(roster$hypothesis_id, rows$hypothesis_id), , drop = FALSE]
    s01_require(identical(rows$endpoint, roster$endpoint) && identical(rows$member, roster$member) && identical(rows$family, 
        roster$family) && all(rows$family_size == roster$family_size), "S12 roster semantics")
    closure <- list()
    for (f in unique(roster$family)) {
        j <- which(roster$family == f)
        s01_require(length(j) == unique(roster$family_size[j]), "S12 family shrinkage")
        v <- pe_close_family(setNames(rows$p_value[j], rows$hypothesis_id[j]), rows$state[j], rows$hypothesis_id[j])
        rows$q_value[j] <- unname(v$q)
        closure[[f]] <- data.frame(family = f, frozen_size = length(j), estimable = sum(rows$state[j] == "ESTIMABLE"), 
            unavailable = sum(rows$state[j] == "NA_UNAVAILABLE"), nonestimable = sum(rows$state[j] == "NA_NONESTIMABLE"), 
            pending = sum(rows$state[j] == "PENDING"), closed = v$closed)
    }
    list(results = rows, families = do.call(rbind, closure))
}

s12_load_bound <- function() {
    p <- context$read_rds(file.path(s12_b, "prepared_populations.rds"))
    raw <- context$read_rds(file.path(s12_b, "original_pseudobulks.rds"))
    setNames(lapply(names(p), function(e) s12_bind(p[[e]], raw[[e]])), names(p))
}

s12_program_target <- function(location = s12_p) s12_target("programs", function() {
    p <- s12_load_bound()
    a <- lapply(p, s12_program)
    alt <- s12_program(p$overall_epithelial, "standardized_member_mean")
    primary <- a$overall_epithelial
    boot <- if (!is.null(primary$source)) 
        s01_pair_bootstrap(primary$source$Tumor_score, primary$source$NAT_score, primary$source$donor, s12_seed("S12_B_overall_H4_paired_bootstrap"))
    else list(state = "NA_UNAVAILABLE", reason = "NO_VALID_OVERALL_SCORE", ci_lower = NA_real_, ci_upper = NA_real_, 
        valid_replicates = 0L)
    b <- p$overall_epithelial$member_qualification
    mat <- p$overall_epithelial$score_counts
    rank_source <- do.call(rbind, lapply(seq_len(ncol(mat)), function(i) {
        ranks <- rank(mat[, i], ties.method = "average")
        j <- match(b$feature_id, rownames(mat))
        data.frame(unit_id = colnames(mat)[i], donor = p$overall_epithelial$sample_data$donor[i], source_tissue = p$overall_epithelial$sample_data$tissue[i], 
            member = b$member, feature_id = b$feature_id, source_symbol = b$source_symbol, source_row = b$source_row, 
            raw_count = mat[j, i], average_rank = ranks[j], rank_percentile = (ranks[j] - 1)/(nrow(mat) - 1), background_size = nrow(mat))
    }))
    x <- mat[b$feature_id, , drop = FALSE]
    moments <- data.frame(member = b$member, feature_id = b$feature_id, pooled_mean = rowMeans(x), pooled_sample_SD = apply(x, 
        1, sd), formal_samples = ncol(x), input_scale = "unintegrated_raw_pseudobulk_counts_pooled_formal_units")
    list(programs = a, alternative = alt, bootstrap = boot, rank_source = rank_source, standardization = moments, 
        bootstrap_seed = s12_seed("S12_B_overall_H4_paired_bootstrap"), score_input = "unintegrated_raw_pseudobulk_counts; within-sample TMM/CPM scaling is rank-monotone; standardized diagnostic uses original approved S06 raw-count scale", 
        backgrounds = lapply(p, function(z) z$background), populations = lapply(p, function(z) list(detection_units = z$eligible_units, 
            model_units = z$modeling_units, donors = z$complete_pair_donors)))
}, location)

s12_native_target <- function(endpoint, location = s12_p) s12_target(paste0("members_", endpoint), function() {
    p <- s12_load_bound()
    x <- s12_native(p[[endpoint]])
    if (is.null(x$native) && any(x$results$state == "NA_NONESTIMABLE" & grepl("NATIVE_MODEL_OR_TEST_FAILURE", x$results$reason))) 
        stop("Unexpected native failure requires review before family closure")
    x$endpoint <- endpoint
    x$runtime <- s12_runtime()
    x
}, location)

s12_composition_target <- function(location = s12_p) s12_target("composition", function() {
    c <- context$read_rds(file.path(s12_b, "composition_structural_population.rds"))
    categories <- sub("^composition_", "", s12_roster()$endpoint[s12_roster()$member == "captured_category"])
    x <- setNames(lapply(categories, function(k) s12_composition(c, k, s12_seed(paste0("S12_B_composition_", k, 
        "_raw_proportion_bootstrap")))), categories)
    list(categories = x, donors = c$complete_pair_donors, population_source = file.path(s12_b, "composition_structural_population.rds"))
}, location)

s12_value <- function(x, k, default = NA_real_) if (is.null(x[[k]])) default else x[[k]]

s12_result_row <- function(roster, type, n, background_size = NA_integer_) {
    data.frame(hypothesis_id = roster$hypothesis_id, endpoint = roster$endpoint, member = roster$member, family = roster$family, 
        family_size = roster$family_size, source = "GSE150290", source_role = "additional_epithelial_context", 
        timing = "post_effect_exploratory", role = "post_effect_exploratory_inference", dependence = "same_source_program_member_and_composition_views_dependent", 
        biological_unit = "donor", biological_n = n, contrast = "Tumor_minus_adjacent_nontumor", type = type, population_reference = if (type == 
            "composition") 
            file.path(s12_b, "composition_structural_population.rds")
        else file.path(s12_b, paste0(roster$endpoint, "_complete_expression_pairs.tsv")), background_reference = if (type == 
            "composition") 
            "ALL_RETAINED_EPITHELIAL_CELLS"
        else file.path(s12_b, paste0(roster$endpoint, if (type == "member") 
            "_modeling_background.tsv"
        else "_detection_background.tsv")), background_size = background_size, technical_state = roster$technical_state, 
        state = "PENDING", reason = "RUNTIME_PENDING", scale = if (type == "program") 
            "ascending_rank_percentile_mean"
        else if (type == "member") 
            "native_log2_fold_change"
        else "paired_empirical_logit_difference", method = if (type == "member") 
            "edgeR_4.10.1_TMM_modern_robust_QL_donor_plus_tissue_full_background"
        else "paired_Student_t", effect = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, standard_error = NA_real_, 
        statistic = NA_real_, numerator_df = NA_real_, denominator_df = NA_real_, p_value = NA_real_, q_value = NA_real_, 
        confidence_interval_status = "NA_NONESTIMABLE", standard_error_status = "NA_NONESTIMABLE", uncertainty_reason = "RUNTIME_PENDING", 
        g_av = NA_real_, g_av_variance = NA_real_, g_av_ci_lower = NA_real_, g_av_ci_upper = NA_real_, g_av_state = "NA_NOT_APPLICABLE", 
        g_av_reason = "NA_NOT_APPLICABLE", standardized_point_magnitude_status = "NA_NOT_APPLICABLE", standardized_precision_status = "NA_NOT_APPLICABLE", 
        feature_id = NA_character_, source_symbol = NA_character_, source_row = NA_integer_, raw_paired_proportion_difference = NA_real_, 
        raw_proportion_bootstrap_lower = NA_real_, raw_proportion_bootstrap_upper = NA_real_, raw_proportion_bootstrap_state = "NA_NOT_APPLICABLE", 
        bootstrap_valid_draws = NA_integer_, bootstrap_seed = NA_integer_, both_zero_pairs = NA_integer_, degeneracy_reason = "NA_NOT_APPLICABLE", 
        reference_tissue = "Adjacent_gastritis_atrophy_metaplasia_qualified_reference_not_healthy", actual_processing_batches = "UNKNOWN", 
        participant_independence = "UNKNOWN_ACROSS_ACCESSIONS", provenance_id = file.path(s12_m, "ENTRY_RELEASE.json"), 
        stringsAsFactors = FALSE)
}

s12_apply_estimate <- function(r, a) {
    for (k in intersect(c("effect", "ci_lower", "ci_upper", "standard_error", "statistic", "p_value"), names(a))) r[[k]] <- a[[k]]
    r$state <- a$state
    r$reason <- s12_value(a, "reason", "NONE")
    r$denominator_df <- s12_value(a, "degrees_of_freedom")
    r$confidence_interval_status <- if (a$state == "ESTIMABLE") 
        "AVAILABLE"
    else a$state
    r$standard_error_status <- if (is.finite(r$standard_error)) 
        "AVAILABLE"
    else a$state
    r$uncertainty_reason <- if (a$state == "ESTIMABLE") 
        "NONE"
    else r$reason
    r
}

s12_tables_compute <- function(location = s12_p) {
    roster <- s12_roster()
    prog <- context$read_rds(file.path(location, "programs.rds"))
    comp <- context$read_rds(file.path(location, "composition.rds"))
    natives <- setNames(lapply(c("overall_epithelial", "enteroendocrine"), function(e) context$read_rds(file.path(location, 
        paste0("members_", e, ".rds")))), c("overall_epithelial", "enteroendocrine"))
    rows <- lapply(seq_len(nrow(roster)), function(i) {
        e <- roster[i, ]
        type <- if (e$member == "H_core4") 
            "program"
        else if (e$member == "captured_category") 
            "composition"
        else "member"
        if (type == "program") {
            a <- prog$programs[[e$endpoint]]
            r <- s12_result_row(e, type, a$biological_n, length(prog$backgrounds[[e$endpoint]]))
            r <- s12_apply_estimate(r, if (is.null(a$estimate)) 
                a
            else a$estimate)
            if (!is.null(a$g_av)) {
                g <- a$g_av
                r$g_av <- g$g_av
                r$g_av_variance <- g$variance
                r$g_av_ci_lower <- g$ci_lower
                r$g_av_ci_upper <- g$ci_upper
                r$g_av_state <- g$state
                r$g_av_reason <- s12_value(g, "reason", "NONE")
                r$standardized_point_magnitude_status <- if (g$state != "ESTIMABLE") 
                  "NA_NONESTIMABLE"
                else if (abs(g$g_av) >= 0.2) 
                  "AT_OR_ABOVE_REGISTERED_0.20_BOUNDARY"
                else "BELOW_REGISTERED_0.20_BOUNDARY"
                r$standardized_precision_status <- if (g$state != "ESTIMABLE") 
                  "NA_NONESTIMABLE"
                else if (g$ci_lower > -0.2 && g$ci_upper < 0.2) 
                  "CI_WITHIN_SMALL_EFFECT_REGION"
                else if (g$ci_lower <= 0 && g$ci_upper >= 0) 
                  "CI_SPANS_ZERO_AND_SUBSTANTIVE_BOUNDARY"
                else "CI_EXCLUDES_ZERO_NOT_BOUNDED_SMALL"
            }
        }
        else if (type == "member") {
            x <- natives[[e$endpoint]]
            a <- x$results[match(e$member, x$results$member), , drop = FALSE]
            r <- s12_result_row(e, type, x$biological_n, length(x$background))
            for (k in c("effect", "statistic", "numerator_df", "denominator_df", "p_value", "feature_id", "source_symbol", 
                "source_row", "state", "reason")) r[[k]] <- a[[k]]
            r$confidence_interval_status <- r$standard_error_status <- "NA_UNAVAILABLE"
            r$uncertainty_reason <- a$uncertainty_reason
        }
        else {
            category <- sub("^composition_", "", e$endpoint)
            x <- comp$categories[[category]]
            r <- s12_result_row(e, type, x$biological_n)
            r <- s12_apply_estimate(r, x$estimate)
            r$raw_paired_proportion_difference <- x$raw_mean_delta
            r$both_zero_pairs <- x$guard$both_zero_pairs
            r$degeneracy_reason <- x$guard$reason
            r$raw_proportion_bootstrap_lower <- x$bootstrap$ci_lower
            r$raw_proportion_bootstrap_upper <- x$bootstrap$ci_upper
            r$raw_proportion_bootstrap_state <- x$bootstrap$state
            r$bootstrap_valid_draws <- x$bootstrap$valid_replicates
            r$bootstrap_seed <- x$seed
        }
        r
    })
    out <- s12_join_close(roster, do.call(rbind, rows))
    s01_require(all(out$families$closed) && all(out$results$state %in% c("ESTIMABLE", "NA_UNAVAILABLE", "NA_NONESTIMABLE")), 
        "S12 unfinished family")
    diagnostics <- list()
    drow <- function(id, category, state, reason, n, effect = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_, 
        valid = NA_integer_, seed = NA_integer_, scale = "ascending_rank_percentile_mean") data.frame(analysis_id = id, 
        category = category, state = state, reason = reason, source = "GSE150290", role = "same_data_robustness", 
        timing = "post_effect_exploratory", family = "NO_FORMAL_FAMILY", q_value = NA_real_, nominal_p_unadjusted = p, 
        biological_n = n, effect = effect, ci_lower = lower, ci_upper = upper, scale = scale, valid_draws = valid, 
        seed = seed, draws = if (is.na(seed)) 
            NA_integer_
        else 10000L, quantile_type = if (is.na(seed)) 
            NA_integer_
        else 7L, minimum_valid_draws = if (is.na(seed)) 
            NA_integer_
        else 9500L)
    a <- prog$alternative
    b <- a$estimate
    diagnostics[[1]] <- drow("S12_B_overall_H4_standardized_member_mean", "standardized_member_mean", a$state, 
        a$reason, a$biological_n, s12_value(b, "effect"), s12_value(b, "ci_lower"), s12_value(b, "ci_upper"), s12_value(b, 
            "p_value"), scale = "mean_pooled_member_z_on_original_raw_count_scale")
    a <- prog$programs$overall_epithelial
    b <- prog$bootstrap
    diagnostics[[2]] <- drow("S12_B_overall_H4_paired_bootstrap", "paired_bootstrap", b$state, s12_value(b, "reason", 
        "NONE"), a$biological_n, s12_value(a$estimate, "effect"), b$ci_lower, b$ci_upper, valid = b$valid_replicates, 
        seed = prog$bootstrap_seed)
    for (k in names(comp$categories)) {
        x <- comp$categories[[k]]
        b <- x$bootstrap
        diagnostics[[length(diagnostics) + 1L]] <- drow(paste0("S12_B_composition_", k, "_raw_proportion_bootstrap"), 
            "composition_raw_proportion_bootstrap", b$state, s12_value(b, "reason", "NONE"), x$biological_n, x$raw_mean_delta, 
            b$ci_lower, b$ci_upper, valid = b$valid_replicates, seed = x$seed, scale = "raw_captured_proportion_difference")
    }
    diagnostics[[5]] <- drow("S12_B_lineage_matched_overall_H4", "conditional_lineage_matched_overall_description", 
        "NA_NOT_APPLICABLE", "L_H4_ELIGIBLE_ZERO_NO_ELIGIBLE_WITHIN_LINEAGE_H4_COMPARISON", NA_integer_)
    out$diagnostics <- do.call(rbind, diagnostics)
    out
}

s12_tables_target <- function(location = s12_p) s12_target("tables", function() s12_tables_compute(location), location)

s12_graph <- function(location = s12_p) {
    s12_program_target(location)
    s12_native_target("overall_epithelial", location)
    s12_native_target("enteroendocrine", location)
    s12_composition_target(location)
    s12_tables_target(location)
}

  s12_target <- function(name,compute,location=s12_p) {
    context$guard(); context$resource_check(8,name)
    x <- compute(); x$public_input_manifest_sha256 <- context$manifest_sha256
    s12_rds(x,file.path(location,paste0(name,".rds")))
  }
s12_exports <- function(location = s12_p, result_location = s12_r) {
    x <- context$read_rds(file.path(location, "tables.rds"))
    p <- context$read_rds(file.path(location, "programs.rds"))
    c <- context$read_rds(file.path(location, "composition.rds"))
    out <- c(s12_tsv(x$results, file.path(result_location, "all_results.tsv")), s12_tsv(x$results[x$results$type == 
        "program", ], file.path(result_location, "program_results.tsv")), s12_tsv(x$results[x$results$type == "member", 
        ], file.path(result_location, "member_results.tsv")), s12_tsv(x$results[x$results$type == "composition", 
        ], file.path(result_location, "composition_results.tsv")), s12_tsv(x$families, file.path(result_location, 
        "family_closure.tsv")), s12_tsv(x$diagnostics, file.path(result_location, "diagnostic_index.tsv")), s12_tsv(p$rank_source, 
        file.path(result_location, "overall_H4_member_rank_source.tsv")), s12_tsv(p$standardization, file.path(result_location, 
        "standardized_member_reference.tsv")))
    if (!is.null(p$programs$overall_epithelial$source)) 
        out <- c(out, s12_tsv(p$programs$overall_epithelial$source, file.path(result_location, "overall_H4_donor_scores.tsv")))
    if (!is.null(p$alternative$source)) 
        out <- c(out, s12_tsv(p$alternative$source, file.path(result_location, "standardized_member_donor_scores.tsv")))
    out <- c(out, s12_tsv(do.call(rbind, lapply(c$categories, function(z) z$source)), file.path(result_location, 
        "captured_composition_donor_source.tsv")))
    for (e in c("overall_epithelial", "enteroendocrine")) {
        z <- context$read_rds(file.path(location, paste0("members_", e, ".rds")))
        lib <- data.frame(unit_id = z$sample_data$unit_id, donor = z$sample_data$donor, tissue = z$sample_data$tissue, 
            full_background_library_size = unname(z$library_sizes), TMM_factor = if (is.null(z$native)) 
                NA_real_
            else z$native$dge$samples$norm.factors)
        design <- data.frame(unit_id = z$sample_data$unit_id, z$design, check.names = FALSE)
        out <- c(out, s12_tsv(lib, file.path(result_location, paste0(e, "_native_library_context.tsv"))), s12_tsv(design, 
            file.path(result_location, paste0(e, "_native_design.tsv"))), s12_tsv(z$results, file.path(result_location, 
            paste0(e, "_native_extraction.tsv"))))
    }
    out
}

s12_figure_sources <- function(result_location = s12_r) {
    r <- s12_read(file.path(result_location, "all_results.tsv"))
    old <- s12_read("results/s09_reproducibility/s06_single_cell_effects/all_program_estimates.tsv")
    o <- old[old$endpoint_id == "GSE206785_epithelial_H", ]
    s01_require(nrow(o) == 1L, "Original epithelial context row")
    oldcontext <- data.frame(source = o$source, endpoint = "epithelial", biological_n = o$biological_n, effect = o$effect, 
        ci_lower = o$ci_lower, ci_upper = o$ci_upper, p_value = o$p_value, q_value = o$q_value, state = o$estimability, 
        role = o$role, timing = o$timing, background_reference = "data/processed/s03_single_cell_ingress/preparation_objects.rds:epithelial", 
        capture_protocol = "EPCAM_depleted", source_result = "results/s09_reproducibility/s06_single_cell_effects/all_program_estimates.tsv", 
        hypothesis_id = o$hypothesis_id)
    z <- r[r$type == "program", ]
    newcontext <- data.frame(source = z$source, endpoint = z$endpoint, biological_n = z$biological_n, effect = z$effect, 
        ci_lower = z$ci_lower, ci_upper = z$ci_upper, p_value = z$p_value, q_value = z$q_value, state = z$state, 
        role = z$role, timing = z$timing, background_reference = z$background_reference, capture_protocol = "Reported_10x_3prime_v2_no_reported_EPCAM_depletion_capture_unbiasedness_not_established", 
        source_result = paste0(result_location, "/program_results.tsv"), hypothesis_id = z$hypothesis_id)
    s12_tsv(rbind(newcontext, oldcontext), file.path(result_location, "program_figure_source.tsv"))
    s12_tsv(r[r$type == "member", ], file.path(result_location, "member_figure_source.tsv"))
    s12_tsv(r[r$type == "composition", ], file.path(result_location, "composition_figure_source.tsv"))
}

  run <- function() {
    context$guard(); context$resource_check(); s12_runtime()
    s12_graph(); s12_exports(); s12_figure_sources()
    list(results=list.files(s12_r,full.names=TRUE),processed=list.files(s12_p,full.names=TRUE))
  }
  list(run=run,bind=s12_bind,program=s12_program,native=s12_native,composition=s12_composition)
}
