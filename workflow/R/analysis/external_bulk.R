# Extracted maintained implementation; stable machine names and seeds retained.
external_bulk_module <- function(context) {
  s11_root <- context$root
  s11_sha <- function(p) digest::digest(file=context$input(p),algo="sha256")
  s11_now <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")
  s11_tsv <- function(x,p) context$write_tsv(x,p)
  s11_rds <- function(x,p) context$write_rds(x,p)
  s11_json <- function(x,p) context$write_json(x,p)
  s11_text <- function(x,p) context$write_text(x,p)
  s11_guard <- function(...) context$guard()
  s11_effect_guard <- function(...) context$guard()
  s11_memory <- function(stage,minimum=8,prestart=FALSE) context$resource_check(minimum,stage)
  s11_hashes <- function(paths) data.frame(path=paths,bytes=file.info(paths)$size,sha256=vapply(paths,s11_sha,""))
  s11_m <- context$output("metadata"); s11_p <- context$output("processed")
  s11_r <- context$output("results"); s11_q <- context$output("qc")
  s11_old <- "data/manifests/post_effect_supplement/s10"
s11_read <- function(p) read.delim(context$input(p), stringsAsFactors = FALSE, check.names = FALSE, quote = "", 
    comment.char = "")

s11_runtime <- function(...) context$runtime("supplementary")

s11_roster <- function() {
    x <- s11_read("config/post_effect_hypothesis_estimand_crosswalk.tsv")
    x <- x[x$module == "A", , drop = FALSE]
    s01_require(nrow(x) == 17 && !anyDuplicated(x$hypothesis_id) && all(x$source == "GSE122401"), "A hypothesis roster")
    x
}

s11_diagnostics <- function() s11_read(file.path(s11_old, "A_diagnostic_registry.tsv"))

s11_order_pairs <- function(md, samples) {
    s01_validate_identity_vector(samples, "S11 matrix sample")
    s01_require(identical(md$sample_id, samples), "S11 sample order mismatch")
    s01_require(all(md$tissue %in% c("NAT", "Tumor")) && !anyDuplicated(paste(md$patient_id, md$tissue)) && all(table(md$patient_id) == 
        2), "S11 complete pair identity")
    ids <- sort(unique(md$patient_id), method = "radix")
    n <- match(paste(ids, "NAT"), paste(md$patient_id, md$tissue))
    t <- match(paste(ids, "Tumor"), paste(md$patient_id, md$tissue))
    s01_require(!anyNA(c(n, t)), "S11 unresolved pair")
    data.frame(patient_id = ids, adjacent_sample = md$sample_id[n], tumor_sample = md$sample_id[t])
}

s11_map_results <- function(roster, rows) {
    s01_require(!anyDuplicated(rows$endpoint_id) && setequal(roster$endpoint_id, rows$endpoint_id) && nrow(rows) == 
        nrow(roster), "S11 complete-roster result join")
    rows <- rows[match(roster$endpoint_id, rows$endpoint_id), , drop = FALSE]
    cbind(roster[, c("endpoint_id", "hypothesis_id", "construct", "family_hypothesis", "multiplicity_family", "frozen_family_size")], 
        rows[, setdiff(names(rows), "endpoint_id"), drop = FALSE])
}

s11_scores <- function() {
    s11_effect_guard()
    s11_memory("score_start")
    a <- context$read_rds("data/processed/post_effect_supplement/s10/GSE122401_adapter.rds")
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    ranks <- apply(a$abundance, 2, rank, ties.method = "average")
    members <- unlist(cfg$programs, use.names = FALSE)
    bridge <- a$member_mapping
    u <- (ranks[match(bridge$ensembl_id[match(members, bridge$gene)], rownames(a$abundance)), , drop = FALSE] - 
        1)/(length(a$background) - 1)
    rownames(u) <- members
    colnames(u) <- a$metadata$sample_id
    scores <- vapply(cfg$programs, function(g) colMeans(u[g, , drop = FALSE]), numeric(ncol(u)))
    scores <- t(scores)
    dimnames(scores) <- list(names(cfg$programs), colnames(u))
    z <- list(scores = scores, member_percentiles = u, metadata = a$metadata, pairs = s11_order_pairs(a$metadata, 
        colnames(scores)), background_sha256 = s11_sha(file.path(s11_old, "GSE122401_background.tsv")), member_mapping = bridge)
    path <- s11_rds(z, file.path(s11_p, "rank_scores.rds"))
    long <- do.call(rbind, lapply(rownames(scores), function(p) data.frame(program = p, sample_id = colnames(scores), 
        score = unname(scores[p, ]))))
    s11_memory("score_complete", 6)
    c(path, s11_tsv(long, file.path(s11_p, "rank_scores.tsv")))
}

s11_native <- function() {
    s11_effect_guard()
    s11_memory("native_start")
    a <- context$read_rds("data/processed/post_effect_supplement/s10/GSE122401_adapter.rds")
    z <- pe_fit_abundance(a$abundance, a$metadata)
    path <- s11_rds(z, file.path(s11_p, "native_limma_fit.rds"))
    s11_memory("native_complete", 6)
    path
}

s11_pair_values <- function(scores, program, pairs) list(tumor = unname(scores[program, pairs$tumor_sample]), adjacent = unname(scores[program, 
    pairs$adjacent_sample]), ids = pairs$patient_id)

s11_endpoint <- function(id, t, n, method, scale, g = NULL) {
    value <- function(x, k) if (is.null(x[[k]])) 
        NA_real_
    else x[[k]]
    data.frame(endpoint_id = id, source = "GSE122401", source_role = "external_cohort_post_effect", timing = "post_effect_exploratory", 
        population = "S11_frozen_80_complete_pairs", population_reference = paste0(s11_m, "/pair_flow.tsv"), biological_unit = "patient", 
        biological_n = n, contrast = "Tumor_minus_adjacent_nontumor", scale = scale, method = method, effect = t$effect, 
        direction = s01_direction_label(t$effect), standard_error = value(t, "standard_error"), statistic = value(t, 
            "statistic"), degrees_of_freedom = value(t, "degrees_of_freedom"), ci_lower = t$ci_lower, ci_upper = t$ci_upper, 
        p_value = t$p_value, q_value = NA_real_, state = t$state, reason = if (is.null(t$reason)) 
            "NONE"
        else t$reason, uncertainty_status = if (t$state == "ESTIMABLE") 
            "AVAILABLE"
        else "NA_NONESTIMABLE", g_av = value(g, "g_av"), g_av_variance = value(g, "variance"), g_av_ci_lower = value(g, 
            "ci_lower"), g_av_ci_upper = value(g, "ci_upper"), g_av_state = if (is.null(g)) 
            "NOT_APPLICABLE_NATIVE_MEMBER"
        else g$state, g_av_reason = if (is.null(g)) 
            "NOT_APPLICABLE"
        else if (is.null(g$reason)) 
            "NONE"
        else g$reason, unmoderated_paired_mean_difference = t$effect, members = "", member_count = NA_integer_, 
        background_genes = 23818L, background_reference = paste0(s11_old, "/GSE122401_background.tsv"), dependence = "A_and_PDC000214_share_all_80_patients;cross_cohort_overlap_not_proven_absent", 
        family_closed = FALSE, native_df_prior = NA_real_, native_df_residual = NA_real_, stringsAsFactors = FALSE)
}

s11_inference <- function() {
    s11_effect_guard()
    z <- context$read_rds(file.path(s11_p, "rank_scores.rds"))
    native <- context$read_rds(file.path(s11_p, "native_limma_fit.rds"))
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    roster <- s11_roster()
    rows <- list()
    for (p in names(cfg$programs)) {
        v <- s11_pair_values(z$scores, p, z$pairs)
        t <- s01_paired_delta_t(v$tumor, v$adjacent, v$ids)
        g <- s01_paired_g_av(v$tumor, v$adjacent, v$ids)
        t$standard_error <- if (t$estimable) 
            sd(v$tumor - v$adjacent)/sqrt(length(v$ids))
        else NA_real_
        r <- s11_endpoint(paste0("PE_A_GSE122401_", p), t, length(v$ids), "paired_delta_Student_t;S01_g_av_SMCRP_LS", 
            "rank_percentile_equal_weight_mean", g)
        r$members <- paste(cfg$programs[[p]], collapse = ";")
        r$member_count <- length(cfg$programs[[p]])
        rows[[length(rows) + 1L]] <- r
    }
    ex <- native$extracted
    planned <- unlist(cfg$programs[1:2], use.names = FALSE)
    for (m in planned) {
        ensg <- z$member_mapping$ensembl_id[match(m, z$member_mapping$gene)]
        i <- match(ensg, ex$feature_id)
        s01_require(!is.na(i), "S11 planned native member absent")
        t <- lapply(ex[c("state", "effect", "standard_error", "statistic", "degrees_of_freedom", "p_value", "ci_lower", 
            "ci_upper")], `[`, i)
        t$reason <- if (t$state == "ESTIMABLE") 
            "NONE"
        else "INVALID_NATIVE_LIMMA_OUTPUT"
        r <- s11_endpoint(paste0("PE_A_GSE122401_", m), t, ex$biological_n, "limma_patient_tissue_eBayes_trend_TRUE_robust_TRUE", 
            "log2(abundance+1)")
        v <- native$prepared$values[ensg, ]
        names(v) <- colnames(native$prepared$values)
        r$unmoderated_paired_mean_difference <- mean(v[z$pairs$tumor_sample] - v[z$pairs$adjacent_sample])
        r$members <- m
        r$member_count <- 1L
        r$native_df_prior <- rep_len(native$fit$df.prior, nrow(native$prepared$values))[i]
        r$native_df_residual <- rep_len(native$fit$df.residual, nrow(native$prepared$values))[i]
        rows[[length(rows) + 1L]] <- r
    }
    all <- s11_map_results(roster, do.call(rbind, rows))
    closure <- list()
    for (f in names(cfg$A$families)) {
        j <- which(all$multiplicity_family == f)
        ids <- all$hypothesis_id[j]
        s01_require(length(j) == cfg$A$families[[f]] && all(all$frozen_family_size[j] == length(j)), "S11 family size")
        closed <- pe_close_family(setNames(all$p_value[j], ids), all$state[j], ids)
        s01_require(closed$closed, "S11 unresolved family slots")
        all$q_value[j] <- unname(closed$q)
        all$family_closed[j] <- TRUE
        closure[[length(closure) + 1L]] <- data.frame(family = f, frozen_size = length(j), terminal_slots = length(j), 
            estimable = sum(all$state[j] == "ESTIMABLE"), nonestimable = sum(all$state[j] == "NA_NONESTIMABLE"), 
            unavailable = sum(all$state[j] == "NA_UNAVAILABLE"), pending = sum(all$state[j] == "PENDING"), closed = TRUE, 
            method = "BH_on_raw_P_with_full_frozen_size")
    }
    prog <- all$family_hypothesis %in% names(cfg$programs)
    c(s11_tsv(all, file.path(s11_r, "all_inferential_results.tsv")), s11_tsv(all[prog, ], file.path(s11_r, "all_program_estimates.tsv")), 
        s11_tsv(all[!prog, ], file.path(s11_r, "all_member_estimates.tsv")), s11_tsv(do.call(rbind, closure), file.path(s11_r, 
            "family_closure.tsv")))
}

s11_diagnostic_run <- function() {
    s11_effect_guard()
    z <- context$read_rds(file.path(s11_p, "rank_scores.rds"))
    a <- context$read_rds("data/processed/post_effect_supplement/s10/GSE122401_adapter.rds")
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    reg <- s11_diagnostics()
    seeds <- s11_read("config/implementation/post_effect_seed_map.tsv")
    members <- unlist(cfg$programs[1:2], use.names = FALSE)
    logx <- log2(a$abundance[match(a$member_mapping$ensembl_id[match(members, a$member_mapping$gene)], rownames(a$abundance)), 
        , drop = FALSE] + 1)
    rownames(logx) <- members
    rows <- list()
    scores <- list()
    draws <- list()
    constants <- data.frame(member = members, pooled_mean = rowMeans(logx), pooled_sample_sd = apply(logx, 1, sd), 
        N_samples = ncol(logx), sd_denominator = ncol(logx) - 1)
    for (id in reg$analysis_id) {
        omit <- if (startsWith(id, "PEA.LOO.")) 
            sub("PEA.LOO.", "", id, fixed = TRUE)
        else "NONE"
        p <- if (omit != "NONE") 
            if (omit %in% cfg$programs$G_core7) 
                "G_core7"
            else "H_core4"
        else sub("^.*[.]", "", id)
        g <- cfg$programs[[p]]
        v <- s11_pair_values(z$scores, p, z$pairs)
        r <- data.frame(analysis_id = id, source = "GSE122401", role = "same_data_robustness", timing = "post_effect_exploratory", 
            program = p, omitted_member = omit, members = paste(g, collapse = ";"), member_count = length(g), population = "S11_frozen_80_complete_pairs", 
            biological_unit = "patient", biological_n = 80L, population_reference = paste0(s11_m, "/pair_flow.tsv"), 
            background_reference = paste0(s11_old, "/GSE122401_background.tsv"), background_genes = 23818L, contrast = "Tumor_minus_adjacent_nontumor", 
            scale = "rank_percentile_equal_weight_mean", method = "", state = "PENDING", reason = "NONE", primary_raw_difference = mean(v$tumor - 
                v$adjacent), effect = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, statistic = NA_real_, 
            degrees_of_freedom = NA_real_, nominal_p_value = NA_real_, p_status = "NOT_APPLICABLE", q_value = NA_real_, 
            multiplicity = "NO_DISCOVERY_FAMILY", uncertainty_status = "NOT_APPLICABLE", direction = "NA_NONESTIMABLE", 
            loo_direction_stable = "NOT_APPLICABLE", score_method_sensitive = "NOT_APPLICABLE", score_method_nonestimable = FALSE, 
            exact_zero_direction_flag = FALSE, g_av = NA_real_, g_av_ci_lower = NA_real_, g_av_ci_upper = NA_real_, 
            g_av_state = "NOT_APPLICABLE", g_av_reason = "NOT_APPLICABLE", concordant_members = NA_integer_, fixed_denominator = NA_integer_, 
            concordance = NA_real_, operation_name = "NOT_APPLICABLE", seed = NA_integer_, draws = NA_integer_, 
            valid_draws = NA_integer_, minimum_valid = NA_integer_, quantile_type = NA_integer_, rng = "NOT_APPLICABLE", 
            dependence = "same_data_robustness_not_another_cohort", stringsAsFactors = FALSE)
        if (startsWith(id, "PEA.LOO.") || startsWith(id, "PEA.standardized_member_mean.")) {
            if (omit != "NONE") {
                g <- setdiff(g, omit)
                s <- colMeans(z$member_percentiles[g, , drop = FALSE])
                sc <- list(estimable = TRUE, score = s)
                r$method <- "leave_one_member_out_same_frozen_rank_background"
            }
            else {
                sc <- s01_standardized_member_mean(logx, g)
                s <- sc$score
                r$method <- "pooled_member_z_equal_weight_mean_then_paired_Student_t"
                r$scale <- "mean_member_z_log2(abundance+1)"
            }
            if (sc$estimable) {
                names(s) <- colnames(logx)
                t <- s01_paired_delta_t(s[z$pairs$tumor_sample], s[z$pairs$adjacent_sample], z$pairs$patient_id)
                gav <- s01_paired_g_av(s[z$pairs$tumor_sample], s[z$pairs$adjacent_sample], z$pairs$patient_id)
                r$g_av <- gav$g_av
                r$g_av_ci_lower <- gav$ci_lower
                r$g_av_ci_upper <- gav$ci_upper
                r$g_av_state <- gav$state
                r$g_av_reason <- if (is.null(gav$reason)) 
                  "NONE"
                else gav$reason
                scores[[length(scores) + 1L]] <- data.frame(analysis_id = id, sample_id = names(s), score = unname(s))
            }
            else t <- list(state = "NA_NONESTIMABLE", estimable = FALSE, reason = sc$reason, effect = NA_real_, 
                ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_, statistic = NA_real_)
            r$members <- paste(g, collapse = ";")
            r$member_count <- length(g)
            r$state <- t$state
            r$reason <- if (is.null(t$reason)) 
                "NONE"
            else t$reason
            r$effect <- t$effect
            r$ci_lower <- t$ci_lower
            r$ci_upper <- t$ci_upper
            r$nominal_p_value <- t$p_value
            r$p_status <- "UNADJUSTED_DIAGNOSTIC"
            r$statistic <- t$statistic
            r$degrees_of_freedom <- if (t$estimable) 
                79
            else NA_real_
            r$uncertainty_status <- if (t$estimable) 
                "AVAILABLE"
            else "NA_NONESTIMABLE"
            if (omit != "NONE") 
                r$loo_direction_stable <- s01_loo_direction_status(r$primary_raw_difference, r$effect, t$estimable)
            else {
                sens <- s01_score_method_sensitivity(r$primary_raw_difference, r$effect, t$estimable)
                r$score_method_sensitive <- sens$score_method_sensitive
                r$score_method_nonestimable <- sens$score_method_nonestimable
                r$exact_zero_direction_flag <- sens$exact_zero_direction_flag
            }
        }
        else if (startsWith(id, "PEA.Wilcoxon.")) {
            if (all(v$tumor - v$adjacent == 0)) 
                t <- list(state = "NA_NONESTIMABLE", reason = "ALL_ZERO_PAIRED_DIFFERENCES", statistic = NA_real_, 
                  p.value = NA_real_)
            else {
                t <- wilcox.test(v$tumor, v$adjacent, paired = TRUE, alternative = "two.sided", exact = FALSE, 
                  correct = FALSE)
                t$state <- if (all(is.finite(c(t$statistic, t$p.value)))) 
                  "ESTIMABLE"
                else "NA_NONESTIMABLE"
            }
            r$method <- "paired_Wilcoxon_two_sided_exact_FALSE_correct_FALSE"
            r$state <- t$state
            r$reason <- if (is.null(t$reason)) 
                "NONE"
            else t$reason
            r$effect <- r$primary_raw_difference
            r$statistic <- unname(t$statistic)
            r$nominal_p_value <- t$p.value
            r$p_status <- "UNADJUSTED_DIAGNOSTIC"
            r$uncertainty_status <- "NOT_PROVIDED_BY_PLANNED_WILCOXON_DIAGNOSTIC"
        }
        else if (startsWith(id, "PEA.bootstrap.")) {
            op <- paste0("post_effect_A.GSE122401.bootstrap.", p)
            seed <- seeds$seed[match(op, seeds$operation_name)]
            boot <- s01_with_production_rng(seed, vapply(seq_len(cfg$common$bootstrap_draws), function(i) mean((v$tumor - 
                v$adjacent)[sample.int(length(v$ids), length(v$ids), replace = TRUE)]), numeric(1)))
            valid <- sum(is.finite(boot))
            ok <- valid >= cfg$common$bootstrap_minimum_valid
            r$method <- "paired_patient_bootstrap_mean_percentile_type7"
            r$state <- if (ok) 
                "ESTIMABLE"
            else "NA_NONESTIMABLE"
            r$reason <- if (ok) 
                "NONE"
            else "FEWER_THAN_9500_VALID_DRAWS"
            r$effect <- r$primary_raw_difference
            if (ok) {
                ci <- quantile(boot[is.finite(boot)], c(0.025, 0.975), type = cfg$common$bootstrap_quantile_type, 
                  names = FALSE)
                r$ci_lower <- ci[1]
                r$ci_upper <- ci[2]
            }
            r$uncertainty_status <- if (ok) 
                "AVAILABLE"
            else "NA_NONESTIMABLE"
            r$operation_name <- op
            r$seed <- seed
            r$draws <- length(boot)
            r$valid_draws <- valid
            r$minimum_valid <- cfg$common$bootstrap_minimum_valid
            r$quantile_type <- cfg$common$bootstrap_quantile_type
            r$rng <- cfg$rng
            draws[[length(draws) + 1L]] <- data.frame(analysis_id = id, draw = seq_along(boot), mean_difference = boot, 
                valid = is.finite(boot))
        }
        else {
            native_delta <- rowMeans(logx[g, z$pairs$tumor_sample, drop = FALSE] - logx[g, z$pairs$adjacent_sample, 
                drop = FALSE])
            c <- s01_member_direction_concordance(r$primary_raw_difference, native_delta, length(g))
            r$method <- "unrounded_native_member_directions_fixed_denominator"
            r$state <- c$state
            r$reason <- if (c$state == "ESTIMABLE") 
                "NONE"
            else c$state
            r$effect <- r$primary_raw_difference
            r$concordant_members <- if (is.null(c$concordant_members)) 
                NA_integer_
            else c$concordant_members
            r$fixed_denominator <- length(g)
            r$concordance <- c$proportion
        }
        r$direction <- s01_direction_label(r$effect)
        rows[[length(rows) + 1L]] <- r
    }
    d <- do.call(rbind, rows)
    s01_require(identical(d$analysis_id, reg$analysis_id), "S11 diagnostic identity order")
    c(s11_tsv(d, file.path(s11_r, "robustness_diagnostics.tsv")), s11_tsv(constants, file.path(s11_p, "standardized_member_constants.tsv")), 
        s11_tsv(do.call(rbind, scores), file.path(s11_p, "diagnostic_scores.tsv")), s11_tsv(do.call(rbind, draws), 
            file.path(s11_p, "bootstrap_draws.tsv")))
}

  run <- function() {
    context$guard(); context$resource_check(); s11_runtime()
    a <- context$read_rds("data/processed/post_effect_supplement/s10/GSE122401_adapter.rds")
    s01_require(identical(dim(a$abundance),c(23818L,160L)),"Frozen A dimensions")
    md <- s11_read(file.path(s11_old,"GSE122401_sample_pair_ledger.tsv"))
    bg <- s11_read(file.path(s11_old,"GSE122401_background.tsv"))$ensembl_gene_id
    bridge <- s11_read(file.path(s11_old,"GSE122401_feature_bridge.tsv"))
    s01_require(identical(a$metadata,md)&&identical(a$background,bg)&&identical(a$member_mapping,bridge),"Frozen A identity")
    s01_require(identical(colnames(a$abundance),md$sample_id)&&identical(rownames(a$abundance),bg),"Frozen A order")
    s11_tsv(s11_order_pairs(md,colnames(a$abundance)),file.path(s11_m,"pair_flow.tsv"))
    s11_scores(); s11_native(); s11_inference(); s11_diagnostic_run()
    list(results=list.files(s11_r,full.names=TRUE),processed=list.files(s11_p,full.names=TRUE))
  }
  list(run=run,roster=s11_roster,order_pairs=s11_order_pairs)
}
