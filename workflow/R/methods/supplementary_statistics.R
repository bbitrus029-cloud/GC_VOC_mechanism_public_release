# Extracted maintained implementation; stable machine names and seeds retained.
pe_extract_limma <- function(fit, prepared) {
    n <- nrow(prepared$values)
    fields <- c("coefficients", "stdev.unscaled", "s2.post", "t", "p.value", "df.total", "df.prior")
    s01_require(all(fields %in% names(fit)), "missing native limma output")
    s01_require(identical(rownames(fit$coefficients), rownames(prepared$values)), "limma gene ordering mismatch")
    for (v in c("s2.post", "df.total", "df.prior")) s01_require(length(fit[[v]]) %in% c(1L, n), paste("invalid limma vector", 
        v))
    for (v in c("coefficients", "stdev.unscaled", "t", "p.value")) s01_require(is.matrix(fit[[v]]) && identical(dim(fit[[v]]), 
        dim(fit$coefficients)) && identical(rownames(fit[[v]]), rownames(prepared$values)), paste("invalid limma matrix", 
        v))
    s01_require(prepared$coefficient <= ncol(fit$coefficients), "missing tissue coefficient")
    s01_extract_limma_contract(fit, prepared, "log2(abundance+1)~tissue+patient", "upstream_UQ_no_renormalization_trend_robust")
}

pe_fit_abundance <- function(abundance, metadata) {
    s01_validate_numeric_matrix(abundance, "continuous abundance")
    s01_require(all(abundance >= 0), "negative abundance")
    p <- s01_prepare_paired_model(log2(abundance + 1), metadata, "patient_id", "patient")
    f <- limma::eBayes(limma::lmFit(p$values, p$design), trend = TRUE, robust = TRUE)
    list(extracted = pe_extract_limma(f, p), fit = f, prepared = p)
}

pe_close_family <- function(p, state, ids) {
    s01_require(length(p) == length(ids) && length(state) == length(ids) && identical(names(p), ids) && !anyDuplicated(ids), 
        "family roster mismatch")
    s01_require(all(state %in% c("PENDING", "ESTIMABLE", "NA_NONESTIMABLE", "NA_UNAVAILABLE")), "unknown family state")
    s01_require(all(is.na(p[state != "ESTIMABLE"])) && all(is.finite(p[state == "ESTIMABLE"])), "family P/state mismatch")
    if (any(state == "PENDING")) 
        return(list(closed = FALSE, q = setNames(rep(NA_real_, length(ids)), ids)))
    list(closed = TRUE, q = s01_frozen_bh_values(p, length(ids), ids))
}

pe_d7_reference <- function(expression, metadata, members) {
    s01_require(is.matrix(expression) && is.numeric(expression) && !is.null(dimnames(expression)), "invalid D7 matrix")
    s01_validate_identity_vector(rownames(expression), "D7 feature")
    s01_validate_identity_vector(colnames(expression), "D7 sample")
    s01_require(identical(metadata$sample_id, colnames(expression)), "D7 metadata ordering mismatch")
    s01_require(!anyDuplicated(paste(metadata$patient_id, metadata$tissue)), "duplicate D7 pair identity")
    s01_require(all(metadata$tissue %in% c("NAT", "Tumor")) && all(table(metadata$patient_id) == 2L), "incomplete D7 identity pairs")
    missing <- setdiff(members, rownames(expression))
    available <- intersect(members, rownames(expression))
    ids <- sort(unique(metadata$patient_id), method = "radix")
    good <- vapply(ids, function(id) all(is.finite(expression[available, metadata$patient_id == id, drop = FALSE])), 
        TRUE)
    md <- metadata[metadata$patient_id %in% ids[good], , drop = FALSE]
    md <- md[order(md$patient_id, match(md$tissue, c("NAT", "Tumor")), method = "radix"), , drop = FALSE]
    x <- matrix(NA_real_, length(members), nrow(md), dimnames = list(members, md$sample_id))
    x[available, ] <- expression[available, md$sample_id, drop = FALSE]
    if (ncol(x) < 2L) 
        return(list(state = "NA_NONESTIMABLE", reason = "NO_COMMON_FINITE_REFERENCE", members = members, excluded = ids[!good]))
    mu <- rowMeans(x)
    sd <- apply(x, 1L, stats::sd)
    constant <- vapply(seq_len(nrow(x)), function(i) if (members[i] %in% missing) 
        FALSE
    else s01_effective_zero(x[i, ]), TRUE)
    list(state = if (length(missing) || any(constant)) "NA_NONESTIMABLE" else "READY", reason = if (length(missing)) "MISSING_OR_AMBIGUOUS_REQUIRED_MEMBER" else if (any(constant)) "CONSTANT_REQUIRED_MEMBER" else "NONE", 
        members = members, missing = missing, mean = mu, sd = sd, constant = constant, metadata = md, patient_ids = ids[good], 
        excluded = ids[!good], reference_sha256 = digest::digest(list(x = x, metadata = md), algo = "sha256"), 
        N_samples = ncol(x), sd_denominator = ncol(x) - 1L)
}

pe_d7_score <- function(expression, reference, samples = reference$metadata$sample_id) {
    s01_require(reference$state == "READY", "D7 reference unavailable")
    s01_require(all(samples %in% reference$metadata$sample_id), "D7 scores must select the frozen reference population")
    x <- expression[reference$members, samples, drop = FALSE]
    s01_require(all(is.finite(x)), "nonfinite D7 measurement")
    colMeans(sweep(sweep(x, 1L, reference$mean, "-"), 1L, reference$sd, "/"))
}

pe_hc3 <- function(fit) {
    X <- model.matrix(fit)
    n <- nrow(X)
    p <- ncol(X)
    fail <- function(reason) list(state = "NA_UNAVAILABLE", reason = reason, covariance = NULL, se = rep(NA_real_, 
        p), lower = rep(NA_real_, p), upper = rep(NA_real_, p))
    if (fit$rank != p || n <= p || any(!is.finite(coef(fit)))) 
        return(fail("SINGULAR_DESIGN_OR_NO_RESIDUAL_DF"))
    h <- hatvalues(fit)
    if (any(!is.finite(h)) || any(1 - h <= sqrt(.Machine$double.eps))) 
        return(fail("NUMERICALLY_ONE_LEVERAGE"))
    a <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(a)) 
        return(fail("HC3_INVERSION_FAILED"))
    v <- a %*% crossprod(X, X * (residuals(fit)/(1 - h))^2) %*% a
    if (any(!is.finite(v)) || any(diag(v) < 0) || !isTRUE(all.equal(v, t(v), tolerance = 1e-10))) 
        return(fail("INVALID_HC3_COVARIANCE"))
    se <- sqrt(diag(v))
    critical <- qt(0.975, df = n - p)
    list(state = "AVAILABLE", reason = "NONE", covariance = v, se = se, lower = coef(fit) - critical * se, upper = coef(fit) + 
        critical * se, df = n - p)
}

pe_zero_support <- function(d, covariates) {
    x <- as.matrix(d[, covariates, drop = FALSE])
    if (ncol(x) == 1L) 
        return(min(x) <= 0 && max(x) >= 0)
    s01_require(ncol(x) == 2L && nrow(x) > 0L && is.numeric(x) && all(is.finite(x)), "origin support requires finite two-covariate points")
    axis_scale <- apply(abs(x), 2L, max)
    axis_scale[axis_scale == 0] <- 1
    x <- unique(sweep(x, 2L, axis_scale, "/"))
    hull <- chull(x)
    if (length(hull) == 1L) 
        return(all(x[hull, ] == 0))
    a <- x[c(hull, hull[1]), , drop = FALSE]
    left <- a[-nrow(a), 1] * a[-1, 2]
    right <- a[-nrow(a), 2] * a[-1, 1]
    cross <- left - right
    tol <- 8 * .Machine$double.eps * (abs(left) + abs(right))
    if (length(hull) == 2L) 
        return(all(abs(cross) <= tol) && all(apply(x[hull, , drop = FALSE], 2L, min) <= 0) && all(apply(x[hull, 
            , drop = FALSE], 2L, max) >= 0))
    all(cross >= -tol) || all(cross <= tol)
}

pe_regression <- function(data, model = c("M1", "M2", "M0"), minimum_n = 20L, influence = TRUE) {
    model <- match.arg(model)
    cov <- switch(model, M0 = character(), M1 = "delta_D", M2 = c("delta_D", "delta_E"))
    required <- c("patient_id", "delta_H", cov)
    s01_require(all(required %in% names(data)), "regression input schema")
    s01_validate_identity_vector(data$patient_id, "regression patient")
    good <- apply(data[, c("delta_H", cov), drop = FALSE], 1L, function(x) all(is.finite(x)))
    d <- data[good, required, drop = FALSE]
    d <- d[order(d$patient_id, method = "radix"), , drop = FALSE]
    fail <- function(reason) list(state = "NA_NONESTIMABLE", reason = reason, biological_n = nrow(d), excluded = data$patient_id[!good], 
        model = model)
    if (nrow(d) < minimum_n) 
        return(fail("INSUFFICIENT_COMPLETE_PAIRS"))
    if (length(cov) && any(vapply(d[cov], s01_effective_zero, TRUE))) 
        return(fail("CONSTANT_COVARIATE"))
    formula <- reformulate(cov, response = "delta_H")
    X <- model.matrix(formula, d)
    if (qr(X)$rank != ncol(X) || nrow(X) <= ncol(X)) 
        return(fail("SINGULAR_DESIGN_OR_NO_RESIDUAL_DF"))
    f <- lm(formula, d)
    tab <- coef(summary(f))
    if (any(!is.finite(tab))) 
        return(fail("INVALID_PRIMARY_OLS_OUTPUT"))
    ci <- confint(f, level = 0.95)
    coefficients <- data.frame(coefficient = rownames(tab), effect = tab[, 1], standard_error = tab[, 2], statistic = tab[, 
        3], residual_df = df.residual(f), p_value = tab[, 4], ci_lower = ci[, 1], ci_upper = ci[, 2], row.names = NULL)
    loo <- NULL
    if (influence && length(cov)) 
        loo <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
            z <- pe_regression(d[-i, , drop = FALSE], model, minimum_n = 1L, influence = FALSE)
            do.call(rbind, lapply(cov, function(v) data.frame(omitted_patient = d$patient_id[i], coefficient = v, 
                state = z$state, reason = z$reason, effect = if (z$state == "ESTIMABLE") 
                  z$coefficients$effect[match(v, z$coefficients$coefficient)]
                else NA_real_)))
        }))
    diagnostics <- list(zero_covariate_supported = if (length(cov)) pe_zero_support(d, cov) else TRUE, ranges = if (length(cov)) lapply(d[cov], 
        range) else list(), intercept_label = if (length(cov) && !pe_zero_support(d, cov)) "EXTRAPOLATIVE" else "WITHIN_OBSERVED_SUPPORT")
    if (length(cov) == 2L) {
        z <- scale(as.matrix(d[cov]))
        r <- cor(z)[1, 2]
        diagnostics$condition_number <- kappa(z, exact = TRUE)
        diagnostics$VIF <- setNames(rep(1/(1 - r^2), 2), cov)
        diagnostics$unstable_conditional_attribution <- any(diagnostics$VIF > 10)
    }
    if (!is.null(loo)) 
        diagnostics$influence_summary <- do.call(rbind, lapply(cov, function(v) {
            a <- loo[loo$coefficient == v, , drop = FALSE]
            good <- a$state == "ESTIMABLE" & is.finite(a$effect)
            main <- coefficients$effect[match(v, coefficients$coefficient)]
            data.frame(coefficient = v, omissions = nrow(a), valid_omissions = sum(good), invalid_omissions = sum(!good), 
                minimum = if (any(good)) 
                  min(a$effect[good])
                else NA_real_, maximum = if (any(good)) 
                  max(a$effect[good])
                else NA_real_, direction_changes = sum(sign(a$effect[good]) != sign(main)))
        }))
    list(state = "ESTIMABLE", reason = "NONE", biological_n = nrow(d), excluded = data$patient_id[!good], model = model, 
        patient_ids = d$patient_id, coefficients = coefficients, fit = f, hc3 = pe_hc3(f), loo = loo, diagnostics = diagnostics, 
        uncertainty = "OLS_and_HC3_condition_on_fixed_scores_and_estimated_E; no_deconvolution_reference_uncertainty")
}

pe_composition_guard <- function(k_T, N_T, k_N, N_N) {
    s01_empirical_logit(k_T, N_T)
    s01_empirical_logit(k_N, N_N)
    s01_require(length(k_T) == length(k_N), "composition pair dimension")
    s01_require(all(c(N_T, N_N) <= floor(sqrt(2^53))), "composition denominator exceeds exact integer-product range")
    reason <- if (all(c(k_T, k_N) == 0)) 
        "ALL_ZERO_CATEGORY"
    else if (all(k_T == N_T) && all(k_N == N_N)) 
        "ALL_FULL_CATEGORY"
    else if (all(k_T * N_N == k_N * N_T)) 
        "UNCHANGED_RAW_PROPORTIONS_IN_EVERY_PAIR"
    else "NONE"
    list(degenerate = reason != "NONE", state = if (reason != "NONE") "NA_NONESTIMABLE" else "RUNTIME_INFERENCE_PENDING", 
        reason = reason, both_zero_pairs = sum(k_T == 0 & k_N == 0), raw_delta = k_T/N_T - k_N/N_N, family_slot_retained = TRUE)
}

pe_match_overall <- function(scores, donors, overall_donors, background, expected_background) {
    s01_require(identical(background, expected_background), "overall background changed")
    s01_require(all(donors %in% overall_donors), "lineage donors not contained in overall population")
    s01_require(!anyDuplicated(donors), "duplicated matched donors")
    scores[scores$donor %in% donors, , drop = FALSE]
}

pe_sc_members <- function(p, members, structural_reason = NULL) {
    result <- data.frame(member = members, state = "NA_NONESTIMABLE", reason = "MEMBER_DETECTION_OR_MAPPING_FAILED", 
        effect = NA_real_, statistic = NA_real_, numerator_df = NA_real_, denominator_df = NA_real_, p_value = NA_real_, 
        standard_error = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, uncertainty_status = "NA_UNAVAILABLE", 
        uncertainty_reason = "native_modern_edgeR_does_not_provide_coefficient_SE_or_95_percent_CI")
    if (!is.null(structural_reason)) {
        result$reason <- structural_reason
        return(list(results = result, program_gate = FALSE))
    }
    if (length(p$complete_pair_donors) < 8L) {
        result$reason <- "INSUFFICIENT_COMPLETE_DONOR_PAIRS"
        return(list(results = result, program_gate = FALSE))
    }
    md <- p$sample_data
    md$sample_id <- md$unit_id
    md$donor_id <- md$donor
    md$tissue <- ifelse(md$tissue == "Normal", "NAT", md$tissue)
    prepared <- s01_prepare_paired_model(p$counts, md, "donor_id", "donor", TRUE)
    detection <- p$detection$available[match(members, p$detection$gene_id)]
    eligible <- members %in% rownames(p$counts) & !is.na(detection) & detection
    native <- tryCatch(withCallingHandlers({
        dge <- edgeR::normLibSizes(edgeR::DGEList(counts = p$counts), method = "TMM")
        s01_require(all(dge$samples$lib.size == colSums(p$counts)) && all(is.finite(dge$samples$norm.factors) & 
            dge$samples$norm.factors > 0), "invalid native library context")
        fit <- edgeR::glmQLFit(dge, prepared$design, legacy = FALSE, robust = TRUE, dispersion = NULL, top.proportion = NULL, 
            abundance.trend = TRUE, covariate.trend = NULL, prior.count = 0.125)
        list(test = edgeR::glmQLFTest(fit, coef = prepared$coefficient), dge = dge, fit = fit)
    }, warning = function(w) stop(paste("UNEXPECTED_EDGER_WARNING_REQUIRES_REVIEW", conditionMessage(w)))), error = function(e) e)
    if (inherits(native, "error")) {
        if (grepl("UNEXPECTED_EDGER_WARNING", conditionMessage(native))) 
            stop(native)
        result$reason[eligible] <- paste("NATIVE_MODEL_OR_TEST_FAILURE", conditionMessage(native))
        return(list(results = result, program_gate = all(eligible)))
    }
    t <- native$test
    n <- nrow(t$table)
    s01_require(length(t$df.test) %in% c(1L, n) && length(t$df.total) %in% c(1L, n), "invalid native edgeR degrees of freedom")
    for (i in which(eligible)) {
        j <- match(members[i], rownames(t$table))
        v <- c(t$table$logFC[j], t$table$F[j], rep_len(t$df.test, n)[j], rep_len(t$df.total, n)[j], t$table$PValue[j])
        if (all(is.finite(v)) && v[2] >= 0 && all(v[3:4] > 0) && v[5] >= 0 && v[5] <= 1) {
            result[i, c("effect", "statistic", "numerator_df", "denominator_df", "p_value")] <- as.list(v)
            result$state[i] <- "ESTIMABLE"
            result$reason[i] <- "NONE"
        }
        else result$reason[i] <- "INVALID_NATIVE_TEST_OUTPUT"
    }
    list(results = result, program_gate = all(eligible), native = native)
}

pe_protein_ratio <- function(ratios, patient_ids, centered = FALSE) {
    s01_validate_identity_vector(patient_ids, "protein patient")
    s01_require(length(ratios) == length(patient_ids), "protein ratio identity length")
    if (length(ratios) < 8L || any(!is.finite(ratios))) 
        return(list(state = "NA_NONESTIMABLE", reason = "INVALID_COMMON_H4_POPULATION"))
    if (s01_effective_zero(ratios)) 
        return(list(state = "NA_NONESTIMABLE", reason = "CONSTANT_PAIRED_RATIO"))
    f <- t.test(ratios, mu = 0, conf.level = 0.95)
    list(state = "ESTIMABLE", effect = mean(ratios), ci_lower = unname(f$conf.int[1]), ci_upper = unname(f$conf.int[2]), 
        p_value = f$p.value, biological_n = length(ratios), g_av = NA_real_, H4_score_state = "NA_UNAVAILABLE", 
        estimand = if (centered) "mean_centered_log2_tumor_adjacent_ratio" else "mean_log2_tumor_adjacent_ratio")
}

