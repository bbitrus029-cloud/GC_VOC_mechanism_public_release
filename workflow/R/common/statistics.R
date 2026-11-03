# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Numerical methods retained from the study implementation.

s01_nonestimable <- function(extra = list()) {
  c(list(estimable = FALSE, state = "NA_NONESTIMABLE"), extra)
}

s01_effective_zero <- function(x) {
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(x)))
  stats::sd(x) <= tolerance
}

s01_validate_numeric_matrix <- function(x, label, require_counts = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) ||
      is.null(colnames(x))) {
    stop(paste0(label, " must be a numeric features-by-samples matrix with dimnames"),
         call. = FALSE)
  }
  s01_validate_identity_vector(rownames(x), paste0(label, " gene row names"))
  s01_validate_identity_vector(colnames(x), paste0(label, " sample column names"))
  if (any(!is.finite(x))) {
    stop(paste0(label, " contains nonfinite values"), call. = FALSE)
  }
  if (isTRUE(require_counts) &&
      (any(x < 0) || any(x != round(x)))) {
    stop(paste0(label, " must contain finite nonnegative integer counts"),
         call. = FALSE)
  }
  invisible(TRUE)
}

s01_validate_program_members <- function(members) {
  s01_validate_identity_vector(members, "program membership")
}

s01_validate_paired_inputs <- function(tumor, reference, pair_id) {
  if (!is.numeric(tumor) || !is.numeric(reference)) {
    stop("paired tissue vectors must be numeric", call. = FALSE)
  }
  if (length(tumor) != length(reference) || length(pair_id) != length(tumor)) {
    stop("paired tissue vectors and pair identities must have equal length",
         call. = FALSE)
  }
  if (length(pair_id)) s01_validate_identity_vector(pair_id, "pair")
  list(tumor = as.numeric(tumor), reference = as.numeric(reference),
       pair_id = as.character(pair_id))
}

s01_rank_program_score <- function(expression, members, background) {
  if (!is.matrix(expression) || !is.numeric(expression) ||
      is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must be a numeric genes-by-samples matrix with dimnames",
         call. = FALSE)
  }
  s01_validate_identity_vector(rownames(expression), "expression gene row names")
  s01_validate_identity_vector(colnames(expression), "expression sample column names")
  members <- s01_validate_program_members(members)
  background <- s01_validate_identity_vector(background, "rank background gene")
  if (length(background) < 2L) {
    return(s01_nonestimable(list(reason = "background_size_below_2")))
  }
  if (!all(members %in% background) || !all(background %in% rownames(expression))) {
    return(s01_nonestimable(list(reason = "missing_duplicated_or_unmapped_member")))
  }
  values <- expression[background, , drop = FALSE]
  if (any(!is.finite(values)) || any(!is.finite(expression[members, , drop = FALSE]))) {
    return(s01_nonestimable(list(reason = "nonfinite_background_or_required_member")))
  }
  ranks <- apply(values, 2L, rank, ties.method = "average")
  if (is.null(dim(ranks))) ranks <- matrix(ranks, ncol = 1L)
  rownames(ranks) <- background
  colnames(ranks) <- colnames(expression)
  scaled <- (ranks - 1) / (length(background) - 1)
  score <- colMeans(scaled[members, , drop = FALSE])
  list(
    estimable = TRUE,
    state = "ESTIMABLE",
    score = unname(score),
    sample_id = colnames(expression),
    member_count = as.integer(length(members)),
    background_size = as.integer(length(background)),
    tie_method = "average"
  )
}

s01_standardized_member_mean <- function(expression, members) {
  if (!is.matrix(expression) || !is.numeric(expression) ||
      is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must be a numeric genes-by-samples matrix with dimnames",
         call. = FALSE)
  }
  s01_validate_identity_vector(rownames(expression), "expression gene row names")
  s01_validate_identity_vector(colnames(expression), "expression sample column names")
  members <- s01_validate_program_members(members)
  if (!all(members %in% rownames(expression))) {
    return(s01_nonestimable(list(reason = "missing_or_duplicated_member")))
  }
  x <- expression[members, , drop = FALSE]
  if (any(!is.finite(x))) {
    return(s01_nonestimable(list(reason = "nonfinite_required_member")))
  }
  member_sd <- apply(x, 1L, stats::sd)
  tolerance <- apply(x, 1L, function(value) {
    sqrt(.Machine$double.eps) * max(1, max(abs(value)))
  })
  if (any(!is.finite(member_sd)) || any(member_sd <= tolerance)) {
    return(s01_nonestimable(list(reason = "effective_zero_member_variance")))
  }
  standardized <- sweep(x, 1L, rowMeans(x), "-")
  standardized <- sweep(standardized, 1L, member_sd, "/")
  list(
    estimable = TRUE,
    state = "ESTIMABLE",
    score = unname(colMeans(standardized)),
    sample_id = colnames(expression),
    member_count = as.integer(length(members))
  )
}

s01_paired_delta_t <- function(tumor, reference, pair_id) {
  paired <- s01_validate_paired_inputs(tumor, reference, pair_id)
  tumor <- paired$tumor
  reference <- paired$reference
  n <- length(tumor)
  delta <- tumor - reference
  raw_mean <- if (length(delta) && all(is.finite(delta))) mean(delta) else NA_real_
  if (n < 3L) {
    return(s01_nonestimable(list(
      reason = "n_less_than_3", biological_n = as.integer(n), effect = raw_mean,
      statistic = NA_real_, p_value = NA_real_, ci_lower = NA_real_,
      ci_upper = NA_real_
    )))
  }
  if (any(!is.finite(delta))) {
    return(s01_nonestimable(list(
      reason = "nonfinite_formal_pair_difference", biological_n = as.integer(n),
      effect = NA_real_, statistic = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(delta)))
  delta_sd <- stats::sd(delta)
  if (!is.finite(delta_sd) || delta_sd^2 <= tolerance^2) {
    return(s01_nonestimable(list(
      reason = "effective_zero_SD_of_pair_differences",
      biological_n = as.integer(n), effect = raw_mean, statistic = NA_real_,
      p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      sd_tolerance = tolerance
    )))
  }
  standard_error <- delta_sd / sqrt(n)
  statistic <- raw_mean / standard_error
  critical <- stats::qt(0.975, df = n - 1L)
  list(
    estimable = TRUE,
    state = "ESTIMABLE",
    biological_n = as.integer(n),
    effect = raw_mean,
    statistic = statistic,
    degrees_of_freedom = as.integer(n - 1L),
    p_value = 2 * stats::pt(-abs(statistic), df = n - 1L),
    ci_lower = raw_mean - critical * standard_error,
    ci_upper = raw_mean + critical * standard_error,
    sd_tolerance = tolerance
  )
}

s01_paired_g_av <- function(tumor, reference, pair_id) {
  paired <- s01_validate_paired_inputs(tumor, reference, pair_id)
  tumor <- paired$tumor
  reference <- paired$reference
  n <- length(tumor)
  native_difference <- if (n && all(is.finite(c(tumor, reference)))) {
    mean(tumor - reference)
  } else {
    NA_real_
  }
  if (n < 3L || any(!is.finite(c(tumor, reference)))) {
    return(s01_nonestimable(list(
      reason = if (n < 3L) "n_less_than_3" else "nonfinite_required_input",
      biological_n = as.integer(n), native_difference = native_difference,
      g_av = NA_real_, variance = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  sd_tumor <- stats::sd(tumor)
  sd_reference <- stats::sd(reference)
  tumor_tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(tumor)))
  reference_tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(reference)))
  if (sd_tumor <= tumor_tolerance || sd_reference <= reference_tolerance) {
    return(s01_nonestimable(list(
      reason = "nonfinite_correlation_or_effective_zero_tissue_SD",
      biological_n = as.integer(n), native_difference = native_difference,
      g_av = NA_real_, variance = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  correlation <- stats::cor(tumor, reference, method = "pearson")
  if (!is.finite(correlation)) {
    return(s01_nonestimable(list(
      reason = "nonfinite_correlation_or_effective_zero_tissue_SD",
      biological_n = as.integer(n), native_difference = native_difference,
      g_av = NA_real_, variance = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  s_av <- sqrt((sd_tumor^2 + sd_reference^2) / 2)
  d_av <- (mean(tumor) - mean(reference)) / s_av
  m <- 2 * (n - 1) / (1 + correlation^2)
  correction <- exp(
    lgamma(m / 2) - 0.5 * log(m / 2) - lgamma((m - 1) / 2)
  )
  g_av <- correction * d_av
  variance <- 2 * (1 - correlation) / n +
    g_av^2 * (1 + correlation^2) / (4 * n)
  half_width <- stats::qnorm(0.975) * sqrt(variance)
  list(
    estimable = TRUE,
    state = "ESTIMABLE",
    biological_n = as.integer(n),
    native_difference = native_difference,
    correlation = correlation,
    s_av = s_av,
    d_av = d_av,
    m = m,
    correction = correction,
    g_av = g_av,
    variance = variance,
    ci_lower = g_av - half_width,
    ci_upper = g_av + half_width
  )
}

s01_frozen_bh_values <- function(p_values, frozen_family_size, family_hypotheses) {
  if (!is.numeric(p_values) || is.factor(p_values) || is.logical(p_values)) {
    stop("P-value vector must be an actual numeric vector", call. = FALSE)
  }
  if (!is.numeric(frozen_family_size) || length(frozen_family_size) != 1L ||
      is.na(frozen_family_size) || !is.finite(frozen_family_size) ||
      frozen_family_size != round(frozen_family_size) || frozen_family_size < 1) {
    stop("frozen family size must be one positive integer-valued number",
         call. = FALSE)
  }
  frozen_family_size <- as.integer(frozen_family_size)
  family_hypotheses <- s01_validate_identity_vector(
    family_hypotheses, "frozen family hypothesis"
  )
  if (length(family_hypotheses) != frozen_family_size ||
      length(p_values) != frozen_family_size) {
    stop("P-value vector length and identities must exactly equal the frozen family size; missing or extra hypotheses are prohibited",
         call. = FALSE)
  }
  if (is.null(names(p_values))) {
    stop("P-value vector must be named in frozen family order", call. = FALSE)
  }
  observed_names <- s01_validate_identity_vector(names(p_values), "P-value hypothesis")
  if (!setequal(observed_names, family_hypotheses)) {
    stop("P-value hypothesis identity is missing, extra, or substituted",
         call. = FALSE)
  }
  if (!identical(observed_names, family_hypotheses)) {
    stop("P-value hypothesis order differs from the frozen family order",
         call. = FALSE)
  }
  if (any(is.infinite(p_values) | is.nan(p_values))) {
    stop("P values must be finite or structural NA", call. = FALSE)
  }
  finite <- !is.na(p_values)
  if (any(p_values[finite] < 0 | p_values[finite] > 1)) {
    stop("P-value range is outside the closed interval 0 to 1", call. = FALSE)
  }
  result <- rep(NA_real_, length(p_values))
  result[finite] <- stats::p.adjust(
    p_values[finite], method = "BH", n = frozen_family_size
  )
  names(result) <- family_hypotheses
  result
}

s01_frozen_bh <- function(p_values, family_name, config = NULL) {
  if (!is.numeric(p_values) || is.factor(p_values) || is.logical(p_values)) {
    stop("P-value vector must be an actual numeric vector", call. = FALSE)
  }
  if (!is.character(family_name) || length(family_name) != 1L ||
      is.na(family_name) || !nzchar(family_name)) {
    stop("frozen family name must be one nonempty character identity",
         call. = FALSE)
  }
  if (!is.list(config) || !is.list(config$multiplicity$families) ||
      !family_name %in% names(config$multiplicity$families)) {
    stop("frozen family is absent from the authoritative family registry",
         call. = FALSE)
  }
  family <- config$multiplicity$families[[family_name]]
  s01_frozen_bh_values(
    p_values, family$frozen_size, as.character(family$hypotheses)
  )
}

s01_direction_label <- function(value) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    return("NA_NONESTIMABLE")
  }
  if (value > 0) "positive" else if (value < 0) "negative" else "zero"
}

s01_dimension_a <- function(delta_1, delta_2, q_1, q_2,
                            valid_1 = TRUE, valid_2 = TRUE) {
  valid <- isTRUE(valid_1) && isTRUE(valid_2) &&
    all(is.finite(c(delta_1, delta_2, q_1, q_2))) &&
    all(c(q_1, q_2) >= 0 & c(q_1, q_2) <= 1)
  if (!valid) {
    return(list(state = "non_estimable", zero_point_estimate_flag = NA))
  }
  zero <- delta_1 == 0 || delta_2 == 0
  if (zero) {
    return(list(state = "uncertain", zero_point_estimate_flag = TRUE))
  }
  same <- delta_1 * delta_2 > 0
  opposite <- delta_1 * delta_2 < 0
  supported <- c(q_1, q_2) < 0.05
  state <- if (same && all(supported)) {
    "replicated_program_change"
  } else if (same && sum(supported) == 1L) {
    "directionally_consistent_partial_support"
  } else if (same && !any(supported)) {
    "concordant_but_inconclusive"
  } else if (opposite && all(supported)) {
    "opposite_direction_statistical_conflict"
  } else if (opposite && sum(supported) < 2L) {
    "directionally_discordant_inconclusive"
  } else {
    "uncertain"
  }
  list(state = state, zero_point_estimate_flag = FALSE)
}

s01_dimension_b <- function(g_av, ci_lower, ci_upper) {
  bound <- s01_registered_resampling_contract()$absolute_g_av_boundary
  if (length(g_av) != 2L || length(ci_lower) != 2L || length(ci_upper) != 2L) {
    stop("Dimension B requires exactly two active-source estimates", call. = FALSE)
  }
  valid <- all(is.finite(c(g_av, ci_lower, ci_upper))) &&
    all(ci_lower <= ci_upper)
  if (!valid) {
    return(list(
      magnitude_class = "magnitude_nonestimable_g_av",
      both_core_CI_bounded_within_small_effect_region = NA,
      precision_uncertain = NA
    ))
  }
  substantive <- abs(g_av) >= bound
  magnitude <- if (all(substantive)) {
    "substantive_in_both"
  } else if (sum(substantive) == 1L) {
    "substantive_in_one"
  } else {
    "small_point_estimates_in_both"
  }
  bounded <- all(ci_lower > -bound & ci_upper < bound)
  spans_zero <- ci_lower <= 0 & ci_upper >= 0
  reaches_boundary <- ci_lower <= -bound | ci_upper >= bound
  precision_uncertain <- !bounded && any(spans_zero & reaches_boundary)
  list(
    magnitude_class = magnitude,
    both_core_CI_bounded_within_small_effect_region = bounded,
    precision_uncertain = precision_uncertain
  )
}

s01_loo_direction_status <- function(full_effect, loo_effect, loo_estimable) {
  if (!isTRUE(loo_estimable) || !is.finite(loo_effect)) return("NA_NONESTIMABLE")
  if (!is.finite(full_effect)) return("NA_NONESTIMABLE")
  if (full_effect == 0) return("NA_NOT_DIRECTIONAL")
  if (loo_effect == 0 || full_effect * loo_effect < 0) return("FALSE")
  if (full_effect * loo_effect > 0) return("TRUE")
  "NA_NOT_DIRECTIONAL"
}

s01_aggregate_loo_status <- function(status) {
  permitted <- c("TRUE", "FALSE", "NA_NONESTIMABLE", "NA_NOT_DIRECTIONAL")
  if (!length(status) || any(!status %in% permitted)) {
    stop("invalid LOO status population", call. = FALSE)
  }
  for (candidate in c("FALSE", "NA_NONESTIMABLE", "NA_NOT_DIRECTIONAL", "TRUE")) {
    if (candidate %in% status) return(candidate)
  }
  stop("unreachable LOO aggregation state", call. = FALSE)
}

s01_loo_replication_sensitivity <- function(label_changed, estimable) {
  if (length(label_changed) != length(estimable) || !length(label_changed) ||
      anyNA(label_changed) || anyNA(estimable)) {
    stop("LOO replication inputs are malformed", call. = FALSE)
  }
  if (any(label_changed & estimable)) return("TRUE")
  if (any(!estimable)) return("NA_NONESTIMABLE")
  "FALSE"
}

s01_member_direction_concordance <- function(full_effect, member_effects,
                                              fixed_denominator) {
  fixed_denominator <- as.integer(fixed_denominator)
  if (length(member_effects) != fixed_denominator || any(!is.finite(member_effects)) ||
      !is.finite(full_effect)) {
    return(list(state = "NA_NONESTIMABLE", proportion = NA_real_))
  }
  if (full_effect == 0) {
    return(list(state = "NA_NOT_DIRECTIONAL", proportion = NA_real_))
  }
  concordant <- sum(member_effects != 0 & full_effect * member_effects > 0)
  list(
    state = "ESTIMABLE",
    concordant_members = as.integer(concordant),
    denominator = fixed_denominator,
    proportion = concordant / fixed_denominator
  )
}

s01_score_method_sensitivity <- function(primary_effect, alternative_effect,
                                          alternative_estimable) {
  if (length(primary_effect) != length(alternative_effect) ||
      length(primary_effect) != length(alternative_estimable) ||
      !length(primary_effect)) {
    stop("score-method inputs must have equal nonzero length", call. = FALSE)
  }
  nonestimable <- !alternative_estimable | !is.finite(alternative_effect) |
    !is.finite(primary_effect)
  estimable <- !nonestimable
  opposite <- estimable & primary_effect != 0 & alternative_effect != 0 &
    primary_effect * alternative_effect < 0
  exact_zero <- any(estimable & (primary_effect == 0 | alternative_effect == 0))
  status <- if (any(opposite)) {
    "TRUE"
  } else if (any(nonestimable)) {
    "NA_NONESTIMABLE"
  } else {
    "FALSE"
  }
  list(
    primary_rank_unrounded_raw_paired_mean_difference = primary_effect,
    standardized_member_mean_unrounded_raw_paired_mean_difference = alternative_effect,
    primary_rank_direction = vapply(primary_effect, s01_direction_label, character(1)),
    standardized_member_mean_direction = vapply(
      alternative_effect, s01_direction_label, character(1)
    ),
    score_method_sensitive = status,
    score_method_nonestimable = any(nonestimable),
    exact_zero_direction_flag = exact_zero
  )
}

s01_paired_wilcoxon <- function(tumor, reference, pair_id) {
  paired <- s01_validate_paired_inputs(tumor, reference, pair_id)
  tumor <- paired$tumor
  reference <- paired$reference
  if (length(tumor) < 1L || any(!is.finite(c(tumor, reference)))) {
    return(s01_nonestimable(list(
      reason = "nonfinite_required_input", statistic = NA_real_, p_value = NA_real_
    )))
  }
  if (all(tumor - reference == 0)) {
    return(s01_nonestimable(list(
      reason = "all_zero_paired_differences_nonfinite_statistic_or_p_value",
      statistic = NA_real_, p_value = NA_real_
    )))
  }
  result <- suppressWarnings(stats::wilcox.test(
    tumor, reference, paired = TRUE, exact = FALSE, correct = FALSE,
    alternative = "two.sided"
  ))
  statistic <- unname(result$statistic)
  p_value <- result$p.value
  if (!is.finite(statistic) || !is.finite(p_value)) {
    return(s01_nonestimable(list(
      reason = "nonfinite_wilcoxon_statistic_or_p_value",
      statistic = NA_real_, p_value = NA_real_
    )))
  }
  list(
    estimable = TRUE, state = "ESTIMABLE",
    statistic = statistic, p_value = p_value,
    multiplicity_status = "NO_FORMAL_FAMILY"
  )
}

s01_with_production_rng <- function(seed, expression) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != round(seed) || seed < 1 ||
      seed > .Machine$integer.max) {
    stop("seed must be one valid non-missing R integer", call. = FALSE)
  }
  seed <- as.integer(seed)
  previous_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  force(expression)
}

s01_registered_resampling_contract <- function() {
  list(
    rng_kind = "L'Ecuyer-CMRG", replicates = 10000L,
    minimum_valid_replicates = 9500L, quantile_type = 7L,
    absolute_g_av_boundary = 0.20
  )
}

s01_pair_bootstrap_internal <- function(tumor, reference, pair_id, seed,
                                         replicates, minimum_valid) {
  paired <- s01_validate_paired_inputs(tumor, reference, pair_id)
  tumor <- paired$tumor
  reference <- paired$reference
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      !is.finite(replicates) || replicates != round(replicates) || replicates < 1 ||
      !is.numeric(minimum_valid) || length(minimum_valid) != 1L ||
      !is.finite(minimum_valid) || minimum_valid != round(minimum_valid) ||
      minimum_valid < 1 || minimum_valid > replicates) {
    stop("bootstrap replicate contract is invalid", call. = FALSE)
  }
  replicates <- as.integer(replicates)
  minimum_valid <- as.integer(minimum_valid)
  if (length(tumor) < 3L || any(!is.finite(c(tumor, reference)))) {
    return(s01_nonestimable(list(
      valid_replicates = 0L, ci_lower = NA_real_, ci_upper = NA_real_,
      rng_kind = "L'Ecuyer-CMRG"
    )))
  }
  if (replicates < 1L || minimum_valid < 1L || minimum_valid > replicates) {
    stop("bootstrap replicate contract is invalid", call. = FALSE)
  }
  delta <- tumor - reference
  boot <- s01_with_production_rng(seed, {
    vapply(seq_len(replicates), function(i) {
      mean(delta[sample.int(length(delta), length(delta), replace = TRUE)])
    }, numeric(1))
  })
  finite <- boot[is.finite(boot)]
  valid <- length(finite)
  if (valid < minimum_valid) {
    return(s01_nonestimable(list(
      valid_replicates = as.integer(valid), ci_lower = NA_real_,
      ci_upper = NA_real_, rng_kind = "L'Ecuyer-CMRG"
    )))
  }
  interval <- stats::quantile(
    finite, probs = c(0.025, 0.975), type = 7L, names = FALSE
  )
  list(
    estimable = TRUE, state = "ESTIMABLE",
    valid_replicates = as.integer(valid), ci_lower = interval[[1L]],
    ci_upper = interval[[2L]], quantile_type = 7L,
    rng_kind = "L'Ecuyer-CMRG"
  )
}

s01_pair_bootstrap <- function(tumor, reference, pair_id, seed) {
  contract <- s01_registered_resampling_contract()
  s01_pair_bootstrap_internal(
    tumor, reference, pair_id, seed, contract$replicates,
    contract$minimum_valid_replicates
  )
}

s01_marker_audit <- function(scores, expected_label) {
  panels <- c(
    "epithelial", "fibroblast", "endothelial", "myeloid",
    "T_NK", "B_plasma", "mast_other_immune"
  )
  required <- c(
    "author_label", "donor_id", "cell_count", "total_umi",
    "marker_mapping_complete", panels
  )
  if (!is.data.frame(scores) || !all(required %in% names(scores)) ||
      !expected_label %in% panels) {
    stop("marker-audit input schema is invalid", call. = FALSE)
  }
  if (anyNA(scores$donor_id)) stop("donor identity contains missing/NA values", call. = FALSE)
  if (any(!nzchar(trimws(as.character(scores$donor_id))))) {
    stop("donor identity contains blank values", call. = FALSE)
  }
  s01_assert_unique_key(scores, c("donor_id", "author_label"), "donor by author label")
  author_label <- as.character(scores$author_label)
  if (anyNA(author_label) || any(!nzchar(trimws(author_label)))) {
    stop("author-label identity contains missing or blank values", call. = FALSE)
  }
  if (any(!author_label %in% panels)) {
    stop("author-label state is outside the exact seven-label vocabulary",
         call. = FALSE)
  }
  integer_count <- function(x) {
    is.numeric(x) && !is.factor(x) && !is.logical(x) && !anyNA(x) &&
      all(is.finite(x)) && all(x >= 0) && all(x == round(x))
  }
  if (!integer_count(scores$cell_count)) {
    stop("cell count must contain nonnegative integer values", call. = FALSE)
  }
  if (!integer_count(scores$total_umi)) {
    stop("total UMI count must contain nonnegative integer values", call. = FALSE)
  }
  if (!is.logical(scores$marker_mapping_complete) ||
      anyNA(scores$marker_mapping_complete)) {
    stop("marker mapping flags must be exact non-missing logical values",
         call. = FALSE)
  }
  if (!all(vapply(scores[, panels, drop = FALSE], function(x) {
    is.numeric(x) && !is.factor(x) && !is.logical(x)
  }, logical(1)))) {
    stop("marker panel scores must be numeric", call. = FALSE)
  }
  eligible <- scores$author_label == expected_label &
    scores$cell_count >= 20L & scores$total_umi >= 10000L
  x <- scores[eligible, , drop = FALSE]
  if (nrow(x)) s01_validate_identity_vector(x$donor_id, "eligible donor")
  if (nrow(x) < 8L || any(!x$marker_mapping_complete) ||
      any(!is.finite(as.matrix(x[, panels, drop = FALSE])))) {
    return(list(
      state = "NA_NONESTIMABLE", eligible_donor_units = as.integer(nrow(x)),
      expected_selection_fraction = NA_real_, unit_adjudication = character(),
      median_rule_pass = NA
    ))
  }
  values <- as.matrix(x[, panels, drop = FALSE])
  selection <- apply(values, 1L, function(row) {
    winners <- panels[row == max(row)]
    if (length(winners) == 1L) winners else "AMBIGUOUS"
  })
  fraction <- sum(selection == expected_label) / nrow(x)
  medians <- apply(values, 2L, stats::median)
  median_pass <- all(medians[[expected_label]] > medians[panels != expected_label])
  list(
    state = if (fraction >= 0.60 && median_pass) "PASS" else "FAIL",
    eligible_donor_units = as.integer(nrow(x)),
    expected_selection_fraction = fraction,
    unit_adjudication = unname(selection),
    donor_balanced_median = unname(medians),
    median_rule_pass = median_pass
  )
}

s01_propagate_marker_states <- function(label_states) {
  required <- c(
    "epithelial", "fibroblast", "endothelial", "myeloid",
    "T_NK", "B_plasma", "mast_other_immune"
  )
  if (is.null(names(label_states)) || !identical(names(label_states), required) ||
      anyNA(label_states) ||
      any(!label_states %in% c("PASS", "FAIL", "NA_NONESTIMABLE"))) {
    stop("marker label-state population is invalid", call. = FALSE)
  }
  eligible <- function(label) if (label_states[[label]] == "PASS") {
    "ELIGIBLE"
  } else {
    "NA_NONESTIMABLE"
  }
  c(
    myeloid_G = eligible("myeloid"),
    myeloid_H = eligible("myeloid"),
    fibroblast_H = eligible("fibroblast"),
    B_plasma_H = eligible("B_plasma"),
    epithelial_H = eligible("epithelial")
  )
}

s01_composition_gate <- function(label_states) {
  required <- c(
    "epithelial", "fibroblast", "endothelial", "myeloid",
    "T_NK", "B_plasma", "mast_other_immune"
  )
  valid <- !is.null(names(label_states)) && identical(names(label_states), required) &&
    !anyNA(label_states) &&
    all(label_states %in% c("PASS", "FAIL", "NA_NONESTIMABLE"))
  if (!valid) stop("composition gate requires exactly seven label states", call. = FALSE)
  pass <- all(label_states[required] == "PASS")
  list(
    estimable = pass,
    state = if (pass) "ELIGIBLE" else "NA_NONESTIMABLE",
    family_states = setNames(rep(
      if (pass) "ELIGIBLE" else "NA_NONESTIMABLE", 4L
    ), c(
      "SC_COMPOSITION_PRIMARY", "SC_COMPOSITION_STRICT200",
      "SC_COMPOSITION_CLR", "SC_COMPOSITION_QUASIBINOMIAL"
    )),
    descriptive_author_label_composition_only = !pass
  )
}

s01_empirical_logit <- function(k, n) {
  if (!is.numeric(k) || !is.numeric(n) || is.factor(k) || is.factor(n) ||
      is.logical(k) || is.logical(n) || length(k) != length(n) ||
      any(!is.finite(c(k, n))) || any(k != round(k) | n != round(n))) {
    stop("empirical-logit counts must be finite integer-valued numeric vectors",
         call. = FALSE)
  }
  if (any(k < 0 | n <= 0 | k > n)) {
    stop("empirical-logit totals must be positive with nonnegative counts",
         call. = FALSE)
  }
  log((k + 0.5) / (n - k + 0.5))
}

s01_clr <- function(counts) {
  if (!is.numeric(counts) || is.factor(counts) || is.logical(counts) ||
      length(counts) != 7L || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != round(counts))) {
    stop("CLR requires exactly seven finite nonnegative integer compartment counts",
         call. = FALSE)
  }
  if (sum(counts) <= 0) stop("CLR requires a positive total", call. = FALSE)
  logged <- log(counts + 0.5)
  logged - mean(logged)
}

s01_adjudicate_quasibinomial <- function(facts, warnings) {
  recognized <- c(
    "glm.fit: algorithm did not converge",
    "glm.fit: algorithm stopped at boundary value"
  )
  if (!is.character(warnings) || anyNA(warnings)) {
    stop("quasibinomial warning audit must contain exact condition messages",
         call. = FALSE)
  }
  unexpected <- setdiff(unique(warnings), recognized)
  if (length(unexpected)) {
    stop(
      paste0("unexpected quasibinomial warning: ", paste(unexpected, collapse = " | ")),
      call. = FALSE
    )
  }
  required <- c(
    "converged", "boundary", "rank", "columns", "fitted",
    "coefficient", "se", "dispersion", "df_residual"
  )
  if (!is.list(facts) || !all(required %in% names(facts))) {
    stop("quasibinomial adjudication facts are incomplete", call. = FALSE)
  }
  if (!is.logical(facts$converged) || length(facts$converged) != 1L ||
      is.na(facts$converged) || !is.logical(facts$boundary) ||
      length(facts$boundary) != 1L || is.na(facts$boundary)) {
    stop("quasibinomial convergence and boundary facts must be exact logical scalars",
         call. = FALSE)
  }
  integer_fact <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      !is.nan(value) && value == round(value)
  }
  if (!integer_fact(facts$rank) || !integer_fact(facts$columns) ||
      !integer_fact(facts$df_residual) || facts$rank < 0 || facts$columns < 1) {
    stop("quasibinomial rank, model columns, and residual df facts are malformed",
         call. = FALSE)
  }
  if (!is.numeric(facts$fitted) || !length(facts$fitted) ||
      is.factor(facts$fitted) || is.logical(facts$fitted) ||
      !is.numeric(facts$coefficient) || length(facts$coefficient) != 1L ||
      !is.numeric(facts$se) || length(facts$se) != 1L ||
      !is.numeric(facts$dispersion) || length(facts$dispersion) != 1L) {
    stop("quasibinomial fitted values or coefficient facts are malformed",
         call. = FALSE)
  }
  epsilon <- 10 * .Machine$double.eps
  flags <- c(
    recognized_warning = any(warnings %in% recognized),
    convergence_failure = !isTRUE(facts$converged),
    boundary_fit = isTRUE(facts$boundary),
    rank_deficient = !is.finite(facts$rank) || !is.finite(facts$columns) ||
      facts$rank < facts$columns,
    fitted_probability_boundary = any(!is.finite(facts$fitted)) ||
      any(facts$fitted < epsilon | facts$fitted > 1 - epsilon),
    nonfinite_coefficient = !is.finite(facts$coefficient),
    nonfinite_se = !is.finite(facts$se),
    nonfinite_dispersion = !is.finite(facts$dispersion),
    nonpositive_df_residual = !is.finite(facts$df_residual) ||
      facts$df_residual <= 0
  )
  estimable <- !any(flags)
  list(
    estimable = estimable,
    state = if (estimable) "ESTIMABLE" else "NA_NONESTIMABLE",
    warnings = warnings,
    flags = flags,
    audit_trail = list(
      warning_messages = warnings,
      recognized_warning_present = any(warnings %in% recognized),
      unexpected_warning_present = FALSE,
      converged = facts$converged, boundary = facts$boundary,
      model_rank = as.integer(facts$rank),
      model_columns = as.integer(facts$columns),
      fitted_probability_boundary = unname(flags[["fitted_probability_boundary"]]),
      coefficient = facts$coefficient, standard_error = facts$se,
      dispersion = facts$dispersion, residual_df = as.integer(facts$df_residual)
    )
  )
}

s01_fit_quasibinomial <- function(data) {
  required <- c("donor", "tissue", "compartment_count", "total_count")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("quasibinomial input schema is invalid", call. = FALSE)
  }
  numeric_counts <- is.numeric(data$compartment_count) &&
    is.numeric(data$total_count) && !is.factor(data$compartment_count) &&
    !is.factor(data$total_count) && !is.logical(data$compartment_count) &&
    !is.logical(data$total_count)
  if (!numeric_counts || any(!is.finite(data$compartment_count)) ||
      any(!is.finite(data$total_count)) ||
      any(data$compartment_count != round(data$compartment_count)) ||
      any(data$total_count != round(data$total_count)) ||
      any(data$compartment_count < 0 | data$total_count <= 0 |
          data$total_count < data$compartment_count)) {
    stop("quasibinomial counts are invalid", call. = FALSE)
  }
  donor <- as.character(data$donor)
  tissue <- as.character(data$tissue)
  if (anyNA(donor)) stop("donor identity contains missing/NA values", call. = FALSE)
  if (any(!nzchar(trimws(donor)))) stop("donor identity contains blank values", call. = FALSE)
  if (anyNA(tissue) || any(!nzchar(trimws(tissue)))) {
    stop("tissue identity contains missing or blank values", call. = FALSE)
  }
  if (!setequal(unique(tissue), c("NAT", "Tumor"))) {
    stop("quasibinomial tissue levels must be exactly NAT and Tumor", call. = FALSE)
  }
  key <- paste(data$donor, data$tissue, sep = "\r")
  if (anyDuplicated(key) || any(table(data$donor) != 2L)) {
    stop("quasibinomial donor pairing is invalid", call. = FALSE)
  }
  data <- data[order(data$donor, match(data$tissue, c("NAT", "Tumor")),
                     method = "radix"), , drop = FALSE]
  data$donor <- factor(data$donor, levels = sort(unique(data$donor), method = "radix"))
  data$tissue <- factor(data$tissue, levels = c("NAT", "Tumor"))
  warnings <- character()
  fit <- withCallingHandlers(
    stats::glm(
      cbind(compartment_count, total_count - compartment_count) ~ tissue + donor,
      data = data,
      family = stats::quasibinomial(link = "logit"),
      control = stats::glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  fit_summary <- summary(fit)
  coefficient_name <- "tissueTumor"
  coefficient <- fit_summary$coefficients[coefficient_name, "Estimate"]
  standard_error <- fit_summary$coefficients[coefficient_name, "Std. Error"]
  facts <- list(
    converged = fit$converged,
    boundary = fit$boundary,
    rank = fit$rank,
    columns = ncol(stats::model.matrix(fit)),
    fitted = stats::fitted(fit),
    coefficient = coefficient,
    se = standard_error,
    dispersion = fit_summary$dispersion,
    df_residual = fit$df.residual
  )
  adjudication <- s01_adjudicate_quasibinomial(facts, warnings)
  if (!adjudication$estimable) {
    return(c(adjudication, list(
      coefficient_name = coefficient_name, effect = NA_real_,
      standard_error = NA_real_, statistic = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, df_residual = fit$df.residual
    )))
  }
  statistic <- coefficient / standard_error
  critical <- stats::qt(0.975, df = fit$df.residual)
  c(adjudication, list(
    coefficient_name = coefficient_name,
    effect = coefficient,
    standard_error = standard_error,
    statistic = statistic,
    p_value = 2 * stats::pt(-abs(statistic), df = fit$df.residual),
    ci_lower = coefficient - critical * standard_error,
    ci_upper = coefficient + critical * standard_error,
    df_residual = as.integer(fit$df.residual)
  ))
}

s01_common_complete_pairs <- function(data, members, expected_n = NULL) {
  members <- s01_validate_identity_vector(members, "protein program member")
  required <- c("patient_id", "tissue", members)
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("protein common-pair input schema is invalid", call. = FALSE)
  }
  patient_id <- as.character(data$patient_id)
  tissue <- as.character(data$tissue)
  if (anyNA(patient_id)) {
    stop("patient identity contains missing/NA values", call. = FALSE)
  }
  if (any(!nzchar(trimws(patient_id)))) {
    stop("patient identity contains blank values", call. = FALSE)
  }
  if (anyNA(tissue) || any(!nzchar(trimws(tissue))) ||
      !setequal(unique(tissue), c("NAT", "Tumor"))) {
    stop("protein tissues must be exactly NAT and Tumor", call. = FALSE)
  }
  key <- paste(data$patient_id, data$tissue, sep = "\r")
  if (anyDuplicated(key)) stop("duplicate patient-tissue row", call. = FALSE)
  split_rows <- split(data, data$patient_id)
  complete <- vapply(split_rows, function(x) {
    nrow(x) == 2L && setequal(x$tissue, c("NAT", "Tumor")) &&
      all(is.finite(as.matrix(x[, members, drop = FALSE])))
  }, logical(1))
  ids <- sort(names(split_rows)[complete], method = "radix")
  if (!is.null(expected_n) && length(ids) != as.integer(expected_n)) {
    stop("common complete-pair count does not match the frozen expectation",
         call. = FALSE)
  }
  ids
}

s01_validate_scalar_number <- function(value, label, integer = FALSE,
                                        lower = -Inf, upper = Inf,
                                        allow_na = FALSE) {
  valid_type <- is.numeric(value) && !is.factor(value) && !is.logical(value) &&
    length(value) == 1L
  if (!valid_type || is.nan(value) || is.infinite(value) ||
      (!allow_na && is.na(value))) {
    stop(paste0(label, " must be one valid numeric scalar"), call. = FALSE)
  }
  if (is.na(value)) return(value)
  if (integer && value != round(value)) {
    stop(paste0(label, " must be integer-valued"), call. = FALSE)
  }
  if (value < lower || value > upper) {
    stop(paste0(label, " is outside its registered domain"), call. = FALSE)
  }
  value
}

s01_bridge_state_internal <- function(n, rho, p_value, ci_lower, ci_upper, q_value,
                                       valid_bootstrap_replicates,
                                       minimum_valid, maximum_replicates,
                                       required_input_valid = TRUE) {
  if (!is.logical(required_input_valid) || length(required_input_valid) != 1L ||
      is.na(required_input_valid)) {
    stop("required-input validity flag must be one exact logical value",
         call. = FALSE)
  }
  n <- s01_validate_scalar_number(n, "bridge n", TRUE, 0)
  valid_bootstrap_replicates <- s01_validate_scalar_number(
    valid_bootstrap_replicates, "valid bootstrap replicate count", TRUE,
    0, maximum_replicates
  )
  if (!is.numeric(rho) || is.factor(rho) || length(rho) != 1L) {
    stop("bridge rho must be one numeric scalar", call. = FALSE)
  }
  if (is.finite(rho)) s01_validate_scalar_number(rho, "bridge rho", FALSE, -1, 1)
  nonestimable <- n < 15L || !required_input_valid || !is.finite(rho) ||
    valid_bootstrap_replicates < minimum_valid
  allow_na <- nonestimable
  structural_missing <- function(value) {
    if (nonestimable && is.numeric(value) && !is.factor(value) &&
        length(value) == 1L && !is.finite(value)) NA_real_ else value
  }
  p_value <- structural_missing(p_value)
  q_value <- structural_missing(q_value)
  ci_lower <- structural_missing(ci_lower)
  ci_upper <- structural_missing(ci_upper)
  p_value <- s01_validate_scalar_number(
    p_value, "bridge nominal P value", FALSE, 0, 1, allow_na
  )
  q_value <- s01_validate_scalar_number(
    q_value, "bridge q value", FALSE, 0, 1, allow_na
  )
  ci_lower <- s01_validate_scalar_number(
    ci_lower, "bridge primary CI lower", FALSE, -1, 1, allow_na
  )
  ci_upper <- s01_validate_scalar_number(
    ci_upper, "bridge primary CI upper", FALSE, -1, 1, allow_na
  )
  if (is.finite(ci_lower) && is.finite(ci_upper) && ci_lower > ci_upper) {
    stop("bridge primary interval is not ordered", call. = FALSE)
  }
  n <- as.integer(n)
  valid_bootstrap_replicates <- as.integer(valid_bootstrap_replicates)
  if (nonestimable) {
    return(list(
      state = "bridge_nonestimable",
      bridge_q_ci_inconsistent = NA,
      bridge_nominal_p_ci_inconsistent = NA
    ))
  }
  excludes_zero <- ci_lower > 0 || ci_upper < 0
  state <- if (q_value < 0.05 && excludes_zero) {
    "bridge_supported"
  } else {
    "bridge_uncertain"
  }
  list(
    state = state,
    bridge_q_ci_inconsistent = xor(q_value < 0.05, excludes_zero),
    bridge_nominal_p_ci_inconsistent = xor(p_value < 0.05, excludes_zero)
  )
}

s01_bridge_state <- function(n, rho, p_value, ci_lower, ci_upper, q_value,
                              valid_bootstrap_replicates,
                              required_input_valid = TRUE) {
  contract <- s01_registered_resampling_contract()
  s01_bridge_state_internal(
    n, rho, p_value, ci_lower, ci_upper, q_value,
    valid_bootstrap_replicates, contract$minimum_valid_replicates,
    contract$replicates, required_input_valid
  )
}

s01_bridge_inference_internal <- function(rna_delta, protein_delta, patient_id,
                                           q_value, seed, replicates,
                                           minimum_valid) {
  if (!is.numeric(rna_delta) || !is.numeric(protein_delta) ||
      length(patient_id) != length(rna_delta)) {
    stop("bridge vectors and patient identities must have equal valid lengths",
         call. = FALSE)
  }
  patient_id <- s01_validate_identity_vector(patient_id, "bridge patient")
  q_value <- s01_validate_scalar_number(
    q_value, "bridge q value", FALSE, 0, 1, allow_na = TRUE
  )
  rna_delta <- as.numeric(rna_delta)
  protein_delta <- as.numeric(protein_delta)
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      !is.finite(replicates) || replicates != round(replicates) || replicates < 1 ||
      !is.numeric(minimum_valid) || length(minimum_valid) != 1L ||
      !is.finite(minimum_valid) || minimum_valid != round(minimum_valid) ||
      minimum_valid < 1 || minimum_valid > replicates) {
    stop("bridge bootstrap replicate contract is invalid", call. = FALSE)
  }
  replicates <- as.integer(replicates)
  minimum_valid <- as.integer(minimum_valid)
  n <- length(rna_delta)
  input_valid <- length(protein_delta) == n && n >= 15L &&
    all(is.finite(c(rna_delta, protein_delta)))
  if (!input_valid) {
    state <- s01_bridge_state_internal(
      n, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 0L,
      minimum_valid, replicates, FALSE
    )
    return(c(state, list(
      biological_n = as.integer(n), rho = NA_real_, p_value = NA_real_,
      q_value = q_value, ci_lower = NA_real_, ci_upper = NA_real_,
      bca_ci_lower = NA_real_, bca_ci_upper = NA_real_,
      valid_bootstrap_replicates = 0L, quantile_type = 7L,
      rng_kind = "L'Ecuyer-CMRG"
    )))
  }
  if (replicates < 1L || minimum_valid < 1L || minimum_valid > replicates) {
    stop("bridge bootstrap replicate contract is invalid", call. = FALSE)
  }
  test <- suppressWarnings(stats::cor.test(
    rna_delta, protein_delta, method = "spearman", exact = FALSE,
    alternative = "two.sided"
  ))
  if (!is.finite(unname(test$estimate))) {
    state <- s01_bridge_state_internal(
      n, unname(test$estimate), NA_real_, NA_real_, NA_real_, q_value, 0L,
      minimum_valid, replicates, FALSE
    )
    return(c(state, list(
      biological_n = as.integer(n), rho = NA_real_, p_value = NA_real_,
      q_value = q_value, ci_lower = NA_real_, ci_upper = NA_real_,
      bca_ci_lower = NA_real_, bca_ci_upper = NA_real_,
      valid_bootstrap_replicates = 0L, quantile_type = 7L,
      rng_kind = "L'Ecuyer-CMRG"
    )))
  }
  data <- data.frame(rna = rna_delta, protein = protein_delta)
  boot_object <- s01_with_production_rng(seed, {
    boot::boot(
      data,
      statistic = function(x, index) {
        suppressWarnings(stats::cor(
          x$rna[index], x$protein[index], method = "spearman"
        ))
      },
      R = replicates,
      sim = "ordinary",
      stype = "i"
    )
  })
  finite <- as.numeric(boot_object$t)
  finite <- finite[is.finite(finite)]
  valid <- length(finite)
  interval <- if (valid >= minimum_valid) {
    stats::quantile(finite, c(0.025, 0.975), type = 7L, names = FALSE)
  } else {
    c(NA_real_, NA_real_)
  }
  bca <- c(NA_real_, NA_real_)
  if (valid >= minimum_valid) {
    bca_attempt <- suppressWarnings(try(
      boot::boot.ci(boot_object, conf = 0.95, type = "bca"), silent = TRUE
    ))
    if (!inherits(bca_attempt, "try-error") && !is.null(bca_attempt$bca) &&
        length(bca_attempt$bca) >= 5L && all(is.finite(bca_attempt$bca[4:5]))) {
      bca <- as.numeric(bca_attempt$bca[4:5])
    }
  }
  state <- s01_bridge_state_internal(
    n = n,
    rho = unname(test$estimate),
    p_value = test$p.value,
    ci_lower = interval[[1L]],
    ci_upper = interval[[2L]],
    q_value = q_value,
    valid_bootstrap_replicates = valid,
    minimum_valid = minimum_valid, maximum_replicates = replicates,
    required_input_valid = input_valid
  )
  c(state, list(
    biological_n = as.integer(n),
    rho = unname(test$estimate),
    p_value = test$p.value,
    q_value = q_value,
    ci_lower = interval[[1L]],
    ci_upper = interval[[2L]],
    bca_ci_lower = bca[[1L]],
    bca_ci_upper = bca[[2L]],
    valid_bootstrap_replicates = as.integer(valid),
    quantile_type = 7L,
    rng_kind = "L'Ecuyer-CMRG"
  ))
}

s01_bridge_inference <- function(rna_delta, protein_delta, patient_id, q_value,
                                  seed) {
  contract <- s01_registered_resampling_contract()
  s01_bridge_inference_internal(
    rna_delta, protein_delta, patient_id, q_value, seed,
    contract$replicates, contract$minimum_valid_replicates
  )
}

s01_cross_layer_state <- function(bulk_dimension_a, bulk_consensus_direction,
                                   protein_estimable, protein_q,
                                   protein_effect) {
  states <- c(
    "non_estimable", "uncertain", "replicated_program_change",
    "directionally_consistent_partial_support", "concordant_but_inconclusive",
    "opposite_direction_statistical_conflict",
    "directionally_discordant_inconclusive"
  )
  if (!is.character(bulk_dimension_a) || length(bulk_dimension_a) != 1L ||
      is.na(bulk_dimension_a) || !bulk_dimension_a %in% states) {
    stop("bulk Dimension A state is outside the registered vocabulary",
         call. = FALSE)
  }
  if (!is.logical(protein_estimable) || length(protein_estimable) != 1L ||
      is.na(protein_estimable)) {
    stop("protein estimability must be one exact logical value", call. = FALSE)
  }
  if (!is.numeric(bulk_consensus_direction) ||
      length(bulk_consensus_direction) != 1L || is.nan(bulk_consensus_direction) ||
      is.infinite(bulk_consensus_direction) ||
      (!is.na(bulk_consensus_direction) &&
       !bulk_consensus_direction %in% c(-1, 0, 1))) {
    stop("bulk consensus direction is outside its registered domain",
         call. = FALSE)
  }
  if (protein_estimable) {
    s01_validate_scalar_number(protein_q, "protein q probability", FALSE, 0, 1)
    s01_validate_scalar_number(protein_effect, "protein effect", FALSE)
  } else {
    s01_validate_scalar_number(protein_q, "protein q probability", FALSE, 0, 1, TRUE)
    s01_validate_scalar_number(protein_effect, "protein effect", FALSE,
                               -Inf, Inf, TRUE)
  }
  if (!isTRUE(protein_estimable)) return("non_estimable")
  consensus <- is.numeric(bulk_consensus_direction) &&
    length(bulk_consensus_direction) == 1L &&
    is.finite(bulk_consensus_direction) && bulk_consensus_direction != 0
  protein_directional <- is.numeric(protein_effect) &&
    length(protein_effect) == 1L && is.finite(protein_effect) && protein_effect != 0
  same <- consensus && protein_directional &&
    bulk_consensus_direction * protein_effect > 0
  opposite <- consensus && protein_directional &&
    bulk_consensus_direction * protein_effect < 0
  supported_bulk <- bulk_dimension_a %in% c(
    "replicated_program_change", "directionally_consistent_partial_support"
  )
  if (supported_bulk && same && is.finite(protein_q) && protein_q < 0.05) {
    return("orthogonally_supported")
  }
  if (same) return("directionally_convergent")
  if (opposite) return("cross_layer_conflict")
  "indeterminate"
}

s01_compensation_compatible <- function(g_state, g_direction,
                                         h_state, h_direction) {
  allowed <- c(
    "replicated_program_change", "directionally_consistent_partial_support"
  )
  all_states <- c(
    "non_estimable", "uncertain", "replicated_program_change",
    "directionally_consistent_partial_support", "concordant_but_inconclusive",
    "opposite_direction_statistical_conflict",
    "directionally_discordant_inconclusive"
  )
  if (!is.character(g_state) || length(g_state) != 1L || is.na(g_state) ||
      !g_state %in% all_states || !is.character(h_state) || length(h_state) != 1L ||
      is.na(h_state) || !h_state %in% all_states) {
    stop("compensation input state is outside the registered vocabulary",
         call. = FALSE)
  }
  for (direction in list(g_direction, h_direction)) {
    if (!is.numeric(direction) || length(direction) != 1L || is.nan(direction) ||
        is.infinite(direction)) {
      stop("compensation direction is outside its registered domain",
           call. = FALSE)
    }
  }
  g_state %in% allowed && h_state %in% allowed &&
    is.finite(g_direction) && is.finite(h_direction) &&
    g_direction > 0 && h_direction > 0
}

s01_read_rank_fixture <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("synthetic rank fixture is missing", call. = FALSE)
  }
  fixture <- utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = character()
  )
  if (!is.data.frame(fixture) || ncol(fixture) < 3L ||
      !identical(names(fixture)[[1L]], "gene_id")) {
    stop("synthetic rank fixture schema is invalid", call. = FALSE)
  }
  gene_id <- s01_validate_identity_vector(fixture$gene_id, "fixture gene")
  sample_id <- s01_validate_identity_vector(names(fixture)[-1L], "fixture sample")
  values <- lapply(fixture[-1L], function(column) {
    converted <- suppressWarnings(as.numeric(column))
    if (length(converted) != length(column) || any(!is.finite(converted))) {
      stop("synthetic rank fixture contains nonfinite or nonnumeric values",
           call. = FALSE)
    }
    converted
  })
  result <- as.matrix(as.data.frame(values, check.names = FALSE))
  rownames(result) <- gene_id
  colnames(result) <- sample_id
  storage.mode(result) <- "double"
  result
}

s01_prepare_paired_model <- function(values, sample_data, unit_column,
                                     unit_label, require_counts = FALSE) {
  s01_validate_numeric_matrix(values, "model matrix", require_counts)
  required <- c("sample_id", unit_column, "tissue")
  if (!is.data.frame(sample_data) || !all(required %in% names(sample_data)) ||
      nrow(sample_data) != ncol(values)) {
    stop("model sample metadata schema or length is invalid", call. = FALSE)
  }
  sample_id <- s01_validate_identity_vector(sample_data$sample_id, "sample")
  if (!identical(sample_id, colnames(values))) {
    stop("model sample identities or order do not match matrix columns", call. = FALSE)
  }
  unit <- as.character(sample_data[[unit_column]])
  tissue <- as.character(sample_data$tissue)
  if (anyNA(unit)) stop(paste0(unit_label, " identity contains missing/NA values"),
                        call. = FALSE)
  if (any(!nzchar(trimws(unit)))) {
    stop(paste0(unit_label, " identity contains blank values"), call. = FALSE)
  }
  if (anyNA(tissue) || any(!nzchar(trimws(tissue)))) {
    stop("tissue identity contains missing or blank values", call. = FALSE)
  }
  if (!setequal(unique(tissue), c("NAT", "Tumor"))) {
    stop("model tissue levels must be exactly NAT and Tumor", call. = FALSE)
  }
  key <- paste(unit, tissue, sep = "\r")
  if (anyDuplicated(key)) {
    stop(paste0(unit_label, "-tissue identity is duplicated"), call. = FALSE)
  }
  unit_table <- table(unit)
  if (any(unit_table != 2L)) {
    stop(paste0(unit_label, " pairing must contain exactly one NAT and one Tumor"),
         call. = FALSE)
  }
  model_data <- data.frame(
    tissue = factor(tissue, levels = c("NAT", "Tumor")),
    unit = factor(unit, levels = sort(unique(unit), method = "radix"))
  )
  design <- stats::model.matrix(~ tissue + unit, data = model_data)
  if (qr(design)$rank != ncol(design) || nrow(design) - ncol(design) <= 0L ||
      !"tissueTumor" %in% colnames(design)) {
    stop("paired model design is rank deficient or has no positive residual df",
         call. = FALSE)
  }
  list(
    values = values, design = design, coefficient = match("tissueTumor", colnames(design)),
    biological_n = as.integer(length(unit_table)), unit_label = unit_label
  )
}

s01_extract_limma_contract <- function(fit, prepared, formula, normalization) {
  coefficient <- prepared$coefficient
  effect <- as.numeric(fit$coefficients[, coefficient])
  standard_error <- as.numeric(
    fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  )
  statistic <- as.numeric(fit$t[, coefficient])
  p_value <- as.numeric(fit$p.value[, coefficient])
  df_total <- rep_len(as.numeric(fit$df.total), length(effect))
  critical <- stats::qt(0.975, df = df_total)
  ci_lower <- effect - critical * standard_error
  ci_upper <- effect + critical * standard_error
  estimable <- is.finite(effect) & is.finite(standard_error) &
    standard_error > 0 & is.finite(statistic) & is.finite(p_value) &
    p_value >= 0 & p_value <= 1 & is.finite(df_total) & df_total > 0 &
    is.finite(ci_lower) & is.finite(ci_upper)
  list(
    valid = all(estimable), state = ifelse(estimable, "ESTIMABLE", "NA_NONESTIMABLE"),
    feature_id = rownames(prepared$values), biological_unit = prepared$unit_label,
    biological_n = prepared$biological_n, formula = formula,
    coefficient_name = "tissueTumor", effect_direction = "Tumor_minus_NAT",
    normalization = normalization, effect = ifelse(estimable, effect, NA_real_),
    standard_error = ifelse(estimable, standard_error, NA_real_),
    statistic = ifelse(estimable, statistic, NA_real_),
    degrees_of_freedom = ifelse(estimable, df_total, NA_real_),
    p_value = ifelse(estimable, p_value, NA_real_),
    ci_lower = ifelse(estimable, ci_lower, NA_real_),
    ci_upper = ifelse(estimable, ci_upper, NA_real_)
  )
}

s01_fit_limma_paired <- function(expression, sample_data) {
  prepared <- s01_prepare_paired_model(
    expression, sample_data, "patient_id", "patient", require_counts = FALSE
  )
  fit <- limma::eBayes(limma::lmFit(prepared$values, prepared$design))
  s01_extract_limma_contract(
    fit, prepared, "expression~tissue+patient", "upstream_frozen_scale_no_voom"
  )
}

s01_fit_voom_limma_paired <- function(counts, sample_data) {
  prepared <- s01_prepare_paired_model(
    counts, sample_data, "patient_id", "patient", require_counts = TRUE
  )
  dge <- edgeR::DGEList(counts = prepared$values)
  dge <- edgeR::normLibSizes(dge, method = "TMM")
  transformed <- limma::voom(dge, prepared$design, plot = FALSE,
                             normalize.method = "none")
  fit <- limma::eBayes(limma::lmFit(transformed, prepared$design))
  s01_extract_limma_contract(
    fit, prepared, "expression~tissue+patient", "TMM_then_voom"
  )
}

s01_sc_member_rows <- function(compartment, config) {
  families <- config$multiplicity$families[c("SC_G_MEMBERS", "SC_H_MEMBERS")]
  rows <- do.call(rbind, lapply(names(families), function(family) {
    ids <- families[[family]]$hypotheses
    data.frame(family_hypothesis = ids, multiplicity_family = family,
      frozen_family_size = as.integer(families[[family]]$frozen_size),
      compartment = sub("__.*$", "", ids), member = sub("^.*__", "", ids),
      stringsAsFactors = FALSE)
  }))
  rownames(rows) <- NULL
  rows[rows$compartment == compartment, , drop = FALSE]
}

s01_adjust_sc_member_families <- function(results, config) {
  s01_validate_sc_member_output(results, "unadjusted", "all", config)
  expected <- do.call(rbind, lapply(c("myeloid", "fibroblast", "B_plasma"),
    s01_sc_member_rows, config = config))
  s01_assert_exact_order(results$family_hypothesis, expected$family_hypothesis,
    "single-cell member reporting")
  s01_require(identical(results$multiplicity_family, expected$multiplicity_family) &&
    identical(results$frozen_family_size, expected$frozen_family_size),
    "single-cell member family identities changed")
  for (family in c("SC_G_MEMBERS", "SC_H_MEMBERS")) {
    rows <- which(results$multiplicity_family == family)
    spec <- config$multiplicity$families[[family]]
    idx <- match(spec$hypotheses, results$family_hypothesis[rows])
    values <- setNames(results$p_value[rows][idx], spec$hypotheses)
    adjusted <- s01_frozen_bh_values(values, as.integer(spec$frozen_size), spec$hypotheses)
    results$q_value[rows[idx]] <- adjusted
  }
  s01_validate_sc_member_output(results, "adjusted", "all", config)
  results
}

s01_synthetic_edger_fixture <- function() {
  donor <- rep(sprintf("D%02d", 1:8), each = 2L)
  tissue <- rep(c("NAT", "Tumor"), 8L)
  samples <- paste(donor, tissue, sep = "_")
  counts <- outer(seq_len(500), seq_len(16), function(g, s) {
    as.integer(30 + (g %% 97) + (s %% 11) + ((g * s) %% 67) +
      ifelse(g <= 8 & s %% 2 == 0, 40 + s, 0))
  })
  genes <- c("ACSL4", "LPCAT3", "ALOX5", "ALOX15B", "POR", "ALDH2", "ALDH3A2", "ALDH7A1", sprintf("BG%04d", 9:500))
  dimnames(counts) <- list(genes, samples)
  list(counts = counts, sample_data = data.frame(sample_id = samples, donor_id = donor,
    tissue = tissue, retained_cells = 30L), gene_mapping = data.frame(feature_id = genes, gene_id = genes))
}

s01_fit_site_interaction <- function(delta, patient_id, site,
                                     reference_level, comparison_level) {
  if (!is.numeric(delta) || length(delta) != length(patient_id) ||
      length(delta) != length(site)) {
    stop("site-interaction vectors must have equal valid lengths", call. = FALSE)
  }
  patient_id <- s01_validate_identity_vector(patient_id, "site-interaction patient")
  site <- as.character(site)
  if (anyNA(site)) stop("site identity contains missing/NA values", call. = FALSE)
  if (any(!nzchar(trimws(site)))) stop("site identity contains blank values", call. = FALSE)
  if (length(reference_level) != 1L || length(comparison_level) != 1L ||
      is.na(reference_level) || is.na(comparison_level) ||
      !nzchar(reference_level) || !nzchar(comparison_level) ||
      identical(reference_level, comparison_level) ||
      !setequal(unique(site), c(reference_level, comparison_level))) {
    stop("site levels and explicit contrast direction are invalid", call. = FALSE)
  }
  if (any(!is.finite(delta))) {
    return(s01_nonestimable(list(
      valid = FALSE, reason = "nonfinite_patient_level_program_delta",
      effect = NA_real_, statistic = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  model_data <- data.frame(
    delta_program_score = as.numeric(delta),
    site = factor(site, levels = c(reference_level, comparison_level))
  )
  design <- stats::model.matrix(~ site, model_data)
  if (qr(design)$rank != ncol(design) || nrow(design) - ncol(design) <= 0L) {
    return(s01_nonestimable(list(
      valid = FALSE, reason = "rank_deficient_or_nonpositive_residual_df",
      effect = NA_real_, statistic = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  fit <- stats::lm(delta_program_score ~ site, data = model_data)
  summary_fit <- summary(fit)
  coefficient_name <- colnames(design)[[2L]]
  coefficient <- summary_fit$coefficients[coefficient_name, , drop = TRUE]
  interval <- stats::confint(fit, parm = coefficient_name, level = 0.95)
  valid <- all(is.finite(c(coefficient, interval))) && fit$df.residual > 0L
  if (!valid) {
    return(s01_nonestimable(list(
      valid = FALSE, reason = "nonfinite_site_model_result",
      effect = NA_real_, statistic = NA_real_, p_value = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_
    )))
  }
  list(
    valid = TRUE, estimable = TRUE, state = "ESTIMABLE",
    biological_unit = "patient", biological_n = as.integer(length(delta)),
    formula = "delta_program_score~site", coefficient_name = coefficient_name,
    effect_direction = paste0(comparison_level, "_minus_", reference_level),
    effect = unname(coefficient[["Estimate"]]),
    standard_error = unname(coefficient[["Std. Error"]]),
    statistic = unname(coefficient[["t value"]]),
    p_value = unname(coefficient[["Pr(>|t|)"]]),
    ci_lower = unname(interval[[1L]]), ci_upper = unname(interval[[2L]]),
    degrees_of_freedom = as.integer(fit$df.residual)
  )
}

s01_validate_sc_member_output <- function(rows, stage = "adjusted", population = "all",
    config) {
  s01_require(length(stage) == 1L && stage %in% c("unadjusted", "adjusted") &&
    length(population) == 1L && population %in% c("all", "myeloid", "fibroblast", "B_plasma") &&
    (stage == "unadjusted" || population == "all"), "invalid member output stage or population")
  compartments <- if (population == "all") c("myeloid", "fibroblast", "B_plasma") else population
  expected <- do.call(rbind, lapply(compartments, s01_sc_member_rows, config = config))
  character_fields <- c("family_hypothesis", "multiplicity_family", "compartment", "member",
    "estimability", "missing_reason", "standard_error_status", "confidence_interval_status",
    "uncertainty_reason", "effect_direction", "biological_unit", "ql_df_status")
  required <- c(character_fields, "frozen_family_size", "biological_n", "effect", "statistic",
    "numerator_df", "denominator_df", "p_value", "q_value", "standard_error", "ci_lower", "ci_upper")
  s01_require(is.data.frame(rows) && !anyDuplicated(names(rows)) &&
    all(required %in% names(rows)) && nrow(rows) == nrow(expected) && nrow(rows) > 0L,
    "member output required columns or population cardinality mismatch")
  s01_require(all(vapply(rows[character_fields], function(v) is.character(v) &&
    length(v) == nrow(rows) && !anyNA(v) && all(nzchar(trimws(v))), logical(1))),
    "member identity/status type or missingness mismatch")
  for (field in names(expected)) s01_require(identical(rows[[field]], expected[[field]]),
    "single-cell member reporting identity or population mismatch")
  s01_require(is.integer(rows$biological_n) && length(rows$biological_n) == nrow(rows) && !anyNA(rows$biological_n) &&
    all(rows$biological_n >= 0L) && all(rows$biological_unit == "donor") &&
    all(rows$effect_direction == "Tumor_minus_NAT"), "invalid member biological unit or n")
  numeric <- c("effect", "statistic", "numerator_df", "denominator_df",
    "p_value", "q_value", "standard_error", "ci_lower", "ci_upper")
  s01_require(is.data.frame(rows) && all(numeric %in% names(rows)) &&
    all(vapply(rows[numeric], is.numeric, logical(1))), "member numeric output schema mismatch")
  s01_require(is.double(rows$numerator_df) && is.double(rows$denominator_df),
    "native member degrees of freedom require double precision")
  s01_require(all(rows$estimability %in% c("ESTIMABLE", "NA_NONESTIMABLE")),
    "member estimability vocabulary mismatch")
  estimable <- rows$estimability == "ESTIMABLE"
  s01_require(identical(rows$ql_df_status, ifelse(estimable, "AVAILABLE", "NA_NONESTIMABLE")),
    "native QL status mismatch")
  s01_require(all(rows$biological_n[estimable] >= config$single_cell_member_model$minimum_complete_pairs) &&
    all(vapply(split(rows$biological_n, rows$compartment), function(n) length(unique(n)) == 1L, logical(1))),
    "invalid member biological n")
  s01_require(all(vapply(rows[numeric], function(v) length(v) == nrow(rows) &&
    !any(is.nan(v) | is.infinite(v)), logical(1))), "invalid member numeric missingness")
  populated_q <- !is.na(rows$q_value)
  s01_require(all(rows$q_value[populated_q] >= 0 & rows$q_value[populated_q] <= 1) &&
    (stage == "unadjusted" || all(is.finite(rows$q_value[estimable]))), "invalid member q value or adjustment stage")
  if (stage == "adjusted") for (family in c("SC_G_MEMBERS", "SC_H_MEMBERS")) {
    spec <- config$multiplicity$families[[family]]
    i <- which(rows$multiplicity_family == family)
    bh <- s01_frozen_bh_values(setNames(rows$p_value[i], rows$family_hypothesis[i]),
      as.integer(spec$frozen_size), spec$hypotheses)
    s01_require(identical(unname(rows$q_value[i]), unname(bh)), "member q values disagree with frozen family BH")
  }
  for (field in c("effect", "statistic", "numerator_df", "denominator_df", "p_value"))
    s01_require(all(is.finite(rows[[field]][estimable])), "invalid estimable native member output")
  s01_require(all(rows$statistic[estimable] >= 0) &&
    all(rows$numerator_df[estimable] > 0 & rows$denominator_df[estimable] > 0) &&
    all(rows$p_value[estimable] >= 0 & rows$p_value[estimable] <= 1),
    "invalid native member statistic probability or degrees of freedom")
  s01_require(all(is.na(rows$p_value[!estimable])) && all(is.na(rows$q_value[!estimable])),
    "non-estimable member cannot have a P or q value")
  s01_require(all(is.na(rows$standard_error)) && all(is.na(rows$ci_lower)) &&
    all(is.na(rows$ci_upper)) &&
    all(rows$standard_error_status == "NA_UNAVAILABLE") &&
    all(rows$confidence_interval_status == "NA_UNAVAILABLE") &&
    all(rows$uncertainty_reason == "native_modern_edgeR_does_not_provide_coefficient_SE_or_95_percent_CI"),
    "native member uncertainty must be explicitly unavailable without fabricated estimates")
  s01_require(all(nzchar(rows$missing_reason)) &&
    all(rows$missing_reason[estimable] == "NA_NOT_APPLICABLE") &&
    all(rows$missing_reason[!estimable] != "NA_NOT_APPLICABLE"),
    "non-estimable member must retain its reason")
  invisible(TRUE)
}
