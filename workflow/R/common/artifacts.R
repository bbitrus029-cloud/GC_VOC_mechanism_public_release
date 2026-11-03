# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Public adaptation: explicit inputs and outputs; original scientific rules retained.

s01_schema_field_contracts <- function(required_fields) {
  integer_fields <- c(
    "biological_n", "frozen_family_size", "donor_count",
    "valid_bootstrap_replicates", "model_rank", "model_columns", "residual_df"
  )
  double_fields <- c(
    "effect", "ci_lower", "ci_upper", "p_value", "q_value", "g_av",
    "g_av_ci_lower", "g_av_ci_upper", "bootstrap_ci_lower",
    "bootstrap_ci_upper", "bca_ci_lower", "bca_ci_upper", "rho",
    "primary_rank_unrounded_raw_paired_mean_difference",
    "standardized_member_mean_unrounded_raw_paired_mean_difference",
    "coefficient", "standard_error", "dispersion", "statistic", "numerator_df", "denominator_df"
  )
  logical_fields <- c(
    "both_core_CI_bounded_within_small_effect_region", "precision_uncertain",
    "score_method_nonestimable", "exact_zero_direction_flag",
    "zero_point_estimate_flag", "bridge_q_ci_inconsistent",
    "bridge_nominal_p_ci_inconsistent", "compensation_compatible",
    "same_cohort_dependence_flag", "recognized_warning_present",
    "unexpected_warning_present", "converged", "boundary",
    "fitted_probability_boundary"
  )
  vector_fields <- "warning_messages"
  types <- setNames(rep("character_nonempty", length(required_fields)), required_fields)
  types[integer_fields] <- "integer"
  types[double_fields] <- "double"
  types[logical_fields] <- "logical"
  types[vector_fields] <- "character_vector"
  s01_require(!anyNA(types), "schema field type assignment is incomplete")

  domains <- setNames(rep("nonempty_UTF8", length(required_fields)), required_fields)
  domains[integer_fields] <- "integer_nonnegative"
  domains["residual_df"] <- "integer_positive_when_estimable"
  domains[double_fields] <- "finite_numeric_when_estimable"
  domains[c("p_value", "q_value")] <- "closed_interval_0_1_when_estimable"
  domains["rho"] <- "closed_interval_minus1_plus1_when_estimable"
  domains[logical_fields] <- "TRUE_or_FALSE"
  domains["warning_messages"] <- "exact_conditionMessage_strings_in_capture_order"
  domains["artifact_sha256"] <- "lowercase_SHA256_exactly_64_hex"
  domains["artifact_path"] <- "nonempty_project_relative_path_no_absolute_or_traversal"
  domains["dimension_a_state"] <- "registered_dimension_A_vocabulary"
  domains["dimension_b_state"] <- "registered_dimension_B_vocabulary"
  domains["bridge_state"] <- "registered_bridge_state_vocabulary"
  domains["cross_layer_state"] <- "registered_cross_layer_state_vocabulary"
  domains["estimability"] <- "registered_estimability_vocabulary"
  domains["marker_audit_state"] <- "PASS_FAIL_or_NA_NONESTIMABLE"
  domains["reconciliation_state"] <- "VERIFIED_or_NA_UNAVAILABLE"
  domains["timing"] <- "registered_timing_vocabulary"
  domains["multiplicity_status"] <- "registered_multiplicity_vocabulary"
  domains[c("primary_rank_direction", "standardized_member_mean_direction",
            "native_direction")] <- "positive_negative_zero_or_NA_NOT_DIRECTIONAL"
  domains["score_method_sensitive"] <- "TRUE_FALSE_or_NA_NONESTIMABLE"
  domains["diagnostic_state"] <- "registered_diagnostic_terminal_state"
  domains["terminal_state"] <- "registered_terminal_state_vocabulary"
  domains["missing_representation"] <- "registered_missing_representation_vocabulary"

  result_dependent <- unique(c(
    double_fields, logical_fields, "member", "missing_reason", "native_direction",
    "primary_rank_direction", "standardized_member_mean_direction",
    "score_method_sensitive", "warning_messages", "diagnostic_state",
    "terminal_state", "missing_representation", "reconciliation_state"
  ))
  nullability <- setNames(rep("never_null_or_blank", length(required_fields)),
                          required_fields)
  nullability[intersect(result_dependent, required_fields)] <-
    "state_dependent_only_as_declared_by_estimability_or_activation_state"
  missing <- setNames(rep("missing_prohibited", length(required_fields)),
                      required_fields)
  missing[intersect(result_dependent, required_fields)] <-
    "use_only_registered_explicit_missing_state_never_silent_drop"
  result <- lapply(required_fields, function(field) list(
    type = unname(types[[field]]), domain = unname(domains[[field]]),
    nullability = unname(nullability[[field]]),
    missing_state_behavior = unname(missing[[field]])
  ))
  names(result) <- required_fields
  for (field in intersect(c("standard_error", "ci_lower", "ci_upper"), required_fields)) {
    result[[field]]$domain <- "finite_numeric_when_estimable_except_14_SC_members_native_uncertainty_unavailable"
    result[[field]]$nullability <- "NA_UNAVAILABLE_for_SC_G_MEMBERS_and_SC_H_MEMBERS_native_uncertainty_only_otherwise_estimability_dependent"
    result[[field]]$missing_state_behavior <- "numeric_NA_with_explicit_uncertainty_status_and_reason_for_14_SC_members_no_fabricated_interval"
  }
  for (field in intersect(c("numerator_df", "denominator_df"), required_fields)) {
    result[[field]]$domain <- "positive_finite_native_QL_double_when_estimable_SC_member_only"
    result[[field]]$nullability <- "numeric_NA_when_QL_not_applicable_or_nonestimable"
    result[[field]]$missing_state_behavior <- "ql_df_status_NA_NOT_APPLICABLE_for_other_methods_NA_NONESTIMABLE_for_failed_QL"
  }
  if ("ql_df_status" %in% required_fields) result$ql_df_status$domain <- "AVAILABLE_NA_NONESTIMABLE_or_NA_NOT_APPLICABLE"
  if ("member" %in% required_fields) result$member$domain <- "registered_member_identity_or_NA_NOT_APPLICABLE_for_program_row"
  for (field in intersect(c("standard_error_status", "confidence_interval_status"), required_fields))
    result[[field]]$domain <- "AVAILABLE_or_NA_NONESTIMABLE_with_NA_UNAVAILABLE_exclusively_for_14_SC_members"
  result
}

s01_output_schema_contract <- function() {
  common <- c(
    "endpoint_id", "hypothesis_id", "source", "biological_unit",
    "biological_n", "estimand", "effect", "ci_lower", "ci_upper",
    "p_value", "q_value", "multiplicity_status", "multiplicity_family",
    "family_hypothesis", "frozen_family_size", "role", "timing",
    "dependence", "activation_condition", "estimability", "missing_reason",
    "claim_ceiling", "provenance_id"
  )
  output <- function(output_name, key, extra,
                     order = c("endpoint_id", "source")) {
    list(
      primary_key = key,
      required_columns = unique(c(common, extra)),
      deterministic_order = order,
      row_population = paste0(
        "all_registered_crosswalk_rows_mapped_to_", output_name,
        "_in_registered_order_including_inactive_and_nonestimable"
      ),
      favorable_row_filter_permitted = FALSE,
      reconciliation = paste0(
        "one_terminal_row_per_crosswalk_endpoint_mapped_to_", output_name
      )
    )
  }
  aggregate_output <- function(key, required, population, reconciliation) {
    list(
      primary_key = key,
      required_columns = required,
      deterministic_order = key,
      row_population = population,
      favorable_row_filter_permitted = FALSE,
      reconciliation = reconciliation
    )
  }
  diagnostic <- function(key, required, population, reconciliation,
                         order = key) {
    list(
      primary_key = key, required_columns = required,
      deterministic_order = order, row_population = population,
      favorable_row_filter_permitted = FALSE, reconciliation = reconciliation
    )
  }
  schema <- list(
    schema_version = "s01-output-schema-1.3.1",
    authority = c(
      "docs/plans/APPROVED_ANALYSIS_PLAN.md#21.3",
      "config/analysis_config.yml",
      "config/hypothesis_estimand_crosswalk.tsv"
    ),
    common_estimate_columns = common,
    missing_representations = c(
      "NA_NONESTIMABLE", "NA_UNAVAILABLE", "NA_NOT_APPLICABLE",
      "NA_NOT_DIRECTIONAL", "NO_FORMAL_FAMILY"
    ),
    terminal_state_values = list(
      dimension_a = c(
        "non_estimable", "uncertain", "replicated_program_change",
        "directionally_consistent_partial_support",
        "concordant_but_inconclusive",
        "opposite_direction_statistical_conflict",
        "directionally_discordant_inconclusive"
      ),
      dimension_b = c(
        "magnitude_nonestimable_g_av", "substantive_in_both",
        "substantive_in_one", "small_point_estimates_in_both"
      ),
      bridge = c("bridge_nonestimable", "bridge_supported", "bridge_uncertain"),
      cross_layer = c(
        "non_estimable", "orthogonally_supported", "directionally_convergent",
        "cross_layer_conflict", "indeterminate"
      ),
      generic_estimability = c(
        "ESTIMABLE", "NA_NONESTIMABLE", "NA_UNAVAILABLE", "NA_NOT_APPLICABLE"
      )
    ),
    multiplicity = list(
      family_population_source = "config.analysis_config.yml#multiplicity.families",
      endpoint_population_source = "config.hypothesis_estimand_crosswalk.tsv",
      q_rule = "BH_on_finite_P_with_complete_vector_equal_frozen_family_size_then_restore_structural_NA_in_frozen_order",
      inactive_GSE122401 = "no_effect_or_P_value_unless_pre_effect_backup_activation_occurs",
      nonestimable_no_test = "nominal_P_NA_and_q_NO_FORMAL_FAMILY",
      exploratory_nominal = "nominal_P_if_estimable_and_q_NO_FORMAL_FAMILY"
    ),
    reporting = list(
      retain_every_registered_hypothesis = TRUE,
      significant_results_only_permitted = FALSE,
      retain_conditional_nonestimable_uncertain_conflicting_unfavorable = TRUE
    ),
    registered_method_contract_resolutions = list(
      edgeR_QL_member_uncertainty = list(
        state = "APPROVED_METHOD_IMPLEMENTED",
        registered_method = "edgeR_QL_tissue_plus_donor",
        approved_reporting_exception = "native_SE_and_CI_NA_UNAVAILABLE_does_not_invalidate_native_QL_test",
        decision_dossier = s01_registered_method_resolution(),
        native_field_mapping = list(logFC = "effect", F = "statistic",
          df.test = "numerator_df", df.total = "denominator_df"),
        invented_rule_permitted = FALSE
      )
    ),
    outputs = list(
      all_program_estimates = output(
        "all_program_estimates",
        c("endpoint_id"),
        c("construct", "native_direction", "g_av", "g_av_ci_lower",
          "g_av_ci_upper", "frozen_family_size")
      ),
      all_member_estimates = output(
        "all_member_estimates",
        c("endpoint_id"), c("construct", "member", "frozen_family_size",
          "standard_error", "statistic", "numerator_df", "denominator_df",
          "standard_error_status", "confidence_interval_status", "uncertainty_reason", "ql_df_status")
      ),
      bulk_reproducibility_states = aggregate_output(
        "construct",
        c(
          "construct", "active_source_1_endpoint_id", "active_source_2_endpoint_id",
          "dimension_a_state", "zero_point_estimate_flag", "dimension_b_state",
          "both_core_CI_bounded_within_small_effect_region", "precision_uncertain",
          "primary_rank_unrounded_raw_paired_mean_difference",
          "standardized_member_mean_unrounded_raw_paired_mean_difference",
          "primary_rank_direction", "standardized_member_mean_direction",
          "score_method_sensitive", "score_method_nonestimable",
          "exact_zero_direction_flag", "claim_ceiling", "provenance_id"
        ),
        "one_row_per_axis_for_active_core_source_pair",
        "exactly_one_state_row_for_each_of_G_and_H"
      ),
      single_cell_state_estimates = output(
        "single_cell_state_estimates",
        c("endpoint_id"), c("broad_label", "donor_count", "marker_audit_state",
          "member", "standard_error", "statistic", "numerator_df", "denominator_df",
          "standard_error_status", "confidence_interval_status", "uncertainty_reason", "ql_df_status")
      ),
      captured_composition_estimates = output(
        "captured_composition_estimates",
        c("endpoint_id"), c("compartment", "composition_method", "donor_count")
      ),
      protein_estimates = output(
        "protein_estimates",
        c("endpoint_id"), c("construct", "member", "common_pair_set_id")
      ),
      cptac_bridge_estimates = output(
        "cptac_bridge_estimates",
        c("endpoint_id"),
        c("construct", "rho", "bootstrap_ci_lower", "bootstrap_ci_upper",
          "bca_ci_lower", "bca_ci_upper", "valid_bootstrap_replicates",
          "bridge_state", "bridge_q_ci_inconsistent",
          "bridge_nominal_p_ci_inconsistent")
      ),
      cross_layer_evidence_matrix = aggregate_output(
        c("construct", "layer"),
        c("construct", "layer", "evidence_identifier", "cross_layer_state",
          "compensation_compatible", "same_cohort_dependence_flag",
          "principal_limitation", "claim_ceiling", "provenance_id"),
        "one_row_per_axis_and_evidence_layer",
        "complete_dependence_aware_layer_inventory_for_each_axis"
      ),
      claim_boundary_ledger = aggregate_output(
        "claim_id",
        c("claim_id", "evidence_identifier", "permitted_claim",
          "prohibited_inference", "uncertainty_statement", "provenance_id"),
        "one_row_per_registered_or_derived_claim_boundary",
        "each_evidence_identifier_has_an_explicit_claim_ceiling"
      ),
      figure_source_data_index = aggregate_output(
        c("figure_id", "source_data_id"),
        c("figure_id", "source_data_id", "artifact_path", "artifact_sha256",
          "reconciliation_state"),
        "one_row_per_figure_source_data_artifact",
        "every_figure_panel_reconciles_to_one_verified_source_data_artifact"
      )
    ),
    indispensable_diagnostics = list(
      exclusions_and_deviations = diagnostic(
        c("source", "record_id"),
        c("source", "record_id", "rule_id", "terminal_state", "provenance_id"),
        "one_row_per_prespecified_exclusion_or_deviation_record",
        "complete_admitted_source_exclusion_and_deviation_population"
      ),
      robustness_diagnostics = diagnostic(
        c("endpoint_id", "diagnostic_id"),
        c("endpoint_id", "diagnostic_id", "diagnostic_state",
          "multiplicity_status", "provenance_id"),
        "one_row_per_registered_endpoint_and_applicable_diagnostic",
        "complete_prespecified_diagnostic_population_with_terminal_states"
      ),
      result_state_ledger = diagnostic(
        "endpoint_id",
        c("endpoint_id", "activation_state", "estimability", "terminal_state",
          "missing_representation", "provenance_id"),
        "exactly_one_row_per_registered_crosswalk_endpoint",
        "exactly_146_ordered_endpoint_terminal_state_rows"
      ),
      quasibinomial_model_audit = diagnostic(
        "endpoint_id",
        c(
          "endpoint_id", "warning_messages", "recognized_warning_present",
          "unexpected_warning_present", "converged", "boundary", "model_rank",
          "model_columns", "fitted_probability_boundary", "coefficient",
          "standard_error", "dispersion", "residual_df", "diagnostic_state",
          "provenance_id"
        ),
        "one_row_per_SC_COMPOSITION_QUASIBINOMIAL_endpoint",
        "exactly_seven_frozen_family_rows_with_exact_warning_and_model_facts"
      )
    )
  )
  required <- unique(unlist(lapply(
    c(schema$outputs, schema$indispensable_diagnostics), `[[`, "required_columns"
  )))
  schema$method_field_applicability <- list(
    authority = "config/hypothesis_estimand_crosswalk.tsv",
    QL_scope = c("SC_G_MEMBERS", "SC_H_MEMBERS"),
    native_QL = "native_F_df.test_df.total_required_when_estimable_no_rounding",
    other_methods = "retain_registered_native_statistic_SE_and_CI_numeric_QL_df_NA_with_ql_df_status_NA_NOT_APPLICABLE",
    uncertainty_exception = "only_14_SC_members_numeric_SE_CI_NA_with_NA_UNAVAILABLE_and_native_unavailability_reason",
    program_member = "member_NA_NOT_APPLICABLE_for_nonmember_program_rows",
    row_validator = "s01_validate_method_output_rows")
  schema$field_contracts <- s01_schema_field_contracts(required)
  schema$column_types <- lapply(schema$field_contracts, `[[`, "type")
  schema
}

s01_validate_output_schemas <- function(schema, rows = NULL) {
  if (!is.null(rows)) s01_validate_method_output_rows(rows)
  expected_outputs <- c(
    "all_program_estimates", "all_member_estimates",
    "bulk_reproducibility_states", "single_cell_state_estimates",
    "captured_composition_estimates", "protein_estimates",
    "cptac_bridge_estimates", "cross_layer_evidence_matrix",
    "claim_boundary_ledger", "figure_source_data_index"
  )
  s01_require(is.list(schema) && identical(names(schema$outputs), expected_outputs),
              "complete output-schema population or deterministic order mismatch")
  s01_require(identical(schema$schema_version, "s01-output-schema-1.3.1"),
              "output schema version mismatch")
  required_common <- c(
    "biological_unit", "biological_n", "estimand", "effect", "ci_lower",
    "ci_upper", "p_value", "q_value", "multiplicity_status", "role", "timing",
    "dependence", "estimability", "claim_ceiling"
  )
  s01_require(all(required_common %in% schema$common_estimate_columns),
              "common estimate schema is incomplete")
  s01_require(isTRUE(schema$reporting$retain_every_registered_hypothesis) &&
                identical(schema$reporting$significant_results_only_permitted, FALSE),
              "schema permits favorable-result filtering")
  conflict <- schema$registered_method_contract_resolutions$edgeR_QL_member_uncertainty
  s01_require(
    identical(conflict$state, "APPROVED_METHOD_IMPLEMENTED") &&
      identical(conflict$registered_method, "edgeR_QL_tissue_plus_donor") &&
      identical(conflict$approved_reporting_exception,
                "native_SE_and_CI_NA_UNAVAILABLE_does_not_invalidate_native_QL_test") &&
      identical(conflict$decision_dossier, s01_registered_method_resolution()) &&
      identical(conflict$native_field_mapping, list(logFC = "effect", F = "statistic",
        df.test = "numerator_df", df.total = "denominator_df")) &&
      identical(conflict$invented_rule_permitted, FALSE),
    "registered edgeR QL uncertainty contract conflict is not preserved"
  )
  s01_require(all(vapply(schema$outputs, function(x) {
    isFALSE(x$favorable_row_filter_permitted) &&
      length(x$row_population) == 1L && nzchar(x$row_population) &&
      length(x$primary_key) > 0L && all(x$primary_key %in% x$required_columns) &&
      length(x$required_columns) > 0L
  }, logical(1))), "one or more output schemas is not complete or fail-closed")
  expected_output_keys <- list(
    all_program_estimates = "endpoint_id",
    all_member_estimates = "endpoint_id",
    bulk_reproducibility_states = "construct",
    single_cell_state_estimates = "endpoint_id",
    captured_composition_estimates = "endpoint_id",
    protein_estimates = "endpoint_id",
    cptac_bridge_estimates = "endpoint_id",
    cross_layer_evidence_matrix = c("construct", "layer"),
    claim_boundary_ledger = "claim_id",
    figure_source_data_index = c("figure_id", "source_data_id")
  )
  expected_populations <- c(
    all_program_estimates = "all_registered_crosswalk_rows_mapped_to_all_program_estimates_in_registered_order_including_inactive_and_nonestimable",
    all_member_estimates = "all_registered_crosswalk_rows_mapped_to_all_member_estimates_in_registered_order_including_inactive_and_nonestimable",
    bulk_reproducibility_states = "one_row_per_axis_for_active_core_source_pair",
    single_cell_state_estimates = "all_registered_crosswalk_rows_mapped_to_single_cell_state_estimates_in_registered_order_including_inactive_and_nonestimable",
    captured_composition_estimates = "all_registered_crosswalk_rows_mapped_to_captured_composition_estimates_in_registered_order_including_inactive_and_nonestimable",
    protein_estimates = "all_registered_crosswalk_rows_mapped_to_protein_estimates_in_registered_order_including_inactive_and_nonestimable",
    cptac_bridge_estimates = "all_registered_crosswalk_rows_mapped_to_cptac_bridge_estimates_in_registered_order_including_inactive_and_nonestimable",
    cross_layer_evidence_matrix = "one_row_per_axis_and_evidence_layer",
    claim_boundary_ledger = "one_row_per_registered_or_derived_claim_boundary",
    figure_source_data_index = "one_row_per_figure_source_data_artifact"
  )
  for (name in expected_outputs) {
    s01_require(
      identical(schema$outputs[[name]]$primary_key, expected_output_keys[[name]]) &&
        identical(schema$outputs[[name]]$row_population,
                  unname(expected_populations[[name]])),
      paste0("output-specific natural key or row population mismatch: ", name)
    )
  }
  expected_diagnostics <- c(
    "exclusions_and_deviations", "robustness_diagnostics",
    "result_state_ledger", "quasibinomial_model_audit"
  )
  s01_require(identical(names(schema$indispensable_diagnostics), expected_diagnostics),
              "indispensable diagnostic schema population mismatch")
  s01_require(all(vapply(schema$indispensable_diagnostics, function(x) {
    isFALSE(x$favorable_row_filter_permitted) && length(x$primary_key) > 0L &&
      all(x$primary_key %in% x$required_columns) &&
      length(x$deterministic_order) > 0L &&
      all(x$deterministic_order %in% x$required_columns) &&
      length(x$row_population) == 1L && nzchar(x$row_population) &&
      length(x$reconciliation) == 1L && nzchar(x$reconciliation)
  }, logical(1))), "indispensable diagnostic natural key or population is invalid")
  expected_diagnostic_keys <- list(
    exclusions_and_deviations = c("source", "record_id"),
    robustness_diagnostics = c("endpoint_id", "diagnostic_id"),
    result_state_ledger = "endpoint_id",
    quasibinomial_model_audit = "endpoint_id"
  )
  for (name in expected_diagnostics) {
    s01_require(identical(
      schema$indispensable_diagnostics[[name]]$primary_key,
      expected_diagnostic_keys[[name]]
    ), paste0("diagnostic-specific natural key mismatch: ", name))
  }
  endpoint_outputs <- c(
    "all_program_estimates", "all_member_estimates", "single_cell_state_estimates",
    "captured_composition_estimates", "protein_estimates", "cptac_bridge_estimates"
  )
  s01_require(all(vapply(schema$outputs[endpoint_outputs], function(x) {
    identical(x$primary_key, "endpoint_id") &&
      all(schema$common_estimate_columns %in% x$required_columns)
  }, logical(1))), "endpoint output schema population or common fields mismatch")
  non_endpoint_outputs <- c(
    "bulk_reproducibility_states", "cross_layer_evidence_matrix",
    "claim_boundary_ledger", "figure_source_data_index"
  )
  s01_require(all(vapply(schema$outputs[non_endpoint_outputs], function(x) {
    !"endpoint_id" %in% x$required_columns
  }, logical(1))), "aggregate or index schema manufactures endpoint rows")
  score_fields <- c(
    "primary_rank_unrounded_raw_paired_mean_difference",
    "standardized_member_mean_unrounded_raw_paired_mean_difference",
    "primary_rank_direction", "standardized_member_mean_direction",
    "score_method_sensitive", "score_method_nonestimable",
    "exact_zero_direction_flag"
  )
  bulk_columns <- schema$outputs$bulk_reproducibility_states$required_columns
  s01_require("zero_point_estimate_flag" %in% bulk_columns &&
                all(score_fields %in% bulk_columns),
              "bulk state schema lacks zero or score-method sensitivity fields")
  s01_require(all(c("NA_NONESTIMABLE", "NA_UNAVAILABLE", "NA_NOT_APPLICABLE",
                    "NA_NOT_DIRECTIONAL", "NO_FORMAL_FAMILY") %in%
                    schema$missing_representations),
              "missing-value representation vocabulary is incomplete")
  s01_require(identical(length(schema$terminal_state_values$cross_layer), 5L),
              "cross-layer terminal states are not exhaustive")
  required_fields <- unique(unlist(lapply(
    c(schema$outputs, schema$indispensable_diagnostics), `[[`, "required_columns"
  )))
  s01_require(is.list(schema$field_contracts) &&
                identical(names(schema$field_contracts), required_fields),
              "field contract population has an omitted, extra, or misordered field")
  field_parts <- c("type", "domain", "nullability", "missing_state_behavior")
  s01_require(all(vapply(schema$field_contracts, function(contract) {
    is.list(contract) && identical(sort(names(contract)), sort(field_parts)) &&
      all(vapply(contract, function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) &&
          nzchar(trimws(value))
      }, logical(1)))
  }, logical(1))), "one or more field contracts lacks type, domain, nullability, or missing behavior")
  allowed_types <- c("character_nonempty", "character_vector", "integer",
                     "double", "logical")
  s01_require(all(vapply(schema$field_contracts, function(contract) {
    contract$type %in% allowed_types &&
      contract$nullability %in% c(
        "never_null_or_blank",
        "numeric_NA_when_QL_not_applicable_or_nonestimable",
        "NA_UNAVAILABLE_for_SC_G_MEMBERS_and_SC_H_MEMBERS_native_uncertainty_only_otherwise_estimability_dependent",
        "state_dependent_only_as_declared_by_estimability_or_activation_state"
      ) && contract$missing_state_behavior %in% c(
        "missing_prohibited",
        "ql_df_status_NA_NOT_APPLICABLE_for_other_methods_NA_NONESTIMABLE_for_failed_QL",
        "numeric_NA_with_explicit_uncertainty_status_and_reason_for_14_SC_members_no_fabricated_interval",
        "use_only_registered_explicit_missing_state_never_silent_drop"
      )
  }, logical(1))), "field contract type, nullability, or missing behavior vocabulary is invalid")
  for (field in c("statistic", "numerator_df", "denominator_df")) {
    s01_require(identical(schema$field_contracts[[field]]$type, "double"),
      "native QL statistics and degrees of freedom must preserve double precision")
  }
  for (field in c("standard_error", "ci_lower", "ci_upper")) {
    s01_require(identical(schema$field_contracts[[field]]$nullability,
      "NA_UNAVAILABLE_for_SC_G_MEMBERS_and_SC_H_MEMBERS_native_uncertainty_only_otherwise_estimability_dependent"),
      "native uncertainty reporting exception scope mismatch")
  }
  expected_scope <- s01_output_schema_contract()
  s01_require(identical(schema$method_field_applicability, expected_scope$method_field_applicability),
    "method field applicability scope mismatch")
  for (field in c("numerator_df", "denominator_df", "member", "ql_df_status",
      "standard_error_status", "confidence_interval_status", "standard_error", "ci_lower", "ci_upper"))
    s01_require(identical(schema$field_contracts[[field]], expected_scope$field_contracts[[field]]),
      "method specific field contract mismatch")
  probability_fields <- c("p_value", "q_value")
  integer_fields <- c("biological_n", "frozen_family_size", "donor_count",
                      "valid_bootstrap_replicates", "model_rank", "model_columns",
                      "residual_df")
  logical_fields <- c(
    "both_core_CI_bounded_within_small_effect_region", "precision_uncertain",
    "score_method_nonestimable", "exact_zero_direction_flag",
    "zero_point_estimate_flag", "bridge_q_ci_inconsistent",
    "bridge_nominal_p_ci_inconsistent", "compensation_compatible",
    "same_cohort_dependence_flag", "recognized_warning_present",
    "unexpected_warning_present", "converged", "boundary",
    "fitted_probability_boundary"
  )
  s01_require(all(vapply(probability_fields, function(field) {
    identical(schema$field_contracts[[field]]$type, "double") &&
      identical(schema$field_contracts[[field]]$domain,
                "closed_interval_0_1_when_estimable")
  }, logical(1))), "probability field type or range contract is invalid")
  s01_require(all(vapply(integer_fields, function(field) {
    identical(schema$field_contracts[[field]]$type, "integer")
  }, logical(1))), "integer count field type contract is invalid")
  s01_require(all(vapply(logical_fields, function(field) {
    identical(schema$field_contracts[[field]]$type, "logical") &&
      identical(schema$field_contracts[[field]]$domain, "TRUE_or_FALSE")
  }, logical(1))), "Boolean field type or vocabulary contract is invalid")
  exact_state_domains <- c(
    dimension_a_state = "registered_dimension_A_vocabulary",
    dimension_b_state = "registered_dimension_B_vocabulary",
    bridge_state = "registered_bridge_state_vocabulary",
    cross_layer_state = "registered_cross_layer_state_vocabulary",
    estimability = "registered_estimability_vocabulary",
    marker_audit_state = "PASS_FAIL_or_NA_NONESTIMABLE",
    reconciliation_state = "VERIFIED_or_NA_UNAVAILABLE"
  )
  s01_require(all(vapply(names(exact_state_domains), function(field) {
    identical(schema$field_contracts[[field]]$domain,
              unname(exact_state_domains[[field]]))
  }, logical(1))), "terminal-state field vocabulary contract is invalid")
  s01_require(identical(schema$field_contracts$artifact_sha256$domain,
                        "lowercase_SHA256_exactly_64_hex"),
              "artifact hash field domain is malformed")
  s01_require(identical(schema$field_contracts$artifact_path$domain,
                        "nonempty_project_relative_path_no_absolute_or_traversal"),
              "artifact path field domain is malformed")
  s01_require(identical(
    names(schema$column_types), required_fields
  ) && identical(unname(unlist(schema$column_types)),
                 unname(vapply(schema$field_contracts, `[[`, character(1), "type"))),
              "legacy column type projection is incomplete")
  qb_required <- c(
    "warning_messages", "converged", "boundary", "model_rank", "model_columns",
    "fitted_probability_boundary", "coefficient", "standard_error", "dispersion",
    "residual_df"
  )
  s01_require(all(qb_required %in%
                    schema$indispensable_diagnostics$quasibinomial_model_audit$required_columns),
              "registered quasibinomial audit trail is incomplete")
  invisible(TRUE)
}

s01_validate_method_output_rows <- function(rows, crosswalk) {
  numeric <- c("standard_error", "ci_lower", "ci_upper", "statistic", "numerator_df", "denominator_df")
  text <- c("endpoint_id", "member", "multiplicity_family", "family_hypothesis", "estimability",
    "standard_error_status", "confidence_interval_status", "uncertainty_reason", "ql_df_status")
  s01_require(is.data.frame(rows) && nrow(rows) > 0L && !anyDuplicated(names(rows)) &&
    all(c(numeric, text) %in% names(rows)), "shared method output required fields missing")
  s01_require(all(vapply(rows[numeric], function(v) is.double(v) &&
    length(v) == nrow(rows) && !any(is.nan(v) | is.infinite(v)), logical(1))) &&
    all(vapply(rows[text], function(v) is.character(v) && !anyNA(v) &&
      all(nzchar(trimws(v))), logical(1))), "shared method output types or missingness invalid")
  i <- match(rows$endpoint_id, crosswalk$endpoint_id)
  s01_require(!anyNA(i) && !anyDuplicated(rows$endpoint_id), "shared output endpoint identity invalid")
  ref <- crosswalk[i, , drop = FALSE]
  s01_require(identical(rows$multiplicity_family, ref$multiplicity_family) &&
    identical(rows$family_hypothesis, ref$family_hypothesis), "shared output family identity mismatch")
  ql <- ref$multiplicity_family %in% c("SC_G_MEMBERS", "SC_H_MEMBERS")
  limma <- grepl("limma", ref$primary_method, fixed = TRUE)
  program <- ref$primary_method == "paired_delta_t_test_Student_t_CI_plus_g_av" |
    (ref$primary_method == "not_run_under_frozen_rules" &
      rows$estimability == "NA_NONESTIMABLE")
  s01_require(all(ql | limma | program), "unsupported shared method output population")
  expected_member <- ifelse(ql, sub("^.*__", "", ref$family_hypothesis),
    ifelse(limma, ref$construct, "NA_NOT_APPLICABLE"))
  s01_require(identical(rows$member, expected_member), "shared output member identity mismatch")
  s01_require(all(rows$estimability %in% c("ESTIMABLE", "NA_NONESTIMABLE")), "shared output estimability invalid")
  ok <- rows$estimability == "ESTIMABLE"
  s01_require(all(is.finite(rows$statistic[ok])), "registered native statistic required")
  s01_require(all(is.na(rows$numerator_df[!ql])) && all(is.na(rows$denominator_df[!ql])) &&
    all(rows$ql_df_status[!ql] == "NA_NOT_APPLICABLE") &&
    all(rows$ql_df_status[ql & ok] == "AVAILABLE") &&
    all(rows$ql_df_status[ql & !ok] == "NA_NONESTIMABLE") &&
    all(is.na(rows$numerator_df[ql & !ok])) && all(is.na(rows$denominator_df[ql & !ok])) &&
    all(is.finite(rows$numerator_df[ql & ok]) & rows$numerator_df[ql & ok] > 0) &&
    all(is.finite(rows$denominator_df[ql & ok]) & rows$denominator_df[ql & ok] > 0) &&
    all(rows$statistic[ql & ok] >= 0), "QL degrees of freedom applicability mismatch")
  s01_require(all(is.na(rows$standard_error[ql])) && all(is.na(rows$ci_lower[ql])) &&
    all(is.na(rows$ci_upper[ql])) && all(rows$standard_error_status[ql] == "NA_UNAVAILABLE") &&
    all(rows$confidence_interval_status[ql] == "NA_UNAVAILABLE") &&
    all(rows$uncertainty_reason[ql] == "native_modern_edgeR_does_not_provide_coefficient_SE_or_95_percent_CI"),
    "native QL uncertainty exception mismatch")
  native <- !ql & ok
  s01_require(all(is.finite(rows$standard_error[native]) & rows$standard_error[native] >= 0) &&
    all(is.finite(rows$ci_lower[native]) & is.finite(rows$ci_upper[native]) &
      rows$ci_lower[native] <= rows$ci_upper[native]) &&
    all(rows$standard_error_status[native] == "AVAILABLE") &&
    all(rows$confidence_interval_status[native] == "AVAILABLE") &&
    all(rows$uncertainty_reason[native] == "NA_NOT_APPLICABLE") &&
    all(rows$standard_error_status[!ql & !ok] == "NA_NONESTIMABLE") &&
    all(rows$confidence_interval_status[!ql & !ok] == "NA_NONESTIMABLE"),
    "registered non-QL uncertainty required; member exception cannot extend to other endpoints")
  invisible(TRUE)
}
