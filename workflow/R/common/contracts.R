# Extracted maintained functions. Stable scientific identifiers and numerical order are retained.
# Numerical methods retained from the study implementation.

s01_fail <- function(message) stop(message, call. = FALSE)

s01_require <- function(condition, message) {
  if (!isTRUE(condition)) s01_fail(message)
  invisible(TRUE)
}

s01_sha256 <- function(path) {
  s01_require(file.exists(path) && !dir.exists(path),
              paste0("missing required file: ", path))
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

s01_validate_fixed_authority_file <- function(path, expected_sha256, label) {
  s01_require(is.character(expected_sha256) && length(expected_sha256) == 1L &&
                grepl("^[0-9a-f]{64}$", expected_sha256),
              paste0(label, " expected SHA-256 is malformed"))
  s01_require(identical(s01_sha256(path), expected_sha256),
              paste0("protected authority hash mismatch: ", label))
  invisible(TRUE)
}

s01_read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = character()
  )
}

s01_validate_identity_vector <- function(x, label, expected = NULL,
                                         exact_order = FALSE) {
  s01_require(length(x) > 0L, paste0(label, " identity population is empty"))
  x <- as.character(x)
  s01_require(!anyNA(x), paste0(label, " identity contains missing/NA values"))
  s01_require(!any(!nzchar(trimws(x))), paste0(label, " identity contains blank values"))
  s01_require(!anyDuplicated(x), paste0(label, " identity contains duplicated values"))
  if (!is.null(expected)) {
    expected <- s01_validate_identity_vector(expected, paste0("expected ", label))
    s01_require(setequal(x, expected),
                paste0(label, " identity contains unexpected or missing values"))
    if (isTRUE(exact_order)) {
      s01_require(identical(x, expected), paste0(label, " identity order mismatch"))
    }
  }
  x
}

s01_assert_exact_order <- function(observed, expected, label) {
  observed <- s01_validate_identity_vector(observed, label)
  expected <- s01_validate_identity_vector(expected, paste0("expected ", label))
  if (length(observed) != length(expected)) {
    s01_fail(paste0(label, " cardinality mismatch"))
  }
  if (anyDuplicated(observed)) s01_fail(paste0(label, " contains duplicated identities"))
  if (!setequal(observed, expected)) {
    s01_fail(paste0(label, " contains a substituted or missing identity"))
  }
  if (!identical(observed, expected)) s01_fail(paste0(label, " order mismatch"))
  invisible(TRUE)
}

s01_assert_unique_key <- function(data, key, label = deparse(substitute(data))) {
  s01_require(is.data.frame(data) && all(key %in% names(data)),
              paste0(label, " lacks its required key"))
  value <- data[, key, drop = FALSE]
  encoded_values <- lapply(value, as.character)
  s01_require(!anyNA(value) && !any(vapply(encoded_values, function(x) {
    any(!nzchar(trimws(x)))
  }, logical(1))),
              paste0(label, " key contains missing values"))
  encoded <- do.call(paste, c(encoded_values, sep = "\r"))
  s01_require(!anyDuplicated(encoded), paste0(label, " key is duplicated"))
  invisible(TRUE)
}

s01_checked_join <- function(left, right, by, relationship = "one-to-one") {
  s01_require(is.data.frame(left) && is.data.frame(right),
              "join inputs must be data frames")
  s01_require(all(by %in% names(left)) && all(by %in% names(right)),
              "join key is absent")
  relationship <- match.arg(relationship, c("one-to-one", "many-to-one"))
  if (identical(relationship, "one-to-one")) s01_assert_unique_key(left, by, "left")
  s01_assert_unique_key(right, by, "right")
  left_key <- do.call(paste, c(left[, by, drop = FALSE], sep = "\r"))
  right_key <- do.call(paste, c(right[, by, drop = FALSE], sep = "\r"))
  index <- match(left_key, right_key)
  s01_require(!anyNA(index), "join contains unmatched left identities")
  appended <- setdiff(names(right), by)
  collision <- intersect(appended, names(left))
  s01_require(!length(collision), "join has non-key column collisions")
  cbind(left, right[index, appended, drop = FALSE], stringsAsFactors = FALSE)
}

s01_coverage_decision <- function(observed, required, mode,
                                   available = observed) {
  observed <- s01_validate_identity_vector(observed, "observed gene")
  required <- s01_validate_identity_vector(required, "required program member")
  available <- s01_validate_identity_vector(available, "available gene")
  s01_require(all(observed %in% required),
              "observed gene identity contains an unexpected program member")
  present <- intersect(required, intersect(observed, available))
  if (identical(mode, "exact")) {
    estimable <- !anyDuplicated(required) && identical(present, required)
  } else if (identical(mode, "single_cell_G")) {
    substrate <- c("ACSL4", "LPCAT3")
    oxidation <- c("ALOX5", "ALOX12", "ALOX15", "ALOX15B", "POR")
    estimable <- length(present) >= 5L && sum(present %in% substrate) >= 1L &&
      sum(present %in% oxidation) >= 2L
  } else if (identical(mode, "single_cell_H")) {
    estimable <- length(present) >= 3L &&
      any(present %in% c("ALDH3A1", "ALDH3A2")) &&
      any(present %in% c("ALDH2", "ALDH7A1"))
  } else {
    s01_fail("unknown coverage decision mode")
  }
  list(
    estimable = estimable,
    state = if (estimable) "ESTIMABLE" else "NA_NONESTIMABLE",
    observed_members = present,
    required_members = required
  )
}

s01_resolve_bulk_roles <- function(core_eligibility, backup_eligible) {
  expected <- c("GSE29272", "GSE27342")
  s01_require(is.logical(core_eligibility) &&
                identical(names(core_eligibility), expected) && !anyNA(core_eligibility),
              "both core eligibility decisions must be named and fixed")
  failures <- sum(!core_eligibility)
  if (failures == 0L) {
    return(list(branch = "normal", first = "GSE29272", replication = "GSE27342"))
  }
  if (failures == 2L) s01_fail("both core sources failed; the backup is not a rescue")
  s01_require(isTRUE(backup_eligible),
              "the conditional backup is unavailable under its technical trigger")
  survivor <- names(core_eligibility)[core_eligibility]
  list(branch = "conditional_backup", first = survivor, replication = "GSE122401")
}

s01_registered_method_resolution <- function() {
  list(
    outcome = "APPROVED_METHOD_IMPLEMENTED",
    method = "edgeR_QL_tissue_plus_donor",
    issue = "RESOLVED_BY_PUBLISHED_S01_AMENDMENT_2026_09_07",
    fixed = list(
      edgeR_version = "4.10.1",
      version_authority = "accepted_receipt.repository.scientific_contract_commit:renv.lock",
      scientific_contract_commit = "b8a42905fc610de81f52d80bb29deeb4e7e06f68",
      renv_lock_git_blob = "92322f52b4d6291e0c3baeffedf6a8cc876951e9",
      renv_lock_sha256 = "ad12c52f4d112c067ef9705d288eabac42b49398514edc4607114e66f383a881",
      scientific_rules = paste0(
        "plan_15.3-15.5: marker_gate; eligible_raw_UMI_donor_by_tissue_by_compartment;",
        "20_cells_and_10000_UMI_per_unit;8_paired_donors;",
        "member_detection_nonzero_in_at_least_50_percent_of_units_and_both_tissues_in_3_donors;",
        "TMM;edgeR_QL;tissue_plus_donor;no_redundant_chemistry_term;",
        "crosswalk_Tumor_minus_NAT;SC_G_MEMBERS_5_and_SC_H_MEMBERS_9"
      )
    ),
    amendment = list(
      receipt = "data/manifests/registration/OSF_S01_METHODS_AMENDMENT_2026-09-07.json",
      revision = "6a9e4b27f508eba67ee25805",
      plan_section = "15.5.1",
      configuration = "single_cell_member_model"
    ),
    production_execution_permitted = FALSE,
    synthetic_check_scope = "ACTUAL_APPROVED_HELPER_WHOLE_SYNTHETIC_BACKGROUND_ONLY",
    new_scientific_rule_invented = FALSE
  )
}
