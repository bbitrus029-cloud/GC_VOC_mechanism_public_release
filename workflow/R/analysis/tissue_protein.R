# Extracted maintained implementation; stable machine names and seeds retained.
tissue_protein_module <- function(context, external_bulk_processed) {
  s13_root <- context$root
  s13_sha <- function(p) digest::digest(file=context$input(p),algo="sha256")
  s13_now <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")
  s13_tsv <- function(x,p) context$write_tsv(x,p)
  s13_rds <- function(x,p) context$write_rds(x,p)
  s13_json <- function(x,p) context$write_json(x,p)
  s13_text <- function(x,p) context$write_text(x,p)
  s13_guard <- function(...) context$guard()
  s13_effect_guard <- function(...) context$guard()
  s13_memory <- function(stage,minimum=8,prestart=FALSE) context$resource_check(minimum,stage)
  s13_hashes <- function(paths) data.frame(path=paths,bytes=file.info(paths)$size,sha256=vapply(paths,s13_sha,""))
  s13_m <- context$output("metadata"); s13_p <- context$output("processed")
  s13_r <- context$output("results"); s13_q <- context$output("qc")
  s13_old <- "data/manifests/post_effect_supplement/s10"
s13_read <- function(p) read.delim(context$input(p), stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", 
    comment.char = "")

s13_runtime <- function(...) context$runtime("supplementary")

s13_roster <- function() {
    x <- s13_read("config/post_effect_hypothesis_estimand_crosswalk.tsv")
    x <- x[x$module %in% c("C-D7", "C-composition", "D"), , drop = FALSE]
    if (nrow(x) != 13L) 
        x <- s13_read("config/post_effect_hypothesis_estimand_crosswalk.tsv")
    x <- x[x$multiplicity_family %in% c("PEC_DIFFERENTIATION", "PEC_COMPOSITION", "PED_H4", "PED_H_MEMBERS"), , 
        drop = FALSE]
    s01_require(nrow(x) == 13L && !anyDuplicated(x$hypothesis_id) && identical(as.integer(table(factor(x$multiplicity_family, 
        levels = c("PEC_DIFFERENTIATION", "PEC_COMPOSITION", "PED_H4", "PED_H_MEMBERS")))), c(4L, 4L, 1L, 4L)), 
        "S13 exact family roster")
    x
}

s14_read <- function(p) read.delim(context$input(p), stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", 
    comment.char = "")

s13_pairs <- function(md, samples) {
    s01_validate_identity_vector(samples, "S13 sample")
    s01_require(identical(md$sample_id, samples) && !anyNA(md[, c("sample_id", "patient_id", "tissue")]), "S13 metadata ordering/identity")
    s01_require(all(md$tissue %in% c("NAT", "Tumor")) && !anyDuplicated(paste(md$patient_id, md$tissue)) && all(table(md$patient_id) == 
        2L), "S13 incomplete/duplicate pair")
    ids <- sort(unique(md$patient_id), method = "radix")
    k <- paste(md$patient_id, md$tissue)
    n <- match(paste(ids, "NAT"), k)
    t <- match(paste(ids, "Tumor"), k)
    s01_require(!anyNA(c(n, t)), "S13 unresolved paired tissue")
    data.frame(patient_id = ids, adjacent_sample = md$sample_id[n], tumor_sample = md$sample_id[t])
}

s13_select_frozen <- function(x, md, reference) {
    s01_validate_identity_vector(rownames(x), "S13 expression feature")
    s13_pairs(md, colnames(x))
    s01_require(all(reference$metadata$sample_id %in% colnames(x)), "S13 reference samples absent")
    j <- match(reference$metadata$sample_id, md$sample_id)
    s01_require(identical(md$patient_id[j], reference$metadata$patient_id) && identical(md$tissue[j], reference$metadata$tissue), 
        "S13 reference identity mismatch")
    list(expression = x[, reference$metadata$sample_id, drop = FALSE], metadata = reference$metadata, pairs = s13_pairs(reference$metadata, 
        reference$metadata$sample_id))
}

s13_expression <- function(source, reference) {
    a <- context$read_rds(reference$input_path)
    if (source == "GSE122401") {
        x <- log2(a$abundance + 1)
        map <- a$member_mapping
        j <- match(reference$members, map$gene)
        s01_require(!anyNA(j) && !anyDuplicated(map$ensembl_id[j]), "S13 D7 feature bridge")
        rownames(x)[match(map$ensembl_id[j], rownames(x))] <- reference$members
        md <- a$metadata
    }
    else if (source == "TCGA_STAD") {
        x <- a$model$voom$E
        md <- a$metadata
    }
    else {
        x <- a
        name <- if (source == "GSE27342_core") 
            "GSE27342"
        else source
        d <- s13_read(paste0("reports/qc/s02_bulk_ingress/", name, "_pair_ledger.tsv"))
        d <- d[!d$pair_excluded & d$core_population, ]
        md <- data.frame(sample_id = d$neutral_id, patient_id = d$patient_id, tissue = d$tissue)
        md <- md[match(colnames(x), md$sample_id), , drop = FALSE]
    }
    s13_select_frozen(x, md, reference)
}

s13_check_reference <- function(input, r, moments, population) {
    fields <- c("sample_id", "patient_id", "tissue")
    s01_require(identical(lapply(input$metadata[fields], as.character), lapply(population[fields], as.character)), 
        "S13 frozen ordered population")
    s01_require(identical(moments$member, r$members) && length(input$pairs$patient_id) == length(r$patient_ids) && 
        identical(input$pairs$patient_id, r$patient_ids), "S13 D7 reference roster")
    x <- matrix(NA_real_, length(r$members), ncol(input$expression), dimnames = list(r$members, colnames(input$expression)))
    available <- intersect(r$members, rownames(input$expression))
    x[available, ] <- input$expression[available, , drop = FALSE]
    s01_require(identical(setdiff(r$members, available), r$missing) && all(is.finite(x[available, , drop = FALSE])), 
        "S13 frozen member coverage")
    s01_require(digest::digest(list(x = x, metadata = input$metadata), algo = "sha256") == r$reference_sha256, 
        "S13 D7 internal identity")
    s01_require(isTRUE(all.equal(unname(r$mean), moments$pooled_mean, tolerance = 1e-13)) && isTRUE(all.equal(unname(r$sd), 
        moments$sample_SD, tolerance = 1e-13)) && all(moments$sd_denominator == ncol(x) - 1) && identical(unname(r$constant), 
        moments$constant), "S13 pooled constants")
    invisible(TRUE)
}

s13_join_h <- function(h, pairs) {
    s01_require(all(pairs$patient_id %in% h$pairs$patient_id), "S13 H patient absent")
    j <- match(pairs$patient_id, h$pairs$patient_id)
    s01_require(identical(pairs$adjacent_sample, h$pairs$adjacent_sample[j]) && identical(pairs$tumor_sample, h$pairs$tumor_sample[j]), 
        "S13 H/D7 bridge conflict")
    n <- match(pairs$adjacent_sample, h$scores$sample_id)
    t <- match(pairs$tumor_sample, h$scores$sample_id)
    data.frame(patient_id = pairs$patient_id, delta_H = h$scores$rank_H[t] - h$scores$rank_H[n], delta_H_standardized = h$scores$standardized_H[t] - 
        h$scores$standardized_H[n])
}

s13_protein_input <- function() {
    x <- read.delim(context$input("data/raw/post_effect_supplement/protein/Gastric_Cancer_Korea_Proteome.itraq.tsv"), 
        check.names = FALSE, stringsAsFactors = FALSE)
    ids <- s13_read(file.path(s13_old, "D_ratio_column_identity.tsv"))
    bridge <- s13_read(file.path(s13_old, "D_sample_RNA_overlap_bridge.tsv"))
    cols <- names(x)[grepl("^Unshared Log T", names(x))]
    actual <- pe_protein_ratio_identities(cols, bridge)
    s01_require(identical(actual, ids), "S13 protein declared ratio identities/order")
    map <- s13_read(file.path(s13_old, "D_H4_feature_mapping.tsv"))
    j <- match(map$Gene, x$Gene)
    s01_require(!anyDuplicated(x$Gene) && !anyNA(j) && all(as.character(x$NCBIGeneID[j]) == as.character(map$NCBIGeneID)) && 
        identical(x$Authority[j], map$Authority) && identical(x$Organism[j], map$Organism), "S13 protein unique H mapping")
    y <- as.matrix(x[j, cols, drop = FALSE])
    rownames(y) <- map$Gene
    s01_require(is.numeric(y), "S13 ratio numeric scale")
    pop <- s13_read(file.path(s13_old, "D_common_H4_population.tsv"))
    patients <- sort(unique(ids$patient_id), method = "radix")
    finite <- vapply(patients, function(p) all(is.finite(y[, ids$patient_id == p, drop = FALSE])), TRUE)
    counts <- as.integer(table(factor(ids$patient_id, levels = patients)))
    s01_require(identical(patients, pop$patient_id) && identical(unname(finite), pop$common_H4_complete) && identical(counts, 
        pop$declared_ratio_replicates) && all(pop$RNA_patient_id == pop$patient_id), "S13 protein frozen common population")
    list(ratios = y, identities = ids, population = pop, bridge = bridge)
}

s13_aggregate <- function(x, ids, pop, members) {
    s01_require(identical(rownames(x), members) && length(members) == 4 && !anyDuplicated(members), "S13 ratio-only exact H4 member roster")
    s01_validate_identity_vector(ids$ratio_column, "S13 declared ratio")
    s01_require(identical(colnames(x), ids$ratio_column) && all(ids$patient_id %in% pop$patient_id), "S13 ratio order/patient identity")
    z <- pe_protein_patient_ratios(x, ids)
    s01_require(identical(z$patient_ids, pop$patient_id[pop$common_H4_complete]) && identical(z$excluded, pop$patient_id[!pop$common_H4_complete]), 
        "S13 common all-repeat all-member population")
    z
}

s13_ratio_test <- function(y, ids) {
    z <- pe_protein_ratio(y, ids, centered = TRUE)
    out <- list(state = z$state, reason = if (is.null(z$reason)) "NONE" else z$reason, effect = if (length(y) && 
        all(is.finite(y))) mean(y) else NA_real_, biological_n = length(ids), standard_error = NA_real_, residual_df = NA_real_, 
        statistic = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_)
    if (z$state == "ESTIMABLE") {
        t <- t.test(y, mu = 0, alternative = "two.sided", conf.level = 0.95)
        for (k in c("effect", "ci_lower", "ci_upper", "p_value")) out[[k]] <- z[[k]]
        out$standard_error <- sd(y)/sqrt(length(y))
        out$residual_df <- unname(t$parameter)
        out$statistic <- unname(t$statistic)
        s01_require(all(is.finite(unlist(out[c("effect", "standard_error", "residual_df", "statistic", "ci_lower", 
            "ci_upper", "p_value")]))), "S13 invalid native t output")
    }
    out
}

s13_preflight <- function() {
    s13_guard()
    s13_memory("preflight")
    refs <- context$read_rds("data/processed/post_effect_supplement/s10/D7_reference_freeze.rds")
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    s01_require(identical(names(refs), cfg$C_D7$cohorts), "S13 cohort roster")
    context$guard() # All reference inputs are explicitly bound by the current input manifest.

    moments <- s13_read(file.path(s13_old, "D7_pooled_reference_moments.tsv"))
    pop <- s13_read(file.path(s13_old, "D7_reference_population.tsv"))
    rows <- list()
    bridges <- list()
    for (s in names(refs)) {
        r <- refs[[s]]
        inp <- s13_expression(s, r)
        s13_check_reference(inp, r, moments[moments$source == s, ], pop[pop$source == s, ])
        h <- s13_h(s)
        j <- match(inp$pairs$patient_id, h$pairs$patient_id)
        s01_require(!anyNA(j) && identical(inp$pairs$adjacent_sample, h$pairs$adjacent_sample[j]) && identical(inp$pairs$tumor_sample, 
            h$pairs$tumor_sample[j]), "S13 frozen H-D7 one-to-one bridge")
        bridges[[s]] <- cbind(source = s, inp$pairs, H_parent_population_n = nrow(h$pairs), D7_reference_n = nrow(inp$pairs), 
            included_M1 = r$state == "READY", reason = if (r$state == "READY") 
                "NONE"
            else paste0(r$reason, ":", paste(r$missing, collapse = ";")))
        rows[[s]] <- data.frame(branch = "C_D7", source = s, biological_n = nrow(inp$pairs), state = if (r$state == 
            "READY") 
            "QUALIFIED_AWAITING_EFFECT_RELEASE"
        else "NA_NONESTIMABLE", reason = if (r$state == "READY") 
            "NONE"
        else paste0(r$reason, ":", paste(r$missing, collapse = ";")), scale = r$scale)
    }
    d <- s13_protein_input()
    assay <- list(estimand=cfg$D$estimand)
    s01_require(cfg$D$selected_branch == "paired_tumor_adjacent_log2_ratio_only" &&
        cfg$D$log_base == 2 && cfg$D$ratio_field == "Unshared_Log_T_over_N" &&
        cfg$D$reference_pool == "none_patient_specific_adjacent_denominator" &&
        cfg$D$normalization == "retain_CDAP_per_ratio_median_centering_no_additional_normalization" &&
        cfg$D$repeated_measurement_rule == "equal_arithmetic_mean_of_all_declared_log2_ratio_replicates_at_later_effect_stage",
        "S13 ratio-only assay semantics")

    s01_require(ncol(d$ratios) == 130L && nrow(d$population) == 80L && sum(d$population$declared_ratio_replicates == 
        2L) == 50L && sum(d$population$declared_ratio_replicates == 1L) == 30L && all(d$population$common_H4_complete), 
        "S13 protein admission counts")
    a <- s13_h("GSE122401")
    s01_require(setequal(d$population$patient_id, a$pairs$patient_id), "S13 full RNA overlap")
    rows$D <- data.frame(branch = "D", source = "PDC000214_Mun2019", biological_n = 80, state = "QUALIFIED_RATIO_ONLY_AWAITING_EFFECT_RELEASE", 
        reason = "NONE", scale = assay$estimand)
    s13_tsv(do.call(rbind, bridges), file.path(s13_m, "C_H_D7_population_bridge.tsv"))
    s13_tsv(do.call(rbind, rows), file.path(s13_m, "branch_readiness.tsv"))
    roster <- s13_roster()
    s13_tsv(roster, file.path(s13_m, "hypothesis_roster.tsv"))

    s13_json(list(status = "PASS", at_utc = s13_now(), D7_reference_identity = "PASS", ordered_population_moments = "PASS", 
        H_D7_bridges = "PASS", protein_assay_mapping_repeat_population_overlap = "PASS", roster_slots = 13, real_D7_scores_or_contrasts = FALSE, 
        real_protein_averages_or_tests = FALSE, C_composition = "NOT_COMPUTED_NOT_RELEASED", runtime = s13_runtime()), 
        file.path(s13_q, "PREFLIGHT.json"))
}

s13_h <- function(source) {
    refs <- context$read_rds("data/processed/post_effect_supplement/s10/D7_reference_freeze.rds")
    r <- refs[[source]]
    a <- s13_expression(source, r)
    x <- a$expression
    md <- a$metadata
    if (source == "GSE122401") {
        z <- context$read_rds(file.path(external_bulk_processed, "rank_scores.rds"))
        alt <- s14_read(file.path(external_bulk_processed, "diagnostic_scores.tsv"))
        alt <- alt[alt$analysis_id == "PEA.standardized_member_mean.H_core4", ]
        out <- data.frame(sample_id = colnames(z$scores), rank_H = unname(z$scores["H_core4", ]), standardized_H = alt$score[match(colnames(z$scores), 
            alt$sample_id)])
        md <- z$metadata
    }
    else {
        members <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))$programs$H_core4
        h <- s01_rank_program_score(x, members, rownames(x))
        z <- s01_standardized_member_mean(x, members)
        s01_require(h$estimable && z$estimable, "S14 frozen H scoring")
        out <- data.frame(sample_id = h$sample_id, rank_H = h$score, standardized_H = z$score)
    }
    list(scores = out, metadata = md, pairs = s13_pairs(md, out$sample_id))
}
s13_C <- function(release) {
    s01_require(file.exists(release), "S13 C release dependency")
    s13_effect_guard("C_D7")
    s13_memory("C_start")
    started <- s13_now()
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    refs <- context$read_rds("data/processed/post_effect_supplement/s10/D7_reference_freeze.rds")
    descriptions <- member_source <- models <- hc3 <- support <- loo <- summary <- paired <- full_H <- list()
    fits <- list()
    states <- list()
    for (s in names(refs)) {
        r <- refs[[s]]
        a <- s13_expression(s, r)
        pp <- a$pairs
        h <- s13_h(s)
        dh <- s13_join_h(h, pp)
        source_rows <- list()
        for (g in r$members) {
            measurable <- g %in% rownames(a$expression)
            sd <- r$sd[match(g, r$members)]
            d <- if (measurable) 
                a$expression[g, pp$tumor_sample] - a$expression[g, pp$adjacent_sample]
            else rep(NA_real_, nrow(pp))
            state <- if (!measurable) 
                "NA_UNAVAILABLE_REQUIRED_MEMBER"
            else if (r$constant[match(g, r$members)]) 
                "CONSTANT_MEASURED_MEMBER"
            else "MEASURED"
            mu <- if (measurable) 
                mean(d)
            else NA_real_
            descriptions[[length(descriptions) + 1]] <- data.frame(source = s, member = g, measurement_state = state, 
                D7_state = r$state, D7_reason = r$reason, scale = r$scale, reference_pairs = nrow(pp), reference_samples = r$N_samples, 
                pooled_mean = r$mean[match(g, r$members)], pooled_variance = sd^2, pooled_sample_SD = sd, paired_n = if (measurable) 
                  nrow(pp)
                else 0L, native_mean_paired_change = mu, direction = if (measurable) 
                  s01_direction_label(mu)
                else "NA_UNAVAILABLE", standardized_mean_paired_change = if (measurable && !r$constant[match(g, 
                  r$members)]) 
                  mu/sd
                else NA_real_, standardization = "unchanged_S10_pooled_gene_mean_and_sample_SD", p_value = NA_real_, 
                q_value = NA_real_, inference_role = "descriptive_no_discovery_family", timing = "post_effect_exploratory", 
                row.names = NULL)
            if (measurable) 
                member_source[[length(member_source) + 1]] <- data.frame(source = s, member = g, patient_id = as.character(pp$patient_id), 
                  adjacent_value = unname(a$expression[g, pp$adjacent_sample]), tumor_value = unname(a$expression[g, 
                    pp$tumor_sample]), native_delta = unname(d), standardized_delta = if (!r$constant[match(g, 
                    r$members)]) 
                    unname(d/sd)
                  else NA_real_, scale = r$scale)
        }
        delta_D <- if (r$state == "READY") {
            sc <- pe_d7_score(a$expression, r)
            unname(sc[pp$tumor_sample] - sc[pp$adjacent_sample])
        }
        else rep(NA_real_, nrow(pp))
        dd <- data.frame(source = s, patient_id = as.character(pp$patient_id), adjacent_sample = pp$adjacent_sample, 
            tumor_sample = pp$tumor_sample, delta_H = dh$delta_H, delta_H_standardized = dh$delta_H_standardized, 
            delta_D = delta_D, D7_state = r$state, model_population = if (r$state == "READY") 
                "M1_COMMON_COMPLETE"
            else "D7_TERMINAL_PLATFORM_MEMBER_FAILURE", missingness_reason = if (r$state == "READY") 
                "NONE"
            else paste(r$missing, collapse = ";"))
        paired[[s]] <- dd
        for (spec in c("rank_H", "standardized_H")) {
            hp <- h$pairs
            v <- h$scores[[spec]]
            n <- match(hp$adjacent_sample, h$scores$sample_id)
            t <- match(hp$tumor_sample, h$scores$sample_id)
            te <- s01_paired_delta_t(v[t], v[n], as.character(hp$patient_id))
            full_H[[length(full_H) + 1]] <- data.frame(source = s, view = "FULL_ACCEPTED_H_POPULATION", H_definition = spec, 
                biological_n = nrow(hp), effect = te$effect, ci_lower = te$ci_lower, ci_upper = te$ci_upper, state = te$state, 
                role = "descriptive", p_value = NA_real_, q_value = NA_real_, timing = "post_effect_exploratory")
        }
        if (r$state != "READY") {
            states[[s]] <- list(state = "NA_NONESTIMABLE", reason = paste0(r$reason, ":", paste(r$missing, collapse = ";")), 
                biological_n = nrow(pp))
            next
        }
        primary <- pe_regression(dd, "M1", minimum_n = cfg$C_D7$minimum_pairs)
        fits[[s]] <- list(primary = primary)
        states[[s]] <- list(state = primary$state, reason = primary$reason, biological_n = primary$biological_n)
        if (primary$state != "ESTIMABLE") 
            next
        keep <- match(primary$patient_ids, dd$patient_id)
        d <- dd[keep, ]
        d_alt <- d
        d_alt$delta_H <- d$delta_H_standardized
        alt <- pe_regression(d_alt, "M1", minimum_n = cfg$C_D7$minimum_pairs, influence = FALSE)
        m0 <- pe_regression(d, "M0", minimum_n = cfg$C_D7$minimum_pairs, influence = FALSE)
        fits[[s]]$alternative <- alt
        fits[[s]]$M0 <- m0
        for (name in c("primary", "alternative", "M0")) {
            fit <- fits[[s]][[name]]
            if (fit$state == "ESTIMABLE") {
                co <- fit$coefficients
                co$source <- s
                co$model <- if (name == "M0") 
                  "M0"
                else "M1"
                co$H_definition <- if (name == "alternative") 
                  "standardized_member_mean_H"
                else "accepted_rank_H"
                co$biological_n <- fit$biological_n
                co$state <- fit$state
                co$reason <- fit$reason
                co$role <- if (name == "primary") 
                  "post_effect_exploratory"
                else if (name == "alternative") 
                  "same_data_robustness"
                else "descriptive"
                co$timing <- "post_effect_exploratory"
                co$units <- if (name == "M0") 
                  "unadjusted_mean_paired_H_score_change"
                else ifelse(co$coefficient == "delta_D", "H_score_units_per_one_D7_paired_change_score_unit", "H_score_units_at_zero_D7_paired_change")
                co$p_value[co$coefficient != "delta_D" | name != "primary"] <- NA_real_
                co$q_value <- NA_real_
                co$multiplicity_family <- ifelse(co$coefficient == "delta_D" & name == "primary", "PEC_DIFFERENTIATION", 
                  "NONE_DESCRIPTIVE")
                models[[length(models) + 1]] <- co
            }
            else {
                co <- data.frame(coefficient = c("(Intercept)", "delta_D"), effect = NA_real_, standard_error = NA_real_, 
                  statistic = NA_real_, residual_df = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, 
                  source = s, model = "M1", H_definition = "standardized_member_mean_H", biological_n = fit$biological_n, 
                  state = fit$state, reason = fit$reason, role = "same_data_robustness", timing = "post_effect_exploratory", 
                  units = c("H_score_units_at_zero_D7_paired_change", "H_score_units_per_one_D7_paired_change_score_unit"), 
                  q_value = NA_real_, multiplicity_family = "NONE_DESCRIPTIVE")
                models[[length(models) + 1]] <- co
            }
        }
        hc <- primary$hc3
        names_coef <- primary$coefficients$coefficient
        hc3[[s]] <- data.frame(source = s, model = "M1", H_definition = "accepted_rank_H", coefficient = names_coef, 
            effect = primary$coefficients$effect, standard_error = unname(hc$se), ci_lower = unname(hc$lower), 
            ci_upper = unname(hc$upper), residual_df = df.residual(primary$fit), state = hc$state, reason = hc$reason, 
            p_value = NA_real_, q_value = NA_real_, role = "same_data_robustness", timing = "post_effect_exploratory", 
            inference = "HC3_Student_t_n_minus_p_finite_sample_approximation")
        ss <- primary$diagnostics
        support[[s]] <- data.frame(source = s, biological_n = primary$biological_n, delta_D_min = min(d$delta_D), 
            delta_D_max = max(d$delta_D), zero_supported = ss$zero_covariate_supported, intercept_label = ss$intercept_label, 
            residual_df = df.residual(primary$fit), R_squared = summary(primary$fit)$r.squared, max_leverage = max(hatvalues(primary$fit)), 
            HC3_state = hc$state, HC3_reason = hc$reason, exclusion_count = length(primary$excluded), excluded_patients = paste(primary$excluded, 
                collapse = ";"), interpretation = "between_patient_paired_change_association_no_percent_of_mean_biologically_explained")
        lo <- primary$loo
        lo$source <- s
        lo$role <- "same_data_robustness"
        lo$timing <- "post_effect_exploratory"
        lo$refit_n <- nrow(d) - 1L
        loo[[s]] <- lo
        su <- ss$influence_summary
        su$source <- s
        su$role <- "same_data_robustness"
        su$timing <- "post_effect_exploratory"
        summary[[s]] <- su
    }
    out <- list(member_descriptions = do.call(rbind, descriptions), member_source = do.call(rbind, member_source), 
        coefficients = do.call(rbind, models), HC3 = do.call(rbind, hc3), support = do.call(rbind, support), LOO = do.call(rbind, 
            loo), influence_summary = do.call(rbind, summary), paired_source = do.call(rbind, paired), full_H = do.call(rbind, 
            full_H), fits = fits, states = states)
    p <- s13_rds(out, file.path(s13_p, "C_D7_analysis.rds"))
    s13_json(list(status = "EXECUTED", started_at_utc = started, completed_at_utc = s13_now(), branch = "C_D7", 
        release = release, release_sha256 = s13_sha(release), output = s13_hashes(p)), file.path(s13_m, "C_D7_EXECUTION.json"))
    s13_memory("C_complete", 6)
    p
}

s13_D <- function(release) {
    s01_require(file.exists(release), "S13 D release dependency")
    s13_effect_guard("D")
    s13_memory("D_start")
    started <- s13_now()
    x <- s13_protein_input()
    cfg <- yaml::read_yaml(context$input("config/post_effect_supplement.yml"))
    z <- s13_aggregate(x$ratios, x$identities, x$population, cfg$D$required_members)
    rows <- lapply(rownames(z$values), function(g) c(list(member = g), s13_ratio_test(unname(z$values[g, ]), z$patient_ids)))
    rows <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
    patients <- do.call(rbind, lapply(rownames(z$values), function(g) data.frame(member = g, patient_id = z$patient_ids, 
        declared_ratio_repeats = x$population$declared_ratio_replicates[match(z$patient_ids, x$population$patient_id)], 
        centered_log2_T_over_adjacent_ratio = unname(z$values[g, ]), scale = "source_per_ratio_median_centered_log2_T_over_N", 
        biological_unit = "patient")))
    out <- list(results = rows, patient_values = z, patient_source = patients, assay = "ratio_only_Unshared_Log_T_over_N_no_renormalization", 
        H4_score = "NA_UNAVAILABLE", g_av = "NA_UNAVAILABLE")
    p <- s13_rds(out, file.path(s13_p, "D_ratio_analysis.rds"))
    s13_json(list(status = "EXECUTED", started_at_utc = started, completed_at_utc = s13_now(), branch = "D", release = release, 
        release_sha256 = s13_sha(release), output = s13_hashes(p)), file.path(s13_m, "D_EXECUTION.json"))
    s13_memory("D_complete", 6)
    p
}

s13_assemble <- function(C, D) {
    s13_effect_guard("C_D7")
    s13_effect_guard("D")
    c <- context$read_rds(C)
    d <- context$read_rds(D)
    r <- s13_roster()
    names(r)[names(r) == "effect_release"] <- "crosswalk_effect_release_at_generation"
    names(r)[names(r) == "frozen_branch_state"] <- "crosswalk_state_at_generation"
    r$state <- "NA_UNAVAILABLE"
    r$reason <- "NONE"
    r$biological_n <- NA_integer_
    for (k in c("effect", "standard_error", "residual_df", "statistic", "ci_lower", "ci_upper", "p_value", "q_value")) r[[k]] <- NA_real_
    r$measurement_scale <- "NA_UNAVAILABLE"
    r$dependence <- "same_existing_dataset_no_new_independent_replication"
    r$current_effect_access <- "NOT_RELEASED"
    r$uncertainty_status <- "NA_UNAVAILABLE"
    r$units <- "NA_UNAVAILABLE"
    for (i in seq_len(nrow(r))) {
        fam <- r$multiplicity_family[i]
        s <- r$source[i]
        if (fam == "PEC_DIFFERENTIATION") {
            st <- c$states[[s]]
            r$state[i] <- st$state
            r$reason[i] <- st$reason
            r$biological_n[i] <- st$biological_n
            r$current_effect_access[i] <- "RELEASED_C_D7"
            r$measurement_scale[i] <- "accepted_rank_H_change_vs_fixed_S10_D7_change"
            r$units[i] <- "rank_H_score_units_per_one_D7_paired_change_score_unit"
            if (st$state == "ESTIMABLE") {
                z <- c$coefficients[c$coefficients$source == s & c$coefficients$H_definition == "accepted_rank_H" & 
                  c$coefficients$model == "M1" & c$coefficients$coefficient == "delta_D", ]
                s01_require(nrow(z) == 1, "S13 primary coefficient identity")
                for (k in c("effect", "standard_error", "residual_df", "statistic", "ci_lower", "ci_upper", "p_value")) r[[k]][i] <- z[[k]]
                r$uncertainty_status[i] <- "AVAILABLE_OLS_Student_t_residual_df"
            }
        }
        else if (fam == "PEC_COMPOSITION") {
            r$state[i] <- if (s == "TCGA_STAD") 
                "PENDING"
            else "NA_UNAVAILABLE"
            r$reason[i] <- if (s == "TCGA_STAD") 
                "PENDING_REFERENCE_QUALIFICATION_NOT_AUTHORIZED_S13"
            else if (s == "GSE122401") 
                "UNRESOLVED_LINEAR_MIXTURE_SCALE"
            else "ARRAY_M1_ONLY"
        }
        else if (fam == "PED_H4") {
            r$reason[i] <- "RATIO_ONLY_NO_SEPARATE_TISSUE_MEASUREMENTS_NO_H4_SCORE"
            r$biological_n[i] <- 80L
            r$current_effect_access[i] <- "RELEASED_D_RATIO_MEMBERS_ONLY"
            r$dependence[i] <- "same_80_patients_as_Module_A_dependent_cross_layer_evidence"
        }
        else {
            z <- d$results[d$results$member == r$family_hypothesis[i], ]
            s01_require(nrow(z) == 1, "S13 D member identity")
            for (k in c("state", "reason", "biological_n", "effect", "standard_error", "residual_df", "statistic", 
                "ci_lower", "ci_upper", "p_value")) r[[k]][i] <- z[[k]]
            r$measurement_scale[i] <- "source_per_ratio_median_centered_log2_T_over_N_equal_patient_repeat_mean"
            r$units[i] <- "mean_centered_log2_tumor_adjacent_ratio"
            r$current_effect_access[i] <- "RELEASED_D"
            r$dependence[i] <- "same_80_patients_as_Module_A_dependent_cross_layer_evidence"
            if (z$state == "ESTIMABLE") 
                r$uncertainty_status[i] <- "AVAILABLE_one_sample_Student_t_95_percent_interval"
        }
    }
    closure <- list()
    for (f in unique(r$multiplicity_family)) {
        j <- which(r$multiplicity_family == f)
        p <- setNames(r$p_value[j], r$hypothesis_id[j])
        cl <- pe_close_family(p, r$state[j], names(p))
        r$q_value[j] <- unname(cl$q)
        closure[[f]] <- data.frame(family = f, frozen_size = length(j), estimable = sum(r$state[j] == "ESTIMABLE"), 
            nonestimable = sum(r$state[j] == "NA_NONESTIMABLE"), unavailable = sum(r$state[j] == "NA_UNAVAILABLE"), 
            pending = sum(r$state[j] == "PENDING"), closed = cl$closed)
    }
    co <- c$coefficients
    j <- which(co$multiplicity_family == "PEC_DIFFERENTIATION")
    rr <- r[r$multiplicity_family == "PEC_DIFFERENTIATION", ]
    co$q_value[j] <- rr$q_value[match(co$source[j], rr$source)]
    tables <- list(all_hypothesis_results = r, family_closure = do.call(rbind, closure), D7_member_descriptions = c$member_descriptions, 
        C_model_coefficients = co, C_HC3 = c$HC3, C_support = c$support, C_leave_one_patient_out = c$LOO, C_influence_summary = c$influence_summary, 
        C_full_unadjusted_H = c$full_H, C_patient_source = c$paired_source, C_member_pair_source = c$member_source, 
        D_patient_ratios = d$patient_source, D_member_results = r[r$multiplicity_family == "PED_H_MEMBERS", ])
    paths <- vapply(names(tables), function(n) s13_tsv(tables[[n]], file.path(s13_r, paste0(n, ".tsv"))), "")
    s13_json(list(status = "COMPLETE_AUTHORIZED_ENDPOINTS", at_utc = s13_now(), outputs = s13_hashes(paths), hypotheses = 13, 
        pending_family = "PEC_COMPOSITION", C_composition_effect_access = FALSE), file.path(s13_m, "ANALYTICAL_OUTPUTS.json"))
    paths
}

  close_terminal_state <- function() {
    p <- file.path(s13_r,"all_hypothesis_results.tsv"); r <- s13_read(p)
    disposition <- "workflow/config/reconstruction_terminal_states.json"
    s01_require(context$config$public_terminal_states$TCGA_composition == "NA_UNAVAILABLE", "Accepted terminal composition state")
    j <- which(r$endpoint_id=="PE_C_TCGA_STAD_beta_E")
    s01_require(length(j)==1L&&r$state[j]=="PENDING","Accepted closeout predecessor")
    r$state[j] <- "NA_UNAVAILABLE"
    r$reason[j] <- paste0("NO_QUALIFIED_WHOLE_TISSUE_REFERENCE_WITHIN_APPROVED_TWO_CANDIDATE_SCOPE;OWNER_DISPOSITION=",disposition)
    s13_tsv(r,p)
    closure <- do.call(rbind,lapply(unique(r$multiplicity_family),function(f) {
      z<-r[r$multiplicity_family==f,];cl<-pe_close_family(setNames(z$p_value,z$hypothesis_id),z$state,z$hypothesis_id)
      data.frame(family=f,frozen_size=nrow(z),estimable=sum(z$state=="ESTIMABLE"),nonestimable=sum(z$state=="NA_NONESTIMABLE"),unavailable=sum(z$state=="NA_UNAVAILABLE"),pending=sum(z$state=="PENDING"),closed=cl$closed)
    }))
    s13_tsv(closure,file.path(s13_r,"family_closure.tsv"))
  }


  run <- function() {
    context$guard(); context$resource_check(); s13_runtime(); s13_preflight()
    C <- s13_C(context$manifest_path)
    D <- s13_D(context$manifest_path)
    s13_assemble(C,D); close_terminal_state()
    list(results=list.files(s13_r,full.names=TRUE),processed=list.files(s13_p,full.names=TRUE))
  }
  list(run=run,pairs=s13_pairs,ratio_test=s13_ratio_test,select_frozen=s13_select_frozen)
}
