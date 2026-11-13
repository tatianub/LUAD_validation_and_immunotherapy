find_survival_pairs <- function(clin_df) {
  cn <- colnames(clin_df)

  # Find all event columns - including _event, _event_5y, etc.
  event_cols <- cn[grepl("_event(_.*)?$", cn, ignore.case = TRUE)]

  if (length(event_cols) == 0) {
    stop("No event columns found in clinical data. 
          Expected columns like OS_event, RFS_event, OS_event_5y, etc.")
  }

  pairs <- list()

  for (event_col in event_cols) {
    # Extract the full pattern including suffixes like _5y
    # For "OS_event_5y" -> base="OS", suffix="_5y"
    # For "OS_event" -> base="OS", suffix=""

    if (grepl("_event_", event_col, ignore.case = TRUE)) {
      # Has suffix: OS_event_5y -> OS and _5y
      base_pattern <- sub("_event_.*$", "", event_col, ignore.case = TRUE)
      suffix_pattern <- sub("^.*_event", "", event_col, ignore.case = TRUE)
      endpoint_name <- paste0(base_pattern, suffix_pattern)
    } else {
      # No suffix: OS_event -> OS
      base_pattern <- sub("_event$", "", event_col, ignore.case = TRUE)
      suffix_pattern <- ""
      endpoint_name <- base_pattern
    }

    # Look for matching time column with same base and suffix
    candidate_time_cols <- c(
      paste0(base_pattern, "_time", suffix_pattern),
      paste0(base_pattern, "_Time", suffix_pattern),
      paste0(base_pattern, "_TIME", suffix_pattern),
      paste0(base_pattern, "_days", suffix_pattern),
      paste0(base_pattern, "_Days", suffix_pattern),
      paste0(base_pattern, "_DAYS", suffix_pattern),
      paste0("time_", base_pattern, suffix_pattern),
      paste0("Time_", base_pattern, suffix_pattern),
      paste0("TIME_", base_pattern, suffix_pattern),
      paste0("ToF_", base_pattern, suffix_pattern),
      paste0("tof_", base_pattern, suffix_pattern)
    )

    matched_time <- candidate_time_cols[candidate_time_cols %in% cn]

    # Fallback: look for any time column with same base and suffix pattern
    if (length(matched_time) == 0) {
      fallback <- cn[
        grepl(paste0("^", base_pattern), cn, ignore.case = TRUE) &
          grepl("time|days|month|surv|tof", cn, ignore.case = TRUE) &
          grepl(paste0(suffix_pattern, "$"), cn, ignore.case = TRUE) &
          !grepl("event", cn, ignore.case = TRUE)
      ]
      matched_time <- fallback
    }

    if (length(matched_time) > 0) {
      pairs[[length(pairs) + 1]] <- data.frame(
        endpoint = endpoint_name,
        event_col = event_col,
        time_col = matched_time[1],
        stringsAsFactors = FALSE
      )
      cat("Found survival pair:", endpoint_name, "->",
          matched_time[1], "/", event_col, "\n")
    } else {
      cat("Warning: no matching time column found for event column:",
          event_col, "\n")
    }
  }

  if (length(pairs) == 0) {
    stop("Could not identify any valid time/event survival column pairs.")
  }

  do.call(rbind, pairs)
}



safe_numeric_event <- function(x) {
  # handles numeric, logical, factor, character
  if (is.logical(x)) return(as.numeric(x))
  if (is.numeric(x)) return(as.numeric(x))
  if (is.factor(x)) x <- as.character(x)

  x <- trimws(tolower(as.character(x)))
  x[x %in% c("dead", "deceased", "yes", "true")] <- "1"
  x[x %in% c("alive", "no", "false")] <- "0"
  as.numeric(x)
}


extract_cox_results <- function(model, model_name, endpoint_name) {
  s <- summary(model)
  coeff_table <- s$coefficients
  conf_int_table <- s$conf.int

  out <- vector("list", nrow(coeff_table))
  for (i in seq_len(nrow(coeff_table))) {
    out[[i]] <- data.frame(
      Endpoint = endpoint_name,
      Model = model_name,
      Variable = rownames(coeff_table)[i],
      HR = coeff_table[i, "exp(coef)"],
      CI_Lower = conf_int_table[i, "lower .95"],
      CI_Upper = conf_int_table[i, "upper .95"],
      P_Value = coeff_table[i, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# Function to run CoxPH for a single gene with specified outcome and
# covariates
run_cox_ranking <- function(
  gene_name,
  exp_data,
  clin_data,
  time_col,
  event_col,
  covariates = NULL
) {
  if (!(gene_name %in% rownames(exp_data))) {
    return(list(status = "skip", reason = "gene_not_found"))
  }

  valid_samples <- intersect(clin_data$sample_id, colnames(exp_data))
  if (length(valid_samples) < 5) {
    return(list(status = "skip", reason = "insufficient_samples"))
  }

  clin_use <- clin_data %>%
    filter(.data$sample_id %in% valid_samples) %>%
    select(
      sample_id,
      all_of(time_col),
      all_of(event_col),
      any_of(covariates)
    ) %>%
    rename(
      time_var = all_of(time_col),
      event_var = all_of(event_col)
    )

  gene_exp <- as.numeric(exp_data[gene_name, clin_use$sample_id, drop = TRUE])
  df <- clin_use %>%
    mutate(gene_exp = gene_exp) %>%
    filter(
      !is.na(gene_exp),
      is.finite(gene_exp),
      !is.na(time_var),
      !is.na(event_var)
    )

  n_events <- sum(df$event_var == 1, na.rm = TRUE)
  if (nrow(df) < 5 || n_events < 2) {
    return(
      list(status = "skip", reason = "insufficient_samples_or_events")
    )
  }

  # Build formula: Surv(time, event) ~ gene_exp [+ covariates]
  formula_str <- "survival::Surv(time_var, event_var) ~ gene_exp"
  if (length(covariates) > 0) {
    # Filter covariates to those actually present in df
    valid_covs <- covariates[covariates %in% colnames(df)]
    if (length(valid_covs) > 0) {
      formula_str <- paste(
        formula_str,
        "+",
        paste(valid_covs, collapse = " + ")
      )
    }
  }

  # Fit Cox Proportional Hazards Model
  fit <- tryCatch(
    eval(
      parse(text = paste("survival::coxph(", formula_str, ", data = df)"))
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(status = "error", reason = conditionMessage(fit)))
  }

  sfit <- summary(fit)
  if (is.null(sfit$coefficients) || nrow(sfit$coefficients) == 0) {
    return(list(status = "skip", reason = "no_coefficients"))
  }

  # Extract first row (gene coefficient, not covariates)
  coef_row <- sfit$coefficients[1, , drop = FALSE]

  return(list(
    status = "ok",
    result = data.frame(
      gene = gene_name,
      z_score = as.numeric(coef_row[1, "z"]),
      p_value = as.numeric(coef_row[1, "Pr(>|z|)"]),
      hazard_ratio = as.numeric(exp(coef_row[1, "coef"]))
    )
  ))
}


  run_one_gene <- function(g) {
    out <- run_cox_ranking(
      g,
      expression_matrix,
      clinical_data,
      time_col,
      event_col,
      all_covariates
    )
    if (!is.list(out) || is.null(out$status)) {
      return(data.frame(
        outcome = outcome_name,
        gene = g,
        z_score = NA_real_,
        p_value = NA_real_,
        hazard_ratio = NA_real_,
        status = "error",
        reason = "invalid_return"
      ))
    }
    if (out$status == "ok") {
      out$result$outcome <- outcome_name
      out$result$status <- "ok"
      out$result$reason <- ""
      return(out$result)
    }
    data.frame(
      outcome = outcome_name,
      gene = g,
      z_score = NA_real_,
      p_value = NA_real_,
      hazard_ratio = NA_real_,
      status = out$status,
      reason = out$reason
    )
  }