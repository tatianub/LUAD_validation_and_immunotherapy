# ========================================
# Utility functions
# ========================================
normalize_sample_id <- function(x) {
  toupper(gsub(
    "[^A-Za-z0-9]", ".",
    as.character(x)
  ))
}
prepare_myc_corr_dt <- function(dt, value_cols, clinical_dt) {
  dt <- merge(
    dt,
    clinical_dt[, .(sample_id, response)],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  non_missing_idx <- !is.na(dt$response)
  for (col_name in value_cols) {
    non_missing_idx <- non_missing_idx & !is.na(dt[[col_name]])
  }
  dt <- dt[non_missing_idx]
  dt[, response := factor(
    response,
    levels = c("resistance", "response")
  )]
  dt
}
style_myc_corr_plot <- function(plot_obj) {
  plot_obj +
    geom_point(
      size = 2.5,
      alpha = 0.6
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      alpha = 0.15
    ) +
    scale_color_manual(values = RESPONSE_COLORS) +
    scale_fill_manual(values = RESPONSE_COLORS) +
    theme_bw() +
    theme(legend.position = "bottom")
}
safe_cor_test <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3) {
    return(list(
      estimate = NA_real_,
      p.value = NA_real_
    ))
  }
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(
      estimate = NA_real_,
      p.value = NA_real_
    ))
  }
  tryCatch(
    cor.test(
      x,
      y,
      method = method,
      exact = FALSE
    ),
    error = function(e) {
      list(
        estimate = NA_real_,
        p.value = NA_real_
      )
    }
  )
}
format_cor_label <- function(cor_test, prefix = "rho") {
  estimate <- as.numeric(cor_test$estimate)
  paste0(
    prefix,
    " = ",
    ifelse(
      is.na(estimate),
      "NA",
      format(estimate, digits = 3)
    ),
    ", p = ",
    ifelse(
      is.na(cor_test$p.value),
      "NA",
      format.pval(
        cor_test$p.value,
        digits = 2,
        eps = 1e-3
      )
    )
  )
}
group_correlations <- function(dt, x_col, y_col, method = "spearman") {
  result <- lapply(levels(dt$response), function(group_name) {
    group_dt <- dt[response == group_name]
    cor_test <- safe_cor_test(
      group_dt[[x_col]],
      group_dt[[y_col]],
      method = method
    )
    data.table(
      response = group_name,
      n = nrow(group_dt),
      estimate = as.numeric(cor_test$estimate),
      p_value = cor_test$p.value
    )
  })
  rbindlist(result)
}
format_group_correlations <- function(results, prefix = "rho") {
  paste(
    vapply(seq_len(nrow(results)), function(i) {
      group_label <- ifelse(
        as.character(results$response[i]) == "resistance",
        "Resistance",
        "Response"
      )
      paste0(
        group_label,
        ": ",
        prefix,
        " = ",
        ifelse(
          is.na(results$estimate[i]),
          "NA",
          format(results$estimate[i], digits = 3)
        ),
        ", p = ",
        ifelse(
          is.na(results$p_value[i]),
          "NA",
          format.pval(
            results$p_value[i],
            digits = 2,
            eps = 1e-3
          )
        ),
        ", n = ",
        results$n[i]
      )
    }, character(1)),
    collapse = "\n"
  )
}

fmt <- function(p) {
  ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-3))
}

get_cor_label <- function(v) {
  paste(
    vapply(seq_len(nrow(cors[variable == v])), function(i) {
      r <- cors[variable == v][i]
      paste0(
        ifelse(r$group == "resistance", "Non-responders", "Responders"),
        ": rho=", round(r$rho, 2),
        ", p=", fmt(r$p),
        ", n=", r$n
      )
    }, character(1)),
    collapse = "\n"
  )
}

cor_fun <- function(data, x) {
  rbindlist(lapply(levels(data$response), function(g) {
    d <- data[response == g]
    ok <- length(unique(d[[x]])) > 1 &&
      length(unique(d$myc_le_targeting)) > 1 &&
      nrow(d) >= 3
    test <- if (ok) {
      suppressWarnings(
        cor.test(d[[x]], d$myc_le_targeting, method = "spearman",
                 exact = FALSE)
      )
    } else {
      NULL
    }
    data.table(
      variable = x,
      group = g,
      n = nrow(d),
      rho = if (ok) unname(test$estimate) else NA_real_,
      p = if (ok) test$p.value else NA_real_
    )
  }))
}