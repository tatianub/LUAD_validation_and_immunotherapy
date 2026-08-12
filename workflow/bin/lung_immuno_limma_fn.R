subset_clinical_data <- function(clinical_data,
                                 subtype = "all",
                                 treatment_type = "all") {
  subtype_value <- tolower(subtype)
  treatment_value <- tolower(treatment_type)

  # Filter by Histological Subtype
  if (subtype_value == "adeno") {
    clinical_data <- clinical_data[
      clinical_data$Histology_Harmonized == "Adeno", 
    ]
  } else if (subtype_value != "all") {
    stop("SUBTYPE must be one of: all, Adeno")
  }

  # Filter by Treatment Type
  if (treatment_value == "pd1") {
    clinical_data <- clinical_data[
      clinical_data$Agent_PD1_Category == "PD(L)1", 
    ]
  } else if (treatment_value != "all") {
    stop("TREATMENT_TYPE must be one of: all, PD1")
  }
  return(clinical_data)
}


run_limma_analysis <- function(data, 
                               response_vector,
                               covariates = NULL) {
  # Input validation
  if (ncol(data) != length(response_vector)) {
    stop(paste0("Number of samples in data must ", 
                "match length of response_vector"))
  }
  # Validate response values
  valid_vals <- c("resistance", "response")
  if (!all(response_vector %in% valid_vals)) {
    stop('response_vector must contain only "resistance" and "response"')
  }
  if (length(unique(response_vector)) < 2) {
    stop(paste0('Both "resistance" and "response" ', 
                'must be present in response_vector'))
  }
  data_matrix <- as.matrix(data)
  gene_names <- rownames(data)
  sample_names <- colnames(data_matrix)
  # Create design data frame with fixed factor levels
  design_data <- data.frame(
    Sample = sample_names,
    Response = factor(response_vector, levels = c("resistance", "response")),
    stringsAsFactors = FALSE
  )
  # Add covariates if provided
  if (!is.null(covariates)) {
    design_data <- cbind(design_data, covariates)
    covariate_terms <- paste(colnames(covariates), collapse = " + ")
    formula_string <- paste("~ 0 + Response +", covariate_terms)
  } else {
    formula_string <- "~ 0 + Response"
  }
  # Create design matrix
  design_matrix <- model.matrix(
    as.formula(formula_string), 
    data = design_data
  )
  # Define single specific contrast
  contrast_strings <- c(
     Response_vs_resistance = "Responseresponse - Responseresistance"
   )
  # Create contrast matrix
  contrast_matrix <- makeContrasts(
    contrasts = contrast_strings,
    levels = design_matrix
  )
  # Run limma analysis
  fit <- lmFit(data_matrix, design_matrix)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit3 <- eBayes(fit2)
  # Extract results for the single contrast
  res <- topTable(fit3, coef = 1, sort.by = "none", n = Inf)
  res$contrast <- names(contrast_strings)[1]
  res$gene <- gene_names
  return(res)
}


filter_fgsea_results <- function(fgseaRes_all,
                                 padj_cutoff = 0.05,
                                 es_threshold = 0.5) {
  if (is.null(fgseaRes_all) || nrow(fgseaRes_all) == 0) {
    return(NULL)
  }
  # Subset based on significance and effect size
  filtered_results <- subset(
    fgseaRes_all,
    padj <= padj_cutoff & abs(ES) >= es_threshold
  )
  # Calculate transformed significance
  filtered_results$log10padj <- -log10(filtered_results$padj)
  return(filtered_results)
}


plot_fgsea_results <- function(res_fgsea_filt, 
                               figure_file_path,
                               max_char = 45,
                               width = 10,
                               height = 8,
                               base_size = 11,
                               low_col = "#2c7bb6",
                               mid_col = "#f7f7f7",
                               high_col = "#d7191c",
                               midpoint = 0,
                               alpha = 0.9,
                               plot_title = "FGSEA Analysis Results",
                               x_lab = "-log10(adjusted p-value)",
                               y_lab = "Pathway") {
  if (is.null(res_fgsea_filt) || nrow(res_fgsea_filt) == 0) {
    message("No significant pathways to plot.")
    return(NULL)
  }
  # Clean and truncate pathway names
  res_fgsea_filt$pathway_clean <- gsub("^REACTOME_", "", res_fgsea_filt$pathway)
  res_fgsea_filt$pathway_clean <- gsub("_", " ", res_fgsea_filt$pathway_clean)
  res_fgsea_filt$pathway_clean <- substr(
    res_fgsea_filt$pathway_clean, 1, max_char
  )
  # Reorder factor levels based on significance (log10padj)
  path_order <- order(res_fgsea_filt$log10padj)
  res_fgsea_filt$pathway_clean <- factor(
    res_fgsea_filt$pathway_clean,
    levels = rev(res_fgsea_filt$pathway_clean[path_order])
  )

  p_bubble <- ggplot(res_fgsea_filt, aes(x = log10padj, y = pathway_clean)) +
    geom_point(aes(size = abs(ES), color = ES), alpha = alpha) +
    scale_color_gradient2(
      low = low_col,
      mid = mid_col,
      high = high_col,
      midpoint = midpoint,
      name = "ES"
    ) +
    scale_size_continuous(name = "|ES|") +
    labs(
      x = x_lab,
      y = y_lab,
      title = plot_title
    ) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(figure_file_path, plot = p_bubble, width = width, height = height)
}

run_fgsea_analysis_TFs <- function(
            res_net,
            gmt_file,
            num_cores = 20,
            min_size = 5,
            max_size = 500
          ) {
    tf_names <- unique(res_net$tf)
    if (length(tf_names) == 0) return(NULL)

    requested_cores <- suppressWarnings(as.integer(num_cores))
    if (is.na(requested_cores) || requested_cores < 1) {
      stop("num_cores must be a positive integer")
    }

    workers <- min(requested_cores, length(tf_names))

    # Constrain internal BLAS/OpenMP/data.table 
    # threading to avoid nested parallelism.
    thread_vars <- c(
      "OMP_NUM_THREADS",
      "MKL_NUM_THREADS",
      "OPENBLAS_NUM_THREADS",
      "BLIS_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS",
      "NUMEXPR_NUM_THREADS",
      "R_DATATABLE_NUM_THREADS"
    )
    old_thread_env <- Sys.getenv(thread_vars, unset = "")
    on.exit(
      do.call(
        Sys.setenv,
        as.list(setNames(old_thread_env, thread_vars))
      ),
      add = TRUE
    )
    do.call(
      Sys.setenv,
      as.list(setNames(rep("1", length(thread_vars)), thread_vars))
    )
    data.table::setDTthreads(1)

    message(
      "Running fGSEA with ", workers,
      " independent R worker(s); each worker forced to 1 internal thread"
    )

    results <- parallel::mclapply(
      tf_names,
      function(tf_name) {
        do.call(
          Sys.setenv,
          as.list(setNames(rep("1", length(thread_vars)), thread_vars))
        )
        data.table::setDTthreads(1)

        res_tf <- res_net[res_net$tf == tf_name, ]
        if (nrow(res_tf) < 10) return(NULL)

        tryCatch({
          result <- run_fgsea(
            res_tf,
            gmt_file,
            min_size = min_size,
            max_size = max_size,
            nproc = 1
          )

          if (!is.null(result) && nrow(result) > 0) {
            result$TF <- tf_name
            return(result)
          }

          NULL
        }, error = function(e) {
          message("Error processing TF ", tf_name, ": ", e$message)
          NULL
        })
      },
      mc.cores = workers
    )

    results <- Filter(Negate(is.null), results)

    if (length(results) == 0) return(NULL)

    data.table::rbindlist(results, fill = TRUE)
  }


# Function to prepare network data for analysis
prepare_network_data <- function(net, edges, response_data) {
    net_clean <- net[, colnames(net) %in% response_data$sample_id, 
                     with = FALSE]
    links <- paste(edges$reg, edges$tar, sep = "_")
    rownames(net_clean) <- links
    return(net_clean)
}

run_network_analysis <- function(net_clean, 
                                 response_vector, 
                                 edges, 
                                 covariates) {
  # Run the differential expression analysis using limma
  res_net <- run_limma_analysis(
    net_clean, 
    response_vector, 
    covariates = covariates
  )
  # Validation: Ensure the number of results matches the edges provided
  if (nrow(res_net) != nrow(edges)) {
    stop(paste0("Dimension mismatch: res_net has ", nrow(res_net), 
                " rows, but edges has ", nrow(edges), " rows."))
  }

  # Map TF (regulator) and Target (gene) from the edges file to the results
  res_net$tf <- edges$reg
  res_net$gene <- edges$tar
  return(res_net)
}


run_fgsea <- function(res_all,
                      gmt_file,
                      min_size = 0,
                      max_size = 1000,
                      nproc = 1) {
  pt <- fgsea::gmtPathways(gmt_file)
  comparisons <- unique(res_all$contrast)

  fgsea_results_list <- lapply(comparisons, function(comparison) {
    # Filter data for the specific contrast
    data <- dplyr::filter(
      res_all, 
      grepl(comparison, contrast, fixed = TRUE)
    )

    # Prepare named vector of ranks
    ranks <- setNames(data$t, data$gene)
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)

    # Run multilevel fGSEA
    fgseaRes <- fgsea::fgseaMultilevel(
      pathways = pt,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      nproc = nproc
    )

    if (nrow(fgseaRes) == 0) return(NULL)

    # Select core columns and add comparison metadata
    fgseaRes <- fgseaRes[, 1:7]
    fgseaRes$cmp <- comparison

    return(fgseaRes)
  })

  # Combine results into a single data.table
  fgseaRes_all <- data.table::rbindlist(
    Filter(Negate(is.null), fgsea_results_list),
    fill = TRUE
  )

  if (nrow(fgseaRes_all) == 0) return(NULL)

  return(fgseaRes_all)
}

save_fgsea_results <- function(results, output_file) {
  if (!is.null(results) && nrow(results) > 0) {
    # Write the results table to a tab-separated file
    write.table(
      results, 
      file = output_file, 
      sep = "\t", 
      quote = FALSE, 
      row.names = FALSE
    )
    
    # Log successful save and metadata
    cat("Results saved to:", output_file, "\n")
    cat("Number of results:", nrow(results), "\n")
  } else {
    # Log that no data was available for output
    cat("No results to save for:", output_file, "\n")
  }
}

plot_tf_lollipop <- function(
            tfs_lollipop,
            vline_xintercept = 0,
            vline_linewidth = 0.4,
            vline_color = "grey45",
            segment_x_start = 0,
            segment_linewidth = 0.7,
            segment_alpha = 0.75,
            point_color = "grey20",
            point_stroke = 0.45,
            point_alpha = 0.95,
            responder_label = "Enriched in responders",
            non_responder_label = "Enriched in non-responders",
            non_responder_color = "#2166AC",
            responder_color = "#B2182B",

            other_tf_shape = 21,
            direct_regulator_shape = 24,
            size_min = 2.2,
            size_max = 6.8,
            size_legend_name = expression(-log[10]("adjusted P-value")),
            base_size = 12,
            legend_position = "right",
            axis_text_y_size = 8,
            title_size = 13,
            margin_top = 10,
            margin_right = 20,
            margin_bottom = 10,
            margin_left = 10,
            plot_title = "",
            x_axis_label = "Enrichment score",
            color_legend_title = NULL,
            fill_legend_title = NULL,
            shape_legend_title = "PD-L1 evidence"
            ,
            top_n_responder = 10,
            top_n_non_responder = 10,
            rank_by = c("padj", "abs_es")
            ) {
  rank_by <- match.arg(rank_by)

  if (!"neglog10_padj" %in% names(tfs_lollipop) && "padj" %in% names(tfs_lollipop)) {
    tfs_lollipop$neglog10_padj <- -log10(tfs_lollipop$padj)
  }

  tfs_lollipop$direction <- ifelse(
    tfs_lollipop$ES > 0,
    responder_label,
    non_responder_label
  )

  tfs_lollipop$responder_rank <- if (rank_by == "padj") {
    tfs_lollipop$neglog10_padj
  } else {
    abs(tfs_lollipop$ES)
  }

  tfs_responder <- 
      tfs_lollipop[tfs_lollipop$direction == responder_label, ]
  tfs_non_responder <- 
      tfs_lollipop[tfs_lollipop$direction == non_responder_label, ]

  tfs_responder <- dplyr::slice_max(
    tfs_responder,
    order_by = .data$responder_rank,
    n = top_n_responder,
    with_ties = FALSE
  )
  tfs_non_responder <- dplyr::slice_max(
    tfs_non_responder,
    order_by = .data$responder_rank,
    n = top_n_non_responder,
    with_ties = FALSE
  )

  tfs_lollipop <- dplyr::bind_rows(tfs_responder, tfs_non_responder)
  tfs_lollipop <- tfs_lollipop[order(
    factor(tfs_lollipop$direction, 
      levels = c(responder_label, non_responder_label)),
    -tfs_lollipop$responder_rank,
    -abs(tfs_lollipop$ES)
  ), ]
  tfs_lollipop$TF <- factor(tfs_lollipop$TF, 
    levels = rev(unique(tfs_lollipop$TF)))

  ggplot2::ggplot(
    tfs_lollipop,
    ggplot2::aes(x = .data$ES, y = .data$TF)
  ) +
    ggplot2::geom_vline(
      xintercept = vline_xintercept,
      linewidth = vline_linewidth,
      color = vline_color
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = segment_x_start,
        xend = .data$ES,
        y = .data$TF,
        yend = .data$TF,
        color = .data$direction
      ),
      linewidth = segment_linewidth,
      alpha = segment_alpha
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$neglog10_padj,
        fill = .data$direction,
        shape = .data$direct_regulator
      ),
      color = point_color,
      stroke = point_stroke,
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(
        c(responder_color, non_responder_color),
        c(responder_label, non_responder_label)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        c(responder_color, non_responder_color),
        c(responder_label, non_responder_label)
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Other TF" = other_tf_shape,
        "Direct PD-L1 regulator" = direct_regulator_shape
      )
    ) +
    ggplot2::scale_size_continuous(
      range = c(size_min, size_max),
      name = size_legend_name
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.text.y = ggplot2::element_text(size = axis_text_y_size),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = title_size),
      plot.margin = ggplot2::margin(
        margin_top,
        margin_right,
        margin_bottom,
        margin_left
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = x_axis_label,
      color = color_legend_title,
      fill = fill_legend_title,
      shape = shape_legend_title
    )
}


plot_fgsea_results_pathway <- function(
        res_fgsea_filt,
        empty_message = "No significant pathways to plot.",
        pathway_prefix_pattern = "^REACTOME_",
        pathway_separator_pattern = "_",
        pathway_separator_replacement = " ",
        pathway_max_chars = 45,
        point_alpha = 0.9,
        highlight_pathway_label = "PD-1 signaling",
        highlight_pathway_aliases = c("PD1 signaling", "PD 1 signaling", "PD 1 SIGNALING"),
        highlight_label_size = 3.6,
        highlight_label_face = "bold",
        highlight_label_vjust = 0,
        highlight_label_hjust = 0,
        highlight_label_fill = "white",
        highlight_label_color = "black",
        highlight_label_alpha = 0.95,
        highlight_outline_color = "black",
        highlight_outline_size = 5.6,
        highlight_outline_stroke = 1.5,
        color_low = "#2166AC",
        color_mid = "#f7f7f7",
        color_high = "#B2182B",
        color_midpoint = 0,
        color_legend_name = "ES",
        size_legend_name = "|ES|",
        x_axis_label = "-log10(adjusted p-value)",
        y_axis_label = "Pathway",
        plot_title = "",
        base_size = 11,
        major_y_grid_color = "grey90"
) {
  if (is.null(res_fgsea_filt) || nrow(res_fgsea_filt) == 0) {
    print(empty_message)
    return(NULL)
  }

  res_fgsea_filt$pathway_clean <- gsub(
    pathway_prefix_pattern,
    "",
    res_fgsea_filt$pathway
  )
  res_fgsea_filt$pathway_clean <- gsub(
    pathway_separator_pattern,
    pathway_separator_replacement,
    res_fgsea_filt$pathway_clean
  )
  res_fgsea_filt$pathway_clean <- substr(
    res_fgsea_filt$pathway_clean,
    1,
    pathway_max_chars
  )
  res_fgsea_filt$pathway_clean <- factor(
    res_fgsea_filt$pathway_clean,
    levels = res_fgsea_filt$pathway_clean[
      order(res_fgsea_filt$log10padj, 
          abs(res_fgsea_filt$ES), decreasing = FALSE)
    ]
  )

  p_bubble <- ggplot2::ggplot(
    res_fgsea_filt,
    ggplot2::aes(x = .data$log10padj, y = .data$pathway_clean)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = abs(.data$ES), color = .data$ES),
      alpha = point_alpha
    ) +
    ggplot2::scale_color_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = color_midpoint,
      name = color_legend_name
    ) +
    ggplot2::scale_size_continuous(name = size_legend_name) +
    ggplot2::labs(
      x = x_axis_label,
      y = y_axis_label,
      title = plot_title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        color = major_y_grid_color
      ),
      panel.grid.minor = ggplot2::element_blank()
    )

  highlight_pattern <- paste(
    c(
      highlight_pathway_label,
      highlight_pathway_aliases,
      "PD[- ]?1.*signaling",
      "PD1.*signaling"
    ),
    collapse = "|"
  )
  highlight_df <- res_fgsea_filt[
    grepl(highlight_pattern, res_fgsea_filt$pathway_clean, ignore.case = TRUE),
  ]

  if (nrow(highlight_df) > 0) {
    p_bubble <- p_bubble + ggplot2::geom_point(
      data = highlight_df,
      ggplot2::aes(x = .data$log10padj, y = .data$pathway_clean),
      inherit.aes = FALSE,
      shape = 21,
      fill = NA,
      color = highlight_outline_color,
      size = highlight_outline_size,
      stroke = highlight_outline_stroke
    )
  }

  p_bubble
}


plot_pathway_lollipop_ab <- function(
  res_fgsea_filt,
  tfs_lollipop,
  output_file = NULL,
  width = 10,
  height = 12,
  panel_gap = grid::unit(4, "mm"),
  panel_layout = c("stacked", "row"),
  pathway_row_height = 1,
  lollipop_row_height = 1,
  pathway_panel_width = 0.95,
  lollipop_panel_width = 0.95,
  tf_top_n_responder = 10,
  tf_top_n_non_responder = 10,
  tf_rank_by = c("padj", "abs_es"),
  label_a = "A.",
  label_b = "B.",
  label_x = 0.01,
  label_y = 0.99,
  label_cex = 1.4
) {
  panel_layout <- match.arg(panel_layout)
  tf_rank_by <- match.arg(tf_rank_by)

  p_pathway <- plot_fgsea_results_pathway(
    res_fgsea_filt,
    highlight_pathway_label = "PD-1 signaling"
  )
  p_lollipop <- plot_tf_lollipop(
    tfs_lollipop,
    responder_label = "Enriched in responders",
    non_responder_label = "Enriched in non-responders",
    top_n_responder = tf_top_n_responder,
    top_n_non_responder = tf_top_n_non_responder,
    rank_by = tf_rank_by
  )

  if (is.null(p_pathway) || is.null(p_lollipop)) {
    return(NULL)
  }

  g_pathway <- ggplot2::ggplotGrob(p_pathway)
  g_lollipop <- ggplot2::ggplotGrob(p_lollipop)

  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  grid::grid.newpage()
  if (panel_layout == "row") {
    lay <- grid::grid.layout(
      nrow = 1,
      ncol = 3,
      widths = grid::unit.c(
        grid::unit(pathway_row_height, "null"),
        panel_gap,
        grid::unit(lollipop_row_height, "null")
      )
    )
  } else {
    lay <- grid::grid.layout(
      nrow = 3,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(pathway_row_height, "null"),
        panel_gap,
        grid::unit(lollipop_row_height, "null")
      )
    )
  }
  grid::pushViewport(grid::viewport(layout = lay))

  pathway_layout_pos <- if (panel_layout == "row") {
    list(row = 1, col = 1)
  } else {
    list(row = 1, col = 1)
  }
  lollipop_layout_pos <- if (panel_layout == "row") {
    list(row = 1, col = 3)
  } else {
    list(row = 3, col = 1)
  }

  grid::pushViewport(
    grid::viewport(
      layout.pos.row = pathway_layout_pos$row,
      layout.pos.col = pathway_layout_pos$col
    )
  )
  grid::pushViewport(
    grid::viewport(
      x = 0.5,
      y = 0.5,
      width = grid::unit(pathway_panel_width, "npc"),
      height = grid::unit(1, "npc"),
      just = c("center", "center")
    )
  )
  grid::grid.draw(g_pathway)
  grid::grid.text(
    label_a,
    x = grid::unit(label_x, "npc"),
    y = grid::unit(label_y, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", cex = label_cex)
  )
  grid::upViewport()
  grid::upViewport()

  grid::pushViewport(
    grid::viewport(
      layout.pos.row = lollipop_layout_pos$row,
      layout.pos.col = lollipop_layout_pos$col
    )
  )
  grid::pushViewport(
    grid::viewport(
      x = 0.5,
      y = 0.5,
      width = grid::unit(lollipop_panel_width, "npc"),
      height = grid::unit(1, "npc"),
      just = c("center", "center")
    )
  )
  grid::grid.draw(g_lollipop)
  grid::grid.text(
    label_b,
    x = grid::unit(label_x, "npc"),
    y = grid::unit(label_y, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontface = "bold", cex = label_cex)
  )
  grid::upViewport()
  grid::upViewport()
  grid::upViewport()

  invisible(
    list(
      pathway_plot = p_pathway,
      lollipop_plot = p_lollipop
    )
  )
}

