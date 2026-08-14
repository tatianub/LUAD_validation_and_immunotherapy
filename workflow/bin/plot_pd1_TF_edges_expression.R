#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "optparse",
  "ggpubr",
  "grid",
  "tools",
  "fgsea",
  "readxl"
)

for (lib in required_libraries) {
  suppressPackageStartupMessages(
    library(lib, character.only = TRUE, quietly = TRUE)
  )
}

####################
## Read arguments ##
####################

option_list <- list(
  optparse::make_option(
    "--clinical_file",
    type = "character",
    default = NULL,
    help = "Path to clinical file."
  ),
  optparse::make_option(
    "--gmt_file",
    type = "character",
    default = NULL,
    help = "Path to pathway GMT file (.gmt)."
  ),
  optparse::make_option(
    "--network_file",
    type = "character",
    default = NULL,
    help = "Path to inferred network file (.RData)."
  ),
  optparse::make_option(
    "--edges_file",
    type = "character",
    default = NULL,
    help = "Path to network edges file."
  ),
  optparse::make_option(
    "--expression_file",
    type = "character",
    default = NULL,
    help = "Path to expression matrix file."
  ),
  optparse::make_option(
    "--samples_file",
    type = "character",
    default = NULL,
    help = "Path to samples file matching expression columns."
  ),
  optparse::make_option(
    "--histo_subtype",
    type = "character",
    default = NULL,
    help = "Histological subtype filter."
  ),
  optparse::make_option(
    "--treatment_type",
    type = "character",
    default = NULL,
    help = "Treatment type filter."
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = NULL,
    help = "Random seed for reproducibility."
  ),
  optparse::make_option(
    "--tf_target_edges_pd1_boxplot",
    type = "character",
    default = NULL,
    help = "Output file for TF target-edge boxplots (PD-1 pathway)."
  ),
  optparse::make_option(
    "--tf_expression_boxplot",
    type = "character",
    default = NULL,
    help = "Output file for TF expression boxplots."
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Avoid set.seed(NULL)
if (!is.null(opt$seed)) {
  set.seed(opt$seed)
}

GMT_FILE <- opt$gmt_file
NETWORK_FILE <- opt$network_file
EDGES_FILE <- opt$edges_file
CLINICAL_FILE <- opt$clinical_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
HISTO_SUBTYPE <- opt$histo_subtype
TREATMENT_TYPE <- opt$treatment_type
TF_TARGET_EDGES_PD1_BOXPLOT <- opt$tf_target_edges_pd1_boxplot
TF_EXPRESSION_BOXPLOT <- opt$tf_expression_boxplot

source("workflow/bin/lung_immuno_limma_fn.R")


# ========================================
# Step 1: Load and filter clinical metadata
# ========================================

clinical_data <- fread(CLINICAL_FILE)

# Keep samples with response information
clinical_data <- clinical_data[!is.na(response)]

cat(
  "Filtering clinical data for histological subtype:",
  HISTO_SUBTYPE,
  "and treatment type:",
  TREATMENT_TYPE,
  "\n"
)

clinical_data_filt <- subset_clinical_data(
  clinical_data = clinical_data,
  subtype = HISTO_SUBTYPE,
  treatment_type = TREATMENT_TYPE
)

cat(
  "Number of samples after clinical filtering:",
  nrow(clinical_data_filt),
  "\n"
)

if (nrow(clinical_data_filt) == 0) {
  stop("No samples remain after clinical filtering.")
}


# ========================================
# Step 2: Load inferred regulatory network
# ========================================

net_env <- new.env()
load(NETWORK_FILE, envir = net_env)

obj_names <- ls(net_env)

cat(
  "Available objects in network file:",
  paste(obj_names, collapse = ", "),
  "\n"
)

if ("net" %in% obj_names) {

  net <- net_env$net
  cat("Loaded 'net' object\n")

} else if ("net_norm" %in% obj_names) {

  net <- net_env$net_norm
  cat("Loaded 'net_norm' object\n")

} else {

  stop(
    "No recognized network object found. Available objects: ",
    paste(obj_names, collapse = ", ")
  )
}

net <- as.matrix(net)
colnames(net) <- gsub("-", "\\.", colnames(net))

# ========================================
# Step 3: Load edge metadata and align
#         network with filtered samples
# ========================================

edges <- fread(EDGES_FILE)

# The rows of the edge table must correspond to the rows of net
if (nrow(edges) != nrow(net)) {
  stop(
    "Number of rows in edges (", nrow(edges),
    ") does not match number of network rows (", nrow(net), ")."
  )
}

# IMPORTANT:
# Use clinical_data_filt here, NOT clinical_data
network_samples <- intersect(
  colnames(net),
  clinical_data_filt$sample_id
)

if (length(network_samples) == 0) {
  stop("No overlapping samples between network and filtered clinical data.")
}

# Preserve network column order
net_clean <- net[
  ,
  colnames(net) %in% network_samples,
  drop = FALSE
]

# Match clinical metadata to network columns
clinical_data_ordered <- clinical_data_filt[
  match(colnames(net_clean), sample_id)
]

stopifnot(
  all(colnames(net_clean) == clinical_data_ordered$sample_id)
)

cat(
  "Samples retained in network analysis:",
  ncol(net_clean),
  "\n"
)


# ========================================
# Step 4: Restrict network to PD-1 pathway
# ========================================

pathways <- fgsea::gmtPathways(GMT_FILE)

pd1_pathway_names <- grep(
  "PD_1",
  names(pathways),
  value = TRUE
)

if (length(pd1_pathway_names) == 0) {
  stop("No pathway containing 'PD_1' found in GMT file.")
}

cat(
  "PD-1 pathway(s):",
  paste(pd1_pathway_names, collapse = ", "),
  "\n"
)

pd1_pathways <- pathways[pd1_pathway_names]

pathways_genes <- unique(
  unlist(pd1_pathways, use.names = FALSE)
)

pd1_edge_idx <- edges$tar %in% pathways_genes

net_pathways <- net_clean[
  pd1_edge_idx,
  ,
  drop = FALSE
]

edges_pathways <- edges[pd1_edge_idx]

stopifnot(
  nrow(net_pathways) == nrow(edges_pathways)
)

cat(
  "Number of PD-1 pathway edges:",
  nrow(edges_pathways),
  "\n"
)


# ========================================
# Selected TFs / genes
# ========================================

selected_tfs <- c(
  "ATF6", "FOSL2", "CREB3L2", "MYC", "MYCN", "TFEC",
  "HEY1", "ZKSCAN3", "MAX", "HEY2", "ZNF354A", "ZNF384",
  "FOXP2", "FOXC2", "MEF2A", "STAT2", "MEF2D", "FOXO1",
  "SOX15", "POU3F1"
)

# CD274 is not a TF, but you may want its expression shown
selected_expression_genes <- c(
  selected_tfs,
  "CD274"
)


# ========================================
# Step 5: Load expression matrix
# ========================================

expression_dt <- fread(EXPRESSION_FILE)

# First column contains gene IDs
gene_ids <- expression_dt[[1]]

# Convert remaining columns explicitly to numeric matrix
expression <- data.matrix(expression_dt[, -1])

samples <- fread(SAMPLES_FILE)

if (ncol(expression) != nrow(samples)) {
  stop(
    "Expression matrix has ",
    ncol(expression),
    " sample columns, but samples file contains ",
    nrow(samples),
    " rows."
  )
}

colnames(expression) <- samples$V2
rownames(expression) <- gene_ids
colnames(expression) <- gsub("-", "\\.", colnames(expression))

# ========================================
# Align expression and clinical samples
# ========================================

expression_samples <- intersect(
  clinical_data_ordered$sample_id,
  colnames(expression)
)

if (length(expression_samples) == 0) {
  stop("No overlapping samples between expression and clinical data.")
}

# Keep the clinical/network order
expression_samples <- clinical_data_ordered$sample_id[
  clinical_data_ordered$sample_id %in% expression_samples
]

expression_clean <- expression[
  ,
  expression_samples,
  drop = FALSE
]

clinical_expression_ordered <- clinical_data_ordered[
  match(
    colnames(expression_clean),
    sample_id
  )
]

stopifnot(
  all(
    colnames(expression_clean) ==
      clinical_expression_ordered$sample_id
  )
)

cat(
  "Samples retained for expression analysis:",
  ncol(expression_clean),
  "\n"
)


# ==========================================================
# PLOTTING: TF target-edge boxplots
# ==========================================================

tf_plot_list <- list()

for (selected_tf in selected_tfs) {

  tf_idx <- edges_pathways$reg == selected_tf

  if (!any(tf_idx)) {
    message("No PD-1 pathway edges found for TF: ", selected_tf)
    next
  }

  tf_edges <- edges_pathways[tf_idx]

  tf_edge_matrix <- net_pathways[
    tf_idx,
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(tf_edge_matrix) == nrow(tf_edges)
  )

  tf_edge_dt <- data.table::as.data.table(
    tf_edge_matrix
  )

  tf_edge_dt[, target := tf_edges$tar]
  tf_edge_dt[, edge_id := paste(
    tf_edges$reg,
    tf_edges$tar,
    sep = "_"
  )]

  tf_edge_long <- data.table::melt(
    tf_edge_dt,
    id.vars = c("edge_id", "target"),
    variable.name = "sample_id",
    value.name = "edge_weight"
  )

  tf_edge_long <- merge(
    tf_edge_long,
    clinical_data_ordered[
      ,
      .(sample_id, response)
    ],
    by = "sample_id",
    all.x = TRUE
  )

  tf_edge_long <- tf_edge_long[
    !is.na(response)
  ]

  tf_edge_long[
    ,
    response := factor(
      response,
      levels = c(
        "resistance",
        "response"
      )
    )
  ]

  # Make sure both response groups are present
  if (length(unique(na.omit(tf_edge_long$response))) < 2) {
    message(
      "Only one response group available for TF: ",
      selected_tf
    )
    next
  }

  p_value <- tryCatch(
    wilcox.test(
      edge_weight ~ response,
      data = tf_edge_long
    )$p.value,
    error = function(e) NA_real_
  )

  med_resistance <- median(
    tf_edge_long[
      response == "resistance",
      edge_weight
    ],
    na.rm = TRUE
  )

  med_response <- median(
    tf_edge_long[
      response == "response",
      edge_weight
    ],
    na.rm = TRUE
  )

  direction <- data.table::fifelse(
    med_response > med_resistance,
    "higher in response",
    data.table::fifelse(
      med_response < med_resistance,
      "higher in resistance",
      "same median"
    )
  )

  y_range <- range(
    tf_edge_long$edge_weight,
    na.rm = TRUE
  )

  # More robust than max * 1.05, especially for negative values
  y_pos <- y_range[2] +
    0.08 * diff(y_range)

  if (!is.finite(y_pos)) {
    y_pos <- y_range[2]
  }

  tf_plot_summary <- data.table(
    p_value = p_value,
    med_resistance = med_resistance,
    med_response = med_response,
    direction = direction,
    y_pos = y_pos
  )

  tf_plot_summary[
    ,
    label := paste0(
      "p = ",
      ifelse(
        is.na(p_value),
        "NA",
        format.pval(
          p_value,
          digits = 2,
          eps = 1e-3
        )
      ),
      "\n",
      direction
    )
  ]

  p_tf_edges <- ggplot(
    tf_edge_long,
    aes(
      x = response,
      y = edge_weight,
      fill = response
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.65
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.35,
      size = 0.8
    ) +
    geom_text(
      data = tf_plot_summary,
      aes(
        x = 1.5,
        y = y_pos,
        label = label
      ),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 3
    ) +
    scale_fill_manual(
      values = c(
        resistance = "#2166AC",
        response = "#B2182B"
      )
    ) +
    labs(
      title = paste0(
        selected_tf,
        " target edges"
      ),
      x = NULL,
      y = "Edge weight"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(
        t = 5.5,
        r = 5.5,
        b = 28,
        l = 5.5
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.02, 0.20)
      )
    )

  tf_plot_list[[selected_tf]] <- p_tf_edges
}


# ========================================
# Export TF target-edge plots
# ========================================

if (
  length(tf_plot_list) > 0 &&
  !is.null(TF_TARGET_EDGES_PD1_BOXPLOT)
) {

  tf_plot_grid <- ggpubr::ggarrange(
    plotlist = tf_plot_list,
    ncol = 5,
    nrow = ceiling(
      length(tf_plot_list) / 5
    ),
    align = "hv"
  )

  pdf(
    TF_TARGET_EDGES_PD1_BOXPLOT,
    width = 18,
    height = max(
      6,
      ceiling(length(tf_plot_list) / 5) * 4
    )
  )

  print(tf_plot_grid)
  dev.off()

  cat(
    "TF target-edge boxplots saved to",
    TF_TARGET_EDGES_PD1_BOXPLOT,
    "\n"
  )

} else {

  message(
    "No TF target-edge plots were generated."
  )
}


# ==========================================================
# PLOTTING: TF / CD274 expression boxplots
# ==========================================================

tf_expr_plot_list <- list()

for (selected_gene in selected_expression_genes) {

  if (!(selected_gene %in% rownames(expression_clean))) {
    message(
      "Gene not found in expression matrix: ",
      selected_gene
    )
    next
  }

  tf_expr <- expression_clean[
    selected_gene,
    ,
    drop = TRUE
  ]

  tf_expr_dt <- data.table(
    sample_id = colnames(expression_clean),
    expression = as.numeric(tf_expr)
  )

  tf_expr_dt <- merge(
    tf_expr_dt,
    clinical_expression_ordered[
      ,
      .(sample_id, response)
    ],
    by = "sample_id",
    all.x = TRUE
  )

  tf_expr_dt <- tf_expr_dt[
    !is.na(response) &
      !is.na(expression)
  ]

  tf_expr_dt[
    ,
    response := factor(
      response,
      levels = c(
        "resistance",
        "response"
      )
    )
  ]

  if (length(unique(na.omit(tf_expr_dt$response))) < 2) {
    message(
      "Only one response group available for gene: ",
      selected_gene
    )
    next
  }

  p_value <- tryCatch(
    wilcox.test(
      expression ~ response,
      data = tf_expr_dt
    )$p.value,
    error = function(e) NA_real_
  )

  med_resistance <- median(
    tf_expr_dt[
      response == "resistance",
      expression
    ],
    na.rm = TRUE
  )

  med_response <- median(
    tf_expr_dt[
      response == "response",
      expression
    ],
    na.rm = TRUE
  )

  direction <- data.table::fifelse(
    med_response > med_resistance,
    "higher in response",
    data.table::fifelse(
      med_response < med_resistance,
      "higher in resistance",
      "same median"
    )
  )

  y_range <- range(
    tf_expr_dt$expression,
    na.rm = TRUE
  )

  y_pos <- y_range[2] +
    0.08 * diff(y_range)

  if (!is.finite(y_pos)) {
    y_pos <- y_range[2]
  }

  tf_expr_summary <- data.table(
    p_value = p_value,
    med_resistance = med_resistance,
    med_response = med_response,
    direction = direction,
    y_pos = y_pos
  )

  tf_expr_summary[
    ,
    label := paste0(
      "p = ",
      ifelse(
        is.na(p_value),
        "NA",
        format.pval(
          p_value,
          digits = 2,
          eps = 1e-3
        )
      ),
      "\n",
      direction
    )
  ]

  p_tf_expr <- ggplot(
    tf_expr_dt,
    aes(
      x = response,
      y = expression,
      fill = response
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.65
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.35,
      size = 0.8
    ) +
    geom_text(
      data = tf_expr_summary,
      aes(
        x = 1.5,
        y = y_pos,
        label = label
      ),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 3
    ) +
    scale_fill_manual(
      values = c(
        resistance = "#2166AC",
        response = "#B2182B"
      )
    ) +
    labs(
      title = paste0(
        selected_gene,
        " expression"
      ),
      x = NULL,
      y = "Log2 expression"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.margin = margin(
        t = 5.5,
        r = 5.5,
        b = 28,
        l = 5.5
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.02, 0.20)
      )
    )

  tf_expr_plot_list[[selected_gene]] <- p_tf_expr
}


# ========================================
# Export expression plots
# ========================================

if (
  length(tf_expr_plot_list) > 0 &&
  !is.null(TF_EXPRESSION_BOXPLOT)
) {

  tf_expr_plot_grid <- ggpubr::ggarrange(
    plotlist = tf_expr_plot_list,
    ncol = 5,
    nrow = ceiling(
      length(tf_expr_plot_list) / 5
    ),
    align = "hv"
  )

  pdf(
    TF_EXPRESSION_BOXPLOT,
    width = 18,
    height = max(
      6,
      ceiling(length(tf_expr_plot_list) / 5) * 4
    )
  )

  print(tf_expr_plot_grid)
  dev.off()

  cat(
    "TF expression boxplots saved to",
    TF_EXPRESSION_BOXPLOT,
    "\n"
  )

} else {

  message(
    "No TF expression plots were generated."
  )
}