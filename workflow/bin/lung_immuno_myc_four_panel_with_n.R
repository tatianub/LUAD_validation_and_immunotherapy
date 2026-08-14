#####################
## Load R packages ##
#####################

required_libraries <- c("data.table",
  "dplyr", "ggplot2",
  "optparse", "ggpubr",
  "grid", "tools",
  "fgsea", "readxl"
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

  optparse::make_option("--clinical_file",
    type = "character", default = NULL,
    help = "Path to clinical metadata file."
  ),

  optparse::make_option("--gmt_file",
    type = "character", default = NULL,
    help = "Path to pathway GMT file."
  ),

  optparse::make_option("--network_file",
    type = "character", default = NULL,
    help = "Path to inferred network file (.RData)."
  ),

  optparse::make_option("--edges_file",
    type = "character", default = NULL,
    help = "Path to network edge metadata file."
  ),

  optparse::make_option("--expression_file",
    type = "character", default = NULL,
    help = "Path to expression matrix."
  ),

  optparse::make_option("--samples_file",
    type = "character", default = NULL,
    help = "Path to sample annotation file matching expression columns."
  ),

  optparse::make_option("--clinical_file_cnv",
    type = "character", default = NULL,
    help = "Path to Excel file containing MYC CNV data."
  ),

  optparse::make_option("--cnv_sheet",
    type = "character", default = "Table_S8_Gistic_Gene_Events",
    help = "Excel sheet containing CNV data."
  ),

  optparse::make_option("--histo_subtype",
    type = "character", default = NULL,
    help = "Histological subtype filter."
  ),

  optparse::make_option("--treatment_type",
    type = "character", default = NULL,
    help = "Treatment type filter."
  ),

  optparse::make_option("--seed",
    type = "integer", default = NULL,
    help = "Random seed."
  ),

  optparse::make_option("--myc_four_panel_summary_file",
    type = "character", default = NULL,
    help = "Output PDF for the four-panel MYC summary."
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)

opt <- optparse::parse_args(opt_parser)

if (!is.null(opt$seed)) {
  set.seed(opt$seed)
}

#########################
## Assign input options ##
#########################

CLINICAL_FILE <- opt$clinical_file
GMT_FILE <- opt$gmt_file
NETWORK_FILE <- opt$network_file
EDGES_FILE <- opt$edges_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file

CLINICAL_FILE_CNV <- opt$clinical_file_cnv
CNV_SHEET <- opt$cnv_sheet

HISTO_SUBTYPE <- opt$histo_subtype
TREATMENT_TYPE <- opt$treatment_type

MYC_FOUR_PANEL_SUMMARY_FILE <- opt$myc_four_panel_summary_file

####################
## Input checks   ##
####################

required_files <- list(clinical_file = CLINICAL_FILE,
  gmt_file = GMT_FILE, network_file = NETWORK_FILE,
  edges_file = EDGES_FILE, expression_file = EXPRESSION_FILE,
  samples_file = SAMPLES_FILE, clinical_file_cnv = CLINICAL_FILE_CNV
)

for (input_name in names(required_files)) {

  input_file <- required_files[[input_name]]

  if (is.null(input_file)) {
    stop("Missing required argument --",
      input_name
    )
  }

  if (!file.exists(input_file)) {
    stop("File does not exist: ",
      input_file
    )
  }
}

if (is.null(MYC_FOUR_PANEL_SUMMARY_FILE)) {
  stop("Missing required argument ",
    "--myc_four_panel_summary_file"
  )
}

######################
## Global settings  ##
######################

RESPONSE_COLORS <- c(resistance = "#2166AC",
  response = "#B2182B"
)

source("workflow/bin/lung_immuno_limma_fn.R")

# ========================================
# Utility functions
# ========================================

normalize_sample_id <- function(x) {

  toupper(gsub(
      "[^A-Za-z0-9]", ".",
      as.character(x)
    )
  )
}

prepare_myc_corr_dt <- function(dt,
  value_cols, clinical_dt
) {

  dt <- merge(dt,
    clinical_dt[
      , .(sample_id, response)
    ], by = "sample_id",
    all.x = TRUE, sort = FALSE
  )

  non_missing_idx <- !is.na(dt$response)

  for (col_name in value_cols) {

    non_missing_idx <- non_missing_idx &
      !is.na(dt[[col_name]])

  }

  dt <- dt[non_missing_idx]

  dt[
    , response := factor(
      response, levels = c(
        "resistance", "response"
      )
    )
  ]

  dt
}

style_myc_corr_plot <- function(plot_obj) {

  plot_obj +

    geom_point(size = 2.5,
      alpha = 0.6
    ) +

    geom_smooth(method = "lm",
      se = TRUE, color = "black",
      alpha = 0.15
    ) +

    scale_color_manual(values = RESPONSE_COLORS
    ) +

    scale_fill_manual(values = RESPONSE_COLORS
    ) +

    theme_bw() +

    theme(legend.position = "bottom")
}

safe_cor_test <- function(x,
  y, method = "spearman"
) {

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3) {

    return(list(
        estimate = NA_real_, p.value = NA_real_
      )
    )
  }

  if (length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {

    return(list(
        estimate = NA_real_, p.value = NA_real_
      )
    )
  }

  tryCatch(

    cor.test(x,
      y, method = method,
      exact = FALSE
    ),

    error = function(e) {

      list(estimate = NA_real_,
        p.value = NA_real_
      )

    }
  )
}

format_cor_label <- function(cor_test,
  prefix = "rho"
) {

  estimate <- as.numeric(cor_test$estimate)

  paste0(prefix,
    " = ", ifelse(
      is.na(estimate), "NA",
      format(estimate,
        digits = 3
      )
    ), ", p = ",
    ifelse(is.na(cor_test$p.value),
      "NA", format.pval(
        cor_test$p.value, digits = 2,
        eps = 1e-3
      )
    )
  )
}

# ========================================
# Step 1: Load and filter clinical data
# ========================================

clinical_data <- fread(CLINICAL_FILE)
head(clinical_data)

if (!"sample_id" %in% names(clinical_data)) {
  stop("'sample_id' column is missing ",
    "from clinical data."
  )
}

if (!"response" %in% names(clinical_data)) {
  stop("'response' column is missing ",
    "from clinical data."
  )
}

clinical_data <- clinical_data[
  !is.na(response)
]

cat("Filtering clinical data for subtype:",
  HISTO_SUBTYPE, "and treatment type:",
  TREATMENT_TYPE, "\n"
)

clinical_data_filt <- subset_clinical_data(

  clinical_data = clinical_data,

  subtype = HISTO_SUBTYPE,

  treatment_type = TREATMENT_TYPE

)
dim(clinical_data_filt)
if (nrow(clinical_data_filt) == 0) {
  stop("No samples remain after ",
    "clinical filtering."
  )
}

if (anyDuplicated(clinical_data_filt$sample_id)
) {

  stop("Duplicated sample IDs found ",
    "in filtered clinical data."
  )
}

cat("Clinical samples after filtering:",
  nrow(clinical_data_filt), "\n"
)

# ========================================
# Step 2: Load network
# ========================================

net_env <- new.env()

load(NETWORK_FILE,
  envir = net_env
)

obj_names <- ls(net_env)

cat("Available network objects:",
  paste(obj_names,
    collapse = ", "
  ), "\n"
)

if ("net" %in% obj_names) {

  net <- net_env$net

  cat("Loaded 'net' object\n")

} else if ("net_norm" %in% obj_names) {

  net <- net_env$net_norm

  cat("Loaded 'net_norm' object\n")

} else {

  stop("No recognized network object found. ",
    "Available objects: ", paste(
      obj_names, collapse = ", "
    )
  )
}

net <- as.matrix(net)
colnames(net) <- gsub("[^A-Za-z0-9]", ".",
  colnames(net)
)

if (is.null(colnames(net))) {
  stop("Network matrix has no sample column names.")
}

# ========================================
# Step 3: Load network edge metadata
# ========================================

edges <- fread(EDGES_FILE)

if (!all(c("reg", "tar") %in% names(edges)
)) {

  stop("Edge file must contain columns ",
    "'reg' and 'tar'."
  )
}

if (nrow(edges) != nrow(net)) {

  stop("Number of edge rows (",
    nrow(edges), ") does not match number of network rows (",
    nrow(net), ")."
  )
}

# ========================================
# Step 4: Load expression matrix
# ========================================

expression_dt <- fread(EXPRESSION_FILE)

gene_ids <- expression_dt[[1]]

head(expression_dt)
expression <- data.matrix(expression_dt[
    , -1,
    with = FALSE
  ]
)

samples <- fread(SAMPLES_FILE)
head(samples)
if (ncol(samples) < 2) {

  stop("Samples file must contain at least ",
    "two columns; sample IDs are expected ", "in column 2."
  )
}

if (ncol(expression) != nrow(samples)) {

  stop("Expression matrix contains ",
    ncol(expression), " sample columns, while samples file contains ",
    nrow(samples), " rows."
  )
}

colnames(expression) <- samples[[2]]
rownames(expression) <- gene_ids
colnames(expression) <- gsub("[^A-Za-z0-9]", ".",
  colnames(expression)
)

if (anyDuplicated(colnames(expression))
) {

  stop("Duplicated sample IDs found ",
    "in expression matrix."
  )
}

# ========================================
# Step 5: Identify common samples
# ========================================

common_samples <- colnames(net)[
  colnames(net) %in%
    clinical_data_filt$sample_id 
]

if (length(common_samples) == 0) {

  stop("No common samples among network, ",
    "expression and filtered clinical data."
  )
}

cat("Common network/expression/clinical samples:",
  length(common_samples), "\n"
)

# Keep NETWORK order as canonical sample order

net_clean <- net[
  , common_samples,
  drop = FALSE
]

expression_clean <- expression[
  , common_samples,
  drop = FALSE
]

clinical_data_ordered <- clinical_data_filt[
    match(common_samples,
      sample_id
    )
  ]

stopifnot(all(
    colnames(net_clean) ==
      clinical_data_ordered$sample_id
  )
)

stopifnot(all(
    colnames(expression_clean) ==
      clinical_data_ordered$sample_id
  )
)

# ========================================
# Step 6: Restrict network to PD-1 pathway
# ========================================

pathways <- fgsea::gmtPathways(GMT_FILE)

pd1_pathway_names <- grep("PD_1",
  names(pathways), value = TRUE
)

if (length(pd1_pathway_names) == 0) {

  stop("No pathway containing 'PD_1' ",
    "was found in GMT file."
  )
}

cat("PD-1 pathway(s):",
  paste(pd1_pathway_names,
    collapse = ", "
  ), "\n"
)

pd1_pathways <- pathways[
    pd1_pathway_names
  ]

pathways_genes <- unique(unlist(
    pd1_pathways, use.names = FALSE
  )
)

pd1_edge_idx <- edges$tar %in%
  pathways_genes

net_pathways <- net_clean[
  pd1_edge_idx, ,
  drop = FALSE
]

edges_pathways <- edges[
  pd1_edge_idx
]

stopifnot(nrow(net_pathways) ==
    nrow(edges_pathways)
)

cat("Number of PD-1 pathway edges:",
  nrow(edges_pathways), "\n"
)

# ========================================
# Step 7: Find MYC and CD274 expression
# ========================================

myc_expr_idx <- which(rownames(expression_clean) == "MYC")

cd274_expr_idx <- which(rownames(expression_clean) == "CD274")

if (length(myc_expr_idx) == 0) {

  stop("MYC was not found ",
    "in expression matrix."
  )
}

if (length(cd274_expr_idx) == 0) {

  stop("CD274 was not found ",
    "in expression matrix."
  )
}

if (length(myc_expr_idx) > 1) {

  stop("Multiple MYC rows found ",
    "in expression matrix."
  )
}

if (length(cd274_expr_idx) > 1) {

  stop("Multiple CD274 rows found ",
    "in expression matrix."
  )
}

# ========================================
# Panel A:
# MYC expression vs MYC -> CD274 edge
# ========================================

myc_cd274_edge_idx <- which(edges_pathways$reg == "MYC" &
  edges_pathways$tar == "CD274"
)

if (length(myc_cd274_edge_idx) == 0) {

  stop("MYC -> CD274 edge was not found ",
    "in the PD-1 network."
  )
}

if (length(myc_cd274_edge_idx) > 1) {

  stop("More than one MYC -> CD274 edge ",
    "was found in the network."
  )
}

myc_expression <- as.numeric(expression_clean[
    myc_expr_idx, ,
    drop = TRUE
  ]
)

myc_cd274_edge <- as.numeric(net_pathways[
    myc_cd274_edge_idx, ,
    drop = TRUE
  ]
)

corr_dt <- data.table(sample_id = common_samples,
  myc_expression = myc_expression, myc_cd274_edge = myc_cd274_edge
)

corr_dt <- prepare_myc_corr_dt(

  corr_dt,

  value_cols = c("myc_expression",
    "myc_cd274_edge"
  ),

  clinical_dt =
    clinical_data_ordered
)

cor_test <- safe_cor_test(

  corr_dt$myc_expression,

  corr_dt$myc_cd274_edge,

  method = "spearman"
)

p_corr <- ggplot(

  corr_dt,

  aes(x = myc_expression,
    y = myc_cd274_edge, fill = response,
    color = response
  )

) +

  labs(

    title =
      "MYC expression vs MYC→CD274 edge weight",

    x =
      "MYC log2 expression",

    y =
      "MYC→CD274 edge weight",

    subtitle = paste0(
      format_cor_label(cor_test,
        prefix = "Spearman rho"
      ),
      "; n = ",
      data.table::uniqueN(corr_dt$sample_id)
    )
  )

p_corr <- style_myc_corr_plot(p_corr)

# ========================================
# Panel B:
# MYC expression vs CD274 expression
# ========================================

cd274_expression <- as.numeric(

  expression_clean[
    cd274_expr_idx, ,
    drop = TRUE
  ]

)

expr_corr_dt <- data.table(

  sample_id = common_samples,

  myc_expression =
    myc_expression,

  cd274_expression =
    cd274_expression

)

expr_corr_dt <- prepare_myc_corr_dt(

  expr_corr_dt,

  value_cols = c("myc_expression",
    "cd274_expression"
  ),

  clinical_dt =
    clinical_data_ordered
)

expr_cor_test <- safe_cor_test(

  expr_corr_dt$myc_expression,

  expr_corr_dt$cd274_expression,

  method = "pearson"
)

p_expr_corr <- ggplot(

  expr_corr_dt,

  aes(x = cd274_expression,
    y = myc_expression, fill = response,
    color = response
  )

) +

  labs(

    title =
      "MYC expression vs CD274 expression",

    x =
      "CD274 log2 expression",

    y =
      "MYC log2 expression",

    subtitle = paste0(
      format_cor_label(expr_cor_test,
        prefix = "Pearson r"
      ),
      "; n = ",
      data.table::uniqueN(expr_corr_dt$sample_id)
    )
  )

p_expr_corr <- style_myc_corr_plot(p_expr_corr)

# ========================================
# Step 8: Load MYC CNV
# ========================================

cnv_data <- data.table::as.data.table(

    readxl::read_excel(CLINICAL_FILE_CNV,
      sheet = CNV_SHEET
    )

  )

gene_col_candidates <- c("Hugo_Symbol",
  "Gene", "Gene_Symbol",
  "Symbol"
)

gene_col <- intersect(gene_col_candidates,
  names(cnv_data)
)

if (length(gene_col) > 0) {

  gene_col <- gene_col[1]

} else {

  gene_col <- names(cnv_data)[1]

  message("No recognized gene-symbol column found. ",
    "Using first column: ", gene_col
  )
}

cnv_data[
  , gene_symbol :=
    toupper(trimws(
        as.character(get(gene_col))
      )
    )
]

myc_cnv_row <- cnv_data[
    gene_symbol == "MYC"
  ]

if (nrow(myc_cnv_row) == 0) {

  stop("MYC was not found in CNV sheet: ",
    CNV_SHEET
  )
}

if (nrow(myc_cnv_row) > 1) {

  stop("More than one MYC row was found ",
    "in CNV sheet: ", CNV_SHEET
  )
}

sample_cols <- setdiff(

  names(cnv_data),

  c(gene_col,
    "gene_symbol"
  )

)

myc_cnv_long <- data.table::melt(

  myc_cnv_row[
    , ..sample_cols
  ],

  variable.name =
    "cnv_sample_col",

  value.name =
    "myc_cnv"

)

myc_cnv_long[
  , sample_id_norm :=
    normalize_sample_id(cnv_sample_col)
]

myc_cnv_long[
  , myc_cnv :=
    suppressWarnings(as.numeric(myc_cnv))
]

# myc_cnv_long <- myc_cnv_long[
#     !is.na(myc_cnv)
#   ]

# ========================================
# Step 9: Calculate MYC PD-1 outdegree
# ========================================

myc_pd1_idx <- which(edges_pathways$reg == "MYC")

if (length(myc_pd1_idx) == 0) {

  stop("No MYC -> PD-1 pathway edges ",
    "found in edges_pathways."
  )
}

myc_pd1_edge_matrix <- net_pathways[
    myc_pd1_idx, ,
    drop = FALSE
  ]

# Sum of MYC -> PD-1 pathway edge weights
# for each patient

myc_outdegree <- colSums(myc_pd1_edge_matrix,
  na.rm = TRUE
)

myc_outdegree_dt <- data.table(

  sample_id =
    names(myc_outdegree),

  myc_outdegree =
    as.numeric(myc_outdegree)

)

myc_outdegree_dt[
  , sample_id_norm :=
    normalize_sample_id(sample_id)
]

# ========================================
# Step 10: Match CNV to clinical/network
# ========================================

clinical_sample_map <- unique(

  clinical_data_ordered[
    , .(
      sample_id, sample_id_norm =
        normalize_sample_id(sample_id
        ),
      response
    )
  ],

  by = "sample_id_norm"

)

# Attach canonical sample_id to CNV data
myc_cnv_matched <- merge(

  myc_cnv_long,

  clinical_sample_map,

  by = "sample_id_norm",

  all = FALSE,

  sort = FALSE

)

cat("MYC CNV samples matched to clinical data:",
  nrow(myc_cnv_matched), "\n"
)

# Merge CNV and MYC outdegree BY NORMALIZED ID
cnv_outdegree_dt <- merge(

  myc_outdegree_dt,

  myc_cnv_matched[
    , .(
      sample_id_norm, myc_cnv,
      response
    )
  ],

  by = "sample_id_norm",

  all = FALSE,

  sort = FALSE

)

cnv_outdegree_dt <- cnv_outdegree_dt[
    !is.na(myc_cnv) & !is.na(myc_outdegree)
  ]

if (nrow(cnv_outdegree_dt) < 3) {

  stop("Too few matched samples for ",
    "MYC CNV vs MYC outdegree analysis."
  )
}

cnv_outdegree_dt[
  , response := factor(
    response, levels = c(
      "resistance", "response"
    )
  )
]

cat("Samples used for MYC CNV analysis:",
  nrow(cnv_outdegree_dt), "\n"
)

# ========================================
# Panel C:
# MYC CNV by response
# ========================================

cnv_range <- range(cnv_outdegree_dt$myc_cnv,
  na.rm = TRUE
)

cnv_y_pos <- cnv_range[2] +
  0.08 * diff(cnv_range)

if (!is.finite(cnv_y_pos) ||
  diff(cnv_range) == 0
) {

  cnv_y_pos <- cnv_range[2] + 0.5

}

myc_cnv_box_summary <- data.table(

    p_value = tryCatch(

      wilcox.test(myc_cnv ~ response,
        data = cnv_outdegree_dt
      )$p.value,

      error = function(e)
        NA_real_

    ),

    med_resistance = median(

      cnv_outdegree_dt[
        response == "resistance", myc_cnv
      ],

      na.rm = TRUE
    ),

    med_response = median(

      cnv_outdegree_dt[
        response == "response", myc_cnv
      ],

      na.rm = TRUE
    ),

    y_pos = cnv_y_pos
  )

myc_cnv_box_summary[
  , direction := fifelse(

    med_response >
      med_resistance,

    "higher in response",

    fifelse(

      med_response <
        med_resistance,

      "higher in resistance",

      "same median"
    )
  )
]

myc_cnv_box_summary[
  , label := paste0(

    "p = ",

    ifelse(

      is.na(p_value),

      "NA",

      format.pval(p_value,
        digits = 2, eps = 1e-3
      )
    ),

    "\n",

    direction
  )
]

n_cnv_total <- data.table::uniqueN(
  cnv_outdegree_dt$sample_id
)

n_cnv_resistance <- data.table::uniqueN(
  cnv_outdegree_dt[
    response == "resistance", sample_id
  ]
)

n_cnv_response <- data.table::uniqueN(
  cnv_outdegree_dt[
    response == "response", sample_id
  ]
)

p_myc_cnv_box <- ggplot(

  cnv_outdegree_dt,

  aes(x = response,
    y = myc_cnv, fill = response
  )

) +

  geom_boxplot(outlier.shape = NA,
    alpha = 0.8, width = 0.65
  ) +

  geom_jitter(width = 0.15,
    alpha = 0.35, size = 0.8
  ) +

  geom_text(

    data =
      myc_cnv_box_summary,

    aes(x = 1.5,
      y = y_pos, label = label
    ),

    inherit.aes = FALSE,

    vjust = -0.4,

    size = 3
  ) +

  scale_fill_manual(values = RESPONSE_COLORS
  ) +

  labs(

    title =
      "MYC CNV by response",

    subtitle = paste0(
      "n = ", n_cnv_total,
      " (resistance = ", n_cnv_resistance,
      ", response = ", n_cnv_response, ")"
    ),

    x = NULL,

    y =
      "MYC CNV"
  ) +

  theme_bw() +

  theme(legend.position = "none"
  ) +

  coord_cartesian(clip = "off"
  ) +

  scale_y_continuous(

    expand = expansion(mult = c(
        0.02, 0.20
      )
    )

  )

# ========================================
# Panel D:
# MYC CNV vs MYC PD-1 outdegree
# ========================================

cnv_cor_test <- safe_cor_test(

  cnv_outdegree_dt$myc_cnv,

  cnv_outdegree_dt$myc_outdegree,

  method = "spearman"
)

p_cnv_outdegree <- ggplot(

  cnv_outdegree_dt,

  aes(x = myc_cnv,
    y = myc_outdegree, fill = response,
    color = response
  )

) +

  geom_point(size = 2.5,
    alpha = 0.7
  ) +

  geom_smooth(method = "lm",
    se = TRUE, color = "black",
    alpha = 0.15
  ) +

  scale_color_manual(values = RESPONSE_COLORS
  ) +

  scale_fill_manual(values = RESPONSE_COLORS
  ) +

  labs(

    title =
      "MYC CNV vs MYC outdegree (PD-1 pathway)",

    x =
      "MYC CNV",

    y =
      "MYC outdegree (sum of edge weights)",

    subtitle = paste0(
      format_cor_label(cnv_cor_test,
        prefix = "Spearman rho"
      ),
      "; n = ",
      data.table::uniqueN(cnv_outdegree_dt$sample_id)
    )
  ) +

  theme_bw() +

  theme(legend.position = "bottom")

# ========================================
# Step 11: Combine and export
# ========================================

four_panel_summary <- ggpubr::ggarrange(

    p_corr,

    p_expr_corr,

    p_myc_cnv_box,

    p_cnv_outdegree,

    labels =
      c("A",
        "B", "C",
        "D"
      ),

    ncol = 2,

    nrow = 2,

    align = "hv"
  )

pdf(

  MYC_FOUR_PANEL_SUMMARY_FILE,

  width = 12,

  height = 12

)

print(four_panel_summary)

dev.off()

cat("MYC four-panel summary saved to ",
  MYC_FOUR_PANEL_SUMMARY_FILE, "\n"
)
