#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table", "dplyr",
  "ggplot2", "optparse",
  "fgsea"
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
    c("--clinical_file"), type = "character",
    default = NULL, help = "Path to clinical metadata file."
  ),

  optparse::make_option(
    c("--gmt_file"), type = "character",
    default = NULL, help = "Path to pathway GMT file."
  ),

  optparse::make_option(
    c("--network_file"), type = "character",
    default = NULL, help = "Path to inferred network file (.RData)."
  ),

  optparse::make_option(
    c("--edges_file"), type = "character",
    default = NULL, help = "Path to network edge file."
  ),

  optparse::make_option(
    c("--expression_file"), type = "character",
    default = NULL, help = "Path to expression matrix."
  ),

  optparse::make_option(
    c("--samples_file"), type = "character",
    default = NULL, help = "Path to sample file matching expression columns."
  ),

  optparse::make_option(
    c("--immune_file"), type = "character",
    default = NULL, help = "Path to immune infiltration score file."
  ),

  optparse::make_option(
    c("--outdegree_file"), type = "character",
    default = NULL, help = paste0(
      "Path to full-network TF outdegree file. ",
      "Required only for MYC total-network outdegree analysis."
    )
  ),

  optparse::make_option(
    c("--histo_subtype"), type = "character",
    default = NULL, help = "Histological subtype filter."
  ),

  optparse::make_option(
    c("--treatment_type"), type = "character",
    default = NULL, help = "Treatment type filter."
  ),

  optparse::make_option(
    c("--seed"), type = "integer",
    default = NULL, help = "Random seed."
  ),

  optparse::make_option(
    c("--myc_pd1_immune_plot"), type = "character",
    default = NULL, help = paste0(
      "Output PDF for MYC PD-1-pathway outdegree ",
      "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--myc_pd1_immune_table"), type = "character",
    default = NULL, help = paste0(
      "Output TSV for MYC PD-1-pathway outdegree ",
      "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--myc_total_immune_plot"), type = "character",
    default = NULL, help = paste0(
      "Output PDF for MYC total-network outdegree ",
      "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--myc_total_immune_table"), type = "character",
    default = NULL, help = paste0(
      "Output TSV for MYC total-network outdegree ",
      "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--tf_pd1_immune_heatmap"), type = "character",
    default = NULL, help = paste0(
      "Output PDF for selected TF PD-1-pathway outdegree ",
      "vs immune infiltration heatmap."
    )
  ),

  optparse::make_option(
    c("--tf_pd1_immune_table"), type = "character",
    default = NULL, help = paste0(
      "Output TSV for selected TF PD-1-pathway outdegree ",
      "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--myc_expression_immune_plot"), type = "character",
    default = NULL, help = paste0(
      "Output PDF for MYC expression ", "vs immune infiltration correlations."
    )
  ),

  optparse::make_option(
    c("--myc_expression_immune_table"), type = "character",
    default = NULL, help = paste0(
      "Output TSV for MYC expression ", "vs immune infiltration correlations."
    )
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)

opt <- optparse::parse_args(opt_parser)

if (!is.null(opt$seed)) {
  set.seed(opt$seed)
}

##########################
## Assign input options ##
##########################

CLINICAL_FILE <- opt$clinical_file
GMT_FILE <- opt$gmt_file
NETWORK_FILE <- opt$network_file
EDGES_FILE <- opt$edges_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
IMMUNE_FILE <- opt$immune_file
OUTDEGREE_FILE <- opt$outdegree_file

HISTO_SUBTYPE <- opt$histo_subtype
TREATMENT_TYPE <- opt$treatment_type

MYC_PD1_IMMUNE_PLOT <- opt$myc_pd1_immune_plot
MYC_PD1_IMMUNE_TABLE <- opt$myc_pd1_immune_table

MYC_TOTAL_IMMUNE_PLOT <- opt$myc_total_immune_plot
MYC_TOTAL_IMMUNE_TABLE <- opt$myc_total_immune_table

TF_PD1_IMMUNE_HEATMAP <- opt$tf_pd1_immune_heatmap
TF_PD1_IMMUNE_TABLE <- opt$tf_pd1_immune_table

MYC_EXPRESSION_IMMUNE_PLOT <- opt$myc_expression_immune_plot
MYC_EXPRESSION_IMMUNE_TABLE <- opt$myc_expression_immune_table

##################
## Input checks ##
##################

required_files <- list(
  clinical_file = CLINICAL_FILE, gmt_file = GMT_FILE,
  network_file = NETWORK_FILE, edges_file = EDGES_FILE,
  expression_file = EXPRESSION_FILE, samples_file = SAMPLES_FILE,
  immune_file = IMMUNE_FILE
)

for (input_name in names(required_files)) {

  input_file <- required_files[[input_name]]

  if (is.null(input_file)) {
    stop(
      "Missing required argument --", input_name
    )
  }

  if (!file.exists(input_file)) {
    stop(
      "File does not exist: ", input_file
    )
  }
}

# Full-network MYC outdegree analysis is optional,
# but if one of its outputs is requested, outdegree_file is required.
if (
  !is.null(MYC_TOTAL_IMMUNE_PLOT) || !is.null(MYC_TOTAL_IMMUNE_TABLE)
) {

  if (is.null(OUTDEGREE_FILE)) {
    stop(
      "--outdegree_file is required when requesting ",
      "MYC total-network immune correlation outputs."
    )
  }

  if (!file.exists(OUTDEGREE_FILE)) {
    stop(
      "Outdegree file does not exist: ", OUTDEGREE_FILE
    )
  }
}

########################
## Helper definitions ##
########################

source("workflow/bin/lung_immuno_limma_fn.R")

RESPONSE_COLORS <- c(
  resistance = "#2166AC", response = "#B2182B"
)

normalize_sample_id <- function(x) {

  toupper(
    gsub(
      "[^A-Za-z0-9]", ".",
      as.character(x)
    )
  )
}

safe_spearman <- function(x, y) {

  x <- suppressWarnings(as.numeric(x))

  y <- suppressWarnings(as.numeric(y))

  valid_idx <- which(
    is.finite(x) &
      is.finite(y)
  )

  n_valid <- length(valid_idx)

  if (n_valid < 3) {

    return(
      list(
        estimate = NA_real_, p.value = NA_real_,
        n = n_valid
      )
    )
  }

  # Correlation is undefined if either variable is constant.
  if (
    length(unique(x[valid_idx])) < 2 || length(unique(y[valid_idx])) < 2
  ) {

    return(
      list(
        estimate = NA_real_, p.value = NA_real_,
        n = n_valid
      )
    )
  }

  cor_result <- tryCatch(

    cor.test(
      x[valid_idx], y[valid_idx],
      method = "spearman", exact = FALSE
    ),

    error = function(e) NULL

  )

  if (is.null(cor_result)) {

    return(
      list(
        estimate = NA_real_, p.value = NA_real_,
        n = n_valid
      )
    )
  }

  list(
    estimate = as.numeric(
      cor_result$estimate
    ), p.value = as.numeric(
      cor_result$p.value
    ), n = n_valid
  )
}

format_cor_subtitle <- function(result) {

  rho_text <- if (
    is.na(result$estimate)
  ) {
    "NA"
  } else {
    format(
      result$estimate, digits = 3
    )
  }

  p_text <- if (
    is.na(result$p.value)
  ) {
    "NA"
  } else {
    format.pval(
      result$p.value, digits = 2,
      eps = 1e-3
    )
  }

  paste0(
    "Spearman rho = ", rho_text,
    ", p = ", p_text
  )
}

# ========================================
# Step 1: Load and filter clinical data
# ========================================

clinical_data <- fread(CLINICAL_FILE)

if (
  !"sample_id" %in%
    names(clinical_data)
) {
  stop(
    "'sample_id' column is missing ", "from clinical data."
  )
}

if (
  !"response" %in%
    names(clinical_data)
) {
  stop(
    "'response' column is missing ", "from clinical data."
  )
}

clinical_data <- clinical_data[
  !is.na(response)
]

cat(
  "Filtering clinical data for histological subtype:", HISTO_SUBTYPE,
  "and treatment type:", TREATMENT_TYPE,
  "\n"
)

clinical_data_filt <- subset_clinical_data(

  clinical_data = clinical_data,

  subtype = HISTO_SUBTYPE,

  treatment_type = TREATMENT_TYPE

)

if (nrow(clinical_data_filt) == 0) {

  stop("No samples remain after clinical filtering.")
}

if (
  anyDuplicated(clinical_data_filt$sample_id)
) {

  stop(
    "Duplicated sample IDs found ", "in filtered clinical data."
  )
}

cat(
  "Clinical samples after filtering:", nrow(clinical_data_filt),
  "\n"
)

# ========================================
# Step 2: Load inferred network
# ========================================

net_env <- new.env()

load(
  NETWORK_FILE, envir = net_env
)

obj_names <- ls(net_env)

cat(
  "Available network objects:", paste(
    obj_names, collapse = ", "
  ), "\n"
)

if ("net" %in% obj_names) {

  net <- net_env$net

  cat("Loaded 'net' object\n")

} else if (
  "net_norm" %in% obj_names
) {

  net <- net_env$net_norm

  cat("Loaded 'net_norm' object\n")

} else {

  stop(
    "No recognized network object found. ", "Available objects: ",
    paste(
      obj_names, collapse = ", "
    )
  )
}

net <- as.matrix(net)
colnames(net) <- gsub(
  "[^A-Za-z0-9]", ".",
  colnames(net)
)
if (
  is.null(colnames(net))
) {

  stop("Network matrix has no sample column names.")
}

# ========================================
# Step 3: Load edge metadata
# ========================================

edges <- fread(EDGES_FILE)

if (
  !all(
    c(
      "reg", "tar"
    ) %in%
      names(edges)
  )
) {

  stop(
    "Edge file must contain ", "'reg' and 'tar' columns."
  )
}

if (
  nrow(edges) !=
    nrow(net)
) {

  stop(
    "Number of rows in edge file (", nrow(edges),
    ") does not match number of ", "network rows (",
    nrow(net), ")."
  )
}

# ========================================
# Step 4: Load expression matrix
# ========================================

expression_dt <- fread(EXPRESSION_FILE)

if (ncol(expression_dt) < 2) {

  stop(
    "Expression file must contain ", "a gene column plus sample columns."
  )
}

gene_ids <- expression_dt[[1]]

expression <- data.matrix(
  expression_dt[
    , -1,
    with = FALSE
  ]
)

samples <- fread(SAMPLES_FILE)

if (ncol(samples) < 2) {

  stop(
    "Samples file must contain at least ",
    "two columns. Sample IDs are expected ", "in column 2."
  )
}

if (
  ncol(expression) !=
    nrow(samples)
) {

  stop(
    "Expression matrix has ", ncol(expression),
    " sample columns, whereas samples file ", "contains ",
    nrow(samples), " rows."
  )
}

colnames(expression) <- as.character(samples[[2]])
colnames(expression) <- gsub(
  "[^A-Za-z0-9]", ".",
  colnames(expression)
)
rownames(expression) <- as.character(gene_ids)

if (
  anyDuplicated(colnames(expression))
) {

  stop(
    "Duplicated sample IDs found ", "in expression matrix."
  )
}

# ========================================
# Step 5: Align network, expression,
#         and clinical samples
# ========================================

common_samples <- colnames(net)[

  colnames(net) %in%
    clinical_data_filt$sample_id &

    colnames(net) %in%
    colnames(expression)

]

if (length(common_samples) == 0) {

  stop(
    "No common samples among network, ",
    "expression and filtered clinical data."
  )
}

# Network order is used as canonical sample order.
net_clean <- net[
  , common_samples,
  drop = FALSE
]

expression_clean <- expression[
  , common_samples,
  drop = FALSE
]

clinical_data_ordered <- clinical_data_filt[
    match(
      common_samples, sample_id
    )
  ]

stopifnot(
  all(
    colnames(net_clean) ==
      clinical_data_ordered$sample_id
  )
)

stopifnot(
  all(
    colnames(expression_clean) ==
      clinical_data_ordered$sample_id
  )
)

cat(
  "Common network/expression/clinical samples:", length(common_samples),
  "\n"
)

# ========================================
# Step 6: Restrict network to PD-1 pathway
# ========================================

pathways <- fgsea::gmtPathways(GMT_FILE)

pd1_pathway_names <- grep(
  "PD_1", names(pathways),
  value = TRUE
)

if (length(pd1_pathway_names) == 0) {

  stop(
    "No pathway containing 'PD_1' ", "was found in the GMT file."
  )
}

cat(
  "PD-1 pathway(s):", paste(
    pd1_pathway_names, collapse = ", "
  ), "\n"
)

pd1_pathway_genes <- unique(
  unlist(
    pathways[
      pd1_pathway_names
    ], use.names = FALSE
  )
)

pd1_edge_idx <- edges$tar %in%
  pd1_pathway_genes

net_pathways <- net_clean[
  pd1_edge_idx, ,
  drop = FALSE
]

edges_pathways <- edges[
  pd1_edge_idx
]

stopifnot(
  nrow(net_pathways) ==
    nrow(edges_pathways)
)

cat(
  "PD-1 pathway edges:", nrow(edges_pathways),
  "\n"
)

# ========================================
# Step 7: Clinical sample map
# ========================================

clinical_sample_map <- unique(

  clinical_data_ordered[
    , .(
      sample_id, sample_id_norm =
        normalize_sample_id(
          sample_id
        ),
      response
    )
  ],

  by =
    "sample_id_norm"

)

# ========================================
# Step 8: Load immune infiltration data
# ========================================

immune_cell_dt <- fread(
  IMMUNE_FILE, sep = "\t",
  header = TRUE, check.names = FALSE
)

if (ncol(immune_cell_dt) < 2) {

  stop(
    "Immune infiltration file must contain ",
    "a sample ID column plus at least one ", "immune score column."
  )
}

# Assume first column contains sample IDs.
immune_id_col <- names(
    immune_cell_dt
  )[1]

immune_cell_dt[
  , sample_id_norm :=
    normalize_sample_id(
      get(immune_id_col)
    )
]

immune_cell_dt[
  , (immune_id_col) := NULL
]

if (
  anyDuplicated(immune_cell_dt$sample_id_norm)
) {

  stop(
    "Duplicated normalized sample IDs found ", "in immune infiltration file."
  )
}

immune_cell_cols <- setdiff(
  names(
    immune_cell_dt
  ), "sample_id_norm"
)

cat(
  "Immune infiltration variables:", length(
    immune_cell_cols
  ), "\n"
)

# ==========================================================
# Analysis 1:
# MYC PD-1-pathway outdegree vs immune infiltration
# ==========================================================

if (
  !is.null(MYC_PD1_IMMUNE_PLOT) || !is.null(MYC_PD1_IMMUNE_TABLE)
) {

  cat(
    "Running MYC PD-1-pathway outdegree ",
    "vs immune infiltration analysis...\n"
  )

  myc_pd1_idx <- which(
    edges_pathways$reg ==
      "MYC"
  )

  if (length(myc_pd1_idx) == 0) {

    warning(
      "No MYC -> PD-1 pathway edges found. ",
      "Skipping MYC PD-1 outdegree immune analysis."
    )

  } else {

    myc_pd1_edge_matrix <- net_pathways[
        myc_pd1_idx, ,
        drop = FALSE
      ]

    # Sum MYC -> PD-1-pathway target edge weights
    # for each sample.
    myc_pd1_outdegree <- colSums(
      myc_pd1_edge_matrix, na.rm = TRUE
    )

    myc_pd1_outdegree_dt <- data.table(

      sample_id =
        names(
          myc_pd1_outdegree
        ),

      sample_id_norm =
        normalize_sample_id(
          names(myc_pd1_outdegree)
        ),

      myc_pd1_outdegree =
        as.numeric(myc_pd1_outdegree)

    )

    myc_pd1_outdegree_dt <- merge(

      myc_pd1_outdegree_dt,

      clinical_sample_map[
        , .(
          sample_id_norm, response
        )
      ],

      by =
        "sample_id_norm",

      all.x = TRUE,

      sort = FALSE

    )

    myc_pd1_immune_dt <- merge(

      myc_pd1_outdegree_dt,

      immune_cell_dt,

      by =
        "sample_id_norm",

      all = FALSE,

      sort = FALSE

    )

    myc_pd1_immune_dt <- myc_pd1_immune_dt[
        !is.na(response) &
          !is.na(myc_pd1_outdegree)
      ]

    myc_pd1_immune_dt[
      , response := factor(
        response, levels = c(
          "resistance", "response"
        )
      )
    ]

    cat(
      "Samples in MYC PD-1 immune analysis:", nrow(
        myc_pd1_immune_dt
      ), "\n"
    )

    myc_pd1_cor_results <- vector(
      "list", length(immune_cell_cols)
    )

    names(
      myc_pd1_cor_results
    ) <- immune_cell_cols

    if (
      !is.null(MYC_PD1_IMMUNE_PLOT)
    ) {

      pdf(
        MYC_PD1_IMMUNE_PLOT, width = 8,
        height = 6
      )
    }

    for (
      immune_cell in
      immune_cell_cols
    ) {

      immune_values <- suppressWarnings(
          as.numeric(myc_pd1_immune_dt[[immune_cell]])
        )

      cor_result <- safe_spearman(

        myc_pd1_immune_dt$
          myc_pd1_outdegree,

        immune_values

      )

      myc_pd1_cor_results[[immune_cell]] <- data.table(

        immune_cell =
          immune_cell,

        n =
          cor_result$n,

        rho =
          cor_result$estimate,

        p_value =
          cor_result$p.value

      )

      if (
        is.null(MYC_PD1_IMMUNE_PLOT)
      ) {
        next
      }

      valid_idx <- which(

        is.finite(
          myc_pd1_immune_dt$
            myc_pd1_outdegree
        ) &

          is.finite(immune_values)

      )

      if (length(valid_idx) < 3) {
        next
      }

      plot_dt <- copy(
        myc_pd1_immune_dt[
          valid_idx
        ]
      )

      plot_dt[
        , immune_score :=
          immune_values[
            valid_idx
          ]
      ]

      p_immune_corr <- ggplot(

        plot_dt,

        aes(
          x = myc_pd1_outdegree, y = immune_score,
          color = response, fill = response
        )

      ) +

        geom_point(
          size = 2.4, alpha = 0.75
        ) +

        geom_smooth(
          method = "lm", se = TRUE,
          color = "black", alpha = 0.15
        ) +

        scale_color_manual(
          values =
            RESPONSE_COLORS
        ) +

        scale_fill_manual(
          values =
            RESPONSE_COLORS
        ) +

        labs(

          title = paste0(
            "MYC PD-1 outdegree vs ", immune_cell
          ),

          x =
            "MYC PD-1 pathway outdegree",

          y = paste0(
            immune_cell, " score"
          ),

          subtitle =
            format_cor_subtitle(cor_result)

        ) +

        theme_bw() +

        theme(
          legend.position =
            "bottom"
        )

      print(p_immune_corr)
    }

    if (
      !is.null(MYC_PD1_IMMUNE_PLOT)
    ) {

      dev.off()

      cat(
        "MYC PD-1 outdegree immune plots saved to:", MYC_PD1_IMMUNE_PLOT,
        "\n"
      )
    }

    myc_pd1_cor_results <- rbindlist(
        myc_pd1_cor_results, fill = TRUE
      )

    myc_pd1_cor_results[
      , fdr_bh :=
        p.adjust(
          p_value, method = "BH"
        )
    ]

    if (
      !is.null(MYC_PD1_IMMUNE_TABLE)
    ) {

      fwrite(
        myc_pd1_cor_results, MYC_PD1_IMMUNE_TABLE,
        sep = "\t"
      )

      cat(
        "MYC PD-1 outdegree immune table saved to:", MYC_PD1_IMMUNE_TABLE,
        "\n"
      )
    }
  }
}

# ==========================================================
# Analysis 2:
# MYC total-network outdegree vs immune infiltration
# ==========================================================

if (
  !is.null(MYC_TOTAL_IMMUNE_PLOT) || !is.null(MYC_TOTAL_IMMUNE_TABLE)
) {

  cat(
    "Running MYC total-network outdegree ",
    "vs immune infiltration analysis...\n"
  )

  outdegree_dt <- fread(
    OUTDEGREE_FILE, sep = "\t",
    header = TRUE
  )

  if (
    !"reg" %in%
      names(outdegree_dt)
  ) {

    stop(
      "Outdegree file must contain ", "a 'reg' column."
    )
  }

  myc_outdegree_row <- outdegree_dt[
      reg == "MYC"
    ]

  if (
    nrow(
      myc_outdegree_row
    ) == 0
  ) {

    warning(
      "MYC row not found in outdegree file. ",
      "Skipping total-network MYC immune analysis."
    )

  } else {

    if (
      nrow(
        myc_outdegree_row
      ) > 1
    ) {

      warning(
        "More than one MYC row found in ", "outdegree file. Using first row."
      )
    }

    sample_cols_outdegree <- setdiff(
      names(
        myc_outdegree_row
      ), "reg"
    )

    myc_total_outdegree <- as.numeric(

      myc_outdegree_row[
        1, ..sample_cols_outdegree
      ]

    )

    names(
      myc_total_outdegree
    ) <- sample_cols_outdegree

    myc_total_outdegree_dt <- data.table(

        sample_id =
          names(
            myc_total_outdegree
          ),

        sample_id_norm =
          normalize_sample_id(
            names(myc_total_outdegree)
          ),

        myc_total_outdegree =
          as.numeric(myc_total_outdegree)

      )

    myc_total_outdegree_dt <- merge(

      myc_total_outdegree_dt,

      clinical_sample_map[
        , .(
          sample_id_norm, response
        )
      ],

      by =
        "sample_id_norm",

      all.x = TRUE,

      sort = FALSE

    )

    myc_total_immune_dt <- merge(

      myc_total_outdegree_dt,

      immune_cell_dt,

      by =
        "sample_id_norm",

      all = FALSE,

      sort = FALSE

    )

    myc_total_immune_dt <- myc_total_immune_dt[
        !is.na(response) &
          !is.na(myc_total_outdegree)
      ]

    myc_total_immune_dt[
      , response := factor(
        response, levels = c(
          "resistance", "response"
        )
      )
    ]

    cat(
      "Samples in MYC total-outdegree immune analysis:", nrow(
        myc_total_immune_dt
      ), "\n"
    )

    total_cor_results <- vector(
      "list", length(immune_cell_cols)
    )

    names(
      total_cor_results
    ) <- immune_cell_cols

    if (
      !is.null(MYC_TOTAL_IMMUNE_PLOT)
    ) {

      pdf(
        MYC_TOTAL_IMMUNE_PLOT, width = 8,
        height = 6
      )
    }

    for (
      immune_cell in
      immune_cell_cols
    ) {

      immune_values <- suppressWarnings(
          as.numeric(myc_total_immune_dt[[immune_cell]])
        )

      cor_result <- safe_spearman(

        myc_total_immune_dt$
          myc_total_outdegree,

        immune_values

      )

      total_cor_results[[immune_cell]] <- data.table(

        immune_cell =
          immune_cell,

        n =
          cor_result$n,

        rho =
          cor_result$estimate,

        p_value =
          cor_result$p.value

      )

      if (
        is.null(MYC_TOTAL_IMMUNE_PLOT)
      ) {
        next
      }

      valid_idx <- which(

        is.finite(
          myc_total_immune_dt$
            myc_total_outdegree
        ) &

          is.finite(immune_values)

      )

      if (length(valid_idx) < 3) {
        next
      }

      plot_dt <- copy(
        myc_total_immune_dt[
          valid_idx
        ]
      )

      plot_dt[
        , immune_score :=
          immune_values[
            valid_idx
          ]
      ]

      p_total_immune_corr <- ggplot(

        plot_dt,

        aes(
          x = myc_total_outdegree, y = immune_score,
          color = response, fill = response
        )

      ) +

        geom_point(
          size = 2.4, alpha = 0.75
        ) +

        geom_smooth(
          method = "lm", se = TRUE,
          color = "black", alpha = 0.15
        ) +

        scale_color_manual(
          values =
            RESPONSE_COLORS
        ) +

        scale_fill_manual(
          values =
            RESPONSE_COLORS
        ) +

        labs(

          title = paste0(
            "MYC total outdegree vs ", immune_cell
          ),

          x =
            "MYC total outdegree (full network)",

          y = paste0(
            immune_cell, " score"
          ),

          subtitle =
            format_cor_subtitle(cor_result)

        ) +

        theme_bw() +

        theme(
          legend.position =
            "bottom"
        )

      print(p_total_immune_corr)
    }

    if (
      !is.null(MYC_TOTAL_IMMUNE_PLOT)
    ) {

      dev.off()

      cat(
        "MYC total-outdegree immune plots saved to:", MYC_TOTAL_IMMUNE_PLOT,
        "\n"
      )
    }

    total_cor_results <- rbindlist(
        total_cor_results, fill = TRUE
      )

    total_cor_results[
      , fdr_bh :=
        p.adjust(
          p_value, method = "BH"
        )
    ]

    if (
      !is.null(MYC_TOTAL_IMMUNE_TABLE)
    ) {

      fwrite(
        total_cor_results, MYC_TOTAL_IMMUNE_TABLE,
        sep = "\t"
      )

      cat(
        "MYC total-outdegree immune table saved to:", MYC_TOTAL_IMMUNE_TABLE,
        "\n"
      )
    }
  }
}

# ==========================================================
# Analysis 3:
# Selected TF PD-1 outdegree vs immune infiltration
# ==========================================================

if (
  !is.null(TF_PD1_IMMUNE_HEATMAP) || !is.null(TF_PD1_IMMUNE_TABLE)
) {

  cat(
    "Running selected TF PD-1 outdegree ",
    "vs immune infiltration analysis...\n"
  )

  selected_tfs <- c(
    "ATF6", "FOSL2",
    "CREB3L2", "MYC",
    "MYCN", "TFEC",
    "HEY1", "ZKSCAN3",
    "MAX", "HEY2",
    "ZNF354A", "ZNF384",
    "FOXP2", "FOXC2",
    "MEF2A", "STAT2",
    "MEF2D", "FOXO1",
    "SOX15", "POU3F1"
  )

  tf_outdegree_long <- rbindlist(

    lapply(

      selected_tfs,

      function(tf_name) {

        tf_idx <- which(
          edges_pathways$reg ==
            tf_name
        )

        if (length(tf_idx) == 0) {

          message(
            "No PD-1 pathway edges found for TF: ", tf_name
          )

          return(
            data.table(
              sample_id_norm =
                character(),
              TF =
                character(),
              outdegree =
                numeric()
            )
          )
        }

        tf_outdegree_vals <- colSums(

          net_pathways[
            tf_idx, ,
            drop = FALSE
          ],

          na.rm = TRUE

        )

        data.table(

          sample_id_norm =
            normalize_sample_id(
              names(tf_outdegree_vals)
            ),

          TF =
            tf_name,

          outdegree =
            as.numeric(tf_outdegree_vals)

        )
      }
    ),

    fill = TRUE
  )

  if (nrow(tf_outdegree_long) == 0) {

    warning(
      "No PD-1 pathway outdegree values ",
      "could be calculated for selected TFs."
    )

  } else {

    tf_outdegree_wide <- dcast(

      tf_outdegree_long,

      sample_id_norm ~ TF,

      value.var =
        "outdegree"

    )

    tf_outdegree_immune_dt <- merge(

      tf_outdegree_wide,

      immune_cell_dt,

      by =
        "sample_id_norm",

      all = FALSE,

      sort = FALSE

    )

    tf_cols_present <- intersect(
      selected_tfs, names(tf_outdegree_immune_dt)
    )

    immune_cols_present <- intersect(
      immune_cell_cols, names(tf_outdegree_immune_dt)
    )

    if (
      length(tf_cols_present) == 0 || length(immune_cols_present) == 0
    ) {

      warning(
        "No overlapping TF outdegree and ",
        "immune infiltration columns available."
      )

    } else {

      cor_results_dt <- rbindlist(

        lapply(

          immune_cols_present,

          function(immune_cell) {

            immune_values <- suppressWarnings(
                as.numeric(tf_outdegree_immune_dt[[immune_cell]])
              )

            rbindlist(

              lapply(

                tf_cols_present,

                function(tf_name) {

                  tf_values <- suppressWarnings(
                      as.numeric(tf_outdegree_immune_dt[[tf_name]])
                    )

                  cor_result <- safe_spearman(
                    tf_values, immune_values
                  )

                  data.table(

                    immune_cell =
                      immune_cell,

                    TF =
                      tf_name,

                    n =
                      cor_result$n,

                    rho =
                      cor_result$estimate,

                    p_value =
                      cor_result$p.value

                  )
                }
              ),

              fill = TRUE
            )
          }
        ),

        fill = TRUE
      )

      cor_results_dt[
        , fdr_bh :=
          p.adjust(
            p_value, method = "BH"
          )
      ]

      # Preserve the original plotting criterion:
      # |rho| >= 0.2 and raw p < 0.05.
      cor_results_dt[
        , sig_label :=
          fifelse(
            !is.na(p_value) &
              !is.na(rho) & abs(rho) >= 0.2 &
              p_value < 0.05,
            "*", ""
          )
      ]

      if (
        !is.null(TF_PD1_IMMUNE_TABLE)
      ) {

        fwrite(
          cor_results_dt, TF_PD1_IMMUNE_TABLE,
          sep = "\t"
        )

        cat(
          "TF PD-1 immune correlation table saved to:", TF_PD1_IMMUNE_TABLE,
          "\n"
        )
      }

      if (
        !is.null(TF_PD1_IMMUNE_HEATMAP)
      ) {

        plot_cor_results <- copy(cor_results_dt)

        plot_cor_results[
          , TF := factor(
            TF, levels =
              selected_tfs
          )
        ]

        plot_cor_results[
          , immune_cell := factor(
            immune_cell, levels =
              rev(immune_cols_present)
          )
        ]

        p_tf_immune_heatmap <- ggplot(

          plot_cor_results,

          aes(
            x = TF, y = immune_cell,
            fill = rho
          )

        ) +

          geom_tile(
            color = "white", linewidth = 0.2
          ) +

          geom_text(
            aes(
              label =
                sig_label
            ), size = 4
          ) +

          scale_fill_gradient2(

            low =
              "#2166AC",

            mid =
              "white",

            high =
              "#B2182B",

            midpoint =
              0,

            na.value =
              "grey90",

            name =
              "Spearman\nrho"

          ) +

          labs(

            title =
              "PD-1-pathway TF outdegree vs immune cell correlations",

            x =
              "TF",

            y =
              "Immune cell type"

          ) +

          theme_bw() +

          theme(

            axis.text.x =
              element_text(
                angle = 45, hjust = 1,
                vjust = 1
              ),

            axis.text.y =
              element_text(
                size = 8
              ),

            panel.grid =
              element_blank(),

            legend.position =
              "right"

          )

        pdf(
          TF_PD1_IMMUNE_HEATMAP, width = 11,
          height = 10
        )

        print(p_tf_immune_heatmap)

        dev.off()

        cat(
          "TF PD-1 immune correlation heatmap saved to:", TF_PD1_IMMUNE_HEATMAP,
          "\n"
        )
      }
    }
  }
}

# ==========================================================
# Analysis 4:
# MYC expression vs immune infiltration
# ==========================================================

if (
  !is.null(MYC_EXPRESSION_IMMUNE_PLOT) || !is.null(MYC_EXPRESSION_IMMUNE_TABLE)
) {

  cat(
    "Running MYC expression ", "vs immune infiltration analysis...\n"
  )

  if (
    !"MYC" %in%
      rownames(expression_clean)
  ) {

    warning(
      "MYC was not found in expression matrix. ",
      "Skipping MYC expression immune analysis."
    )

  } else {

    if (
      sum(
        rownames(
          expression_clean
        ) == "MYC"
      ) > 1
    ) {

      stop(
        "More than one MYC row found ", "in expression matrix."
      )
    }

    myc_expression_dt <- data.table(

      sample_id =
        colnames(
          expression_clean
        ),

      sample_id_norm =
        normalize_sample_id(
          colnames(expression_clean)
        ),

      myc_expression =
        as.numeric(
          expression_clean[
            "MYC", ,
            drop = TRUE
          ]
        )

    )

    myc_expression_dt <- merge(

      myc_expression_dt,

      clinical_sample_map[
        , .(
          sample_id_norm, response
        )
      ],

      by =
        "sample_id_norm",

      all.x = TRUE,

      sort = FALSE

    )

    myc_expression_immune_dt <- merge(

      myc_expression_dt,

      immune_cell_dt,

      by =
        "sample_id_norm",

      all = FALSE,

      sort = FALSE

    )

    myc_expression_immune_dt <- myc_expression_immune_dt[
        !is.na(response) &
          !is.na(myc_expression)
      ]

    myc_expression_immune_dt[
      , response := factor(
        response, levels = c(
          "resistance", "response"
        )
      )
    ]

    cat(
      "Samples in MYC expression immune analysis:", nrow(
        myc_expression_immune_dt
      ), "\n"
    )

    expression_cor_results <- vector(
      "list", length(immune_cell_cols)
    )

    names(
      expression_cor_results
    ) <- immune_cell_cols

    if (
      !is.null(MYC_EXPRESSION_IMMUNE_PLOT)
    ) {

      pdf(
        MYC_EXPRESSION_IMMUNE_PLOT, width = 8,
        height = 6
      )
    }

    for (
      immune_cell in
      immune_cell_cols
    ) {

      immune_values <- suppressWarnings(
          as.numeric(myc_expression_immune_dt[[immune_cell]])
        )

      cor_result <- safe_spearman(

        myc_expression_immune_dt$
          myc_expression,

        immune_values

      )

      expression_cor_results[[immune_cell]] <- data.table(

        immune_cell =
          immune_cell,

        n =
          cor_result$n,

        rho =
          cor_result$estimate,

        p_value =
          cor_result$p.value

      )

      if (
        is.null(MYC_EXPRESSION_IMMUNE_PLOT)
      ) {
        next
      }

      valid_idx <- which(

        is.finite(
          myc_expression_immune_dt$
            myc_expression
        ) &

          is.finite(immune_values)

      )

      if (length(valid_idx) < 3) {
        next
      }

      plot_dt <- copy(
        myc_expression_immune_dt[
          valid_idx
        ]
      )

      plot_dt[
        , immune_score :=
          immune_values[
            valid_idx
          ]
      ]

      p_expr_immune_corr <- ggplot(

        plot_dt,

        aes(
          x = myc_expression, y = immune_score,
          color = response, fill = response
        )

      ) +

        geom_point(
          size = 2.4, alpha = 0.75
        ) +

        geom_smooth(
          method = "lm", se = TRUE,
          color = "black", alpha = 0.15
        ) +

        scale_color_manual(
          values =
            RESPONSE_COLORS
        ) +

        scale_fill_manual(
          values =
            RESPONSE_COLORS
        ) +

        labs(

          title = paste0(
            "MYC expression vs ", immune_cell
          ),

          x =
            "MYC log2 expression",

          y = paste0(
            immune_cell, " score"
          ),

          subtitle =
            format_cor_subtitle(cor_result)

        ) +

        theme_bw() +

        theme(
          legend.position =
            "bottom"
        )

      print(p_expr_immune_corr)
    }

    if (
      !is.null(MYC_EXPRESSION_IMMUNE_PLOT)
    ) {

      dev.off()

      cat(
        "MYC expression immune plots saved to:", MYC_EXPRESSION_IMMUNE_PLOT,
        "\n"
      )
    }

    expression_cor_results <- rbindlist(
        expression_cor_results, fill = TRUE
      )

    expression_cor_results[
      , fdr_bh :=
        p.adjust(
          p_value, method = "BH"
        )
    ]

    if (
      !is.null(MYC_EXPRESSION_IMMUNE_TABLE)
    ) {

      fwrite(
        expression_cor_results, MYC_EXPRESSION_IMMUNE_TABLE,
        sep = "\t"
      )

      cat(
        "MYC expression immune table saved to:", MYC_EXPRESSION_IMMUNE_TABLE,
        "\n"
      )
    }
  }
}

cat("Immune infiltration correlation analysis completed.\n")
