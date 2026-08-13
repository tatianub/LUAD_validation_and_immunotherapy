#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table",
  "EnsDb.Hsapiens.v86",
  "ggplot2",
  "ggpubr",
  "optparse"
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
    "--tpm_file",
    type = "character",
    default = NULL,
    help = "Path to harmonized RNA-seq TPM GCT file.",
    metavar = "character"
  ),
  optparse::make_option(
    "--clinical_file",
    type = "character",
    default = NULL,
    help = "Path to clinical metadata file.",
    metavar = "character"
  ),
  optparse::make_option(
    "--expression_output_file",
    type = "character",
    default = NULL,
    help = "Output TSV for processed gene expression.",
    metavar = "character"
  ),
  optparse::make_option(
    "--samples_output_file",
    type = "character",
    default = NULL,
    help = "Output TSV for sample IDs.",
    metavar = "character"
  ),
  optparse::make_option(
    "--pca_plot_file",
    type = "character",
    default = NULL,
    help = "Output PDF for PCA plots.",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(
  option_list = option_list
)

opt <- optparse::parse_args(opt_parser)

#######################
## Initialize values ##
#######################

DATA_TPM_FILE <- opt$tpm_file
CLINICAL_FILE <- opt$clinical_file
EXPRESSION_OUTPUT_FILE <- opt$expression_output_file
SAMPLES_OUTPUT_FILE <- opt$samples_output_file
PCA_PLOT_FILE <- opt$pca_plot_file

##################
## Input checks ##
##################

required_files <- list(
  tpm_file = DATA_TPM_FILE,
  clinical_file = CLINICAL_FILE
)

for (input_name in names(required_files)) {
  input_file <- required_files[[input_name]]

  if (is.null(input_file)) {
    stop("Missing required argument --", input_name)
  }

  if (!file.exists(input_file)) {
    stop("File does not exist: ", input_file)
  }
}

required_outputs <- list(
  expression_output_file = EXPRESSION_OUTPUT_FILE,
  samples_output_file = SAMPLES_OUTPUT_FILE,
  pca_plot_file = PCA_PLOT_FILE
)

for (output_name in names(required_outputs)) {
  if (is.null(required_outputs[[output_name]])) {
    stop("Missing required argument --", output_name)
  }
}

# ========================================
# Step 1: Load TPM and clinical data
# ========================================

data_tpm <- fread(
  DATA_TPM_FILE,
  skip = 2,
  sep = "\t",
  header = TRUE
)

clinical_data <- fread(
  CLINICAL_FILE
)

if (!all(c("Name", "Description") %in% names(data_tpm))) {
  stop(
    "TPM GCT file must contain 'Name' and 'Description' columns."
  )
}

if (!"sample_id" %in% names(clinical_data)) {
  stop("'sample_id' column is missing from clinical data.")
}

if (!"response" %in% names(clinical_data)) {
  stop("'response' column is missing from clinical data.")
}

# ========================================
# Step 2: Clean Ensembl gene IDs
# ========================================

original_ids <- as.character(
  data_tpm$Name
)

ensembl_ids_clean <- sub(
  "\\..*$",
  "",
  original_ids
)

# ========================================
# Step 3: Fetch protein-coding annotation
# ========================================

genes_info <- genes(
  EnsDb.Hsapiens.v86,
  filter = GeneIdFilter(ensembl_ids_clean),
  return.type = "data.frame"
)

protein_coding_ids <- genes_info[
  genes_info$gene_biotype == "protein_coding",
  ,
  drop = FALSE
]

if (nrow(protein_coding_ids) == 0) {
  stop(
    "No protein-coding genes were found in the annotation."
  )
}

# ========================================
# Step 4: Prepare TPM matrix
# ========================================

data_tpm_matrix <- data.matrix(
  data_tpm[
    ,
    !c("Name", "Description"),
    with = FALSE
  ]
)

rownames(data_tpm_matrix) <- original_ids

if (anyDuplicated(colnames(data_tpm_matrix))) {
  stop(
    "Duplicated sample IDs found in TPM matrix."
  )
}

# ========================================
# Step 5: Transform and filter expression
# ========================================

data_tpm_log <- log2(
  data_tpm_matrix + 1
)

data_tpm_log <- data_tpm_log[
  rowSums(
    data_tpm_log,
    na.rm = TRUE
  ) > 0,
  ,
  drop = FALSE
]

gene_fraction <- rowMeans(
  data_tpm_log > 0,
  na.rm = TRUE
)

data_tpm_filtered <- data_tpm_log[
  gene_fraction >= 0.1,
  ,
  drop = FALSE
]

if (nrow(data_tpm_filtered) == 0) {
  stop(
    "No genes remain after expression filtering."
  )
}

# ========================================
# Step 6: Clean IDs and collapse duplicates
# ========================================

clean_rows <- sub(
  "\\..*$",
  "",
  rownames(data_tpm_filtered)
)

if (anyDuplicated(clean_rows)) {
  message(
    paste(
      "Duplicates found after stripping versions.",
      "Aggregating by mean."
    )
  )

  data_tpm_filtered <- rowsum(
    data_tpm_filtered,
    group = clean_rows,
    reorder = FALSE
  )

  duplicate_counts <- table(
    clean_rows
  )

  data_tpm_filtered <- data_tpm_filtered /
    as.numeric(
      duplicate_counts[
        rownames(data_tpm_filtered)
      ]
    )
} else {
  rownames(data_tpm_filtered) <- clean_rows
}

# ========================================
# Step 7: Map Ensembl IDs to gene symbols
# ========================================

annot2 <- protein_coding_ids[
  protein_coding_ids$gene_id %in%
    rownames(data_tpm_filtered),
  ,
  drop = FALSE
]

annot2 <- annot2[
  !is.na(annot2$gene_name) &
    annot2$gene_name != "",
  ,
  drop = FALSE
]

annot2 <- annot2[
  !duplicated(
    annot2[
      ,
      c("gene_id", "gene_name"),
      drop = FALSE
    ]
  ),
  ,
  drop = FALSE
]

if (nrow(annot2) == 0) {
  stop(
    paste(
      "No annotated protein-coding genes remain",
      "after filtering."
    )
  )
}

expr2 <- data_tpm_filtered[
  annot2$gene_id,
  ,
  drop = FALSE
]

# ========================================
# Step 8: Collapse expression to gene symbols
# ========================================

# rowsum() uses reorder = TRUE by default.
# This reproduces the alphabetical gene-symbol ordering
# of the original preprocessing script.

expr_gene_sum <- rowsum(
  expr2,
  group = annot2$gene_name
)

gene_counts <- table(
  annot2$gene_name
)

expr_gene_mean <- expr_gene_sum /
  as.numeric(
    gene_counts[
      rownames(expr_gene_sum)
    ]
  )

expr_gene_mean <- round(
  expr_gene_mean,
  3
)

if (anyDuplicated(rownames(expr_gene_mean))) {
  stop(
    "Duplicated gene symbols remain after collapsing expression."
  )
}

# ========================================
# Step 9: Run PCA
# ========================================

pca <- prcomp(
  t(expr_gene_mean),
  scale. = TRUE
)

percent_var <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Sample = rownames(pca$x)
)

pca_df$response <- clinical_data$response[
  match(
    pca_df$Sample,
    clinical_data$sample_id
  )
]

# ========================================
# Step 10: Create PCA plots
# ========================================

pca_subtitle <- paste0(
  "Based on ",
  nrow(expr_gene_mean),
  " protein-coding genes"
)

pca_x_label <- paste0(
  "PC1 (",
  round(percent_var[1], 1),
  "% variance)"
)

pca_y_label <- paste0(
  "PC2 (",
  round(percent_var[2], 1),
  "% variance)"
)

pca_caption <- paste0(
  "Data: log2(TPM + 1) | ",
  "Filter: expressed in >=10% of samples"
)

g1 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    size = 4,
    alpha = 0.8,
    color = "steelblue"
  ) +
  labs(
    title = "Principal Component Analysis",
    subtitle = pca_subtitle,
    x = pca_x_label,
    y = pca_y_label,
    caption = pca_caption
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

g2 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = response
  )
) +
  geom_point(
    size = 4,
    alpha = 0.8
  ) +
  labs(
    title = "Principal Component Analysis by response",
    subtitle = pca_subtitle,
    x = pca_x_label,
    y = pca_y_label,
    caption = pca_caption,
    color = "Response"
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

pdf(
  PCA_PLOT_FILE,
  width = 10,
  height = 10
)

print(
  ggpubr::ggarrange(
    g1,
    g2,
    ncol = 1,
    nrow = 2
  )
)

dev.off()

# ========================================
# Step 11: Write processed expression data
# ========================================

samples <- data.frame(
  samples = colnames(expr_gene_mean)
)

write.table(
  expr_gene_mean,
  file = EXPRESSION_OUTPUT_FILE,
  col.names = FALSE,
  sep = "\t",
  row.names = TRUE,
  quote = FALSE
)

write.table(
  samples,
  file = SAMPLES_OUTPUT_FILE,
  col.names = FALSE,
  sep = "\t",
  row.names = TRUE,
  quote = FALSE
)

cat(
  "Processed expression matrix:",
  nrow(expr_gene_mean),
  "genes x",
  ncol(expr_gene_mean),
  "samples\n"
)

cat(
  "PCA plot saved to:",
  PCA_PLOT_FILE,
  "\n"
)

cat(
  "Expression preprocessing completed.\n"
)