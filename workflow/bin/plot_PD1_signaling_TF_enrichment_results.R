#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table", "dplyr", "ggplot2", "optparse", "ggrepel", "forcats"
)

for (lib in required_libraries) {
  suppressPackageStartupMessages(
    library(lib, character.only = TRUE, quietly = TRUE)
  )
}

## Read arguments ##
####################
option_list <- list(
  make_option(c("--results_enrichment_tfs"), type = "character", 
              default = NULL, help = "Path to TF enrichment results."),
  make_option(c("--results_enrichment_pathways"), type = "character", 
              default = NULL, help = "Path to pathway enrichment results."),
  make_option(c("--motif_prior"), type = "character", 
              default = NULL, help = "Path to motif prior file."),
  make_option(c("--output_file"), type = "character", 
              default = NULL, help = "Path to save the PDF plot."),
  make_option(c("--padj_threshold_pathways"), type = "numeric", 
              default = 0.05, help = "P-adj threshold for pathways."),
  make_option(c("--padj_threshold_tfs"), type = "numeric", 
              default = 0.01, help = "P-adj threshold for TFs."),
  make_option(c("--es_threshold_pathways"), type = "numeric", 
              default = 0.5, help = "ES threshold for pathways."),
  make_option(c("--es_threshold_tfs"), type = "numeric", 
              default = 0.5, help = "ES threshold for TFs.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Source utility functions
source("workflow/bin/lung_immuno_limma_fn.R")

########################
## Initialize values  ##
########################
ENRICHMENT_TFS <- opt$results_enrichment_tfs
ENRICHMENT_PATHWAYS <- opt$results_enrichment_pathways
MOTIF_PRIOR_FILE <- opt$motif_prior
OUTPUT_PLOT <- opt$output_file

# Thresholds from arguments
PADJ_PATH <- opt$padj_threshold_pathways
PADJ_TFS <- opt$padj_threshold_tfs
ES_PATH <- opt$es_threshold_pathways
ES_TFS <- opt$es_threshold_tfs

# 1. Process Motif Prior for CD274 (PD-L1)
motif_prior <- fread(MOTIF_PRIOR_FILE)
colnames(motif_prior) <- c("TF", "target", "direct_regulator")
motif_prior <- motif_prior[grep("^CD274$", motif_prior$target), ]

# 2. Process Pathway Results
res_pathways <- fread(ENRICHMENT_PATHWAYS)
res_pathways <- filter_fgsea_results(
  res_pathways, 
  padj_cutoff = PADJ_PATH, 
  es_threshold = ES_PATH
)
res_pathways <- res_pathways[order(abs(NES), decreasing = TRUE), ]

# 3. Process TF Results
res_tfs <- fread(ENRICHMENT_TFS)
res_tfs <- res_tfs[grep("PD_1", pathway), ]
res_tfs <- res_tfs %>% filter(!is.na(ES), !is.na(padj))

res_tfs_sig <- res_tfs %>%
  filter(padj < PADJ_TFS, abs(ES) > ES_TFS) %>%
  arrange(desc(ES))

# Define direct regulators
res_tfs_sig$direct_regulator <- ifelse(
  res_tfs_sig$TF %in% motif_prior$TF, 
  "Direct PD-L1 regulator",
  "Other TF"
)

# Keep top 10 and bottom 10 TFs
res_tfs_sig <- rbind(head(res_tfs_sig, 10), tail(res_tfs_sig, 10))

# Prepare data for lollipop plot
tfs_lollipop <- res_tfs_sig %>%
  mutate(
    neglog10_padj = -log10(padj),
    direction = ifelse(ES > 0, "Responder-enriched", "Non-responder-enriched"),
    TF = fct_reorder(TF, ES),
    direct_regulator = factor(
      direct_regulator,
      levels = c("Other TF", "Direct PD-L1 regulator")
    )
  )

# 4. Generate combined plot
plot_pathway_lollipop_ab(
  res_fgsea_filt = res_pathways,
  tfs_lollipop = tfs_lollipop,
  output_file = OUTPUT_PLOT,
  panel_layout = "row",
  width = 14,
  height = 6,
  pathway_row_height = 1.1,
  lollipop_row_height = 0.9,
  pathway_panel_width = 0.95,
  lollipop_panel_width = 0.95,
  label_a = "A.",
  label_b = "B.",
  tf_rank_by = "abs_es"
)