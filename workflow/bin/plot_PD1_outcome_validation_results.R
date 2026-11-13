#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggrepel",
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
    c("-g", "--result_files"),
    type = "character",
    default = NULL,
    help = paste0(
      "Comma-separated paths ",
      "to GSEA univariate Cox result files."
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("-o", "--validation_plot_file"),
    type = "character",
    default = NULL,
    help = "Path to output validation plot PDF.",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


#######################
## Initialize values ##
#######################
RESULT_FILES <- opt$result_files
VALIDATION_PLOT_FILE <- opt$validation_plot_file

if (is.null(RESULT_FILES)) stop("--result_files is required")
if (is.null(VALIDATION_PLOT_FILE)) {
  stop("--validation_plot_file is required")
}


#######################
## Read result files ##
#######################
# split the comma-separated string into individual file paths
result_files <- unlist(strsplit(RESULT_FILES, ","))
result_files <- trimws(result_files)

# stop if any of the input files do not exist
missing_files <- result_files[!file.exists(result_files)]
if (length(missing_files) > 0) {
  stop("Missing result files:\n", paste(missing_files, collapse = "\n"))
}

# keep full path and file name (used as dataset label)
result_files <- data.frame(
  original = result_files,
  dataset = basename(result_files),
  stringsAsFactors = FALSE
)

# read each result file and stack them into one table
res_all <- NULL
for (i in seq_len(nrow(result_files))) {
  res <- data.table::fread(result_files$original[i])
  res$dataset <- result_files$dataset[i]
  res_all <- dplyr::bind_rows(res_all, res)
}

# strip the file suffix to get the dataset name
res_all$dataset <- gsub(
  "_gsea_univariate_cox_results_norm\\.txt$",
  "",
  res_all$dataset
)


#############################
## Prepare plotting data   ##
#############################
# create plotting dataframe using raw p-values
plot_df <- res_all %>%
  mutate(
    neglog10_pval = -log10(pval),
    NES_plot = ifelse(pval > 0.1, NA, NES),
    outcome = gsub("5yr", "5y", outcome)
  )

# rename datasets for display
plot_df$dataset <- gsub("LUAD_Chen", "OncoSG", plot_df$dataset)
plot_df$dataset <- gsub("LUAD_Okayama", "GSE31210", plot_df$dataset)
plot_df$dataset <- gsub("LUAD_GSE81089", "GSE81089", plot_df$dataset)
plot_df$dataset <- gsub("LUAD_Dolgalev", "GSE229705", plot_df$dataset)

# prepare data: cap size and order datasets by median significance
plot_df <- plot_df %>%
  mutate(neglog10_pval_cap = pmin(neglog10_pval, 8))

dataset_order <- plot_df %>%
  group_by(dataset) %>%
  summarise(med = median(neglog10_pval, na.rm = TRUE)) %>%
  arrange(desc(med)) %>%
  pull(dataset)

plot_df <- plot_df %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

first_ds <- levels(plot_df$dataset)[1]


##################
## Create plot  ##
##################
# plotting with pale yellow box for first dataset and highlight points
p <- ggplot(plot_df, aes(x = dataset, y = neglog10_pval)) +
  annotate(
    "rect",
    xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
    fill = "#d3d3d3", alpha = 0.65
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.1),
    linetype = "dotted",
    color = "grey70"
  ) +
  annotate(
    "text",
    x = Inf, y = -log10(0.05), label = "p = 0.05",
    hjust = 1.05, vjust = -0.3, size = 3.2, color = "grey50"
  ) +
  annotate(
    "text",
    x = Inf, y = -log10(0.1), label = "p = 0.1",
    hjust = 1.05, vjust = -0.3, size = 3.2, color = "grey70"
  ) +
  # background points (all datasets)
  geom_point(
    aes(fill = NES_plot, size = neglog10_pval_cap),
    colour = "black", shape = 21, stroke = 0.35, alpha = 0.9
  ) +
  # overlay: emphasize the first dataset points (thicker border)
  geom_point(
    data = subset(plot_df, dataset == first_ds),
    aes(fill = NES_plot, size = neglog10_pval_cap),
    colour = "black", shape = 21, stroke = 1.0, alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, pval < 0.1),
    aes(label = outcome),
    size = 3.2,
    max.overlaps = 20
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#f7f7f7", high = "#B2182B",
    midpoint = 0, na.value = "grey80", name = "NES"
  ) +
  scale_size_continuous(
    name = expression(-log[10](p)),
    range = c(2.5, 8)
  ) +
  labs(x = NULL, y = expression(-log[10](p))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


###################################
## Meta-analysis (Stouffer)      ##
###################################
# meta-analysis (Stouffer) across datasets per outcome
meta_stouffer <- function(NES, p, weights = NULL) {
  ok <- is.finite(NES) & is.finite(p) & p > 0 & p <= 1
  NES <- NES[ok]; p <- p[ok]
  z <- sign(NES) * qnorm(1 - p / 2)
  if (is.null(weights)) {
    combined_z <- sum(z) / sqrt(length(z))
  } else {
    weights <- weights[ok]
    combined_z <- sum(weights * z) / sqrt(sum(weights^2))
  }
  combined_p <- 2 * pnorm(-abs(combined_z))
  data.frame(
    n_datasets = length(z),
    combined_z = combined_z,
    combined_p = combined_p,
    direction = ifelse(combined_z > 0, "positive", "negative")
  )
}

# exclude TCGA for meta-analysis
res_all <- res_all[-grep("TCGA_LUAD", res_all$dataset), ]

# compute meta p per outcome and print
meta_by_outcome <- res_all %>%
  group_by(outcome) %>%
  summarise(meta = list(meta_stouffer(NES = NES, p = pval))) %>%
  tidyr::unnest_wider(meta)

print(meta_by_outcome)

# optional: extract OS meta and add to plot later
os_meta <- meta_by_outcome %>% filter(outcome == "OS")
# os_meta$combined_p is the Stouffer p-value for OS


#################################
## Annotate plot and save PDF  ##
#################################
# annotate the plot with OS meta p
p <- p + annotate(
  "text",
  x = Inf, y = Inf,
  label = paste0("OS meta p=", signif(os_meta$combined_p, 2)),
  hjust = 1.1, vjust = 2, size = 3.2
)

ggsave(
  VALIDATION_PLOT_FILE,
  p,
  width = 8,
  height = 6,
  device = cairo_pdf
)