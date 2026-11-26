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
nes_limit <- ceiling(max(abs(plot_df$NES_plot), na.rm = TRUE))
p <- ggplot(plot_df, aes(x = dataset, y = neglog10_pval)) +

  # Highlight TCGA-LUAD with a light gray background
  annotate(
    "rect",
    xmin = 0.5, xmax = 1.5,
    ymin = -Inf, ymax = Inf,
    fill = "#F0F0F0",
    alpha = 1
  ) +

  # Significance threshold: P = 0.05
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.5
  ) +

  # Significance threshold: P = 0.1
  geom_hline(
    yintercept = -log10(0.1),
    linetype = "dotted",
    colour = "grey65",
    linewidth = 0.5
  ) +

  # Annotate significance thresholds
  annotate(
    "text",
    x = Inf,
    y = -log10(0.05),
    label = "italic(P) == 0.05",
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.4,
    size = 3.2,
    colour = "grey45"
  ) +

  annotate(
    "text",
    x = Inf,
    y = -log10(0.1),
    label = "italic(P) == 0.1",
    parse = TRUE,
    hjust = 1.05,
    vjust = -0.4,
    size = 3.2,
    colour = "grey65"
  ) +

  # Points for all datasets
  geom_point(
    aes(
      fill = NES_plot,
      size = neglog10_pval_cap
    ),
    colour = "black",
    shape = 21,
    stroke = 0.35,
    alpha = 0.9
  ) +

  # Emphasize TCGA-LUAD points with thicker borders
  geom_point(
    data = subset(plot_df, dataset == first_ds),
    aes(
      fill = NES_plot,
      size = neglog10_pval_cap
    ),
    colour = "black",
    shape = 21,
    stroke = 1.0,
    alpha = 1
  ) +

# Label only endpoints with P < 0.1
ggrepel::geom_text_repel(
  data = subset(plot_df, pval < 0.1),
  aes(
    label = outcome,
    point.size = neglog10_pval_cap
  ),
  size = 3.2,
  colour = "black",
  seed = 123,
  box.padding = 0.5,
  point.padding = 0.5,
  force = 3,
  force_pull = 0.5,
  max.overlaps = Inf,
  max.iter = 50000,
  max.time = 5,
  min.segment.length = 0,
  segment.colour = "grey50",
  segment.size = 0.25,
  show.legend = FALSE
) +
scale_size_continuous(
  name = expression(-log[10](italic(P))),
  range = c(2.5, 8),
  aesthetics = c("size", "point.size")
) +
  # NES color scale
  # Negative NES: red
  # Positive NES: blue
  # Masked NES: gray
scale_fill_gradient2(
  low = "#2166AC",
  mid = "#F7F7F7",
  high = "#B2182B",
  midpoint = 0,
  limits = c(-nes_limit, nes_limit),
  breaks = c(-nes_limit, 0, nes_limit),
  labels = c(
    paste0(-nes_limit, " (Better\noutcome)"),
    "0",
    paste0(nes_limit, " (Worse\noutcome)")
  ),
  na.value = "grey80",
  name = "NES"
  ) +
  # Additional space for labels and meta-analysis annotation
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.12))
  ) +

  labs(
    x = NULL,
    y = expression(-log[10](italic(P)))
  ) +

  coord_cartesian(clip = "off") +

  # Clean publication-style theme
  theme_bw(base_size = 13) +

  theme(
    # Axis labels
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black",
      size = 11
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 11
    ),

    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 8)
    ),

    # Minimal horizontal gridlines
    panel.grid.major.y = element_line(
      colour = "grey92",
      linewidth = 0.3
    ),

    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    # Remove full panel border
    panel.border = element_blank(),

    # Retain x and y axis lines
    axis.line = element_line(
      colour = "black",
      linewidth = 0.4
    ),

    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),

    # Legends
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.spacing.y = grid::unit(0.3, "cm"),

    # Plot margins
    plot.margin = margin(
      t = 10, r = 15,
      b = 10, l = 10
    )
  )


###################################
## Meta-analysis (Stouffer)      ##
###################################

# Meta-analysis across datasets per outcome
meta_stouffer <- function(NES, p, weights = NULL) {

  ok <- is.finite(NES) &
        is.finite(p) &
        p > 0 &
        p <= 1

  NES <- NES[ok]
  p <- p[ok]

  z <- sign(NES) * qnorm(1 - p / 2)

  if (is.null(weights)) {

    combined_z <- sum(z) / sqrt(length(z))

  } else {

    weights <- weights[ok]

    combined_z <- sum(weights * z) /
                  sqrt(sum(weights^2))

  }

  combined_p <- 2 * pnorm(-abs(combined_z))

  data.frame(
    n_datasets = length(z),
    combined_z = combined_z,
    combined_p = combined_p,
    direction = ifelse(
      combined_z > 0,
      "positive",
      "negative"
    )
  )
}


# Exclude TCGA from the meta-analysis
# without modifying the original results object
res_validation <- res_all %>%
  filter(!grepl("TCGA_LUAD", dataset))


# Compute meta-analysis P-values per outcome
meta_by_outcome <- res_validation %>%
  group_by(outcome) %>%
  summarise(
    meta = list(
      meta_stouffer(
        NES = NES,
        p = pval
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(meta)

print(meta_by_outcome)


# Extract OS meta-analysis results
os_meta <- meta_by_outcome %>%
  filter(outcome == "OS")


#################################
## Annotate plot and save PDF  ##
#################################

# Format meta-analysis P-value in scientific notation
meta_p <- os_meta$combined_p

exponent <- floor(log10(meta_p))
mantissa <- signif(meta_p / 10^exponent, 2)

meta_label <- sprintf(
  "'OS meta-analysis:'~italic(P)==%g%%*%%10^{%d}",
  mantissa,
  exponent
)


# Add OS meta-analysis P-value in the upper-right corner
p <- p + annotate(
  "text",
  x = Inf,
  y = Inf,
  label = meta_label,
  parse = TRUE,
  hjust = 1.05,
  vjust = 1.5,
  size = 3.5,
  colour = "black"
)


# Display plot
print(p)


# Save PDF
ggsave(
  VALIDATION_PLOT_FILE,
  p,
  width = 8,
  height = 6,
  device = cairo_pdf
)