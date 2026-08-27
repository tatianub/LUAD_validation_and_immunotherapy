#####################
## Load R packages ##
#####################
required_libraries <- c(
    "data.table", "dplyr",
    "ggplot2", "optparse",
    "ggrepel"
)
for (lib in required_libraries) {
    suppressPackageStartupMessages(
        library(lib, character.only = TRUE, quietly = TRUE)
    )
}
set.seed(123)

####################
## Read arguments ##
####################
option_list <- list(
    optparse::make_option(
        c("-e", "--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the expression file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-x", "--exp_clean"),
        type = "character",
        default = NULL,
        help = "Output path for cleaned expression file (.tsv).",
        metavar = "character"),
    optparse::make_option(
        c("-p", "--pca_plot"),
        type = "character",
        default = NULL,
        help = "Output path for PCA plot (PDF).",
        metavar = "character"),
    optparse::make_option(
        c("-s", "--samples_file"),
        type = "character",
        default = NULL,
        help = "Output path for samples file",
        metavar = "character")

)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source("workflow/bin/preprocess_expression_fn.R")
EXPRESSION_FILE <- opt$expression_file
OUTPUT_PCA_FILE <- opt$pca_plot
OUTPUT_EXPRESSION_CLEAN_FILE <- opt$exp_clean
SAMPLES_OUTPUT <- opt$samples_file

# loading the expression data
cat("Loading expression data...\n")
expr <- fread(EXPRESSION_FILE)
length(unique(expr$gene))
length(expr$gene)
genes <- expr$gene
expr <- expr[, grep("-T|-N", colnames(expr)), with = FALSE]
expr <- log2(expr + 1)
cat("Expression data dimensions:", dim(expr), "\n")
head(expr)

# Filter genes expressed in at least 5% of samples
n_samples <- ncol(expr)
min_samples <- ceiling(n_samples * 0.05)  # At least 5% of samples
expressed_samples <- apply(expr, 1, function(x) sum(x > 0))  # Count samples with expression > 0
genes_expressed_filter <- expressed_samples >= min_samples
cat("Filtering genes: keeping", sum(genes_expressed_filter), 
    "out of", length(genes_expressed_filter), 
    "genes that are expressed in at least", min_samples, 
    "samples (5%)\n")

genes_clean <- genes[genes_expressed_filter]
expr <- expr[genes_expressed_filter, ]

# Additional filtering: remove genes with zero variance and zero total expression
genes_clean <- genes_clean[apply(expr, 1, var) > 0]
expr <- expr[apply(expr, 1, var) > 0, ]
genes_clean <- genes_clean[rowSums(expr) > 0]
expr <- expr[rowSums(expr) > 0, ]
rownames(expr) <- genes_clean

expr <- round(expr, 3)
rownames(expr) <- genes_clean

pdf(OUTPUT_PCA_FILE, width = 8, height = 8)
create_pca_outlier_plot(expr, 
                            sample_names = colnames(expr),
                            outlier_method = "distance",
                            outlier_threshold = 4,
                            scale_data = TRUE,
                            center_data = TRUE,
                            plot_title = "PCA Plot with Outlier Detection",
                            color_palette = c("FALSE" = "black",
                                "TRUE" = "red"),
                            point_size = 3,
                            point_alpha = 0.7,
                            label_outliers = TRUE)
dev.off()

expr_clean <- expr
genes <- rownames(expr)
genes_clean <- genes[rowSums(expr_clean) > 0]
expr_clean <- expr_clean[rowSums(expr_clean) > 0, ]
rownames(expr_clean) <- genes_clean

# Save complete clean expression data for multiple cohorts too
write.table(expr_clean, OUTPUT_EXPRESSION_CLEAN_FILE,
                col.names = FALSE, sep = "\t", 
                row.names = TRUE, quote = FALSE)
cat("Complete expression data saved to:", OUTPUT_EXPRESSION_CLEAN_FILE, "\n")

# Also save complete samples file for multiple cohorts
samples <- data.frame(samples = colnames(expr_clean))
write.table(samples, SAMPLES_OUTPUT,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)
