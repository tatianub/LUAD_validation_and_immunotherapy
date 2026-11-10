#####################
## Load R packages ##
#####################
required_libraries <- c(
    "data.table", "dplyr",
    "ggplot2", "optparse",
    "ggrepel", "readxl",
    "clusterProfiler",
    "org.Hs.eg.db"
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
        c("-o", "--output_pca_file"),
        type = "character",
        default = NULL,
        help = "Path to output PCA file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-m", "--output_expression_clean_file"),
        type = "character",
        default = NULL,
        help = "Path to output expression matrix file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--output_samples_file"),
        type = "character",
        default = NULL,
        help = "Path to output samples file.",
        metavar = "character"
    )
)


opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

source("workflow/bin/preprocess_expression_fn.R")
EXPRESSION_FILE <- opt$expression_file
OUTPUT_PCA_FILE <- opt$output_pca_file
OUTPUT_EXPRESSION_CLEAN_FILE <- 
    opt$output_expression_clean_file
OUTPUT_SAMPLES_FILE <- opt$output_samples_file


# loading the expression data
cat("Loading expression data...\n")
expr <- readRDS(EXPRESSION_FILE)
expr[1:5, 1:5]
genes <- 
    sapply(strsplit(as.character(expr$Description), 
            split = " - "), function(x) x[1])
genes_clean <- genes[genes != "NA" & !duplicated(genes)]
length(unique(genes_clean))
length(genes_clean)

expr_clean <- expr[genes != "NA"  & !duplicated(genes), -c(1:2)]

log_exp <- log2(expr_clean + 1)
rownames(log_exp) <- genes_clean
data_matrix <- as.matrix(log_exp)
head(data_matrix)

pdf(OUTPUT_PCA_FILE, width = 8, height = 8)
create_pca_outlier_plot(data_matrix, 
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
dim(data_matrix)
write.table(log_exp,  OUTPUT_EXPRESSION_CLEAN_FILE,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)
samples <- data.frame(samples = colnames(log_exp))
write.table(samples,
            OUTPUT_SAMPLES_FILE,
            col.names = FALSE,
            sep = "\t", 
            row.names = TRUE,
            quote = FALSE)

