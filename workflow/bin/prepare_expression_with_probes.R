#####################
## Load R packages ##
#####################
required_libraries <- c(
    "data.table", "dplyr",
    "ggplot2", "optparse",
    "ggrepel", "readxl",
    "clusterProfiler",
    "org.Hs.eg.db",
    "hgu133plus2.db",
    "AnnotationDbi"
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
        c("--GSE_id"),
        type = "character",
        default = NULL,
        help = "GSE ID for the dataset",
        metavar = "character"
    ),
    optparse::make_option(
        c("--expression_file"),
        type = "character",
        default = "expression.RData",
        help = "Path to saved expression data file (.RData) ",
        metavar = "character"
    ),
    optparse::make_option(
        c("--probe_file"),
        type = "character",
        default = NULL,
        help = "Path to probe annotation file (.txt).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-x", "--exp_clean"),
        type = "character",
        default = NULL,
        help = "Output path for cleaned expression file (.tsv).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-p", "--pca_plot"),
        type = "character",
        default = NULL,
        help = "Output path for PCA plot (PDF).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Output directory for cohort-specific files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-s", "--samples_file"),
        type = "character",
        default = NULL,
        help = "Output path for samples file",
        metavar = "character"
    )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

EXPRESSION_FILE <- opt$expression_file
PROBE_FILE <- opt$probe_file
OUTPUT_EXPRESSION_CLEAN_FILE <- opt$exp_clean
OUTPUT_PCA_FILE <- opt$pca_plot
OUTPUT_DIR <- opt$output_dir
SAMPLES_OUTPUT <- opt$samples_file
GSE_ID <- opt$GSE_id

source("workflow/bin/preprocess_expression_fn.R")
source("workflow/bin/check_log_transformation_fn.R")

############################
## Load expression data   ##
############################
cat("Loading expression file...\n")
exp_env <- new.env()
load(EXPRESSION_FILE, envir = exp_env)
exp <- exp_env[["expr"]]
exp <- exp[complete.cases(exp), ]

# # Check current data
# exp <- check_log_transformation(exp)
# dim(exp)

annot <- fread(PROBE_FILE)

# keep only probes present in the expression matrix
annot2 <- annot[annot$PROBEID %in% rownames(exp),
               c("PROBEID", "SYMBOL")]
annot2 <- annot2[annot2$SYMBOL != "" & !is.na(annot2$SYMBOL), ]

# subset expression matrix to annotated probes
expr2 <- exp[annot2$PROBEID, , drop = FALSE]

# make sure order matches
annot2 <- annot2[match(rownames(expr2), annot2$PROBEID), ]
dim(expr2)

expr_gene_mean <-
    rowsum(expr2, group = annot2$SYMBOL) /
    as.vector(table(annot2$SYMBOL))
head(expr_gene_mean)

data_matrix <- as.matrix(expr_gene_mean)

############################
## Process expression data ##
############################

# Remove genes with zero variance
data_matrix <- data_matrix[
    apply(data_matrix, 1, var) > 0, , drop = FALSE]

data_matrix <- round(data_matrix, 3)

cat("Saving expression data to:", OUTPUT_EXPRESSION_CLEAN_FILE, "\n")
write.table(data_matrix, OUTPUT_EXPRESSION_CLEAN_FILE,
            col.names = FALSE, sep = "\t",
            row.names = TRUE, quote = FALSE)

samples <- data.frame(samples = colnames(data_matrix))
write.table(samples, SAMPLES_OUTPUT,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)

# Create PCA plot
pdf(OUTPUT_PCA_FILE, width = 10, height = 8)
cat("Creating PCA plot...\n")
create_pca_outlier_plot(data_matrix)
dev.off()
cat("PCA plot saved to:", OUTPUT_PCA_FILE, "\n")
