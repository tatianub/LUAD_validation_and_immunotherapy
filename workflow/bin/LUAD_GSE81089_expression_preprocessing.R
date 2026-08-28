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
        c("--expression_file"),
        type = "character",
        default = NULL,
        help = "Path to the expression file.",
        metavar = "character"
    ),
    optparse::make_option(
         c("--annotation_file"),
        type = "character",
        default = NULL,
        help = "Path to the annotation file.",
        metavar = "character"
    ),
    optparse::make_option(
         c("--samples_file"),
        type = "character",
        default = NULL,
        help = "Path to the samples file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("--output_pca_file"),
        type = "character",
        default = NULL,
        help = "Path to output PCA file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("--output_expression_clean_file"),
        type = "character",
        default = NULL,
        help = "Path to output expression matrix file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("--output_samples_file"),
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
ANNOTATION_FILE <- opt$annotation_file
OUTPUT_PCA_FILE <- opt$output_pca_file
OUTPUT_EXPRESSION_CLEAN_FILE <- opt$output_expression_clean_file
OUTPUT_SAMPLES_FILE <- opt$output_samples_file
SAMPLE_ANNOTATION_FILE <- opt$samples_file

# loading the expression data
cat("Loading expression data...\n")
expr <- fread(EXPRESSION_FILE)
annotation <- fread(ANNOTATION_FILE)
sample_annotation <- fread(SAMPLE_ANNOTATION_FILE, head = F)
colnames(expr)[-1] <- 
    sample_annotation$V2[match(colnames(expr)[-1], sample_annotation$V1)]
expr <- expr[, -grep("matched", colnames(expr)), with = FALSE]
dim(expr)
head(expr)
min(expr[, -1, with = FALSE])
max(expr[, -1, with = FALSE])
ids <- expr$GeneID
table(ids %in% annotation$GeneID)
genes <- annotation$Symbol[match(ids, annotation$GeneID)]
length(genes)
length(unique(genes))
genes[duplicated(genes) | duplicated(genes, fromLast = TRUE)]
#  "TRNAV-CAC" "TRNAV-CAC" "TRNAV-CAC"

expr <- as.matrix(expr[, -1, with = FALSE])
expr <- log2(expr + 1)
expr <- rowsum(expr, group = genes) / as.vector(table(genes))

cat("Expression data dimensions:", dim(expr), "\n")
genes <- rownames(expr)
genes <- genes[apply(expr, 1, var) > 0]
expr <- expr[apply(expr, 1, var) > 0, ]
genes <- genes[rowSums(expr) > 0]
expr <- expr[rowSums(expr) > 0, ]
rownames(expr) <- genes

data_matrix <- expr
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
dev.off()
log_exp <- expr
save(log_exp, file = OUTPUT_EXPRESSION_CLEAN_FILE)
samples <- data.frame(samples = colnames(log_exp))
write.table(samples, OUTPUT_SAMPLES_FILE,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)

