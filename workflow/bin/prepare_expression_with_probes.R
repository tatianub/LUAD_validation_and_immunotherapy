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
        c("-c", "--cohorts_file"),
        type = "character", 
        default = NULL,
        help = "Path to cohorts file (.txt) with columns GSM and cohort.",
        metavar = "character"),
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
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Output directory for cohort-specific files.",
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

EXPRESSION_FILE <- opt$expression_file
COHORTS_FILE <- opt$cohorts_file
PROBE_FILE <- opt$probe_file
OUTPUT_EXPRESSION_CLEAN_FILE <- opt$exp_clean
OUTPUT_PCA_FILE <- opt$pca_plot
OUTPUT_DIR <- opt$output_dir
SAMPLES_OUTPUT <- opt$samples_file
GSE_ID <- opt$GSE_id

source("workflow/bin/preprocess_expression_fn.R")
source("workflow/bin/check_log_transformation_fn.R")
############################
## Load expression data      ##
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

# Read in the cohorts file
cohorts <- fread(COHORTS_FILE)
head(cohorts)

# Check cohort composition
cat("\n=== Cohort Analysis ===\n")
cohort_summary <- table(cohorts$cohort)
print(cohort_summary)

# Check if single cohort or multiple cohorts
unique_cohorts <- unique(cohorts$cohort)

if (length(unique_cohorts) == 1) {
    cat("Single cohort detected:", unique_cohorts[1], "\n")
    cat("Total samples:", nrow(cohorts), "\n")
    # For single cohort, just save the main expression file
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

    # Create single PCA plot
    pdf(OUTPUT_PCA_FILE, width = 10, height = 8)
    cat("Creating PCA plot for single cohort...\n")
    create_pca_outlier_plot(data_matrix)
    dev.off()
    cat("PCA plot saved to:", OUTPUT_PCA_FILE, "\n")
    
    } else {
    # Multiple cohorts - original processing
    cat("Multiple cohorts detected:", length(unique_cohorts), "\n")
    cat("Total samples:", nrow(cohorts), "\n")
    for (cohort_name in names(cohort_summary)) {
        cat(paste("Cohort", cohort_name, ":", 
                  cohort_summary[cohort_name], "samples\n"))
    }
    covariate_cohort <- cohorts$cohort[match(colnames(data_matrix), 
                                         cohorts$GSM)]


    pdf(OUTPUT_PCA_FILE, width = 10, height = 8)

    # 1. PCA plot colored by cohort
    cat("\nCreating PCA plot colored by cohort...\n")
    create_pca_outlier_plot(data_matrix, covariate = covariate_cohort)

    # Extract cohort-specific data
    unique_cohorts <- unique(cohorts$cohort)
    cat("\nProcessing individual cohorts...\n")
    cat("Found", length(unique_cohorts), 
    "cohorts:", paste(unique_cohorts, collapse = ", "), "\n")

    # Store output files for summary
    cohort_output_files <- c()

    for (cohort_name in unique_cohorts) {
        cat(paste("Processing", cohort_name, "...\n"))
        # Get samples for this cohort
        cohort_samples <- cohorts$GSM[cohorts$cohort == cohort_name]
        # Subset data matrix
        data_matrix_cohort <- data_matrix[, cohort_samples, drop = FALSE]

        # Remove genes with zero variance
        data_matrix_cohort <- data_matrix_cohort[
            apply(data_matrix_cohort, 1, var) > 0, , drop = FALSE]
        cat(paste("  Cohort", cohort_name, ":", 
              ncol(data_matrix_cohort), "samples,", 
              nrow(data_matrix_cohort), "genes\n"))
        # Create PCA plot for this cohort
        create_pca_outlier_plot(data_matrix_cohort)
        samples <- data.frame(samples = colnames(data_matrix_cohort))
        # Generate output filename dynamically for any cohort
        output_file <- paste0(cohort_name, "_", GSE_ID, "_expression.tsv")
        output_file <- file.path(OUTPUT_DIR, output_file)
        samples_file <- paste0(cohort_name, "_", GSE_ID, "_samples.tsv")
        samples_file <- file.path(OUTPUT_DIR, samples_file)

        write.table(data_matrix_cohort, output_file,
                col.names = FALSE, sep = "\t", 
                row.names = TRUE, quote = FALSE)
        write.table(samples, samples_file,
                col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)
        cat(paste("  Saved", cohort_name, "expression to:", 
              output_file, "\n"))
        cat(paste("  Saved", cohort_name, "samples to:", 
              samples_file, "\n"))
    }

    dev.off()
    cat("PCA plots saved to:", OUTPUT_PCA_FILE, "\n")
    
    # Save complete clean expression data for multiple cohorts too
    write.table(data_matrix, OUTPUT_EXPRESSION_CLEAN_FILE,
                col.names = FALSE, sep = "\t", 
                row.names = TRUE, quote = FALSE)
    cat("Complete expression data saved to:", OUTPUT_EXPRESSION_CLEAN_FILE, "\n")
    
    # Also save complete samples file for multiple cohorts
    samples <- data.frame(samples = colnames(data_matrix))
    write.table(samples, SAMPLES_OUTPUT,
            col.names = FALSE, sep = "\t", row.names = TRUE, quote = FALSE)
    cat("Complete samples data saved to:", SAMPLES_OUTPUT, "\n")
}
