#####################
## Load R packages ##
#####################
required_libraries <- c(
    "GEOquery", "SummarizedExperiment",
    "data.table", "dplyr", "optparse"
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
        help = "GEO Series ID (e.g., GSE30219) to download data for.",
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
        c("--probe_annotation_file"),
        type = "character",
        default = "probes.RData",
        help = "Path to probe annotation file (.RData)",
        metavar = "character"
    ),
    optparse::make_option(
        c("--clinical_file"),
        type = "character",
        default = "clinical.RData",
        help = "Path to clinical data file (.RData)",
        metavar = "character"
    )

)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

GSE_ID <- opt$GSE_id

# Set default file names based on GSE_ID if not provided
if (opt$expression_file == "expression.RData") {
    EXPRESSION_FILE <- paste0(GSE_ID, "_expression.RData")
} else {
    EXPRESSION_FILE <- opt$expression_file
}

if (opt$probe_annotation_file == "probes.RData") {
    PROBE_ANNOTATION_FILE <- paste0(GSE_ID, "_probes.RData")
} else {
    PROBE_ANNOTATION_FILE <- opt$probe_annotation_file
}

if (opt$clinical_file == "clinical.RData") {
    CLINICAL_FILE <- paste0(GSE_ID, "_clinical.RData")
} else {
    CLINICAL_FILE <- opt$clinical_file
}

## Download and process GEO data
############################

# Check if files already exist
all_files_exist <- 
    all(file.exists(c(EXPRESSION_FILE, 
                    PROBE_ANNOTATION_FILE, 
                    CLINICAL_FILE)))
if (all_files_exist) {
    print(paste("Data files for", GSE_ID, "already exist. Skipping download."))
    print(paste("Files found:"))
    print(paste("  Expression:", EXPRESSION_FILE))
    print(paste("  Probes:", PROBE_ANNOTATION_FILE))
    print(paste("  Clinical:", CLINICAL_FILE))
    quit(save = "no")
}

print(paste("Downloading data for", GSE_ID))
print(paste("Will save to:"))
print(paste("  Expression:", EXPRESSION_FILE))
print(paste("  Probes:", PROBE_ANNOTATION_FILE))
print(paste("  Clinical:", CLINICAL_FILE))  
# Download and parse the Series Matrix file(s)
gse <- getGEO(GSE_ID, GSEMatrix = TRUE)
eset <- gse[[1]]
expr <- exprs(eset)      # expression matrix
pheno <- pData(eset)     # sample metadata
feature <- fData(eset)   # feature/probe annotations
head(expr)
head(feature)
colnames(feature)
# colnames_to_select <- c("GENE_SYMBOL", "GENE_NAME", "UNIGENE_ID", "ID", "Gene Symbol", "Gene Title")  
# feature <- feature[, colnames(feature) %in% colnames_to_select, drop = FALSE]
head(feature)
head(pheno)
save(expr, file = EXPRESSION_FILE)
save(feature, file = PROBE_ANNOTATION_FILE)
save(pheno, file = CLINICAL_FILE)