#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table",
  "estimate",
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
    c("--expr_file"),
    type = "character",
    default = NULL,
    help = "Path to expression file for ESTIMATE.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--samples_file"),
    type = "character",
    default = NULL,
    help = "Path to samples file for expression.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--output_expression_gct_file"),
    type = "character",
    default = NULL,
    help = "Path to output expression GCT file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--output_expression_gct_file_filtered"),
    type = "character",
    default = NULL,
    help = "Path to output filtered expression GCT file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--output_estimate_scores_gct_file"),
    type = "character",
    default = NULL,
    help = "Path to output ESTIMATE scores GCT file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--output_estimate_scores_tsv_file_transposed"),
    type = "character",
    default = NULL,
    help = "Path to output transposed ESTIMATE scores file.",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

EXPRESSION_FILE <- opt$expr_file
SAMPLES_FILE <- opt$samples_file
OUTPUT_EXPRESSION_GCT_FILE <- opt$output_expression_gct_file
OUTPUT_EXPRESSION_GCT_FILE_FILTERED <- opt$output_expression_gct_file_filtered
OUTPUT_ESTIMATE_SCORES_GCT_FILE <- opt$output_estimate_scores_gct_file
OUTPUT_ESTIMATE_SCORES_TSV_FILE_TRANSPOSED <- opt$output_estimate_scores_tsv_file_transposed


##########################
## Read expression data ##
##########################

expr <- fread(EXPRESSION_FILE)
samples <- fread(SAMPLES_FILE)
genes <- expr[[1L]]
expr <- as.data.frame(expr[, -1, with = FALSE])
rownames(expr) <- genes
head(expr)
colnames(expr) <- samples$V2
write.table(
  cbind(NAME = rownames(expr), Description = rownames(expr), expr),
  file = OUTPUT_EXPRESSION_GCT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
##########################
## Filter common genes ##
##########################
filterCommonGenes(
  input.f = OUTPUT_EXPRESSION_GCT_FILE,
  output.f = OUTPUT_EXPRESSION_GCT_FILE_FILTERED,
  id = "GeneSymbol"
)
##########################
## Estimate scores ##
##########################
estimateScore(
  OUTPUT_EXPRESSION_GCT_FILE_FILTERED,
  OUTPUT_ESTIMATE_SCORES_GCT_FILE,
  platform = "illumina"
)

scores <- read.table(
  OUTPUT_ESTIMATE_SCORES_GCT_FILE,
  skip = 2,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove Description columns
scores <- scores[, !grepl("^Description", colnames(scores))]

# Transpose
scores <- as.data.frame(t(scores))

# Convert all columns to numeric
scores[] <- lapply(scores, as.numeric)
head(scores)
write.table(
  scores,
  file = OUTPUT_ESTIMATE_SCORES_TSV_FILE_TRANSPOSED,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)
