#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table",
  "dplyr",
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
    c("--fgsea_results_file"),
    type = "character",
    default = NULL,
    help = "Path to merged fgsea results file (all models).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--pathway"),
    type = "character",
    default = "REACTOME_PD_1_SIGNALING",
    help = "Pathway name to extract from fgsea results.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--output_table"),
    type = "character",
    default = NULL,
    help = "Path to output summary table (.txt).",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

FGSEA_RESULTS_FILE <- opt$fgsea_results_file
PATHWAY <- opt$pathway
OUTPUT_TABLE <- opt$output_table

if (is.null(FGSEA_RESULTS_FILE)) {
  stop("--fgsea_results_file is required.")
}

if (is.null(OUTPUT_TABLE)) {
  stop("--output_table is required.")
}

####################
## Model metadata ##
####################
model_meta <- data.frame(
  model = c(
    "M0_response_only_all_samples",
    "M1_age_sex_smoking",
    "M2_age_sex_smoking_PDL1",
    "S2a_response_only_41_subset",
    "S2b_response_TMB_purity"
  ),
  label = c(
    "M0_response_only_all_samples",
    "M1_age_sex_smoking",
    "M2_age_sex_smoking_PDL1",
    "S2a_response_only_41_subset",
    "S2b_response_TMB_purity"
  ),
  samples = c(
    "all samples",
    "all with age/sex/smoking",
    "all with age/sex/smoking/PDL1",
    "41 pts (TMB+purity available)",
    "41 pts (TMB+purity available)"
  ),
  formula = c(
    "~ response",
    "~ response + age + sex + smoking",
    "~ response + age + sex + smoking + PDL1_expression",
    "~ response",
    "~ response + TMB + purity"
  ),
  stringsAsFactors = FALSE
)


######################
## Read fgsea table ##
######################
fgsea <- fread(FGSEA_RESULTS_FILE)

pathway_res <- fgsea[fgsea$pathway == PATHWAY, ]

if (nrow(pathway_res) == 0) {
  stop(
    "Pathway '", PATHWAY, "' not found in ", FGSEA_RESULTS_FILE, ".\n",
    "Available pathways (first 20): ",
    paste(head(unique(fgsea$pathway), 20), collapse = ", ")
  )
}

pathway_res <- pathway_res[, .(
  model,
  n_samples,
  n_responders,
  n_non_responders,
  NES,
  padj
)]

summary_table <- merge(
  model_meta,
  pathway_res,
  by = "model",
  all.x = TRUE
)

summary_table <- summary_table[
  match(model_meta$model, summary_table$model), ]

summary_table <- summary_table[, c(
  "label",
  "samples",
  "formula",
  "n_samples",
  "n_responders",
  "n_non_responders",
  "NES",
  "padj"
)]

colnames(summary_table)[1] <- "model"

cat("\nSummary table for pathway:", PATHWAY, "\n\n")
print(summary_table, row.names = FALSE)

write.table(
  summary_table,
  file = OUTPUT_TABLE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSaved to:", OUTPUT_TABLE, "\n")
