#####################
## Load R packages ##
#####################
required_libraries <- c(
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
  optparse::make_option("--indegree_file",
    type = "character", default = NULL,
    help = "Path to indegree/GTS file."
  ),
  optparse::make_option("--clinical_file",
    type = "character", default = NULL,
    help = "Path to clinical file."
  ),
  optparse::make_option("--expression_file",
    type = "character", default = NULL,
    help = "Path to expression file."
  ),
  optparse::make_option("--samples_file",
    type = "character", default = NULL,
    help = "Path to samples file."
  ),
  optparse::make_option("--immune_estimate",
    type = "character", default = NULL,
    help = "Path to ESTIMATE scores."
  ),
  optparse::make_option("--subtype_type",
    type = "character", default = "all",
    help = "Histological subtype filter."
  ),
  optparse::make_option("--treatment_type",
    type = "character", default = "all",
    help = "Treatment type filter."
  ),
  optparse::make_option("--fgsea_results",
    type = "character", default = "all",
    help = "Path to fgsea results file."
  ),

  optparse::make_option("--results_logistic_models",
    type = "character", default = NULL,
    help = "Output Table S3 TSV."
  )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

##########################
## Assign input options ##
##########################
INDEGREE_FILE <- opt$indegree_file
CLINICAL_FILE <- opt$clinical_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
IMMUNE_ESTIMATE_FILE <- opt$immune_estimate
SUBTYPE_TYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
FGSEA_RESULTS <- opt$fgsea_results
RESULTS_LOGISTIC_MODELS <- opt$results_logistic_models

##################
## Input checks ##
##################
required_files <- list(
  indegree_file = INDEGREE_FILE,
  clinical_file = CLINICAL_FILE,
  expression_file = EXPRESSION_FILE,
  samples_file = SAMPLES_FILE,
  immune_estimate = IMMUNE_ESTIMATE_FILE,
  fgsea_results = FGSEA_RESULTS
)
for (input_name in names(required_files)) {
  input_file <- required_files[[input_name]]
  if (is.null(input_file) || !nzchar(input_file)) {
    stop("Missing required argument --", input_name)
  }
  if (!file.exists(input_file)) {
    stop("File not found: ", input_file)
  }
}
if (is.null(RESULTS_LOGISTIC_MODELS) ||
    !nzchar(RESULTS_LOGISTIC_MODELS)) {
  stop("Missing required argument --results_logistic_models")
}
dir.create(
  dirname(RESULTS_LOGISTIC_MODELS),
  recursive = TRUE,
  showWarnings = FALSE
)

#####################
## Global settings ##
#####################
# Retain original clinical filtering function and sample ID
# normalization.
source("workflow/bin/lung_immuno_limma_fn.R")
normalize_sample_id <- function(x) {
  gsub("-", ".", trimws(as.character(x)))
}
# ========================================
# Step 1: Load CD274 expression
# ========================================
expression <- fread(EXPRESSION_FILE)
genes <- expression[[1]]
expression <- as.matrix(
  expression[, -1, with = FALSE]
)
samples <- fread(SAMPLES_FILE)
stopifnot(ncol(expression) == nrow(samples))
colnames(expression) <- normalize_sample_id(samples$V2)
rownames(expression) <- genes
stopifnot(
  sum(rownames(expression) == "CD274") == 1,
  !anyDuplicated(colnames(expression))
)
pdl1_data <- data.frame(
  sample_id = colnames(expression),
  PDL1_expression = as.numeric(expression["CD274", ])
)

# Extract leading-edge genes from FGSEA results
fgsea_results <- fread(FGSEA_RESULTS)
pd1_row <- fgsea_results[pathway == "REACTOME_PD_1_SIGNALING"]
if (nrow(pd1_row) != 1) {
  stop("Could not uniquely identify REACTOME_PD_1_SIGNALING in FGSEA results.")
}
PD1_LEADING_EDGE <- unlist(strsplit(pd1_row$leadingEdge, ";"))
# ========================================
# Step 2: Load ESTIMATE ImmuneScore
# ========================================
immune_estimate <- fread(IMMUNE_ESTIMATE_FILE)
setnames(immune_estimate, 1, "sample_id")
immune_data <- data.frame(
  sample_id = normalize_sample_id(immune_estimate$sample_id),
  ImmuneScore = as.numeric(immune_estimate$ImmuneScore)
)
stopifnot(!anyDuplicated(immune_data$sample_id))

# ========================================
# Step 3: Load clinical cohort and response labels
# ========================================
clinical_data <- fread(CLINICAL_FILE)
clinical_data <- clinical_data[!is.na(response)]
clinical_data <- subset_clinical_data(
  clinical_data = clinical_data,
  subtype = SUBTYPE_TYPE,
  treatment_type = TREATMENT_TYPE
)
clinical_data$sample_id <- normalize_sample_id(
  clinical_data$sample_id
)
stopifnot(!anyDuplicated(clinical_data$sample_id))
if (!all(clinical_data$response %in% c("response", "resistance"))) {
  stop("Unexpected response labels")
}

# ========================================
# Step 4: Sum GTS over PD-1 leading-edge genes
# ========================================
indegree <- fread(INDEGREE_FILE)
gts_genes <- indegree$tar
indegree <- as.matrix(
  indegree[, -1, with = FALSE]
)
rownames(indegree) <- gts_genes
colnames(indegree) <- normalize_sample_id(colnames(indegree))
stopifnot(!anyDuplicated(colnames(indegree)))
if (!all(PD1_LEADING_EDGE %in% rownames(indegree))) {
  stop(
    "Missing leading-edge genes: ",
    paste(
      setdiff(PD1_LEADING_EDGE, rownames(indegree)),
      collapse = ", "
    )
  )
}
gts <- indegree[PD1_LEADING_EDGE, , drop = FALSE]
if (anyNA(gts)) {
  stop("Missing GTS values")
}
gts_sum <- colSums(gts)
gts_data <- data.frame(
  sample_id = names(gts_sum),
  PD1_LE_GTS_SUM = as.numeric(gts_sum)
)

# ========================================
# Step 5: Build identical complete-case cohort
# ========================================
d <- data.frame(
  sample_id = clinical_data$sample_id,
  response = clinical_data$response
) %>%
  left_join(gts_data, by = "sample_id") %>%
  left_join(pdl1_data, by = "sample_id") %>%
  left_join(immune_data, by = "sample_id")
stopifnot(nrow(d) == nrow(clinical_data))
d$response_binary <- as.integer(d$response == "response")
model_vars <- c(
  "response_binary",
  "PD1_LE_GTS_SUM",
  "PDL1_expression",
  "ImmuneScore"
)
d <- d[complete.cases(d[, model_vars]), , drop = FALSE]
if (nrow(d) < 15 || length(unique(d$response_binary)) != 2) {
  stop("Insufficient complete cases or only one response class")
}
scale_checked <- function(x) {
  if (any(!is.finite(x)) || !is.finite(sd(x)) || sd(x) == 0) {
    stop("Non-finite values or zero variance in predictor")
  }
  as.numeric(scale(x))
}
d$gts_z <- scale_checked(d$PD1_LE_GTS_SUM)
d$pdl1_z <- scale_checked(d$PDL1_expression)
d$immune_z <- scale_checked(d$ImmuneScore)
cat(
  "Patients:", nrow(d),
  "| Responders:", sum(d$response_binary == 1),
  "| Non-responders:", sum(d$response_binary == 0),
  "\n"
)

# ========================================
# Step 6: Fit the four models for Table S3
# ========================================
models <- list(
  "PD-1 pathway regulatory score" =
    response_binary ~ gts_z,
  "PD-L1 expression" =
    response_binary ~ gts_z + pdl1_z,
  "ESTIMATE ImmuneScore" =
    response_binary ~ gts_z + immune_z,
  "PD-L1 expression + ESTIMATE ImmuneScore" =
    response_binary ~ gts_z + pdl1_z + immune_z
)
extract_gts <- function(formula, label) {
  fit <- glm(formula, data = d, family = binomial())
  cf <- summary(fit)$coefficients["gts_z", ]
  beta <- unname(cf["Estimate"])
  se <- unname(cf["Std. Error"])
  data.frame(
    Model = label,
    n = nobs(fit),
    OR = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    P = unname(cf["Pr(>|z|)"]),
    check.names = FALSE
  )
}
results <- do.call(
  rbind,
  Map(extract_gts, models, names(models))
)
rownames(results) <- NULL
results$`OR (95% CI)` <- sprintf(
  "%.2f (%.2f–%.2f)",
  results$OR, results$CI_low, results$CI_high
)
results$P_formatted <- sprintf("%.3f", results$P)
print(
  results[, c("Model", "n", "OR (95% CI)", "P_formatted")],
  row.names = FALSE
)

# ========================================
# Step 7: Save results
# ========================================
fwrite(results, RESULTS_LOGISTIC_MODELS, sep = "\t")
cat("Saved:", RESULTS_LOGISTIC_MODELS, "\n")
