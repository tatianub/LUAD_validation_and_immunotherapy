#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "optparse",
  "survival",
  "ggpubr",
  "grid",
  "tools",
  "limma",
  "fgsea"
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
    "--expression_file",
    type = "character",
    default = NULL,
    help = "Path to expression file.",
    metavar = "character"
  ),
  optparse::make_option(
    "--samples_file",
    type = "character",
    default = NULL,
    help = "Path to samples file.",
    metavar = "character"
  ),
  optparse::make_option(
    "--clinical_file",
    type = "character",
    default = NULL,
    help = "Path to clinical file.",
    metavar = "character"
  ),
  optparse::make_option(
    "--subtype_type",
    type = "character",
    default = "all",
    help = "Type of histological subtype.",
    metavar = "character"
  ),
  optparse::make_option(
    "--treatment_type",
    type = "character",
    default = "all",
    help = "Type of treatment.",
    metavar = "character"
  ),
  optparse::make_option(
    "--covariates",
    type = "logical",
    default = TRUE,
    help = "Whether to include covariates in the analysis.",
    metavar = "logical"
  ),
  optparse::make_option(
    "--pathway_gmt_file",
    type = "character",
    default = NULL,
    help = "Path to pathway GMT file (.gmt).",
    metavar = "character"
  ),
  optparse::make_option(
    "--limma_results_file",
    type = "character",
    default = NULL,
    help = "Path to limma results file (.txt).",
    metavar = "character"
  ),
  optparse::make_option(
    "--fgsea_results_file",
    type = "character",
    default = NULL,
    help = "Path to fgsea results file (.txt).",
    metavar = "character"
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = 2026,
    help = "Random seed for reproducibility [default %default].",
    metavar = "integer"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#######################
## Initialize values ##
#######################

CLINICAL_FILE <- opt$clinical_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
GMT_FILE <- opt$pathway_gmt_file
HISTO_SUBTYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
COVARIATES <- opt$covariates
LIMMA_RESULTS_FILE <- opt$limma_results_file
FGSEA_RESULTS_FILE <- opt$fgsea_results_file
SEED <- opt$seed

if (!is.null(SEED) && !is.na(SEED)) {
  set.seed(SEED)
  cat("Using random seed:", SEED, "\n")
}

##################
## Input checks ##
##################

required_files <- list(
  clinical_file = CLINICAL_FILE,
  expression_file = EXPRESSION_FILE,
  samples_file = SAMPLES_FILE,
  pathway_gmt_file = GMT_FILE
)

for (input_name in names(required_files)) {
  input_file <- required_files[[input_name]]

  if (is.null(input_file)) {
    stop("Missing required argument --", input_name)
  }

  if (!file.exists(input_file)) {
    stop("File does not exist: ", input_file)
  }
}

required_outputs <- list(
  limma_results_file = LIMMA_RESULTS_FILE,
  fgsea_results_file = FGSEA_RESULTS_FILE
)

for (output_name in names(required_outputs)) {
  if (is.null(required_outputs[[output_name]])) {
    stop("Missing required argument --", output_name)
  }
}

source("workflow/bin/lung_immuno_limma_fn.R")

# ========================================
# Step 1: Load and filter clinical data
# ========================================

clinical_data <- fread(CLINICAL_FILE)

if (!"sample_id" %in% names(clinical_data)) {
  stop("'sample_id' column is missing from clinical data.")
}

if (!"response" %in% names(clinical_data)) {
  stop("'response' column is missing from clinical data.")
}

clinical_data <- clinical_data[!is.na(response)]

cat(
  "Filtering clinical data for histological subtype:",
  HISTO_SUBTYPE,
  "and treatment type:",
  TREATMENT_TYPE,
  "\n"
)

clinical_data_filt <- subset_clinical_data(
  clinical_data = clinical_data,
  subtype = HISTO_SUBTYPE,
  treatment_type = TREATMENT_TYPE
)

if (nrow(clinical_data_filt) == 0) {
  stop("No samples remain after clinical filtering.")
}

if (anyDuplicated(clinical_data_filt$sample_id)) {
  stop("Duplicated sample IDs found in filtered clinical data.")
}

# ========================================
# Step 2: Load expression matrix
# ========================================

expression_dt <- fread(EXPRESSION_FILE)

if (ncol(expression_dt) < 2) {
  stop(
    "Expression file must contain a gene column plus sample columns."
  )
}

tars <- expression_dt[[1]]

expression <- data.matrix(
  expression_dt[, -1, with = FALSE]
)

samples <- fread(SAMPLES_FILE)

if (ncol(samples) < 2) {
  stop(
    "Samples file must contain at least two columns; sample IDs are ",
    "expected in column 2."
  )
}

if (ncol(expression) != nrow(samples)) {
  stop(
    "Expression matrix contains ",
    ncol(expression),
    " sample columns, while samples file contains ",
    nrow(samples),
    " rows."
  )
}

colnames(expression) <- as.character(samples[[2]])
rownames(expression) <- as.character(tars)

if (anyDuplicated(colnames(expression))) {
  stop("Duplicated sample IDs found in expression matrix.")
}
colnames(expression) <- gsub("-", ".", colnames(expression))
# ========================================
# Step 3: Align expression and clinical data
# ========================================

common_samples <- colnames(expression)[
  colnames(expression) %in% clinical_data_filt$sample_id
]

if (length(common_samples) == 0) {
  stop("No common samples between expression and clinical data.")
}

expression_cl <- expression[
  ,
  common_samples,
  drop = FALSE
]

clinical_data_ordered <- clinical_data_filt[
  match(common_samples, sample_id)
]

if (!all(colnames(expression_cl) == clinical_data_ordered$sample_id)) {
  stop("Expression and clinical samples are not aligned.")
}

response <- clinical_data_ordered$response

cat(
  "Samples used in differential analysis:",
  length(common_samples),
  "\n"
)

# ========================================
# Step 4: Run limma
# ========================================

if (!COVARIATES) {
  cat("Running limma analysis without covariates...\n")

  res_limma <- run_limma_analysis(
    data = expression_cl,
    response_vector = response
  )
} else {
  cat("Running limma analysis with covariates...\n")

  required_covariates <- c(
    "Patient_Age_at_Diagnosis",
    "Patient_Sex",
    "Patient_Smoking_Pack_Years_Harmonized"
  )

  missing_covariates <- setdiff(
    required_covariates,
    names(clinical_data_ordered)
  )

  if (length(missing_covariates) > 0) {
    stop(
      "Missing clinical covariate column(s): ",
      paste(missing_covariates, collapse = ", ")
    )
  }

  covariates <- data.frame(
    age = suppressWarnings(
      as.numeric(clinical_data_ordered$Patient_Age_at_Diagnosis)
    ),
    sex = factor(clinical_data_ordered$Patient_Sex),
    smoking = suppressWarnings(
      as.numeric(
        clinical_data_ordered$Patient_Smoking_Pack_Years_Harmonized
      )
    ),
    row.names = clinical_data_ordered$sample_id
  )

  complete_covariates <- complete.cases(covariates)

  if (sum(complete_covariates) == 0) {
    stop("No samples have complete covariate information.")
  }

  covariates <- covariates[
    complete_covariates,
    ,
    drop = FALSE
  ]

  expression_cl_covariates <- expression_cl[
    ,
    rownames(covariates),
    drop = FALSE
  ]

  if (!all(
    colnames(expression_cl_covariates) == rownames(covariates)
  )) {
    stop(
      "Column names of expression data do not match row names of ",
      "covariates."
    )
  }

  response_updated <- clinical_data_ordered$response[
    match(
      colnames(expression_cl_covariates),
      clinical_data_ordered$sample_id
    )
  ]

  if (anyNA(response_updated)) {
    stop("Missing response values after covariate filtering.")
  }

  cat(
    "Samples with complete covariates:",
    ncol(expression_cl_covariates),
    "\n"
  )

  res_limma <- run_limma_analysis(
    data = expression_cl_covariates,
    response_vector = response_updated,
    covariates = covariates
  )
}

# ========================================
# Step 5: Run fgsea
# ========================================

res_fgsea <- run_fgsea(
  res_all = res_limma,
  gmt_file = GMT_FILE
)

# ========================================
# Step 6: Write results
# ========================================

write.table(
  res_limma,
  file = LIMMA_RESULTS_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

save_fgsea_results(res_fgsea, FGSEA_RESULTS_FILE)

cat("Limma and fgsea analyses completed.\n")
