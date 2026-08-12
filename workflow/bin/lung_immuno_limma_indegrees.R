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
    "--indegree_file",
    type = "character",
    default = NULL,
    help = "Path to indegree file.",
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
INDEGREE_FILE <- opt$indegree_file
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
  indegree_file = INDEGREE_FILE,
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
# Step 2: Load indegree matrix
# ========================================

indegree_dt <- fread(INDEGREE_FILE)

if (!"tar" %in% names(indegree_dt)) {
  stop("'tar' column is missing from indegree file.")
}

if (ncol(indegree_dt) < 2) {
  stop(
    "Indegree file must contain a target column plus sample columns."
  )
}

tars <- indegree_dt$tar

indegree <- data.matrix(
  indegree_dt[, -1, with = FALSE]
)

rownames(indegree) <- as.character(tars)

if (is.null(colnames(indegree))) {
  stop("Indegree matrix has no sample column names.")
}

if (anyDuplicated(colnames(indegree))) {
  stop("Duplicated sample IDs found in indegree matrix.")
}

# ========================================
# Step 3: Align indegree and clinical data
# ========================================
normalize_sample_id <- function(x) {
  gsub("[^A-Za-z0-9]", ".", as.character(x))
}

indegree_sample_ids_norm <- normalize_sample_id(colnames(indegree))
clinical_sample_ids_norm <- normalize_sample_id(clinical_data_filt$sample_id)

common_samples_norm <- indegree_sample_ids_norm[
  indegree_sample_ids_norm %in% clinical_sample_ids_norm
]

if (length(common_samples_norm) == 0) {
  stop("No common samples between indegree and clinical data.")
}

keep_indegree_cols <- indegree_sample_ids_norm %in% common_samples_norm

indegree_cl <- indegree[
  ,
  keep_indegree_cols,
  drop = FALSE
]

colnames(indegree_cl) <- indegree_sample_ids_norm[keep_indegree_cols]

clinical_data_filt[, sample_id_norm := clinical_sample_ids_norm]

clinical_data_ordered <- clinical_data_filt[
  match(colnames(indegree_cl), sample_id_norm)
]

if (!all(colnames(indegree_cl) == clinical_data_ordered$sample_id_norm)) {
  stop("Indegree and clinical samples are not aligned.")
}

response <- clinical_data_ordered$response

cat(
  "Samples used in differential analysis:",
  ncol(indegree_cl),
  "\n"
)

# ========================================
# Step 4: Run limma
# ========================================

if (!COVARIATES) {
  cat("Running limma analysis without covariates...\n")

  res_limma <- run_limma_analysis(
    data = indegree_cl,
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
    row.names = clinical_data_ordered$sample_id_norm
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

  indegree_cl_covariates <- indegree_cl[
    ,
    rownames(covariates),
    drop = FALSE
  ]

  if (!all(
    colnames(indegree_cl_covariates) == rownames(covariates)
  )) {
    stop(
      "Column names of indegree data do not match row names of ",
      "covariates."
    )
  }

  response_updated <- clinical_data_ordered$response[
    match(
      colnames(indegree_cl_covariates),
      clinical_data_ordered$sample_id_norm
    )
  ]

  if (anyNA(response_updated)) {
    stop("Missing response values after covariate filtering.")
  }

  cat(
    "Samples with complete covariates:",
    ncol(indegree_cl_covariates),
    "\n"
  )

  res_limma <- run_limma_analysis(
    data = indegree_cl_covariates,
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

write.table(
  res_fgsea,
  file = FGSEA_RESULTS_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

cat("Limma and fgsea analyses completed.\n")
