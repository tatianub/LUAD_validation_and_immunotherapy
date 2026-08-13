#####################
## Load R packages ##
#####################
required_libraries <- c(
  "readxl",
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
    c("--indegree_file"),
    type = "character",
    default = NULL,
    help = "Path to indegree file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--clinical_file"),
    type = "character",
    default = NULL,
    help = "Path to clinical file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--subtype_type"),
    type = "character",
    default = "all",
    help = "Type of subtype.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--treatment_type"),
    type = "character",
    default = "all",
    help = "Type of treatment.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--pathway_gmt_file"),
    type = "character",
    default = NULL,
    help = "Path to pathway GMT file (.gmt).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--limma_results_file"),
    type = "character",
    default = NULL,
    help = "Path to limma results file (.txt).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--fgsea_results_file"),
    type = "character",
    default = NULL,
    help = "Path to fgsea results file (.txt).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--expression_file"),
    type = "character",
    default = NULL,
    help = "Path to expression file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--samples_file"),
    type = "character",
    default = NULL,
    help = "Path to samples file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--clinical_file_extended"),
    type = "character",
    default = NULL,
    help = "Path to extended clinical file (xlsx).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--seed"),
    type = "integer",
    default = 2026,
    help = "Random seed for reproducibility [default %default].",
    metavar = "integer"
  )
  )

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


########################
## Initialize values  ##
########################

CLINICAL_FILE <- opt$clinical_file
INDEGREE_FILE <- opt$indegree_file
GMT_FILE <- opt$pathway_gmt_file
HISTO_SUBTYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
LIMMA_RESULTS_FILE <- opt$limma_results_file
FGSEA_RESULTS_FILE <- opt$fgsea_results_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
CLINICAL_FILE_EXTENDED <- opt$clinical_file_extended
SEED <- opt$seed

if (!is.null(SEED) && !is.na(SEED)) {
  set.seed(SEED)
  cat("Using random seed:", SEED, "\n")
}

source("workflow/bin/lung_immuno_limma_fn.R")


# read in the expression file
expression <- fread(EXPRESSION_FILE)
head(expression)
tars <- expression$V1
expression <- as.matrix(expression[, -1])

samples <- fread(SAMPLES_FILE)
colnames(expression) <- samples$V2
rownames(expression) <- tars

exp_pdl1 <- expression[grep("CD274", rownames(expression), value = TRUE), ]
exp_pdl1 <- data.table("sample_id" = names(exp_pdl1),
            "PDL1_expression" = as.numeric(exp_pdl1))
exp_pdl1$sample_id <- gsub("-", ".", exp_pdl1$sample_id)
# mutat
mutation_burden <- read_excel(
  CLINICAL_FILE_EXTENDED,
  sheet = "Table_S5_Mutation_Burden"
)
mutation_burden$sample_id <- 
    mutation_burden$Harmonized_SU2C_WES_Tumor_Sample_ID_v2
mutation_burden$sample_id <- gsub("-", ".", mutation_burden$sample_id) 


purity_data <- read_excel(
  CLINICAL_FILE_EXTENDED,
  sheet = "Table_S4_Purity_and_Ploidy"
)
purity_data$sample_id <- 
    purity_data$Harmonized_SU2C_WES_Tumor_Sample_ID_v2
purity_data$sample_id <- gsub("-", ".", purity_data$sample_id)

# read in clinical file and filter to only samples with response data

clinical_data <- fread(CLINICAL_FILE)
clinical_data <- clinical_data[!is.na(clinical_data$response)]

# read in the indegree file
indegree <- fread(INDEGREE_FILE)
tars <- indegree$tar
indegree <- as.matrix(indegree[, -1])

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

colnames(indegree) <- gsub("-", ".", colnames(indegree))
indegree_cl <- indegree[, colnames(indegree) %in% clinical_data_filt$sample_id]
dim(indegree_cl)
rownames(indegree_cl) <- tars
response <- clinical_data_filt$response[
  match(colnames(indegree_cl), clinical_data_filt$sample_id)]

# Baseline model on all available samples: response only (no covariates)
cat("Model 0: ~ response (all samples, no covariates)\n")
res_limma_m0 <- run_limma_analysis(
  data = indegree_cl,
  response_vector = response
)
res_limma_m0$model <- "M0_response_only_all_samples"
res_fgsea_m0 <- run_fgsea(res_all = res_limma_m0, gmt_file = GMT_FILE)
if (!is.null(res_fgsea_m0)) {
  res_fgsea_m0$model            <- "M0_response_only_all_samples"
  res_fgsea_m0$n_samples        <- length(response)
  res_fgsea_m0$n_responders     <- sum(response == "response")
  res_fgsea_m0$n_non_responders <- sum(response == "resistance")
}

cat("Running limma analysis with covariates...\n")

  # Build full covariates table aligned to indegree_cl samples
  clinical_data_ordered <- clinical_data_filt[
    match(colnames(indegree_cl), clinical_data_filt$sample_id), ]

  covariates_full <- data.frame(
    age             = as.numeric(clinical_data_ordered$Patient_Age_at_Diagnosis),
    sex             = factor(clinical_data_ordered$Patient_Sex),
    smoking         = as.numeric(
      clinical_data_ordered$Patient_Smoking_Pack_Years_Harmonized),
    pdl1_expression = exp_pdl1$PDL1_expression[
      match(clinical_data_ordered$sample_id, exp_pdl1$sample_id)],
    tmb_status      = mutation_burden$TMB[
      match(clinical_data_ordered$sample_id, mutation_burden$sample_id)],
    purity          = purity_data$Purity[
      match(clinical_data_ordered$sample_id, purity_data$sample_id)],
    row.names       = colnames(indegree_cl)
  )

  ## ---- Model 1: ~ response + age + sex + smoking -------------------------
  cat("Model 1: ~ response + age + sex + smoking\n")
  cov_m1 <- covariates_full[
    !is.na(covariates_full$age) &
    !is.na(covariates_full$sex) &
    !is.na(covariates_full$smoking),
    c("age", "sex", "smoking")]
  indegree_m1 <- indegree_cl[, rownames(cov_m1), drop = FALSE]
  rownames(indegree_m1) <- tars
  dim(indegree_m1)
  resp_m1 <- clinical_data_filt$response[
    match(rownames(cov_m1), clinical_data_filt$sample_id)]
  res_limma_m1 <- run_limma_analysis(
    data = indegree_m1,
    response_vector = resp_m1,
    covariates = cov_m1
  )
  res_limma_m1$model <- "M1_age_sex_smoking"
  res_fgsea_m1 <- run_fgsea(res_all = res_limma_m1, gmt_file = GMT_FILE)
  if (!is.null(res_fgsea_m1)) {
    res_fgsea_m1$model          <- "M1_age_sex_smoking"
    res_fgsea_m1$n_samples      <- length(resp_m1)
    res_fgsea_m1$n_responders   <- sum(resp_m1 == "response")
    res_fgsea_m1$n_non_responders <- sum(resp_m1 == "resistance")
  }

  ## ---- Model 2: ~ response + age + sex + smoking + PDL1_expression -------
  cat("Model 2: ~ response + age + sex + smoking + PDL1_expression\n")
  cov_m2 <- covariates_full[
    !is.na(covariates_full$age) &
    !is.na(covariates_full$sex) &
    !is.na(covariates_full$smoking) &
    !is.na(covariates_full$pdl1_expression),
    c("age", "sex", "smoking", "pdl1_expression")]
  indegree_m2 <- indegree_cl[, rownames(cov_m2), drop = FALSE]
  rownames(indegree_m2) <- tars
  resp_m2 <- clinical_data_filt$response[
    match(rownames(cov_m2), clinical_data_filt$sample_id)]
  res_limma_m2 <- run_limma_analysis(
    data = indegree_m2,
    response_vector = resp_m2,
    covariates = cov_m2
  )
  res_limma_m2$model <- "M2_age_sex_smoking_PDL1"
  res_fgsea_m2 <- run_fgsea(res_all = res_limma_m2, gmt_file = GMT_FILE)
  if (!is.null(res_fgsea_m2)) {
    res_fgsea_m2$model          <- "M2_age_sex_smoking_PDL1"
    res_fgsea_m2$n_samples      <- length(resp_m2)
    res_fgsea_m2$n_responders   <- sum(resp_m2 == "response")
    res_fgsea_m2$n_non_responders <- sum(resp_m2 == "resistance")
  }

  ## ---- Sensitivity 2: 41 patients (TMB + purity available) ---------------
  cov_s2_full <- covariates_full[
    !is.na(covariates_full$tmb_status) &
    !is.na(covariates_full$purity), ]
  indegree_s2 <- indegree_cl[, rownames(cov_s2_full), drop = FALSE]
  rownames(indegree_s2) <- tars
  resp_s2 <- clinical_data_filt$response[
    match(rownames(cov_s2_full), clinical_data_filt$sample_id)]

  # Sensitivity 2a (41-sample subset): ~ response only
  cat("Sensitivity 2a (41-sample subset): ~ response (no covariates)\n")
  res_limma_s2a <- run_limma_analysis(
    data = indegree_s2,
    response_vector = resp_s2
  )
  res_limma_s2a$model <- "S2a_response_only_41_subset"
  res_fgsea_s2a <- run_fgsea(res_all = res_limma_s2a, gmt_file = GMT_FILE)
  if (!is.null(res_fgsea_s2a)) {
    res_fgsea_s2a$model          <- "S2a_response_only_41_subset"
    res_fgsea_s2a$n_samples      <- length(resp_s2)
    res_fgsea_s2a$n_responders   <- sum(resp_s2 == "response")
    res_fgsea_s2a$n_non_responders <- sum(resp_s2 == "resistance")
  }

  # Sensitivity 2b: ~ response + TMB + purity
  cat("Sensitivity 2b: ~ response + TMB + purity\n")
  cov_s2b <- cov_s2_full[, c("tmb_status", "purity")]
  res_limma_s2b <- run_limma_analysis(
    data = indegree_s2,
    response_vector = resp_s2,
    covariates = cov_s2b
  )
  res_limma_s2b$model <- "S2b_response_TMB_purity"
  res_fgsea_s2b <- run_fgsea(res_all = res_limma_s2b, gmt_file = GMT_FILE)
  if (!is.null(res_fgsea_s2b)) {
    res_fgsea_s2b$model          <- "S2b_response_TMB_purity"
    res_fgsea_s2b$n_samples      <- length(resp_s2)
    res_fgsea_s2b$n_responders   <- sum(resp_s2 == "response")
    res_fgsea_s2b$n_non_responders <- sum(resp_s2 == "resistance")
  }

## ---- Merge all results --------------------------------------------------
res_limma <- rbind(
  res_limma_m0,
  res_limma_m1,
  res_limma_m2,
  res_limma_s2a,
  res_limma_s2b
)
res_fgsea <- data.table::rbindlist(
  Filter(Negate(is.null),
    list(
      res_fgsea_m0,
      res_fgsea_m1,
      res_fgsea_m2,
      res_fgsea_s2a,
      res_fgsea_s2b
    )),
  fill = TRUE
)


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
