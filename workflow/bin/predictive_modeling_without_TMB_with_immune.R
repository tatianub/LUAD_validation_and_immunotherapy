#####################
## Load R packages ##
#####################

required_libraries <- c(
  "data.table",
  "dplyr",
  "optparse",
  "pROC",
  "ggplot2",
  "ggpubr",
  "patchwork",
  "pROC",
  "dplyr"
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
    help = "Path to indegree/GTS file."
  ),

  optparse::make_option(
    c("--clinical_file"),
    type = "character",
    default = NULL,
    help = "Path to clinical file."
  ),

  optparse::make_option(
    c("--expression_file"),
    type = "character",
    default = NULL,
    help = "Path to expression file."
  ),

  optparse::make_option(
    c("--samples_file"),
    type = "character",
    default = NULL,
    help = "Path to samples file."
  ),

  optparse::make_option(
    c("--pathway_gmt_file"),
    type = "character",
    default = NULL,
    help = "Path to pathway GMT file."
  ),

  optparse::make_option(
    c("--immune_estimate"),
    type = "character",
    default = NULL,
    help = "Path to ESTIMATE immune score file."
  ),

  optparse::make_option(
    c("--subtype_type"),
    type = "character",
    default = "all"
  ),

  optparse::make_option(
    c("--treatment_type"),
    type = "character",
    default = "all"
  ),

  optparse::make_option(
    c("--results_logistic_models"),
    type = "character",
    default = NULL,
    help = "Path to results logistics file."
  ),
  optparse::make_option(
    c("--results_roc_comparisons"),
    type = "character",
    default = NULL,
    help = "Path to results ROC comparisons file."
  ),

  optparse::make_option(
    c("--seed"),
    type = "integer",
    default = 2026
  )
)

opt_parser <- optparse::OptionParser(
  option_list = option_list
)

opt <- optparse::parse_args(opt_parser)


#######################
## Initialize values ##
#######################

INDEGREE_FILE <- opt$indegree_file
CLINICAL_FILE <- opt$clinical_file
EXPRESSION_FILE <- opt$expression_file
SAMPLES_FILE <- opt$samples_file
GMT_FILE <- opt$pathway_gmt_file
IMMUNE_ESTIMATE <- opt$immune_estimate

HISTO_SUBTYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
OUTPUT_PREFIX <- opt$output_prefix
SEED <- opt$seed

RESULTS_LOGISTIC_MODELS <- opt$results_logistic_models
RESULTS_ROC_COMPARISONS <- opt$results_roc_comparisons


set.seed(SEED)


# Load required functions

source("workflow/bin/lung_immuno_limma_fn.R")
source("workflow/bin/predictive_modeling_fn.R")


##########################
## Read expression data ##
##########################

expression <- fread(
  EXPRESSION_FILE
)

genes_exp <- expression$V1

expression <- as.matrix(
  expression[, -1]
)

samples <- fread(
  SAMPLES_FILE
)

colnames(expression) <- samples$V2
rownames(expression) <- genes_exp

colnames(expression) <- gsub(
  "-",
  ".",
  colnames(expression)
)


########################
## PD-L1 expression   ##
########################

exp_pdl1 <- expression[
  grep(
    "^CD274$",
    rownames(expression)
  ),
  ,
  drop = FALSE
]

if (nrow(exp_pdl1) != 1) {
  stop(
    "Could not uniquely identify CD274 in expression matrix."
  )
}

pdl1_data <- data.frame(
  sample_id = colnames(expression),
  PDL1_expression = as.numeric(
    exp_pdl1[1, ]
  ),
  stringsAsFactors = FALSE
)


#################################
## ESTIMATE immune infiltration ##
#################################

immune_estimate <- fread(
  IMMUNE_ESTIMATE
)

colnames(immune_estimate)[1] <- "sample_id"

immune_estimate$sample_id <- gsub(
  "-",
  ".",
  immune_estimate$sample_id
)

immune_data <- data.frame(
  sample_id = immune_estimate$sample_id,
  ImmuneScore = as.numeric(
    immune_estimate$ImmuneScore
  ),
  stringsAsFactors = FALSE
)

cat(
  "\nESTIMATE immune scores loaded: ",
  nrow(immune_data),
  "\n\n",
  sep = ""
)


########################
## Read clinical data ##
########################

clinical_data <- fread(
  CLINICAL_FILE
)

clinical_data <- clinical_data[
  !is.na(clinical_data$response)
]

clinical_data_filt <- subset_clinical_data(
  clinical_data = clinical_data,
  subtype = HISTO_SUBTYPE,
  treatment_type = TREATMENT_TYPE
)

clinical_data_filt$sample_id <- gsub(
  "-",
  ".",
  clinical_data_filt$sample_id
)


########################
## Read indegree/GTS  ##
########################

indegree <- fread(
  INDEGREE_FILE
)

genes_gts <- indegree$tar

indegree <- as.matrix(
  indegree[, -1]
)

rownames(indegree) <- genes_gts

colnames(indegree) <- gsub(
  "-",
  ".",
  colnames(indegree)
)


##########################
## Read PD-1 gene set   ##
##########################

gmt_lines <- readLines(
  GMT_FILE
)

gmt_split <- strsplit(
  gmt_lines,
  "\t"
)

pathway_names <- sapply(
  gmt_split,
  `[`,
  1
)

pd1_index <- which(
  pathway_names == "REACTOME_PD_1_SIGNALING"
)

if (length(pd1_index) != 1) {
  stop(
    paste0(
      "Could not uniquely identify ",
      "REACTOME_PD_1_SIGNALING in GMT file."
    )
  )
}

# GMT:
# column 1 = pathway
# column 2 = description
# columns 3+ = genes

pd1_genes <- gmt_split[[pd1_index]][
  -c(1, 2)
]

pd1_genes_present <- intersect(
  pd1_genes,
  rownames(indegree)
)

cat(
  "\nPD-1 pathway genes in GMT: ",
  length(pd1_genes),
  "\n",
  sep = ""
)

cat(
  "PD-1 pathway genes present in GTS matrix: ",
  length(pd1_genes_present),
  "\n",
  sep = ""
)

cat(
  "Genes used: ",
  paste(
    pd1_genes_present,
    collapse = ", "
  ),
  "\n\n",
  sep = ""
)

if (length(pd1_genes_present) == 0) {
  stop(
    "No PD-1 pathway genes found in indegree matrix."
  )
}


#####################################
## Calculate patient-level PD1 GTS ##
#####################################

pd1_gts_sum <- colSums(
  indegree[
    pd1_genes_present,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

pd1_gts <- data.frame(
  sample_id = names(pd1_gts_sum),
  PD1_GTS_sum = as.numeric(
    pd1_gts_sum
  )
)


## ---------------------------------
## 2. Whole PD-1 pathway: PCA
## ---------------------------------
pd1_matrix <- indegree[
  pd1_genes_present,
  ,
  drop = FALSE
]

# PCA expects samples in rows and genes in columns
pd1_matrix_pca <- t(pd1_matrix)

# Check for missing values
if (anyNA(pd1_matrix_pca)) {
  stop("Missing values detected in PD-1 GTS matrix before PCA.")
}

# Remove genes with zero variance, if any
nonzero_var <- apply(
  pd1_matrix_pca,
  2,
  sd
) > 0

pd1_matrix_pca <- pd1_matrix_pca[
  ,
  nonzero_var,
  drop = FALSE
]

# PCA
# scale. = TRUE means each gene contributes on a comparable scale
pd1_pca <- prcomp(
  pd1_matrix_pca,
  center = TRUE,
  scale. = TRUE
)

# Patient-level PC1 score
pd1_gts_pc1 <- pd1_pca$x[, "PC1"]

# Variance explained by each PC
pd1_pca_variance <- (
  pd1_pca$sdev^2 /
    sum(pd1_pca$sdev^2)
) * 100

cat(
  "PD-1 GTS PC1 variance explained:",
  round(pd1_pca_variance[1], 2),
  "%\n"
)

cat(
  "PD-1 GTS PC2 variance explained:",
  round(pd1_pca_variance[2], 2),
  "%\n"
)


## ---------------------------------
## 3. PD-1 leading-edge GTS
## ---------------------------------

pd1_leading_edge <- c(
  "CD3D",
  "HLA-DRB1",
  "PDCD1",
  "HLA-DQB1",
  "HLA-DQA1",
  "HLA-DQB2",
  "HLA-DRB5",
  "CD3G",
  "HLA-DPA1",
  "HLA-DRA",
  "CD274",
  "HLA-DPB1",
  "CD3E"
)



# Determine which leading-edge genes are present
pd1_leading_edge_present <- intersect(
  pd1_leading_edge,
  rownames(indegree)
)

cat(
  "PD-1 leading-edge genes present:",
  length(pd1_leading_edge_present),
  "/",
  length(pd1_leading_edge),
  "\n"
)

print(pd1_leading_edge_present)


# Extract leading-edge GTS matrix
pd1_le_matrix <- indegree[
  pd1_leading_edge_present,
  ,
  drop = FALSE
]


# Sum leading-edge GTS
pd1_le_gts_sum <- colSums(
  pd1_le_matrix,
  na.rm = TRUE
)



#########################################
## Extract PDCD1 and CD274 indegrees   ##
#########################################

if (!"PDCD1" %in% rownames(indegree)) {
  stop(
    "PDCD1 not found in indegree matrix."
  )
}

if (!"CD274" %in% rownames(indegree)) {
  stop(
    "CD274 not found in indegree matrix."
  )
}

checkpoint_gts <- data.frame(

  sample_id = colnames(indegree),

  PDCD1_GTS = as.numeric(
    indegree["PDCD1", ]
  ),

  CD274_GTS = as.numeric(
    indegree["CD274", ]
  ),

  stringsAsFactors = FALSE
)

############################################################
## Build patient-level analysis dataset
############################################################

analysis_data <- data.frame(
  sample_id = clinical_data_filt$sample_id,
  response = clinical_data_filt$response,
  stringsAsFactors = FALSE
)

analysis_data <- merge(
  analysis_data,
  pd1_gts,
  by = "sample_id",
  all.x = TRUE
)

analysis_data <- merge(
  analysis_data,
  checkpoint_gts,
  by = "sample_id",
  all.x = TRUE
)

analysis_data <- merge(
  analysis_data,
  pdl1_data,
  by = "sample_id",
  all.x = TRUE
)

analysis_data <- merge(
  analysis_data,
  immune_data,
  by = "sample_id",
  all.x = TRUE
)


analysis_data <- merge(
  analysis_data,
  data.frame(
    sample_id = names(pd1_le_gts_sum),
    PD1_LE_GTS_SUM = as.numeric(pd1_le_gts_sum),
    stringsAsFactors = FALSE
  ),
  by = "sample_id",
  all.x = TRUE
)

analysis_data <- merge(
  analysis_data,
  data.frame(
    sample_id = names(pd1_gts_pc1),
    PD1_GTS_PC1 = as.numeric(pd1_gts_pc1),
    stringsAsFactors = FALSE
  ),
  by = "sample_id",
  all.x = TRUE
)


############################################################
## Encode treatment response
############################################################

# Responders = 1
# Non-responders = 0

analysis_data$response_binary <- ifelse(
  analysis_data$response == "response",
  1,
  0
)



############################################################
## Complete-case data for biomarker comparison
############################################################

# IMPORTANT:
# All biomarkers are evaluated in exactly the same patients.
# This ensures that ROC curves, DeLong tests, logistic models,
# and likelihood-ratio tests are directly comparable.

comparison_data <- analysis_data[
  complete.cases(
    analysis_data[
      ,
      c(
        "response_binary",
        "PD1_GTS_sum",
        "PDCD1_GTS",
        "CD274_GTS",
        "PDL1_expression",
        "PD1_GTS_PC1",
        "PD1_LE_GTS_SUM",
        "PD1_GTS_PC1",
        "ImmuneScore"
      )
    ]
  ),
]

cat(
  "Complete-case patients used: ",
  nrow(comparison_data),
  "\n",
  sep = ""
)

cat(
  "Responders: ",
  sum(comparison_data$response_binary == 1),
  "\n",
  sep = ""
)

cat(
  "Non-responders: ",
  sum(comparison_data$response_binary == 0),
  "\n\n",
  sep = ""
)

############################################################
## Standardize predictors
############################################################

comparison_data$PD1_GTS_z <- as.numeric(
  scale(
    comparison_data$PD1_GTS_sum
  )
)

comparison_data$PDCD1_GTS_z <- as.numeric(
  scale(
    comparison_data$PDCD1_GTS
  )
)

comparison_data$CD274_GTS_z <- as.numeric(
  scale(
    comparison_data$CD274_GTS
  )
)

comparison_data$PDL1_z <- as.numeric(
  scale(
    comparison_data$PDL1_expression
  )
)

comparison_data$ImmuneScore_z <- as.numeric(
  scale(
    comparison_data$ImmuneScore
  )
)

comparison_data$PD1_LE_GTS_SUM_z <- as.numeric(
  scale(
    comparison_data$PD1_LE_GTS_SUM
  )
)

comparison_data$PD1_GTS_PC1_z <- as.numeric(
  scale(
    comparison_data$PD1_GTS_PC1
  )
)


############################################################
## Logistic regression models
############################################################

# ----------------------------------------------------------
# Univariable models
# ----------------------------------------------------------

model_pd1_pathway <- glm(
  response_binary ~ PD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pd1_le_gts_sum <- glm(
  response_binary ~ PD1_LE_GTS_SUM_z,
  data = comparison_data,
  family = binomial
)

model_pd1_gts_pc1 <- glm(
  response_binary ~ PD1_GTS_PC1_z,
  data = comparison_data,
  family = binomial
)

model_pdcd1 <- glm(
  response_binary ~ PDCD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_cd274 <- glm(
  response_binary ~ CD274_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1 <- glm(
  response_binary ~ PDL1_z,
  data = comparison_data,
  family = binomial
)

model_immune <- glm(
  response_binary ~ ImmuneScore_z,
  data = comparison_data,
  family = binomial
)


# ----------------------------------------------------------
# Adjusted for PD-L1 expression
# ----------------------------------------------------------

# contribute beyond PD-L1 expression
model_pdl1_pd1 <- glm(
  response_binary ~
    PDL1_z +
    PD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_pdcd1 <- glm(
  response_binary ~
    PDL1_z +
    PDCD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_cd274 <- glm(
  response_binary ~
    PDL1_z +
    CD274_GTS_z,
  data = comparison_data,
  family = binomial
)
model_pdl1_pd1_le_gts_sum <- glm(
  response_binary ~
    PDL1_z +
    PD1_LE_GTS_SUM_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_pd1_gts_pc1 <- glm(
  response_binary ~
    PDL1_z +
    PD1_GTS_PC1_z,
  data = comparison_data,
  family = binomial
)


# ----------------------------------------------------------
# Adjusted for ESTIMATE ImmuneScore
# ----------------------------------------------------------

model_immune_pd1 <- glm(
  response_binary ~
    ImmuneScore_z +
    PD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_immune_pdcd1 <- glm(
  response_binary ~
    ImmuneScore_z +
    PDCD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_immune_cd274 <- glm(
  response_binary ~
    ImmuneScore_z +
    CD274_GTS_z,
  data = comparison_data,
  family = binomial
)

model_immune_pd1_le_gts_sum <- glm(
  response_binary ~
    ImmuneScore_z +
    PD1_LE_GTS_SUM_z,
  data = comparison_data,
  family = binomial
)

model_immune_pd1_gts_pc1 <- glm(
  response_binary ~
    ImmuneScore_z +
    PD1_GTS_PC1_z,
  data = comparison_data,
  family = binomial
)

# ----------------------------------------------------------
# PD-L1 + ImmuneScore reference model
# ----------------------------------------------------------

model_pdl1_immune <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z,
  data = comparison_data,
  family = binomial
)


# ----------------------------------------------------------
# Adjusted for both PD-L1 and ImmuneScore
# ----------------------------------------------------------

model_pdl1_immune_pd1 <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z +
    PD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_immune_pdcd1 <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z +
    PDCD1_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_immune_cd274 <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z +
    CD274_GTS_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_immune_pd1_le_gts_sum <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z +
    PD1_LE_GTS_SUM_z,
  data = comparison_data,
  family = binomial
)

model_pdl1_immune_pd1_gts_pc1 <- glm(
  response_binary ~
    PDL1_z +
    ImmuneScore_z +
    PD1_GTS_PC1_z,
  data = comparison_data,
  family = binomial
)



############################################################
## Combine logistic regression results
############################################################

logistic_results <- rbind(

  # Univariable models

  extract_model_results(
    model_pd1_pathway,
    "PD1_pathway_GTS"
  ),
  
  extract_model_results(
    model_pdcd1,
    "PDCD1_GTS"
  ),

  extract_model_results(
    model_cd274,
    "CD274_GTS"
  ),

  extract_model_results(
    model_pdl1,
    "PDL1_expression"
  ),

  extract_model_results(
    model_immune,
    "ESTIMATE_ImmuneScore"
  ),

  extract_model_results(
    model_pd1_le_gts_sum,
    "PD1_LE_GTS_SUM"
  ),
  extract_model_results(
    model_pd1_gts_pc1,
    "PD1_GTS_PC1"
  ),

  # PD-L1-adjusted models

  extract_model_results(
    model_pdl1_pd1,
    "PDL1_plus_PD1_pathway_GTS"
  ),

  extract_model_results(
    model_pdl1_pdcd1,
    "PDL1_plus_PDCD1_GTS"
  ),

  extract_model_results(
    model_pdl1_cd274,
    "PDL1_plus_CD274_GTS"
  ),

  extract_model_results(
    model_pdl1_pd1_le_gts_sum,
    "PDL1_plus_PD1_LE_GTS_SUM"
  ),
  extract_model_results(
  model_pdl1_pd1_gts_pc1,
  "PDL1_plus_PD1_GTS_PC1"
  ),

  # ImmuneScore-adjusted models

  extract_model_results(
    model_immune_pd1,
    "ImmuneScore_plus_PD1_pathway_GTS"
  ),

  extract_model_results(
    model_immune_pdcd1,
    "ImmuneScore_plus_PDCD1_GTS"
  ),

  extract_model_results(
    model_immune_cd274,
    "ImmuneScore_plus_CD274_GTS"
  ),
  extract_model_results(
    model_immune_pd1_le_gts_sum,
    "ImmuneScore_plus_PD1_LE_GTS_SUM"
  ),  
  extract_model_results(
  model_immune_pd1_gts_pc1,
    "ImmuneScore_plus_PD1_GTS_PC1"
  ), 
  
  # PD-L1 + ImmuneScore reference model

  extract_model_results(
    model_pdl1_immune,
    "PDL1_plus_ImmuneScore"
  ),

  # Fully adjusted models

  extract_model_results(
    model_pdl1_immune_pd1,
    "PDL1_plus_ImmuneScore_plus_PD1_pathway_GTS"
  ),

  extract_model_results(
    model_pdl1_immune_pdcd1,
    "PDL1_plus_ImmuneScore_plus_PDCD1_GTS"
  ),

  extract_model_results(
    model_pdl1_immune_cd274,
    "PDL1_plus_ImmuneScore_plus_CD274_GTS"
  ),

  extract_model_results(
    model_pdl1_immune_pd1_le_gts_sum,
    "PDL1_plus_ImmuneScore_plus_PD1_LE_GTS_SUM"
  ),

  extract_model_results(
    model_pdl1_immune_pd1_gts_pc1,
    "PDL1_plus_ImmuneScore_plus_PD1_GTS_PC1"
  )
)

cat(
  "\nLogistic regression results:\n\n"
)

print(
  logistic_results
)

logistic_results


############################################################
## Network measures added to PD-L1
############################################################

roc_vs_pdl1 <- bind_rows(

  compare_nested_roc(
    model_pdl1,
    model_pdl1_pd1,
    "PD-L1 + PD1 pathway GTS"
  ),

  compare_nested_roc(
    model_pdl1,
    model_pdl1_pdcd1,
    "PD-L1 + PDCD1 GTS"
  ),

  compare_nested_roc(
    model_pdl1,
    model_pdl1_cd274,
    "PD-L1 + CD274 GTS"
  ),

  compare_nested_roc(
    model_pdl1,
    model_pdl1_pd1_le_gts_sum,
    "PD-L1 + PD1 leading-edge GTS"
  ),

  compare_nested_roc(
    model_pdl1,
    model_pdl1_pd1_gts_pc1,
    "PD-L1 + PD1 GTS PC1"
  )
)

roc_vs_pdl1


############################################################
## Network measures added to ImmuneScore
############################################################

roc_vs_immune <- bind_rows(

  compare_nested_roc(
    model_immune,
    model_immune_pd1,
    "ImmuneScore + PD1 pathway GTS"
  ),

  compare_nested_roc(
    model_immune,
    model_immune_pdcd1,
    "ImmuneScore + PDCD1 GTS"
  ),

  compare_nested_roc(
    model_immune,
    model_immune_cd274,
    "ImmuneScore + CD274 GTS"
  ),

  compare_nested_roc(
    model_immune,
    model_immune_pd1_le_gts_sum,
    "ImmuneScore + PD1 leading-edge GTS"
  ),

  compare_nested_roc(
    model_immune,
    model_immune_pd1_gts_pc1,
    "ImmuneScore + PD1 GTS PC1"
  )
)

roc_vs_immune


############################################################
## Network measures added to PD-L1 + ImmuneScore
############################################################

roc_vs_pdl1_immune <- bind_rows(

  compare_nested_roc(
    model_pdl1_immune,
    model_pdl1_immune_pd1,
    "PD-L1 + ImmuneScore + PD1 pathway GTS"
  ),

  compare_nested_roc(
    model_pdl1_immune,
    model_pdl1_immune_pdcd1,
    "PD-L1 + ImmuneScore + PDCD1 GTS"
  ),

  compare_nested_roc(
    model_pdl1_immune,
    model_pdl1_immune_cd274,
    "PD-L1 + ImmuneScore + CD274 GTS"
  ),

  compare_nested_roc(
    model_pdl1_immune,
    model_pdl1_immune_pd1_le_gts_sum,
    "PD-L1 + ImmuneScore + PD1 leading-edge GTS"
  ),

  compare_nested_roc(
    model_pdl1_immune,
    model_pdl1_immune_pd1_gts_pc1,
    "PD-L1 + ImmuneScore + PD1 GTS PC1"
  )
)

roc_vs_pdl1_immune
############################################################
## Final ROC comparison table
############################################################

roc_model_results <- bind_rows(

  roc_vs_pdl1 %>%
    mutate(reference = "PD-L1"),

  roc_vs_immune %>%
    mutate(reference = "ImmuneScore"),

  roc_vs_pdl1_immune %>%
    mutate(reference = "PD-L1 + ImmuneScore")

) %>%
  select(
    reference,
    comparison,
    AUC_base,
    AUC_extended,
    delta_AUC,
    DeLong_p,
    LRT_p
  )

print(roc_model_results)
############################################################
## Save model and ROC results
############################################################

## Logistic regression results
write.table(
  logistic_results,
  RESULTS_LOGISTIC_MODELS,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## Model-based ROC / nested-model comparisons
write.table(
  roc_model_results,
  RESULTS_ROC_COMPARISONS,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
