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
  "dplyr",
  "readxl"
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
    c("--clinical_file_extended"),
    type = "character",
    default = NULL,
    help = "Path to extended clinical file."
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
    c("--fgsea_results_tmb_subset"),
    type = "character",
    default = NULL,
    help = "Path to FGSEA results for TMB subset."
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

IMMUNE_ESTIMATE <- opt$immune_estimate
FGSEA_RESULTS_FILE <- opt$fgsea_results_tmb_subset
HISTO_SUBTYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
OUTPUT_PREFIX <- opt$output_prefix
GMT_FILE <- opt$pathway_gmt_file
CLINICAL_FILE_EXTENDED <- opt$clinical_file_extended
SEED <- opt$seed

RESULTS_LOGISTIC_MODELS <- opt$results_logistic_models
RESULTS_ROC_COMPARISONS <- opt$results_roc_comparisons

set.seed(SEED)


# Load required functions

source("workflow/bin/lung_immuno_limma_fn.R")
source("workflow/bin/predictive_modeling_fn.R")


#####################
## Mutation burden ##
#####################


mutation_burden <- read_excel(
  CLINICAL_FILE_EXTENDED,
  sheet = "Table_S5_Mutation_Burden"
)

mutation_burden$sample_id <-
  mutation_burden$Harmonized_SU2C_WES_Tumor_Sample_ID_v2

mutation_burden$sample_id <- gsub(
  "-",
  ".",
  mutation_burden$sample_id
)


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
indegree[1:5, 1:5]

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

indegree <- indegree[, colnames(indegree) %in% mutation_burden$sample_id]
indegree <- indegree[, colnames(indegree) %in% clinical_data_filt$sample_id]

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
## 3. PD-1 leading-edge GTS
## ---------------------------------
fgsea_results <- read.table(
  FGSEA_RESULTS_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
pd1_row <- fgsea_results[fgsea_results$pathway == "REACTOME_PD_1_SIGNALING", ]
pd1_row <- pd1_row[pd1_row$model == "S2a_response_only_41_subset", ]
print(pd1_row$leadingEdge)
pd1_leading_edge <- unlist(strsplit(pd1_row$leadingEdge, ";"))

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


analysis_data$response_binary <- ifelse(
  analysis_data$response == "response",
  1,
  0
)

# ############################################################
# ## TMB subset analysis
# ############################################################

# ## Add TMB data
tmb_data <- data.frame(
  sample_id = mutation_burden$sample_id,
  TMB = as.numeric(mutation_burden$TMB),
  stringsAsFactors = FALSE
)

## Make sample IDs consistent with the other datasets
tmb_data$sample_id <- gsub(
  "-",
  ".",
  tmb_data$sample_id
)

analysis_data <- merge(
  analysis_data,
  tmb_data,
  by = "sample_id",
  all.x = TRUE
)


############################################################
## Complete-case dataset for TMB analysis
############################################################

## IMPORTANT:
## All biomarkers are evaluated in exactly the same patients.
## This makes all models and ROC comparisons directly comparable.

tmb_comparison_data <- analysis_data[
  complete.cases(
    analysis_data[
      ,
      c(
        "response_binary",
        "PD1_GTS_sum",
        "PDCD1_GTS",
        "CD274_GTS",
        "PDL1_expression",
        "ImmuneScore",
        "PD1_LE_GTS_SUM",
        "TMB"
      )
    ]
  ),
]


cat(
  "\nTMB complete-case patients used: ",
  nrow(tmb_comparison_data),
  "\n",
  sep = ""
)

cat(
  "Responders: ",
  sum(tmb_comparison_data$response_binary == 1),
  "\n",
  sep = ""
)

cat(
  "Non-responders: ",
  sum(tmb_comparison_data$response_binary == 0),
  "\n\n",
  sep = ""
)


############################################################
## Standardize predictors within the TMB subset
############################################################

tmb_comparison_data$PD1_GTS_z <- as.numeric(
  scale(tmb_comparison_data$PD1_GTS_sum)
)

tmb_comparison_data$PDCD1_GTS_z <- as.numeric(
  scale(tmb_comparison_data$PDCD1_GTS)
)

tmb_comparison_data$CD274_GTS_z <- as.numeric(
  scale(tmb_comparison_data$CD274_GTS)
)

tmb_comparison_data$PDL1_z <- as.numeric(
  scale(tmb_comparison_data$PDL1_expression)
)

tmb_comparison_data$ImmuneScore_z <- as.numeric(
  scale(tmb_comparison_data$ImmuneScore)
)

tmb_comparison_data$PD1_LE_GTS_SUM_z <- as.numeric(
  scale(tmb_comparison_data$PD1_LE_GTS_SUM)
)

# Convert mutation count to mutations/Mb and log10-transform
tmb_comparison_data$TMB_log <- log10(
  (tmb_comparison_data$TMB + 1) / 33
)

# Standardize transformed TMB for logistic regression
tmb_comparison_data$TMB_z <- as.numeric(
  scale(tmb_comparison_data$TMB_log)
)

############################################################
## Univariable models in TMB subset
############################################################

tmb_model_tmb <- glm(
  response_binary ~ TMB_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_pd1 <- glm(
  response_binary ~ PD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_pdcd1 <- glm(
  response_binary ~ PDCD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_cd274 <- glm(
  response_binary ~ CD274_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_pdl1 <- glm(
  response_binary ~ PDL1_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_immune <- glm(
  response_binary ~ ImmuneScore_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_le <- glm(
  response_binary ~ PD1_LE_GTS_SUM_z,
  data = tmb_comparison_data,
  family = binomial
)


############################################################
## Network features added to TMB
############################################################

tmb_model_tmb_pd1 <- glm(
  response_binary ~ TMB_z + PD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_tmb_pdcd1 <- glm(
  response_binary ~ TMB_z + PDCD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_tmb_cd274 <- glm(
  response_binary ~ TMB_z + CD274_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_tmb_le <- glm(
  response_binary ~ TMB_z + PD1_LE_GTS_SUM_z,
  data = tmb_comparison_data,
  family = binomial
)



############################################################
## TMB + PD-L1 + ImmuneScore reference model
############################################################

tmb_model_reference <- glm(
  response_binary ~
    TMB_z +
    PDL1_z +
    ImmuneScore_z,
  data = tmb_comparison_data,
  family = binomial
)


############################################################
## Network features added to TMB + PD-L1 + ImmuneScore
############################################################

tmb_model_full_pd1 <- glm(
  response_binary ~
    TMB_z +
    PDL1_z +
    ImmuneScore_z +
    PD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_full_pdcd1 <- glm(
  response_binary ~
    TMB_z +
    PDL1_z +
    ImmuneScore_z +
    PDCD1_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_full_cd274 <- glm(
  response_binary ~
    TMB_z +
    PDL1_z +
    ImmuneScore_z +
    CD274_GTS_z,
  data = tmb_comparison_data,
  family = binomial
)

tmb_model_full_le <- glm(
  response_binary ~
    TMB_z +
    PDL1_z +
    ImmuneScore_z +
    PD1_LE_GTS_SUM_z,
  data = tmb_comparison_data,
  family = binomial
)



############################################################
## Logistic regression results
############################################################

tmb_logistic_results <- rbind(

  ## Univariable models
  extract_model_results(
    tmb_model_tmb,
    "TMB"
  ),

  extract_model_results(
    tmb_model_pd1,
    "PD1_pathway_GTS"
  ),

  extract_model_results(
    tmb_model_pdcd1,
    "PDCD1_GTS"
  ),

  extract_model_results(
    tmb_model_cd274,
    "CD274_GTS"
  ),

  extract_model_results(
    tmb_model_pdl1,
    "PDL1_expression"
  ),

  extract_model_results(
    tmb_model_immune,
    "ESTIMATE_ImmuneScore"
  ),

  extract_model_results(
    tmb_model_le,
    "PD1_LE_GTS_SUM"
  ),


  ## TMB-adjusted models
  extract_model_results(
    tmb_model_tmb_pd1,
    "TMB_plus_PD1_pathway_GTS"
  ),

  extract_model_results(
    tmb_model_tmb_pdcd1,
    "TMB_plus_PDCD1_GTS"
  ),

  extract_model_results(
    tmb_model_tmb_cd274,
    "TMB_plus_CD274_GTS"
  ),

  extract_model_results(
    tmb_model_tmb_le,
    "TMB_plus_PD1_LE_GTS_SUM"
  ),

  ## TMB + PD-L1 + ImmuneScore reference model
  extract_model_results(
    tmb_model_reference,
    "TMB_plus_PDL1_plus_ImmuneScore"
  ),


  ## Fully adjusted models
  extract_model_results(
    tmb_model_full_pd1,
    "TMB_plus_PDL1_plus_ImmuneScore_plus_PD1_pathway_GTS"
  ),

  extract_model_results(
    tmb_model_full_pdcd1,
    "TMB_plus_PDL1_plus_ImmuneScore_plus_PDCD1_GTS"
  ),

  extract_model_results(
    tmb_model_full_cd274,
    "TMB_plus_PDL1_plus_ImmuneScore_plus_CD274_GTS"
  ),

  extract_model_results(
    tmb_model_full_le,
    "TMB_plus_PDL1_plus_ImmuneScore_plus_PD1_LE_GTS_SUM"
  )
)


cat(
  "\nTMB-subset logistic regression results:\n\n"
)

print(tmb_logistic_results)


############################################################
## ROC: network features added to TMB
############################################################

tmb_roc_vs_tmb <- bind_rows(

  compare_nested_roc(
    tmb_model_tmb,
    tmb_model_tmb_pd1,
    "TMB + PD1 pathway GTS"
  ),

  compare_nested_roc(
    tmb_model_tmb,
    tmb_model_tmb_pdcd1,
    "TMB + PDCD1 GTS"
  ),

  compare_nested_roc(
    tmb_model_tmb,
    tmb_model_tmb_cd274,
    "TMB + CD274 GTS"
  ),

  compare_nested_roc(
    tmb_model_tmb,
    tmb_model_tmb_le,
    "TMB + PD1 leading-edge GTS"
  )
)


############################################################
## ROC: network features added to
## TMB + PD-L1 + ImmuneScore
############################################################

tmb_roc_vs_full <- bind_rows(

  compare_nested_roc(
    tmb_model_reference,
    tmb_model_full_pd1,
    "TMB + PD-L1 + ImmuneScore + PD1 pathway GTS"
  ),

  compare_nested_roc(
    tmb_model_reference,
    tmb_model_full_pdcd1,
    "TMB + PD-L1 + ImmuneScore + PDCD1 GTS"
  ),

  compare_nested_roc(
    tmb_model_reference,
    tmb_model_full_cd274,
    "TMB + PD-L1 + ImmuneScore + CD274 GTS"
  ),

  compare_nested_roc(
    tmb_model_reference,
    tmb_model_full_le,
    "TMB + PD-L1 + ImmuneScore + PD1 leading-edge GTS"
  )
)


############################################################
## Final TMB ROC comparison table
############################################################

tmb_roc_results <- bind_rows(

  tmb_roc_vs_tmb %>%
    mutate(
      reference = "TMB"
    ),

  tmb_roc_vs_full %>%
    mutate(
      reference = "TMB + PD-L1 + ImmuneScore"
    )

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


cat(
  "\nTMB-subset ROC/model comparisons:\n\n"
)

print(tmb_roc_results)


############################################################
## Save TMB subset results
############################################################

write.table(
  tmb_logistic_results,
  RESULTS_LOGISTIC_MODELS,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tmb_roc_results,
  RESULTS_ROC_COMPARISONS,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)