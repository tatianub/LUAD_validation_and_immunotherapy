#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table", "dplyr", "ggplot2", "optparse", "survival",
  "ggpubr", "grid", "tools", "limma", "fgsea", "ggrepel"
)

for (lib in required_libraries) {
  suppressPackageStartupMessages(
    library(lib, character.only = TRUE, quietly = TRUE)
  )
}

## Read arguments ##
####################
option_list <- list(
  make_option(c("--network_file"), type = "character", default = NULL,
              help = "Path to network file.", metavar = "character"),
  make_option(c("--edges_file"), type = "character", default = NULL,
              help = "Path to edges file.", metavar = "character"),
  make_option(c("--gmt_file"), type = "character", default = NULL,
              help = "Path to GMT file.", metavar = "character"),
  make_option(c("--clinical_file"), type = "character", default = NULL,
              help = "Path to clinical file.", metavar = "character"),
  make_option(c("--subtype_type"), type = "character", default = "all",
              help = "Type of subtype.", metavar = "character"),
  make_option(c("--treatment_type"), type = "character", default = "all",
              help = "Type of treatment.", metavar = "character"),
  make_option(c("--limma_results_edges"), type = "character", default = NULL,
              help = "Path to save limma results.", metavar = "character"),
  make_option(c("--fgsea_results"), type = "character", default = NULL,
              help = "Path to save fgsea results.", metavar = "character"),
  make_option(c("--num_cores"), type = "integer", default = 1,
              help = "Number of cores.", metavar = "integer")
)
set.seed(12345)
# Crucial fix: Must parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



########################
## Initialize values  ##
########################
CLINICAL_FILE <- opt$clinical_file
NETWORK_FILE <- opt$network_file
EDGES_FILE <- opt$edges_file
GMT_FILE <- opt$gmt_file
HISTO_SUBTYPE <- opt$subtype_type
TREATMENT_TYPE <- opt$treatment_type
RESULTS_DIFFERENTIAL_EDGES <- opt$limma_results_edges
RESULTS_FGSEA <- opt$fgsea_results
NUM_CORES <- opt$num_cores


source("workflow/bin/lung_immuno_limma_fn.R")

# read in clinical data
clinical_data <- fread(CLINICAL_FILE)
head(clinical_data)
clinical_data <- clinical_data[!is.na(clinical_data$response)]

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

# Load network (handles both 'net' and 'net_norm' objects)
net_env <- new.env()
load(NETWORK_FILE, envir = net_env)

# Check what objects are available
obj_names <- ls(net_env)
cat("Available objects in network file:", 
    paste(obj_names, collapse = ", "), "\n")

# Try to get the network object (prefer 'net', fallback to 'net_norm')
if ("net" %in% obj_names) {
    net <- net_env$net
    cat("Loaded 'net' object\n")
} else if ("net_norm" %in% obj_names) {
    net <- net_env$net_norm
    cat("Loaded 'net_norm' object\n")
} else {
    stop("No recognized network object found. Available objects: ", 
         paste(obj_names, collapse = ", "))
}

# read in edges file

edges <- fread(EDGES_FILE)

net <- as.matrix(net)
colnames(net) <- gsub("-", "\\.", colnames(net))
## clean network
net_clean <- net[, colnames(net) %in% clinical_data$sample_id]
colnames(net_clean)

clinical_data_filt <- 
    clinical_data[match(colnames(net_clean), clinical_data$sample_id)]

stopifnot(all(colnames(net_clean) == clinical_data$sample_id))

clinical_data_ordered <- clinical_data_filt[
    match(colnames(net_clean), clinical_data_filt$sample_id), ]

covariates <- data.frame(
    age = clinical_data_ordered$Patient_Age_at_Diagnosis,
    sex = clinical_data_ordered$Patient_Sex,
    smoking = clinical_data_ordered$Patient_Smoking_Pack_Years_Harmonized)
rownames(covariates) <- clinical_data_ordered$sample_id
covariates$age <- as.numeric(covariates$age)
covariates$sex <- factor(covariates$sex)
covariates$smoking <- as.numeric(covariates$smoking)
covariates <- covariates[
    !is.na(covariates$age) & !is.na(covariates$sex) & !is.na(covariates$smoking), ]
net_clean_covariates <- net_clean[, rownames(covariates)]

if (!all(colnames(net_clean_covariates) == rownames(covariates))) {
    stop(
      "Column names of net_clean_covariates do not match ",
      "row names of covariates"
    )
  }
response_updated <- clinical_data_filt$response[
    match(colnames(net_clean_covariates), clinical_data_filt$sample_id)]


res_net <- run_network_analysis(net_clean_covariates, 
                            response_updated, 
                            edges, 
                            covariates)
head(res_net)
rm(net_clean, net)
write.table(res_net, 
    RESULTS_DIFFERENTIAL_EDGES, 
    sep = "\t", quote = FALSE, row.names = FALSE)

fgseaRes_net <- run_fgsea_analysis_TFs(
    res_net,
    GMT_FILE,
    num_cores = NUM_CORES
    )
save_fgsea_results(fgseaRes_net, RESULTS_FGSEA)


