#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "optparse",
  "GEOquery"
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
    c("-g", "--clinical_file"),
    type = "character",
    default = NULL,
    help = "Path to clinical file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-c", "--clinical_output_file"),
    type = "character",
    default = NULL,
    help = "Path to output clinical data file.",
    metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

CLINICAL_FILE <- opt$clinical_file
CLINICAL_OUTPUT_FILE <- opt$clinical_output_file

if (is.null(CLINICAL_FILE) || is.null(CLINICAL_OUTPUT_FILE)) {
  stop("Both --clinical_file and --clinical_output_file are required.")
}

if (!file.exists(CLINICAL_FILE)) {
  stop("clinical_file not found: ", CLINICAL_FILE)
}

#######################
## Load clinical data ##
#######################
clinical_env <- new.env()
load(CLINICAL_FILE, envir = clinical_env)

if (!"pheno" %in% ls(envir = clinical_env)) {
  stop("Expected object 'pheno' not found in RData file.")
}

clinical <- clinical_env[["pheno"]]

############################
## Process clinical table ##
############################
clin <- clinical[clinical[["Histology:ch1"]] == "AD", ]

clin$OS_event <- ifelse(
  is.na(clin[["Status:ch1"]]),
  NA,
  ifelse(clin[["Status:ch1"]] == "Dead", 1, 0)
)

clin$OS_time <- as.numeric(clin[["Survival (days):ch1"]])

clin$gender <- factor(
  clin[["Sex:ch1"]],
  levels = c("F", "M"),
  labels = c("female", "male")
)

clin$age <- as.numeric(clin[["Age:ch1"]])
clin$geo_accession <- rownames(clin)

clin$stage <- clin[["Stage (Pathological ):ch1"]]
clin$stage <- gsub("A|B", "", clin$stage)
clin$stage <- paste0("Stage_", clin$stage)

clin$smoking_status <- ifelse(
  is.na(clin[["Smoking (BI):ch1"]]),
  NA,
  ifelse(clin[["Smoking (BI):ch1"]] == 0, "Never smoker", "Smoker")
)

#########################
## 5-year OS censoring ##
#########################
# Set to NULL to auto-detect based on the range of OS_time values
time_unit <- NULL
cutoff_years <- 5

cat("OS_time summary (to help verify time unit):\n")
print(summary(clin$OS_time))

if (is.null(time_unit)) {
  med <- median(clin$OS_time, na.rm = TRUE)
  time_unit <- if (med > 365) "days" else if (med > 20) "months" else "years"
  cat("Auto-detected time_unit:", time_unit,
      "(set manually if incorrect)\n")
}

cutoff <- switch(time_unit,
  days   = cutoff_years * 365.25,
  months = cutoff_years * 12,
  years  = cutoff_years,
  stop("Unknown time_unit: must be 'days', 'months', or 'years'")
)
cat("Using cutoff:", cutoff, time_unit, "(=", cutoff_years, "years)\n")

clin$OS_event_5y <- clin$OS_event
clin$OS_time_5y <- as.numeric(clin$OS_time)

clin$OS_event_5y[clin$OS_time_5y > cutoff] <- 0

# clin$stage <- ifelse(
#   clin$stage %in% c("stage_II", "stage_III"),
#   "Stage_II-III",
#   ifelse(clin$stage == "stage_I", "Stage_I", NA)
# )


########################
## Select final fields ##
########################
clin <- clin[, c(
  "geo_accession",
  "OS_time",
  "OS_event",
  "OS_time_5y",
  "OS_event_5y",
  "age",
  "gender",
  "smoking_status",
  "stage"
)]

##################
## Write output ##
##################
write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)