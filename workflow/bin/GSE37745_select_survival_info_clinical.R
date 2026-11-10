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

cat("Loading clinical file...\n")
clin_env <- new.env()
load(CLINICAL_FILE, envir = clin_env)

if (!"pheno" %in% ls(envir = clin_env)) {
  stop("Expected object 'pheno' not found in RData file.")
}

clin <- clin_env[["pheno"]]

required_cols <- c(
  "geo_accession",
  "adjuvant treatment:ch1",
  "dead:ch1",
  "days to determined death status:ch1",
  "recurrence:ch1",
  "days to recurrence / to last visit:ch1",
  "age:ch1",
  "gender:ch1",
  "histology:ch1",
  "tumor stage:ch1"
)

missing_cols <- setdiff(required_cols, colnames(clin))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

clin$treatment <- clin[["adjuvant treatment:ch1"]]

clin$OS_event <- clin[["dead:ch1"]]
clin$OS_event <- ifelse(
  clin$OS_event == "no", 0,
  ifelse(clin$OS_event == "yes", 1, NA)
)

clin$OS_time <- as.numeric(clin[["days to determined death status:ch1"]])

cutoff_years <- 5
time_unit <- NULL

cat("OS_time summary (to help verify time unit):\n")
print(summary(clin$OS_time))

if (is.null(time_unit)) {
  med <- median(clin$OS_time, na.rm = TRUE)
  time_unit <- if (med > 365) {
    "days"
  } else if (med > 20) {
    "months"
  } else {
    "years"
  }
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

clin$OS_event_5y <- ifelse(clin$OS_time > cutoff, 0, clin$OS_event)
clin$OS_time_5y <- clin$OS_time

clin$RFS_event <- clin[["recurrence:ch1"]]
clin$RFS_event <- ifelse(
  clin$RFS_event == "no", 0,
  ifelse(clin$RFS_event == "yes", 1, NA)
)

clin$RFS_time <- as.numeric(clin[["days to recurrence / to last visit:ch1"]])

clin$age <- as.numeric(clin[["age:ch1"]])
clin$gender <- clin[["gender:ch1"]]

clin <- clin[clin[["histology:ch1"]] == "adeno", ]

clin$stage <- clin[["tumor stage:ch1"]]
clin$stage <- gsub("a|b|c", "", clin$stage)
clin$stage <- paste("Stage", clin$stage, sep = "_")

# clin$stage <- ifelse(
#   clin$stage %in% c("Stage_1", "Stage_2"), "Stage_I_II",
#   ifelse(clin$stage %in% c("Stage_3", "Stage_4"), "Stage_III_IV", NA)
# )

clin <- clin[, c(
  "geo_accession",
  "OS_event",
  "OS_time",
  "RFS_event",
  "RFS_time",
  "OS_event_5y",
  "OS_time_5y",
  "age",
  "gender",
  "stage",
  "treatment"
)]

write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

