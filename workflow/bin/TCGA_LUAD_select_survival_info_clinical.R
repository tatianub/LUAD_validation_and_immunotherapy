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
clin <- fread(CLINICAL_FILE)
head(clin)
table(clin$stage)
clin$stage <- gsub("A|B", "", clin$stage)
table(clin$stage)
clin$stage <- gsub("[Not vailable]", NA, clin$stage, fixed = TRUE)
clin$stage <- gsub("[Discrepancy]", NA, clin$stage, fixed = TRUE)
table(clin$stage)

clin$stage[clin$stage == "Stage I"] <- "Stage_I"
clin$stage[clin$stage == "Stage II"] <- "Stage_II"
clin$stage[clin$stage == "Stage III"] <- "Stage_III"
clin$stage[clin$stage == "Stage IV"] <- "Stage_IV"

head(clin)
clin$age <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
clin$gender[clin$gender == "FEMALE"] <- "female"
clin$gender[clin$gender == "MALE"] <- "male"
head(clin)

clin$DSS_event <- clin$DSS
clin$DSS_time <- clin$DSS.time
clin$DFI_event <- clin$DFI  
clin$DFI_time <- clin$DFI.time
clin$PFI_event <- clin$PFI
clin$PFI_time <- clin$PFI.time
clin$geo_accession <- clin$bcr_patient_barcode
clin$OS_time <- as.numeric(clin$OS_time)

#########################
## 5-year OS censoring ##
#########################
# Set to NULL to auto-detect based on the range of OS_time values

cutoff_years <- 5
# Set time_unit to match the unit of OS_time in the data:
# "days", "months", or "years"
# Set to NULL to auto-detect based on the range of OS_time values
time_unit <- NULL

cat("OS_time summary (to help verify time unit):\n")
print(summary(clin$OS_time))
clin$OS_time <- as.numeric(clin$OS_time)
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


clin$OS_event_5y <- clin$OS_event
clin$OS_time_5y <- as.numeric(clin$OS_time)

clin$OS_event_5y[clin$OS_time_5y > cutoff] <- 0
clin <- clin[!is.na(clin$stage), ]
clin <- clin[, c(
  "geo_accession", "OS_event", "OS_time",
  "OS_event_5y", "OS_time_5y", 
  "DSS_event", "DSS_time", "DFI_event", 
  "DFI_time", "PFI_event", "PFI_time",
  "age", "gender", "stage"
)]
head(clin)


write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

