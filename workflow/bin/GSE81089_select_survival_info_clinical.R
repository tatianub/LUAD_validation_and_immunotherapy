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
    c("-g", "--gse_series_matrix_file"),
    type = "character",
    default = NULL,
    help = "Path to GSE series matrix file.",
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

########################
## Initialize values  ##
########################
GSE_SERIES_MATRIX_FILE <- opt$gse_series_matrix_file
CLINICAL_OUTPUT_FILE <- opt$clinical_output_file


if (is.null(GSE_SERIES_MATRIX_FILE) || is.null(CLINICAL_OUTPUT_FILE)) {
  stop("Both --gse_series_matrix_file and --clinical_output_file are required.")
}

if (!file.exists(GSE_SERIES_MATRIX_FILE)) {
  stop("gse_series_matrix_file not found: ", GSE_SERIES_MATRIX_FILE)
}

#######################
## Load clinical data ##
#######################
gset <- getGEO(filename = GSE_SERIES_MATRIX_FILE)
pheno <- pData(gset)

required_cols <- c(
  "histology:ch1",
  "tumor (t) or normal (n):ch1",
  "gender:ch1",
  "age:ch1",
  "smoking:ch1",
  "vital date:ch1",
  "surgery date:ch1",
  "dead:ch1",
  "stage tnm:ch1"
)

missing_cols <- setdiff(required_cols, colnames(pheno))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

clin <- pheno[pheno[["histology:ch1"]] == 2, ]
clin <- clin[grep("T", clin[["tumor (t) or normal (n):ch1"]]), ]

clin$gender <- clin[["gender:ch1"]]
clin$age <- as.numeric(clin[["age:ch1"]])

clin$smoking_status <- clin[["smoking:ch1"]]
clin$smoking_status <- paste0("smoking_", clin$smoking_status)

clin$OS_time <- as.numeric(
  as.Date(clin[["vital date:ch1"]]) -
    as.Date(clin[["surgery date:ch1"]])
)

clin$OS_event <- clin[["dead:ch1"]]
clin$OS_event <- as.numeric(clin$OS_event)

#########################
## 5-year OS censoring ##
#########################
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

################
## Stage data ##
################
clin$stage <- clin[["stage tnm:ch1"]]

clin$stage <- gsub("1", "I", clin$stage)
clin$stage <- gsub("2", "I", clin$stage)
clin$stage <- gsub("3", "II", clin$stage)
clin$stage <- gsub("4", "II", clin$stage)
clin$stage <- gsub("5", "III", clin$stage)
clin$stage <- gsub("6", "III", clin$stage)
clin$stage <- gsub("7", "IV", clin$stage)

clin$stage <- paste0("Stage_", clin$stage)

clin$stage_modified <- clin$stage
# clin$stage_modified[clin$stage %in% c("Stage_I", "Stage_II")] <- "Stage_I-II"
# clin$stage_modified[
#   clin$stage %in% c("Stage_III", "Stage_IV")
# ] <- "Stage_III-IV"

clin$stage <- clin$stage_modified

clin$geo_accession <- clin[["tumor (t) or normal (n):ch1"]]

clin <- clin[, c(
  "geo_accession",
  "OS_event",
  "OS_time",
  "OS_event_5y",
  "OS_time_5y",
  "age",
  "gender",
  "stage",
  "smoking_status"
)]

write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


