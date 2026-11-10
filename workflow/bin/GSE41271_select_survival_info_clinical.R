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
clin <- clinical[clinical[["histology:ch1"]] == "Adenocarcinoma", ]

clin[["date of surgery:ch1"]] <- as.Date(clin[["date of surgery:ch1"]])
clin[["last follow-up survival:ch1"]] <- as.Date(
  clin[["last follow-up survival:ch1"]]
)

clin$OS_event <- ifelse(
  is.na(clin[["vital statistics:ch1"]]),
  NA,
  ifelse(clin[["vital statistics:ch1"]] == "D", 1, 0)
)

clin$OS_time <- as.numeric(
  clin[["last follow-up survival:ch1"]] -
    clin[["date of surgery:ch1"]]
)

clin$smoking_status <- clin[["tobacco history:ch1"]]

clin[["date of birth:ch1"]] <- as.Date(clin[["date of birth:ch1"]])

clin$age <- floor(
  as.numeric(clin[["date of surgery:ch1"]] -
               clin[["date of birth:ch1"]]) / 365.25
)

clin$gender <- factor(
  clin[["gender:ch1"]],
  levels = c("F", "M"),
  labels = c("female", "male")
)

clin[["last follow-up recurrence:ch1"]] <- as.Date(
  clin[["last follow-up recurrence:ch1"]]
)

clin$RFS_event <- ifelse(
  is.na(clin[["recurrence:ch1"]]),
  NA,
  ifelse(clin[["recurrence:ch1"]] == "Y", 1, 0)
)

clin$RFS_time <- as.numeric(
  clin[["last follow-up recurrence:ch1"]] -
    clin[["date of surgery:ch1"]]
)

clin$stage <- trimws(clin[["final patient stage:ch1"]])
clin$stage[clin$stage == "IA Vs IIIB"] <- NA
clin <- clin[!is.na(clin$stage), ]

clin$stage <- gsub("A|B", "", clin$stage)
clin$stage <- paste0("Stage_", clin$stage)

clin$geo_accession <- rownames(clin)


#########################
## 5-year OS censoring ##
#########################
# Set to NULL to auto-detect based on the range of OS_time values
cutoff_years <- 5
time_unit <- NULL

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

#########################
## Recode broad stage  ##
#########################
clin$stage[clin$stage %in% c("stage_III", "stage_IV")] <- "stage_III_IV"
clin$stage[clin$stage %in% c("stage_I", "stage_II")] <- "stage_I_II"

########################
## Select final fields ##
########################
clin <- clin[, c(
  "geo_accession",
  "OS_time",
  "OS_event",
  "RFS_time",
  "RFS_event",
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

