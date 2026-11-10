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
  stop("Expected object 'pheno' not found in RData
  file: ", CLINICAL_FILE)
}

clin <- clin_env[["pheno"]]
clin$OS_event <- clin$"status:ch1"
clin$OS_event <- 
    ifelse(clin$OS_event == "alive", 0, ifelse(clin$OS_event == "dead", 1, NA))
clin$OS_time <- clin$"survival time:ch1"
clin$DFS_time <- clin$"disease-free survival time:ch1"
clin$DFS_event <- clin$"recurrence:ch1"
clin$DFS_event <- 
    ifelse(clin$DFS_event == "Y", 1, ifelse(clin$DFS_event == "N", 0, NA))

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


# Create a new variable for 5-year overall survival event
clin$OS_event_5y <- ifelse(clin$OS_time <= cutoff & clin$OS_event == 1, 1, 0)
clin$OS_time_5y <- clin$OS_time
colnames(clin) <- gsub("[: -]", "_", colnames(clin))
clin$stage <- clin$"Stage_ch1"
clin$stage <- gsub("A|B", "", clin$stage)
clin$stage
clin$stage <- gsub("1", "I", clin$stage)
clin$stage <- gsub("2", "II", clin$stage)
clin$stage <- paste("Stage", clin$stage, sep = "_")
clin <- clin[clin$histology_ch1 == "adenocarcinoma", ]
clin$age <- as.numeric(clin$age_ch1)

clin$gender  <- clin$Sex_ch1
clin$smoking_status <- clin$"smoking_ch1"
clin$gender[clin$gender == "M"] <- "male"
clin$gender[clin$gender == "F"] <- "female"


clin <- clin[, c("geo_accession", "OS_time", 
                "OS_event", "DFS_time", "DFS_event", 
                "OS_time_5y", "OS_event_5y",
                 "age", "gender", "smoking_status", "stage")]

write.table(clin, 
          file = CLINICAL_OUTPUT_FILE, 
          sep = "\t", quote = FALSE, row.names = FALSE)



