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

clin <- data.table::fread(CLINICAL_FILE, data.table = FALSE)

required_cols <- c(
  "type", "OS_status", "OS", "Gender", "Smoking_status",
  "PFS", "PFS_status", "progression_5y", "Age",
  "Cytologic_grade", "Sample_ID", "Stage"
)

missing_cols <- setdiff(required_cols, colnames(clin))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# gender male - 1, female - 2
# select only tumor
clin <- clin[grep("LUAD", clin$type), ]

clin$OS_event <- clin$OS_status
clin$OS_time <- as.numeric(clin$OS)

clin$OS_event_5y <- clin$OS_status
clin$OS_time_5y <- as.numeric(clin$OS)

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

clin$gender <- clin$Gender
clin$gender[clin$gender == "1"] <- "male"
clin$gender[clin$gender == "2"] <- "female"

clin$smoking_status <- clin$Smoking_status
clin$smoking_status[clin$smoking_status == "y"] <- "smoking"
clin$smoking_status[clin$smoking_status == "n"] <- "non-smoking"

clin$PFS_time <- as.numeric(clin$PFS)
clin$PFS_event <- clin$PFS_status

clin$PFS_time_5y <- as.numeric(clin$PFS)
clin$PFS_event_5y <- clin$progression_5y
clin$PFS_event_5y <- ifelse(
  clin$PFS_event_5y == "Progression", 1,
  ifelse(clin$PFS_event_5y == "No Progression (>5y)", 0, NA)
)

clin$age <- clin$Age
clin$age[clin$age == "40-49"] <- 45
clin$age[clin$age == "50-59"] <- 55
clin$age[clin$age == "60-69"] <- 65
clin$age[clin$age == "70-79"] <- 75
clin$age[clin$age == "80-89"] <- 85
clin$age <- as.numeric(clin$age)

clin$grade <- clin$Cytologic_grade
clin$geo_accession <- clin$Sample_ID

clin$PFS_event <- ifelse(
  clin$PFS_event == "y", 1,
  ifelse(clin$PFS_event == "n", 0, NA)
)

clin$geo_accession <- gsub("\\.", "-", clin$geo_accession)

clin$stage <- clin$Stage
clin$stage[clin$stage == "1.0"] <- "Stage_I"
clin$stage[clin$stage == "IA"] <- "Stage_I"
clin$stage[clin$stage == "IB"] <- "Stage_I"


# clin$stage[clin$stage == "IA"] <- "Stage_IA"
# clin$stage[clin$stage == "IB"] <- "Stage_IB"

clin <- clin[, c(
  "geo_accession", "OS_time", "OS_event", "PFS_time",
  "PFS_event", "OS_time_5y", "OS_event_5y",
  "PFS_time_5y", "PFS_event_5y",
  "gender", "smoking_status",
  "age", "stage"
)]

write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

