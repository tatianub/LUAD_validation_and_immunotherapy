#####################
## Load R packages ##
#####################
required_libraries <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "optparse",
  "GEOquery",
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
    c("-g", "--clinical_file"),
    type = "character",
    default = NULL,
    help = "Path to clinical file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("-p", "--clinical_file_from_publication"),
    type = "character",
    default = NULL,
    help = "Path to clinical file from publication.",
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
CLINICAL_FILE <- opt$clinical_file
CLINICAL_FILE_FROM_PUBLICATION <- opt$clinical_file_from_publication
CLINICAL_OUTPUT_FILE <- opt$clinical_output_file

########################
## Validate arguments ##
########################
if (is.null(CLINICAL_FILE)) {
  stop("--clinical_file is required.", call. = FALSE)
}

if (is.null(CLINICAL_FILE_FROM_PUBLICATION)) {
  stop("--clinical_file_from_publication is required.", call. = FALSE)
}

if (is.null(CLINICAL_OUTPUT_FILE)) {
  stop("--clinical_output_file is required.", call. = FALSE)
}

if (!file.exists(CLINICAL_FILE)) {
  stop("clinical_file not found: ", CLINICAL_FILE, call. = FALSE)
}

if (!file.exists(CLINICAL_FILE_FROM_PUBLICATION)) {
  stop(
    "clinical_file_from_publication not found: ",
    CLINICAL_FILE_FROM_PUBLICATION,
    call. = FALSE
  )
}


clin_from_article <- readxl::read_excel(
  CLINICAL_FILE_FROM_PUBLICATION,
  sheet = "TableS1",
  skip = 1
)

head(clin_from_article)
colnames(clin_from_article)
table(clin_from_article$"TKI.treatment")
table(clin_from_article$"Chemo.treatment")

table(clin_from_article$"TKI.treatment", clin_from_article$"Chemo.treatment")
clin_from_article$"Patient.ID"
clin_from_article$"TKI.treatment"

clin_from_article <- clin_from_article[
  -grep("BGI-", clin_from_article$"Patient.ID"),
]

tki_yes <- !is.na(clin_from_article$"TKI.treatment") &
  tolower(trimws(clin_from_article$"TKI.treatment")) == "yes"

chemo_yes <- !is.na(clin_from_article$"Chemo.treatment") &
  tolower(trimws(clin_from_article$"Chemo.treatment")) == "yes"

ids_treatments <- clin_from_article$"Patient.ID"[tki_yes | chemo_yes]
ids_treatments <- ids_treatments[!is.na(ids_treatments)]
length(unique(ids_treatments))
ids_treatments <- paste0(ids_treatments, "T")
ids_treatments


# read clinical file
clin <- fread(CLINICAL_FILE)
head(clin)
str(clin)

# select only tumor
clin <- clin[grep("LUAD", clin$type), ]

clin$treatment <- ifelse(
  clin$Sample_ID %in% ids_treatments,
  "treatment",
  "no treatment"
)

clin$treatment

clin$OS_event <- clin$OS_status
clin$OS_event <- ifelse(
  clin$OS_event == "Dead",
  1,
  ifelse(clin$OS_event == "Alive", 0, NA)
)

clin$OS_time <- as.numeric(clin$OS)

clin$OS_event_5y <- clin$OS_status
clin$OS_time_5y <- as.numeric(clin$OS)

cutoff_years <- 5

# Set time_unit to match the unit of OS_time in the data:
# "days", "months", or "years"
# Set to NULL to auto-detect based on the range of OS_time values
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
  days = cutoff_years * 365.25,
  months = cutoff_years * 12,
  years = cutoff_years,
  stop("Unknown time_unit: must be 'days', 'months', or 'years'")
)

cat("Using cutoff:", cutoff, time_unit, "(=", cutoff_years, "years)\n")

# Censor patients with OS_time > cutoff
clin$OS_event_5y[clin$OS_time_5y > cutoff] <- 0

clin$gender <- clin$Gender
clin$gender[clin$gender == "Male"] <- "male"
clin$gender[clin$gender == "Female"] <- "female"

clin$smoking_status <- clin$Smoking_status
clin$smoking_status[clin$smoking_status == "y"] <- "smoking"
clin$smoking_status[clin$smoking_status == "n"] <- "non-smoking"

clin$age <- clin$Age
clin$stage <- clin$Stage

clin$geo_accession <- clin$Sample_ID
clin$geo_accession <- gsub("\\.", "-", clin$geo_accession)

clin$stage[clin$stage == ""] <- NA
clin$stage

clin$stage <- ifelse(
  is.na(clin$stage),
  NA,
  paste0("Stage_", clin$stage)
)

table(clin$stage)
table(clin$OS_event)

# clin$stage[clin$stage %in% c("Stage_I", "Stage_II")] <- "Stage_I_II"
# clin$stage[clin$stage %in% c("Stage_III", "Stage_IV")] <- "Stage_III_IV"

# table(clin$stage)

clin <- clin[!is.na(clin$stage), ]

clin <- clin[, c(
  "geo_accession", "OS_time", "OS_event",
  "OS_time_5y", "OS_event_5y",
  "gender", "smoking_status", "age",
  "stage", "treatment"
)]

write.table(
  clin,
  file = CLINICAL_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
