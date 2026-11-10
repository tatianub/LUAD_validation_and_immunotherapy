## Load R packages
required_libraries <- c("data.table", "dplyr", "ggplot2",
                        "optparse", "GEOquery")
for (lib in required_libraries) {
  suppressPackageStartupMessages(
    library(lib, character.only = TRUE, quietly = TRUE)
  )
}

## Read arguments
option_list <- list(
  optparse::make_option(c("-g", "--clinical_file"),
                        type = "character", default = NULL,
                        help = "Path to clinical file."),
  optparse::make_option(c("-p", "--clinical_file_from_publication"),
                        type = "character", default = NULL,
                        help = "Path to clinical file from publication."),
  optparse::make_option(c("-c", "--clinical_output_file"),
                        type = "character", default = NULL,
                        help = "Path to output clinical data file.")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

CLINICAL_FILE <- opt$clinical_file
CLINICAL_FILE_FROM_PUBLICATION <- opt$clinical_file_from_publication
CLINICAL_OUTPUT_FILE <- opt$clinical_output_file


if (is.null(CLINICAL_FILE) || is.null(CLINICAL_FILE_FROM_PUBLICATION) ||
    is.null(CLINICAL_OUTPUT_FILE)) {
  stop("All three arguments --clinical_file, ",
       "--clinical_file_from_publication and ",
       "--clinical_output_file are required.")
}
if (!file.exists(CLINICAL_FILE)) {
  stop("clinical_file not found: ", CLINICAL_FILE)
}
if (!file.exists(CLINICAL_FILE_FROM_PUBLICATION)) {
  stop("clinical_file_from_publication not found: ",
       CLINICAL_FILE_FROM_PUBLICATION)
}

clinical_from_publication <- fread(CLINICAL_FILE_FROM_PUBLICATION)
chemo_col <- "Adj. CHEMOTHERAPY (N=no; Y=yes; ND: not known)"
radio_col <- "Adj. RADIOTHERAPY (N=no; Y=yes; ND: not known)"
if (!(chemo_col %in% colnames(clinical_from_publication))) {
  stop("chemo column not found in publication file")
}
if (!(radio_col %in% colnames(clinical_from_publication))) {
  stop("radio column not found in publication file")
}

cat("Loading clinical RData...\n")
clin_env <- new.env()
load(CLINICAL_FILE, envir = clin_env)
if (!"pheno" %in% ls(envir = clin_env)) {
  stop("Expected object 'pheno' not found in RData")
}
clin <- clin_env[["pheno"]]

if ("source_name_ch1" %in% colnames(clin)) {
  clin <- clin[!clin$source_name_ch1 == "Non Tumoral Lung", ]
}

if (!("characteristics_ch1.3" %in% colnames(clin))) {
  stop("characteristics_ch1.3 column missing in clinical data")
}
clin <- clin[grep("ADC|BAS", clin$characteristics_ch1.3), ]

if (!("description" %in% colnames(clin))) {
  stop("description column missing in clinical data for ID matching")
}
id_match <- match(clin$description, clinical_from_publication$Tumor_ID)
clin$adjuvant_chemotherapy <-
  clinical_from_publication[[chemo_col]][id_match]
clin$adjuvant_radiotherapy <-
  clinical_from_publication[[radio_col]][id_match]

clin <- clin %>%
  mutate(treatment = case_when(
    adjuvant_chemotherapy == "Y" |
    adjuvant_radiotherapy == "Y" ~ "treatment",
    adjuvant_chemotherapy == "N" &
    adjuvant_radiotherapy == "N" ~ "None",
    TRUE ~ "Unknown"
  )) %>%
  mutate(treatment = factor(treatment,
                            levels = c("None", "treatment", "Unknown")))

clin$treatment[clin$treatment == "Unknown"] <- NA

if (!("status:ch1" %in% colnames(clin))) {
  stop("'status:ch1' missing in clinical data")
}
clin$OS_event <- ifelse(clin[["status:ch1"]] == "ALIVE", 0,
                        ifelse(clin[["status:ch1"]] == "DEAD",
                               1, NA))
clin$OS_time <- as.numeric(clin[["follow-up time (months):ch1"]])

cutoff_years <- 5
time_unit <- NULL

cat("OS_time summary:\n")
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
  cat("Auto-detected time_unit:", time_unit, "\n")
}

cutoff <- switch(time_unit,
                 days = cutoff_years * 365.25,
                 months = cutoff_years * 12,
                 years = cutoff_years,
                 stop("Unknown time_unit"))

cat("Using cutoff:", cutoff, time_unit, "\n")

clin$OS_event_5y <- clin$OS_event
clin$OS_time_5y <- as.numeric(clin$OS_time)
clin$OS_event_5y[clin$OS_time_5y > cutoff] <- 0

clin$DFS_event <- clin[["relapse (event=1; no event=0):ch1"]]
clin$DFS_time <- clin[["disease free survival in months:ch1"]]

colnames(clin) <- gsub("\\:|\\ ", "_", colnames(clin))

clin$age <- clin$age_at_surgery_ch1
clin$gender <- clin$gender_ch1
clin$gender[clin$gender == "M"] <- "male"
clin$gender[clin$gender == "F"] <- "female"

if (all(c("pm_stage_ch1", "pn_stage_ch1") %in% colnames(clin))) {
  clin$stage <- with(clin,
                     ifelse(pm_stage_ch1 == "M1", "Stage_III-IV",
                       ifelse(pn_stage_ch1 %in% c("N1", "N2"),
                              "Stage_III-IV",
                              ifelse(pn_stage_ch1 == "N0",
                                     "Stage_I-II", NA))))
} else {
  clin$stage <- NA
}

clin <- clin[, c("geo_accession", "OS_event", "OS_time",
                 "OS_event_5y", "OS_time_5y", "DFS_event",
                 "DFS_time", "age", "gender", "stage",
                 "treatment")]

write.table(clin, file = CLINICAL_OUTPUT_FILE, sep = "\t",
            quote = FALSE, row.names = FALSE)



