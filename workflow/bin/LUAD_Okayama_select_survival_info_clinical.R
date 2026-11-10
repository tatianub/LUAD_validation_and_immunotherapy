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

#######################

gset <- getGEO(
  filename = GSE_SERIES_MATRIX_FILE 
)

expr <- exprs(gset)
pheno <- pData(gset)
# filter out the normal lung samples
clin <- pheno[grep("primary", pheno$characteristics_ch1), ]


clin$OS_time <- as.numeric(clin$"days before death/censor:ch1")
clin$OS_event <- ifelse(is.na(clin$"death:ch1"), NA, 
                   ifelse(clin$"death:ch1" == "dead", 1, 0))

clin$gender <- clin$"gender:ch1"
clin$age <- as.numeric(clin$"age (years):ch1")
clin$smoking_status <- clin$"smoking status:ch1"

clin$RFS_time <- clin$"days before relapse/censor:ch1"
clin$RFS_event <- clin$"relapse:ch1"
clin$RFS_event <- ifelse(is.na(clin$RFS_event), NA, 
                   ifelse(clin$RFS_event == "relapsed", 1, 0))

clin$gene_alteration <- clin$"gene alteration status:ch1"
clin$stage <- clin$"pathological stage:ch1"
clin$stage <- paste0("Stage_", clin$stage)
clin$stage[clin$stage == "Stage_IA"] <- "Stage_I"
clin$stage[clin$stage == "Stage_IB"] <- "Stage_I"

# exclude patients with incomplete resection or adjuvant therapy
# (this is from the original paper, see methods section)
exclude_col <- paste0(
  "exclude for prognosis analysis due to ",
  "incomplete resection or adjuvant therapy:ch1"
)
clin <- clin[!clin[[exclude_col]] == "exclude", ]
dim(clin)

cutoff_years <- 5
# Set time_unit to match the unit of OS_time in the data:
# "days", "months", or "years"
# Set to NULL to auto-detect based on the range of OS_time values
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

clin$OS_event_5y <- ifelse(clin$OS_time > cutoff, 0, clin$OS_event)
clin$OS_event <- clin$OS_event_5y
clin$OS_time_5y <- clin$OS_time
dim(clin)
table(clin$stage)
table(clin$OS_event)
table(clin$RFS_event)
table(clin$gender)
median(as.numeric(clin$OS_time))
table(clin$gene_alteration)
clin <- clin[, c("geo_accession", "OS_time", "OS_event",
              "RFS_time", "RFS_event", "OS_time_5y", "OS_event_5y",
              "gender", "age", "smoking_status", "stage", "gene_alteration")]

table(clin$OS_event)

write.table(clin, 
          file = CLINICAL_OUTPUT_FILE, 
          sep = "\t", quote = FALSE, row.names = FALSE)

