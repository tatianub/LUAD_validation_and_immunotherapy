#####################
## Load R packages ##
#####################
suppressPackageStartupMessages(library(optparse))


####################
## Read arguments ##
####################
option_list <- list(
  make_option(
    "--series_matrix",
    type = "character"
  ),
  make_option(
    "--output_file",
    type = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

OUTPUT_FILE <- opt$output_file
########################
## Read series matrix ##
########################
# Read GEO series matrix
lines <- readLines(opt$series_matrix, warn = FALSE)

# Extract sample metadata rows
geo_line <- grep("^!Sample_geo_accession\\t", lines, value = TRUE)
title_line <- grep("^!Sample_title\\t", lines, value = TRUE)

if (length(geo_line) != 1 || length(title_line) != 1) {
  stop("Could not find unique sample accession and title rows.")
}


################################
## Extract sample information ##
################################
# Split tab-separated values and remove quotation marks
extract_values <- function(line) {
  values <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", values)
}

geo_ids <- extract_values(geo_line)
sample_names <- extract_values(title_line)

if (length(geo_ids) != length(sample_names)) {
  stop("Number of GEO IDs does not match number of sample names.")
}

# Create sample annotation table
sample_annotation <- data.frame(
  GEO_ID = geo_ids,
  Sample = sample_names,
  stringsAsFactors = FALSE
)


#####################
## Save annotation ##
#####################
# Write without column names
write.table(
  sample_annotation,
  file = OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Sample annotation generated successfully.\n")
cat("Number of samples:", nrow(sample_annotation), "\n")
cat("Output:", OUTPUT_FILE, "\n")