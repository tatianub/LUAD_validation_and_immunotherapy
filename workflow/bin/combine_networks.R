
#####################
## Load R packages ##
#####################
required_libraries <- c("gtools", "purrr", "tidyverse", "data.table", 
                        "optparse", "plyr")
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
        c("-f", "--network_dir"),
        type = "character",
        default = NULL,
        help = "Path to the directory containing LIONESS network TXT files.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-i", "--lioness_sample_mapping"),
        type = "character",
        default = NULL,
        help = "Path to the LIONESS sample mapping file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to output file with combined network RData file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

NETWORK_DIR <- opt$network_dir
LIONESS_SAMPLE_MAPPING <- opt$lioness_sample_mapping
OUTPUT_FILE <- opt$output_file

filelist <- list.files(NETWORK_DIR, pattern = "^lioness\\.[0-9]+\\.txt$", full.names = TRUE)
cat("Starting LIONESS network combination and saving...\n")
cat(sprintf("Sample mapping file: %s\n", LIONESS_SAMPLE_MAPPING))
cat(sprintf("Output file: %s\n", OUTPUT_FILE))

####################
## Load data      ##
####################

source("workflow/bin/combine_networks_fn.R")

# Parse the comma-separated list of network files
cat("Parsing LIONESS network file list...\n")

if (length(filelist) == 0) {
    stop("Error: No LIONESS network files provided")
}

# Check if files exist
missing_files <- filelist[!file.exists(filelist)]
if (length(missing_files) > 0) {
    stop(sprintf("Error: The following files do not exist:\n%s", 
                 paste(missing_files, collapse = "\n")))
}

cat(sprintf("Found %d LIONESS network files in directory\n", length(filelist)))

# Load LIONESS network manifest
cat("Loading LIONESS sample mapping...\n")
if (!file.exists(LIONESS_SAMPLE_MAPPING)) {
    stop(sprintf("Cannot find LIONESS sample mapping file: %s", 
                 LIONESS_SAMPLE_MAPPING))
}
info_net <- fread(LIONESS_SAMPLE_MAPPING)

# Get all unique cancer types

file_basenames <- basename(filelist)

net <- read_networks_txt(info_net, start = 1, end = nrow(info_net))
save(net, file = OUTPUT_FILE)
cat("Combined network saved to:", OUTPUT_FILE, "\n")


