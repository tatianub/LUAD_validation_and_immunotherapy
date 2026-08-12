#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "optparse")
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
        c("-s", "--samples_panda_file"),
        type = "character",
        default = NULL,
        help = "Path to the PANDA samples file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the output LIONESS sample mapping file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

NETWORK_DIR <- opt$network_dir
SAMPLES_PANDA_FILE <- opt$samples_panda_file
OUTPUT_FILE <- opt$output_file

####################
cat("Creating LIONESS sample mapping file...\n")
cat(sprintf("PANDA samples file: %s\n", SAMPLES_PANDA_FILE))
cat(sprintf("Output file: %s\n", OUTPUT_FILE))


####################
## Process data   ##
####################

# Parse the comma-separated list of network files
cat("Parsing LIONESS network file list...\n")
filelist <- list.files(NETWORK_DIR, pattern = "^lioness\\.[0-9]+\\.txt$", full.names = TRUE)

if (length(filelist) == 0) {
    cat("Error: No LIONESS network files provided\n")
    quit(status = 1)
}

# Check if files exist
missing_files <- filelist[!file.exists(filelist)]
if (length(missing_files) > 0) {
    cat("Error: The following files do not exist:\n")
    cat(paste(missing_files, collapse = "\n"), "\n")
    quit(status = 1)
}

cat(sprintf("Found %d LIONESS network files\n", length(filelist)))

# Create file information table with proper path extraction
cat("Processing file paths and names...\n")
filelist_dat <- data.table(
    "file" = basename(filelist),
    "file_path" = filelist
)

# Extract network numbers and sort by them
numbers <- as.numeric(regmatches(filelist_dat$file,
                                 regexpr("[0-9]+", filelist_dat$file)))

if (any(is.na(numbers))) {
    cat("Warning: Some files don't contain numeric identifiers\n")
    cat("Files without numbers will be placed at the end\n")
    numbers[is.na(numbers)] <- max(numbers, na.rm = TRUE) +
                               seq_along(which(is.na(numbers)))
}

filelist_dat <- filelist_dat[order(numbers), ]
cat(sprintf("Ordered %d files by network number\n", nrow(filelist_dat)))

# Load PANDA samples file
cat("Loading PANDA samples information...\n")
samples <- fread(SAMPLES_PANDA_FILE)
head(samples)
# Validate samples file structure
colnames(samples) <- c("number", "sample_id")
cat(sprintf("Loaded %d samples from PANDA file\n", nrow(samples)))


# Generate expected LIONESS network filenames
samples$lioness_net <- paste0("lioness.", seq_len(nrow(samples)), ".txt")

# Match files to samples
cat("Mapping LIONESS files to samples...\n")
filelist_dat$sample_id <-
    samples$sample_id[match(filelist_dat$file, samples$lioness_net)]


####################
## Save output    ##
####################

cat("Writing LIONESS sample mapping file...\n")
write.table(filelist_dat, OUTPUT_FILE,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            quote = FALSE)
