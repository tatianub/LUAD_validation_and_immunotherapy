#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "tidyverse", "purrr", "optparse",
                        "preprocessCore")
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
        c("-n", "--network_file"),
        type = "character",
        default = NULL,
        help = "Path to a network file (RData).",
        metavar = "character"
    ),

    optparse::make_option(
        c("-o", "--output_file"),
        type = "character",
        default = NULL,
        help = "Path to the quantile normalzed network output file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
NETWORK_FILE <- opt$network_file
OUTPUT_FILE <- opt$output_file

source("workflow/bin/quantile_normalize_networks_fn.R")
####################
## Process networks ##
####################
net_norm <- quantile_normalize_net(NETWORK_FILE)
cat(sprintf("  Saving normalized network: %s\n", OUTPUT_FILE))
save(net_norm, file = OUTPUT_FILE)
rm(net_norm)
