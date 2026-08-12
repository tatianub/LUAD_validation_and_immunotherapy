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
        c("-p", "--panda_network_file"),
        type = "character",
        default = NULL,
        help = "Path to the panda network file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--output_edge_file"),
        type = "character",
        default = NULL,
        help = "Path to the output edge file.",
        metavar = "character"
    )
)

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Assign parsed arguments to variables
PANDA_NETWORK_FILE <- opt$panda_network_file
OUTPUT_EDGE_FILE <- opt$output_edge_file

panda_network <- fread(PANDA_NETWORK_FILE)
edges <- data.table("reg" = panda_network$V1, "tar" = panda_network$V2)
write.table(edges, OUTPUT_EDGE_FILE,
            col.names = T, row.names= F, sep="\t", quote = F)
