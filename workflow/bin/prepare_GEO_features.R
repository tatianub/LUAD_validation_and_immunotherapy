#####################
## Load R packages ##
#####################
library(optparse)


####################
## Read arguments ##
####################
option_list <- list(
  make_option(
    "--dataset",
    type = "character"
  ),
  make_option(
    "--probe_file",
    type = "character"
  ),
  make_option(
    "--output_file",
    type = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#################
## Load probes ##
#################
# load the probe file into a separate environment
probes <- new.env()
load(opt$probe_file, envir = probes)

features <- probes[["feature"]]


###################################
## Map probe IDs to gene symbols ##
###################################
# Same mappings as in the original script
if (opt$dataset %in% c("GSE41271", "GSE42127")) {
  features$PROBEID <- rownames(features)
  features$SYMBOL <- features$ILMN_Gene

} else if (opt$dataset == "GSE13213") {
  features$PROBEID <- rownames(features)
  features$SYMBOL <- features$GENE_SYMBOL

} else {
  stop("Unsupported dataset: ", opt$dataset)
}

# keep probe IDs and gene symbols only
features <- features[, c("PROBEID", "SYMBOL")]


#####################
## Save annotation ##
#####################
write.table(
  features,
  file = opt$output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)