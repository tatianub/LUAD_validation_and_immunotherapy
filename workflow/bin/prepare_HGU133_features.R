#####################
## Load R packages ##
#####################
suppressPackageStartupMessages({
  library(optparse)
  library(AnnotationDbi)
  library(hgu133plus2.db)
})


####################
## Read arguments ##
####################
option_list <- list(
  make_option(
    "--output_file",
    type = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#######################
## Initialize values ##
#######################
OUTPUT_FILE <- opt$output_file


###############################
## Retrieve probe annotation ##
###############################
# Retrieve Affymetrix probe annotations
annot <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = keys(hgu133plus2.db, keytype = "PROBEID"),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

# Keep probe IDs and gene symbols
features <- annot[, c("PROBEID", "SYMBOL")]


#####################
## Save annotation ##
#####################
write.table(
  features,
  file = OUTPUT_FILE,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)