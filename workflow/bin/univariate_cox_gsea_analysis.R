#####################
## Load R packages ##
#####################
required_libraries <- c( 
  "data.table",
  "dplyr",
  "ggplot2",
  "optparse",
  "survival",
  "ggpubr",
  "grid",
  "tools",
  "fgsea"
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
    c("--indegree_file"),
    type = "character",
    default = NULL,
    help = "Path to indegree file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--clinical_file"),
    type = "character",
    default = NULL,
    help = "Path to clinical file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--pathway_gmt_file"),
    type = "character",
    default = NULL,
    help = "Path to pathway GMT file (.gmt).",
    metavar = "character"
  ),
  optparse::make_option(
    c("--covariates"),
    type = "character",
    default = NULL,
    help = paste(
      "Comma-separated list of covariate column names to adjust",
      "for in Cox models (e.g. 'age,gender,stage')."
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("--survival_outcomes"),
    type = "character",
    default = NULL,
    help = paste(
      "Comma-separated list of outcome names to test",
      "(e.g. 'OS,PFS'). If not specified, all detected outcomes",
      "in clinical data are tested."
    ),
    metavar = "character"
  ),
  optparse::make_option(
    c("--cores"),
    type = "integer",
    default = 1,
    help = "Number of CPU cores for parallel processing (default: 1).",
    metavar = "integer"
  ),
  optparse::make_option(
    c("--results_output_file"),
    type = "character",
    default = NULL,
    help = "Path to the results output file.",
    metavar = "character"
  ),
  optparse::make_option(
    c("--gsea_seed"),
    type = "integer",
    default = 12345L,
    help = "Random seed for GSEA.",
    metavar = "integer"
  ),
  optparse::make_option(
    c("--gsea_sample_size"),
    type = "integer",
    default = 1001L,
    help = "Sample size for GSEA.",
    metavar = "integer"
  ),
  optparse::make_option(
    c("--gsea_min_size"),
    type = "integer",
    default = 10L,
    help = "Minimum gene set size for GSEA.",
    metavar = "integer"
  ),
  optparse::make_option(
    c("--gsea_max_size"),
    type = "integer",
    default = 500L,
    help = "Maximum gene set size for GSEA.",
    metavar = "integer"
  )
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

set.seed(GSEA_SEED)

INDEGREE_FILE <- opt$indegree_file
CLINICAL_FILE <- opt$clinical_file
PATHWAY_GMT_FILE <- opt$pathway_gmt_file
COVARIATES <- opt$covariates
OUTCOMES <- opt$survival_outcomes
NCORES <- as.integer(opt$cores)
RESULTS_OUTPUT_FILE <- opt$results_output_file

# GSEA settings
GSEA_SEED <- as.integer(opt$gsea_seed)
GSEA_SAMPLE_SIZE <- as.integer(opt$gsea_sample_size)
GSEA_MIN_SIZE <- as.integer(opt$gsea_min_size)
GSEA_MAX_SIZE <- as.integer(opt$gsea_max_size)

if (is.na(NCORES) || NCORES < 1L) {
  NCORES <- 1L
}

cat("Arguments loaded successfully\n")
cat("INDEGREE_FILE       :", INDEGREE_FILE, "\n")
cat("CLINICAL_FILE       :", CLINICAL_FILE, "\n")
cat("PATHWAY_GMT_FILE    :", PATHWAY_GMT_FILE, "\n")
cat("COVARIATES          :",
    paste(COVARIATES, collapse = ", "), "\n")
cat("OUTCOMES            :", OUTCOMES, "\n")
cat("NCORES              :", NCORES, "\n")
cat("RESULTS_OUTPUT_FILE :", RESULTS_OUTPUT_FILE, "\n")

#########################
## Load helper functions
#########################

source("workflow/bin/helper_fn.R")

outcome_args <- if (!is.null(OUTCOMES)) {
  trimws(strsplit(OUTCOMES, ",")[[1]])
} else {
  NULL
}

#########################
## Read input data
#########################

indegree <- fread(INDEGREE_FILE)
clinical <- fread(CLINICAL_FILE)

colnames(indegree) <- trimws(colnames(indegree))
clinical$geo_accession <- trimws(
  as.character(clinical$geo_accession)
)

if (anyDuplicated(clinical$geo_accession)) {
  stop("Duplicated sample IDs in clinical data.")
}

if (!"tar" %in% colnames(indegree)) {
  stop("Column 'tar' not found in indegree data.")
}

matched_samples <- intersect(
  colnames(indegree),
  clinical$geo_accession
)

if (length(matched_samples) == 0) {
  stop("No overlapping sample IDs.")
}

indegree_clinical <- indegree[, ..matched_samples]

expression_matrix <- as.matrix(indegree_clinical)
rownames(expression_matrix) <- indegree$tar
storage.mode(expression_matrix) <- "numeric"

clinical_indegree <- clinical[
  match(
    colnames(expression_matrix),
    clinical$geo_accession
  ),
]

stopifnot(
  identical(
    colnames(expression_matrix),
    clinical_indegree$geo_accession
  )
)

cat("Matched samples:", ncol(expression_matrix), "\n")
cat("Number of genes:", nrow(expression_matrix), "\n")

#########################
## Detect survival outcomes
#########################

available_outcomes <- find_survival_pairs(
  clinical_indegree
)

if (nrow(available_outcomes) == 0) {
  stop("No survival outcome pairs found.")
}

if (!is.null(outcome_args)) {

  outcomes_to_test <- available_outcomes %>%
    filter(.data$endpoint %in% outcome_args)

  if (nrow(outcomes_to_test) == 0) {
    stop(
      "Requested outcomes not found. Available: ",
      paste(
        unique(available_outcomes$endpoint),
        collapse = ", "
      )
    )
  }

} else {
  outcomes_to_test <- available_outcomes
}

cat(
  "Testing outcomes:",
  paste(
    unique(outcomes_to_test$endpoint),
    collapse = ", "
  ),
  "\n"
)

#########################
## Prepare covariates
#########################

all_covariates <- if (!is.null(COVARIATES)) {
  COVARIATES[
    COVARIATES %in% colnames(clinical_indegree)
  ]
} else {
  character(0)
}

cat(
  "Using covariates:",
  paste(all_covariates, collapse = ", "),
  "\n"
)

clinical_data <- clinical_indegree %>%
  rename(sample_id = geo_accession)

#########################
## Run Cox regression
#########################

all_outcome_results <- list()

gene_names <- rownames(expression_matrix)

for (i in seq_len(nrow(outcomes_to_test))) {

  outcome_name <- outcomes_to_test$endpoint[i]
  time_col <- outcomes_to_test$time_col[i]
  event_col <- outcomes_to_test$event_col[i]

  cat(
    "\n=== Processing outcome:",
    outcome_name,
    "(", time_col, ",", event_col, ") ===\n"
  )

  results_list <- if (NCORES > 1) {

    parallel::mclapply(
      gene_names,
      run_one_gene,
      mc.cores = NCORES
    )

  } else {

    lapply(gene_names, run_one_gene)

  }

  outcome_results <- bind_rows(results_list)

  all_outcome_results[[length(all_outcome_results) + 1]] <-
    outcome_results

  cat("Status breakdown:\n")
  print(table(
    outcome_results$status,
    useNA = "ifany"
  ))

  cat(
    "Successful Cox models:",
    sum(outcome_results$status == "ok", na.rm = TRUE),
    "of",
    nrow(expression_matrix),
    "\n"
  )
}

#########################
## Combine Cox results
#########################

all_results <- bind_rows(all_outcome_results)

final_results <- all_results %>%
  filter(status == "ok") %>%
  select(
    outcome,
    gene,
    z_score,
    p_value,
    hazard_ratio
  )

if (nrow(final_results) == 0) {
  stop("No valid Cox models were fitted.")
}

cat("\n=== FINAL COX RESULTS ===\n")
print(table(final_results$outcome))

#########################
## Prepare ranked lists
#########################

ranked_lists <- final_results %>%
  group_split(outcome) %>%
  setNames(
    sapply(., function(x) unique(x$outcome))
  )

# Save Cox results to enable comparison between runs
cox_output_file <- paste0(
  RESULTS_OUTPUT_FILE,
  ".cox_rankings.rds"
)

saveRDS(ranked_lists, cox_output_file)

cat("Saved Cox rankings:", cox_output_file, "\n")

#########################
## Load pathways
#########################

mypathway_pathways <- gmtPathways(
  PATHWAY_GMT_FILE
)

cat(
  "Loaded",
  length(mypathway_pathways),
  "pathways\n"
)

# Test PD-1 pathway only
mypathway_pathways <- mypathway_pathways[
  grep("PD_1", names(mypathway_pathways))
]

if (length(mypathway_pathways) == 0) {
  stop("No PD-1 pathway found in GMT file.")
}

cat(
  "Selected pathways:",
  paste(names(mypathway_pathways), collapse = ", "),
  "\n"
)

#########################
## GSEA function
#########################

run_one_gsea <- function(
  outcome_name,
  seed = GSEA_SEED
) {

  set.seed(seed)

  cat(
    "\nRunning GSEA for",
    outcome_name,
    "with seed",
    seed,
    "...\n"
  )

  ranked_df <- ranked_lists[[outcome_name]] %>%
    select(gene, z_score) %>%
    filter(
      !is.na(gene),
      gene != "",
      !is.na(z_score),
      is.finite(z_score)
    )

  # Diagnose duplicated gene identifiers
  n_duplicate_genes <- sum(
    duplicated(ranked_df$gene)
  )

  cat(
    "Duplicated gene identifiers:",
    n_duplicate_genes,
    "\n"
  )

  # Preserve original strategy:
  # retain the largest absolute z-score per gene.
  # Add deterministic ordering for exact ties.
  ranked_df <- ranked_df %>%
    arrange(gene, desc(abs(z_score)), desc(z_score)) %>%
    group_by(gene) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(desc(z_score), gene)

  if (nrow(ranked_df) < 2) {
    cat("Too few genes for GSEA\n")
    return(data.frame())
  }

  ranks <- setNames(
    ranked_df$z_score,
    ranked_df$gene
  )

  if (length(unique(ranks)) < 2) {
    cat("Insufficient variation in z-scores\n")
    return(data.frame())
  }

  #########################
  ## Ranking diagnostics
  #########################

  n_ties <- sum(duplicated(ranks))

  pct_ties <- 100 * n_ties / length(ranks)

  cat("Number of ranked genes:", length(ranks), "\n")
  cat("Number of tied statistics:", n_ties, "\n")
  cat(
    "Percentage tied:",
    round(pct_ties, 2),
    "%\n"
  )

  if (pct_ties > 1) {
    warning(
      "More than 1% tied z-scores for ",
      outcome_name
    )
  }

  #########################
  ## Multilevel GSEA
  #########################

  fgseaRes <- tryCatch(

    fgseaMultilevel(
      pathways = mypathway_pathways,
      stats = ranks,
      minSize = GSEA_MIN_SIZE,
      maxSize = GSEA_MAX_SIZE,
      sampleSize = GSEA_SAMPLE_SIZE,
      eps = 0,
      nproc = 1
    ),

    error = function(e) {
      cat(
        "GSEA failed for",
        outcome_name,
        ":",
        e$message,
        "\n"
      )
      return(NULL)
    }
  )

  if (is.null(fgseaRes) || nrow(fgseaRes) == 0) {
    return(data.frame())
  }

  fgseaRes$outcome <- outcome_name

  # Additional diagnostics in output
  fgseaRes$n_ranked_genes <- length(ranks)
  fgseaRes$n_tied_stats <- n_ties
  fgseaRes$pct_tied_stats <- pct_ties
  fgseaRes$gsea_seed <- seed
  fgseaRes$gsea_sample_size <- GSEA_SAMPLE_SIZE

  cat("\nGSEA results:\n")

  print(
    fgseaRes[, .(
      pathway,
      pval,
      padj,
      ES,
      NES,
      log2err
    )]
  )

  as.data.frame(fgseaRes)
}

#########################
## Run GSEA for each outcome
#########################

outcome_names <- names(ranked_lists)

gsea_results <- if (NCORES > 1) {

  parallel::mclapply(
    seq_along(outcome_names),
    function(i) {
      run_one_gsea(
        outcome_name = outcome_names[i],
        seed = GSEA_SEED + i
      )
    },
    mc.cores = min(
      NCORES,
      length(outcome_names)
    ),
    mc.preschedule = FALSE
  )

} else {

  lapply(
    seq_along(outcome_names),
    function(i) {
      run_one_gsea(
        outcome_name = outcome_names[i],
        seed = GSEA_SEED + i
      )
    }
  )

}

names(gsea_results) <- outcome_names

#########################
## Combine GSEA results
#########################

all_gsea <- bind_rows(gsea_results)

if (nrow(all_gsea) == 0) {
  stop("No GSEA results were produced.")
}

# Convert list columns to character
all_gsea <- all_gsea %>%
  mutate(
    across(
      where(is.list),
      ~ vapply(
        .,
        function(x) paste(x, collapse = ";"),
        character(1)
      )
    )
  )

#########################
## Save results
#########################

write.table(
  all_gsea,
  file = RESULTS_OUTPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nGSEA results saved to:",
  RESULTS_OUTPUT_FILE,
  "\n"
)

cat("\n=== GSEA COMPLETE ===\n")
