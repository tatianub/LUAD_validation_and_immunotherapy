# ---------------------------------------------------------------------------
# Load R packages
# ---------------------------------------------------------------------------
required_libraries <- c(
  "data.table", "ggplot2", "optparse", "readxl", "grid", "MASS"
)
for (lib in required_libraries) {
  suppressPackageStartupMessages(
    library(lib, character.only = TRUE, quietly = TRUE)
  )
}
options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Read arguments
# ---------------------------------------------------------------------------
option_list <- list(
  optparse::make_option("--clinical_file", type = "character"),
  optparse::make_option("--network_file", type = "character"),
  optparse::make_option("--edges_file", type = "character"),
  optparse::make_option("--expression_file", type = "character"),
  optparse::make_option("--samples_file", type = "character"),
  optparse::make_option("--clinical_file_cnv", type = "character"),
  optparse::make_option(
    "--cnv_sheet",
    type = "character",
    default = "Table_S8_Gistic_Gene_Events"
  ),
  optparse::make_option("--tf_fgsea_file", type = "character"),
  optparse::make_option(
    "--pathway",
    type = "character",
    default = "REACTOME_PD_1_SIGNALING"
  ),
  optparse::make_option(
    "--comparison",
    type = "character",
    default = "Response_vs_resistance"
  ),
  optparse::make_option(
    "--histo_subtype",
    type = "character",
    default = NULL
  ),
  optparse::make_option(
    "--treatment_type",
    type = "character",
    default = NULL
  ),
  optparse::make_option("--linear_models_file", type = "character"),
  optparse::make_option("--logistic_models_file", type = "character"),
  optparse::make_option("--logistic_lrt_file", type = "character"),
  optparse::make_option("--four_panel_summary_file", type = "character"),
  optparse::make_option("--interaction_results_file", type = "character")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# ---------------------------------------------------------------------------
# Assign input options
# ---------------------------------------------------------------------------
CLINICAL_FILE            <- opt$clinical_file
NETWORK_FILE             <- opt$network_file
EDGES_FILE               <- opt$edges_file
EXPRESSION_FILE          <- opt$expression_file
SAMPLES_FILE             <- opt$samples_file
CLINICAL_FILE_CNV        <- opt$clinical_file_cnv
CNV_SHEET                <- opt$cnv_sheet
TF_FGSEA_FILE            <- opt$tf_fgsea_file
PATHWAY                  <- opt$pathway
COMPARISON               <- opt$comparison
HISTO_SUBTYPE            <- opt$histo_subtype
TREATMENT_TYPE           <- opt$treatment_type
LINEAR_MODELS_FILE       <- opt$linear_models_file
LOGISTIC_MODELS_FILE     <- opt$logistic_models_file
LOGISTIC_LRT_FILE        <- opt$logistic_lrt_file
FOUR_PANEL_SUMMARY_FILE  <- opt$four_panel_summary_file
INTERACTION_RESULTS_FILE <- opt$interaction_results_file

colors <- c(resistance = "#2166AC", response = "#B2182B")

source("workflow/bin/lung_immuno_limma_fn.R")
source("workflow/bin/utility_correlation_fn.R")

# ---------------------------------------------------------------------------
# Load and filter clinical data
# ---------------------------------------------------------------------------
clinical <- fread(CLINICAL_FILE)
stopifnot(all(c("sample_id", "response") %in% names(clinical)))
clinical <- subset_clinical_data(
  clinical_data = clinical[!is.na(response)],
  subtype = HISTO_SUBTYPE,
  treatment_type = TREATMENT_TYPE
)
clinical[, sample_id := gsub("[^A-Za-z0-9]", ".", sample_id)]
clinical <- clinical[response %in% c("resistance", "response")]
if (anyDuplicated(clinical$sample_id)) stop("Duplicate clinical sample IDs")

# ---------------------------------------------------------------------------
# Extract MYC-specific PD-1 leading edge
# ---------------------------------------------------------------------------
fg <- fread(TF_FGSEA_FILE)
stopifnot(all(c("pathway", "leadingEdge", "cmp", "TF") %in% names(fg)))
sel <- fg[TF == "MYC" & pathway == PATHWAY & cmp == COMPARISON]
if (nrow(sel) != 1L) {
  available <- unique(fg[TF == "MYC" & grepl("PD_1", pathway), pathway])
  stop(
    "Expected one MYC/pathway/comparison row, found ", nrow(sel),
    ". Available MYC PD-1 pathways: ", paste(available, collapse = ", ")
  )
}
le <- unique(trimws(strsplit(sel$leadingEdge[1], ";", fixed = TRUE)[[1]]))
le <- le[nzchar(le)]
if (!length(le)) stop("Selected leading edge is empty")

# message("MYC-specific leading edge: ", paste(le, collapse = ", "))

# ---------------------------------------------------------------------------
# Load LIONESS networks
# ---------------------------------------------------------------------------
e <- new.env()
load(NETWORK_FILE, envir = e)
if (exists("net", envir = e, inherits = FALSE)) {
  net <- e$net
} else if (exists("net_norm", envir = e, inherits = FALSE)) {
  net <- e$net_norm
} else {
  stop("No net or net_norm in RData")
}
net <- as.matrix(net)
if (is.null(colnames(net))) stop("Network has no sample column names")
colnames(net) <- gsub("[^A-Za-z0-9]", ".", colnames(net))
if (anyDuplicated(colnames(net))) stop("Duplicate normalized network IDs")

edges <- fread(EDGES_FILE)
stopifnot(all(c("reg", "tar") %in% names(edges)), nrow(edges) == nrow(net))

# ---------------------------------------------------------------------------
# Load expression matrix
# ---------------------------------------------------------------------------
expr_dt <- fread(EXPRESSION_FILE)
expr_genes <- as.character(expr_dt[[1]])
expr <- data.matrix(expr_dt[, -1, with = FALSE])
samples <- fread(SAMPLES_FILE)
rownames(expr) <- expr_genes
colnames(expr) <- gsub("[^A-Za-z0-9]", ".", as.character(samples$V2))
if (anyDuplicated(colnames(expr))) stop("Duplicate normalized expression IDs")

common <- colnames(net)[
  colnames(net) %in% clinical$sample_id &
    colnames(net) %in% colnames(expr)
]
if (!length(common)) stop("No shared network, expression and clinical samples")
if (sum(rownames(expr) == "MYC") != 1L) stop("Expected one MYC expression row")

# ---------------------------------------------------------------------------
# Calculate MYC leading-edge targeting score
# ---------------------------------------------------------------------------
idx <- which(edges$reg == "MYC" & edges$tar %in% le)
if (!length(idx)) stop("No MYC leading-edge edges found")
if (anyDuplicated(edges$tar[idx])) {
  stop("Duplicate MYC-target edges; resolve duplicates before summing")
}
missing <- setdiff(le, edges$tar[idx])
if (length(missing)) {
  warning(
    "Leading-edge genes missing from MYC network: ",
    paste(missing, collapse = ", ")
  )
}
edge_mat <- net[idx, common, drop = FALSE]
if (anyNA(edge_mat)) {
  stop("Missing MYC edge weights; cannot sum incomplete networks")
}
score <- colSums(edge_mat)

dt <- data.table(
  sample_id = common,
  myc_le_targeting = as.numeric(score),
  myc_expression = as.numeric(expr["MYC", common])
)
dt <- merge(dt, clinical[, .(sample_id, response)], by = "sample_id",
            sort = FALSE)
dt[, response := factor(response, levels = c("resistance", "response"))]

# ---------------------------------------------------------------------------
# Load MYC CNV
# ---------------------------------------------------------------------------
cnv <- as.data.table(read_excel(CLINICAL_FILE_CNV, sheet = CNV_SHEET))
gene_col <- intersect(
  c("Hugo_Symbol", "Gene", "Gene_Symbol", "Symbol"),
  names(cnv)
)
if (!length(gene_col)) {
  gene_col <- names(cnv)[1]
} else {
  gene_col <- gene_col[1]
}
myc_row <- cnv[toupper(trimws(as.character(get(gene_col)))) == "MYC"]
if (nrow(myc_row) != 1L) stop("Expected one MYC CNV row, found ", nrow(myc_row))

# Keep only CNV columns corresponding to patients in the analysis.
# The Excel sheet also contains metadata (Gene ID, Cytoband, etc.).

dt[, id_norm := normalize_id(sample_id)]
cnv_cols <- setdiff(names(myc_row), gene_col)
cnv_cols <- cnv_cols[normalize_id(cnv_cols) %in% dt$id_norm]
if (!length(cnv_cols)) {
  stop("No CNV sample columns match the clinical/network sample IDs")
}
cnv_long <- melt(
  myc_row[, ..cnv_cols],
  measure.vars = cnv_cols,
  variable.name = "cnv_sample",
  value.name = "myc_cnv",
  variable.factor = FALSE
)
cnv_long[, id_norm := normalize_id(cnv_sample)]
cnv_long[, myc_cnv := suppressWarnings(as.numeric(myc_cnv))]
if (anyDuplicated(cnv_long$id_norm)) {
  stop("Duplicate normalized CNV sample IDs")
}
message("Matched MYC CNV sample columns: ", length(cnv_cols))

if (anyDuplicated(dt$id_norm)) stop("Duplicate normalized sample IDs")
dt <- merge(dt, cnv_long[, .(id_norm, myc_cnv)], by = "id_norm",
            all.x = TRUE, sort = FALSE)
setcolorder(
  dt,
  c("sample_id", "response", "myc_le_targeting", "myc_expression",
    "myc_cnv", "id_norm")
)

full <- dt[
  is.finite(myc_le_targeting) & is.finite(myc_expression) & !is.na(response)
]
sub <- full[is.finite(myc_cnv)]
if (nrow(full) < 8 || length(unique(full$response)) < 2) {
  stop("Insufficient full-cohort data")
}
if (nrow(sub) < 8 || length(unique(sub$response)) < 2) {
  stop("Insufficient CNV-subset data")
}

# ---------------------------------------------------------------------------
cors <- rbind(cor_fun(full, "myc_expression"), cor_fun(sub, "myc_cnv"))

# ---------------------------------------------------------------------------
# Linear models: targeting outcome
# ---------------------------------------------------------------------------
model_specs <- list(
  list(
    name = "Unadjusted (all)",
    data = full,
    formula = myc_le_targeting ~ response
  ),
  list(
    name = "Expression-adjusted (all)",
    data = full,
    formula = myc_le_targeting ~ myc_expression + response
  ),
  list(
    name = "Unadjusted (CNV subset)",
    data = sub,
    formula = myc_le_targeting ~ response
  ),
  list(
    name = "Expression-adjusted (CNV subset)",
    data = sub,
    formula = myc_le_targeting ~ myc_expression + response
  ),
  list(
    name = "Expression + CNV adjusted",
    data = sub,
    formula = myc_le_targeting ~ myc_expression + myc_cnv + response
  )
)

model_results <- rbindlist(lapply(model_specs, function(s) {
  m <- lm(s$formula, data = s$data)
  b <- summary(m)$coefficients
  ci <- confint(m)
  data.table(
    model = s$name,
    n = nobs(m),
    term = rownames(b),
    estimate = b[, 1],
    SE = b[, 2],
    p_value = b[, 4],
    CI_low = ci[, 1],
    CI_high = ci[, 2],
    r_squared = summary(m)$r.squared
  )
}))
fwrite(model_results, LINEAR_MODELS_FILE, sep = "\t")

# ---------------------------------------------------------------------------
# Interaction analysis: does the slope differ by response group?
# ---------------------------------------------------------------------------

interaction_specs <- list(
  list(
    name = "MYC expression x response (all)",
    data = full,
    formula = myc_le_targeting ~ myc_expression * response,
    term = "myc_expression:responseresponse"
  ),
  list(
    name = "MYC CNV x response (CNV subset)",
    data = sub,
    formula = myc_le_targeting ~ myc_cnv * response,
    term = "myc_cnv:responseresponse"
  )
)

interaction_results <- rbindlist(lapply(interaction_specs, function(s) {
  m <- lm(s$formula, data = s$data)
  b <- summary(m)$coefficients
  if (!s$term %in% rownames(b)) {
    return(data.table(
      model = s$name, n = nobs(m), term = s$term,
      estimate = NA_real_, SE = NA_real_, p_value = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_
    ))
  }
  ci <- confint(m)
  data.table(
    model = s$name,
    n = nobs(m),
    term = s$term,
    estimate = b[s$term, 1],
    SE = b[s$term, 2],
    p_value = b[s$term, 4],
    CI_low = ci[s$term, 1],
    CI_high = ci[s$term, 2]
  )
}))
fwrite(interaction_results, INTERACTION_RESULTS_FILE, sep = "\t")

interaction_p <- interaction_results[
  term == "myc_expression:responseresponse", p_value
]


# ---------------------------------------------------------------------------
# Logistic models: response outcome
# ---------------------------------------------------------------------------
# Nested models use the same patients within each comparison.
logit_specs <- list(
  list(
    name = "MYC LE alone (all)",
    data = full,
    formula = response ~ myc_le_targeting
  ),
  list(
    name = "Expression alone (all)",
    data = full,
    formula = response ~ myc_expression
  ),
  list(
    name = "Expression + LE (all)",
    data = full,
    formula = response ~ myc_expression + myc_le_targeting
  ),
  list(
    name = "Expression + CNV (subset)",
    data = sub,
    formula = response ~ myc_expression + myc_cnv
  ),
  list(
    name = "Expression + CNV + LE (subset)",
    data = sub,
    formula = response ~ myc_expression + myc_cnv + myc_le_targeting
  )
)

logit_results <- rbindlist(lapply(logit_specs, function(s) {
  m <- glm(s$formula, data = s$data, family = binomial())
  b <- summary(m)$coefficients
  data.table(
    model = s$name,
    n = nobs(m),
    term = rownames(b),
    log_odds = b[, 1],
    OR = exp(b[, 1]),
    p_value = b[, 4]
  )
}))
fwrite(logit_results, LOGISTIC_MODELS_FILE, sep = "\t")

logit_lrt <- rbindlist(list(
  data.table(
    comparison = "Add LE to MYC expression (all)",
    n = nrow(full),
    p = anova(
      glm(response ~ myc_expression, data = full, family = binomial()),
      glm(response ~ myc_expression + myc_le_targeting, data = full,
          family = binomial()),
      test = "Chisq"
    )$`Pr(>Chi)`[2]
  ),
  data.table(
    comparison = "Add LE to MYC expression + CNV (subset)",
    n = nrow(sub),
    p = anova(
      glm(response ~ myc_expression + myc_cnv, data = sub,
          family = binomial()),
      glm(response ~ myc_expression + myc_cnv + myc_le_targeting,
          data = sub, family = binomial()),
      test = "Chisq"
    )$`Pr(>Chi)`[2]
  )
))
fwrite(logit_lrt, LOGISTIC_LRT_FILE, sep = "\t")

# ---------------------------------------------------------------------------
# Four-panel figure
# ---------------------------------------------------------------------------
plot_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

y_label <- "MYC leading-edge outdegree"

# Panel A: targeting by response
wilcox_p <- wilcox.test(myc_le_targeting ~ response, data = full)$p.value

p_a <- ggplot(full, aes(response, myc_le_targeting,
                        fill = response, color = response)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .2) +
  geom_jitter(width = .12, alpha = .65, size = 1.8) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(
    title = "A. MYC leading-edge outdegree by response",
    subtitle = paste0("n=", nrow(full), "; Wilcoxon p=", fmt(wilcox_p)),
    x = NULL,
    y = y_label
  ) +
  plot_theme +
  theme(legend.position = "none")

# Panel B: MYC expression vs targeting
p_b <- ggplot(full, aes(myc_expression, myc_le_targeting,
                        color = response, fill = response)) +
  geom_point(alpha = .65, size = 1.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .7, alpha = .13) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    title = "B. MYC expression vs leading-edge outdegree",
    subtitle = get_cor_label("myc_expression"),
    x = "MYC expression",
    y = y_label
  ) +
  plot_theme

# Panel C: MYC CNV vs targeting
p_c <- ggplot(sub, aes(myc_cnv, myc_le_targeting,
                       color = response, fill = response)) +
  geom_point(alpha = .7, size = 2,
             position = position_jitter(width = .06, height = 0)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = .7, alpha = .13) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    title = "C. MYC CNV vs leading-edge outdegree",
    subtitle = get_cor_label("myc_cnv"),
    x = "MYC CNV (GISTIC values)",
    y = y_label
  ) +
  plot_theme

# Panel D: response effect across linear models
adjusted <- model_results[term == "responseresponse"]
adjusted[, model := factor(
  model,
  levels = rev(vapply(model_specs, `[[`, character(1), "name"))
)]

# NOTE: the offsets below (+5, +13, -2) are in data units, so they only
# look right when the response effects are on a similar scale.

label_x <- max(adjusted$CI_high, na.rm = TRUE) + 5
adjusted[, label := paste0("p=", fmt(p_value), "; n=", n)]

p_d <- ggplot(adjusted, aes(estimate, model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                orientation = "y", width = .18) +
  geom_point(size = 2.5) +
  geom_text(aes(x = label_x, label = label), hjust = 0, size = 2.8) +
  scale_x_continuous(
    limits = c(min(0, adjusted$CI_low) - 2, label_x + 13),
    expand = expansion(mult = c(.02, .02))
  ) +
  labs(
    title = "D. Response effect on MYC leading-edge outdegree",
    subtitle = "Responder minus non-responder; 95% confidence intervals",
    x = "Difference in MYC leading-edge outdegree",
    y = NULL
  ) +
  plot_theme

# ---------------------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------------------

pdf(FOUR_PANEL_SUMMARY_FILE, width = 12, height = 12, onefile = TRUE)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
plots <- list(p_a, p_b, p_c, p_d)
for (i in seq_along(plots)) {
  row <- (i - 1L) %/% 2L + 1L
  col <- (i - 1L) %% 2L + 1L
  print(
    plots[[i]],
    newpage = FALSE,
    vp = grid::viewport(layout.pos.row = row, layout.pos.col = col)
  )
}
grid::popViewport()
dev.off()

cat("Saved: ", FOUR_PANEL_SUMMARY_FILE, "\n", sep = "")
cat("Leading-edge genes used: ", length(idx), "/", length(le), "\n", sep = "")
cat("Full cohort: ", nrow(full), "; CNV subset: ", nrow(sub), "\n", sep = "")
print(adjusted[, .(model, n, estimate, CI_low, CI_high, p_value)])
print(logit_lrt)
