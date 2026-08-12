#####################
## Load R packages ##
#####################
required_libraries <- c("data.table", "dplyr", "ggplot2", "optparse")
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
        help = "Path to TCGA network .RData file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-e", "--edges_file"),
        type = "character",
        default = NULL,
        help = "Path to edges file (.txt).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-o", "--target_network"),
        type = "character",
        default = NULL,
        help = "Output path for target gene network (.RData).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-i", "--indegrees"),
        type = "character",
        default = NULL,
        help = "Output path for indegrees file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-l", "--outdegrees"),
        type = "character",
        default = NULL,
        help = "Output path for outdegree file.",
        metavar = "character"
    ),
    optparse::make_option(
        c("-p", "--pca_plot"),
        type = "character",
        default = NULL,
        help = "Output path for PCA plot (.pdf).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-t", "--target_gene"),
        type = "character",
        default = "CD274",
        help = "Target gene for network analysis (default: CD274).",
        metavar = "character"
    ),
    optparse::make_option(
        c("-k", "--seed"),
        type = "integer",
        default = 2025,
        help = "Random seed (default: 2025).",
        metavar = "integer"
    )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

### Initialize variables ###
NETWORK_FILE <- opt$network_file
EDGES_FILE <- opt$edges_file
TARGET_NETWORK_OUTPUT <- opt$target_network
INDEGREES_OUTPUT <- opt$indegrees
OUTDEGREES_OUTPUT <- opt$outdegrees
PCA_PLOT_OUTPUT <- opt$pca_plot
TARGET_GENE <- opt$target_gene
SEED <- opt$seed

# Validate required arguments
required_args <- list(
    "Network file" = NETWORK_FILE,
    "Edges file" = EDGES_FILE,
    "Target network output" = TARGET_NETWORK_OUTPUT,
    "Indegrees output" = INDEGREES_OUTPUT,
    "Outdegrees output" = OUTDEGREES_OUTPUT,
    "PCA plot output" = PCA_PLOT_OUTPUT
)
missing_args <- names(required_args)[sapply(required_args, is.null)]
if (length(missing_args) > 0) {
    stop("Required arguments missing: ", 
         paste(missing_args, collapse = ", "))
}

cat("Starting network analysis...\n")
cat("Network file:", NETWORK_FILE, "\n")
cat("Edges file:", EDGES_FILE, "\n")
cat("Target gene:", TARGET_GENE, "\n")
cat("Target network output:", TARGET_NETWORK_OUTPUT, "\n")
cat("Indegrees output:", INDEGREES_OUTPUT, "\n")
cat("Outdegrees output:", OUTDEGREES_OUTPUT, "\n")
cat("PCA plot output:", PCA_PLOT_OUTPUT, "\n")


# Create output directories
for (output_file in c(TARGET_NETWORK_OUTPUT, INDEGREES_OUTPUT, 
                     PCA_PLOT_OUTPUT)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
}

# Set seed
set.seed(SEED)

############################
## Load network data      ##
############################
cat("Loading network file...\n")
net_env <- new.env()
load(NETWORK_FILE, envir = net_env)

# Try different variable names
net_var_names <- c("net_norm", "network", "net")
net <- NULL
for (var_name in net_var_names) {
    if (var_name %in% names(net_env)) {
        net <- net_env[[var_name]]
        break
    }
}

if (is.null(net)) {
    cat("Available variables in network file:", 
        paste(names(net_env), collapse = ", "), "\n")
    stop("Could not find network data in file: ", NETWORK_FILE)
}

cat("Network dimensions:", dim(net), "\n")

############################
## Load edges file        ##
############################

cat("Loading edges file...\n")
edges <- fread(EDGES_FILE)
cat("Edges dimensions:", dim(edges), "\n")

############################
## Calculate indegrees    ##
############################

cat("Calculating indegrees...\n")
net_edges <- cbind(edges, net)

# Use summarise with modern syntax

indegree_norm <-  net_edges %>%
                select(-c("reg")) %>%
                group_by(tar) %>%
                summarise_all(funs(sum))

outdegree_norm <-  net_edges %>%
                select(-c("tar")) %>%
                group_by(reg) %>%
                summarise_all(funs(sum))


# Save indegrees
head(indegree_norm)
write.table(indegree_norm, INDEGREES_OUTPUT,
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("Saved indegrees to:", INDEGREES_OUTPUT, "\n")

head(outdegree_norm)
write.table(outdegree_norm, OUTDEGREES_OUTPUT,
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
cat("Saved outdegrees to:", OUTDEGREES_OUTPUT, "\n")



############################
## PCA analysis           ##
############################

cat("Performing PCA analysis...\n")
ind_norm_scaled <- t(apply(indegree_norm[, -1], 1, function(x) scale(x)))
pca_result_ind_norm <- prcomp(t(ind_norm_scaled), retx = TRUE, center = TRUE)
PCA_ind_norm <- as.data.frame(pca_result_ind_norm$x[, 1:20])
PCA_ind_norm$sample_id <- colnames(indegree_norm)[-1]
# Create PCA plot
pdf(PCA_PLOT_OUTPUT, width = 10, height = 10)
p1 <- ggplot(PCA_ind_norm, aes(PC1, PC2)) + 
    geom_point(size = 3) + 
    ggtitle("PCA network indegrees") + 
    theme_minimal() +  
    theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
print(p1)
dev.off()
cat("Saved PCA plot to:", PCA_PLOT_OUTPUT, "\n")

############################
## Extract target gene network ##
############################

cat("Extracting", TARGET_GENE, "subnetwork...\n")
idx <- which(edges$tar == TARGET_GENE)
cat("Found", length(idx), "edges for", TARGET_GENE, "\n")

if (length(idx) == 0) {
    stop("No edges found for target gene: ", TARGET_GENE)
}

links <- edges[idx, ]
pd1_net <- net[idx, ]
pd1_net <- as.matrix(pd1_net)
pd1_net <- cbind(links, pd1_net)

# Save target gene network
save(pd1_net, file = TARGET_NETWORK_OUTPUT)
cat("Saved", TARGET_GENE, "network to:", TARGET_NETWORK_OUTPUT, "\n")




