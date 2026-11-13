###############################################################################
# NETWORK POST-PROCESSING — LUAD VALIDATION COHORTS
###############################################################################
# This file contains rules for:
#   1. Creating LIONESS sample mappings
#   2. Combining sample-specific LIONESS networks
#   3. Quantile-normalizing the combined networks
#   4. Creating PANDA network edge files

###############################################################################

rule all_network_postprocessing:
    input:
        expand(
            MERGED_NETWORKS_FILE_DATASET_NORMALIZED,
            dataset=NETWORK_POSTPROCESSING_DATASETS
        ),
        expand(
            NETWORK_EDGE_FILE_DATASET,
            dataset=NETWORK_POSTPROCESSING_DATASETS
        ),
        expand(
            INDEGREES_DATASET,
            dataset=NETWORK_POSTPROCESSING_DATASETS,
            norm=["norm"]
        )

###############################################################################
# COPY PRECOMPUTED TCGA_LUAD INDEGREES
###############################################################################

rule copy_TCGA_LUAD_indegrees:
    input:
        indegrees = TCGA_LUAD_INDEGREES_RESOURCE
    output:
        indegrees = TCGA_LUAD_INDEGREES_FILE
    message:
        "Copying precomputed TCGA_LUAD indegrees from resources"
    shell:
        """
        mkdir -p $(dirname {output.indegrees:q})
        cp {input.indegrees:q} {output.indegrees:q}
        """
        
###############################################################################
# 1. CREATE LIONESS SAMPLE MAPPINGS
###############################################################################
# Match sample IDs from the PANDA input sample list to network filenames.

rule create_lioness_mapping:
    input:
        samples_file = DATASET_SAMPLES_PANDA_DATASET,
        networks_dir = NETWORK_DIR_DATASET
    output:
        mapping = LIONESS_SAMPLE_MAPPING_DATASET
    log:
        "logs/create_lioness_sample_mapping_{dataset}.log"
    message:
        "Creating LIONESS sample mapping file for {wildcards.dataset}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/create_lioness_sample_mapping_file.R \
            --network_dir {input.networks_dir} \
            --samples_panda_file {input.samples_file} \
            --output_file {output.mapping} \
            > {log} 2>&1
        """

###############################################################################
# 2. COMBINE LIONESS NETWORKS
###############################################################################
# Save sample-specific networks together in one dataset-level RData file.

rule combine_lioness_networks:
    """Save LIONESS sample-specific networks to dataset-specific RData files."""
    input:
        network_dir = NETWORK_DIR_DATASET,
        lioness_sample_mapping = LIONESS_SAMPLE_MAPPING_DATASET
    output:
        output_file = MERGED_NETWORKS_FILE_DATASET
    log:
        "logs/save_cancer_specific_lioness_networks_{dataset}.log"
    message:
        "Saving LIONESS networks to cancer-specific Rdata files for {wildcards.dataset}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/combine_networks.R \
            --network_dir {input.network_dir} \
            --lioness_sample_mapping {input.lioness_sample_mapping} \
            --output_file {output.output_file} \
            > {log} 2>&1
        """

###############################################################################
# 3. QUANTILE-NORMALIZE COMBINED NETWORKS
###############################################################################

rule normalize_networks:
    """Apply quantile normalization to LIONESS networks."""
    input:
        network_file = MERGED_NETWORKS_FILE_DATASET
    output:
        output_file = MERGED_NETWORKS_FILE_DATASET_NORMALIZED
    log:
        "logs/quantile_normalize_network_{dataset}.log"
    message:
        "Applying quantile normalization to network"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/quantile_normalize_networks.R \
            --network_file "{input.network_file}" \
            --output_file {output.output_file} \
            > {log} 2>&1
        """

###############################################################################
# 4. CREATE PANDA NETWORK EDGE FILE
###############################################################################
# Convert the PANDA network output to an edge file for downstream analysis.

rule create_network_edge_file:
    """Create a network edge file from the PANDA network file."""
    input:
        panda_input = PANDA_NETWORK_FILE_DATASET
    output:
        edge_file = NETWORK_EDGE_FILE_DATASET
    log:
        "logs/create_network_edge_file_{dataset}.log"
    message:
        "Creating network edge file"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/create_edge_file.R \
            --panda_network_file {input.panda_input} \
            --output_edge_file {output.edge_file} \
            > {log} 2>&1
        """
###############################################################################
# 5. NETWORK ANALYSIS
###############################################################################

rule network_analysis:
    """
    Basic network analysis
    
    Output files:
    - Target gene specific network (.RData)
    - Network indegrees (.tsv)
    - PCA plot (.pdf)

    """
    input:
        network = MERGED_NETWORKS_FILE_DATASET_NORM_RAW,
        edges = NETWORK_EDGE_FILE_DATASET
    output:
        target_network = PD1_NETWORK_DATASET,
        indegrees = INDEGREES_DATASET,
        outdegrees = OUTDEGREES_DATASET,
        pca_plot = PCA_INDEGREES_DATASET
    params:
        bin = config["bin"],
        target_gene = "CD274",
        seed = 2025
    wildcard_constraints:
        dataset = "|".join(re.escape(d) for d in NETWORK_POSTPROCESSING_DATASETS),
        norm = "norm"
    log:
        "logs/network_analysis_{dataset}_{norm}.log"
    message:
        "Performing basic network analysis on {wildcards.dataset} dataset with {wildcards.norm} network"
    shell:
        """        
        Rscript {params.bin}/network_analysis.R \
            --network_file {input.network} \
            --edges_file {input.edges} \
            --target_network {output.target_network} \
            --indegrees {output.indegrees} \
            --outdegrees {output.outdegrees} \
            --pca_plot {output.pca_plot} \
            --target_gene {params.target_gene} \
            --seed {params.seed} \
        2> {log}
        """