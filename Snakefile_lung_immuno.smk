## 2. Load required modules: 
##    - module load snakemake/7.23.1-foss-2022a
##    - module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
##  dry run snakemake --cores 1 -np
## 3. Source conda: source ~/.bashrc  # conda version 24.1.2
## 4. Clean conda cache: rm -rf .snakemake/conda
## 5. Execute pipeline: snakemake --use-conda --conda-frontend conda --cores 1



## Libraries
import os 
import sys
import glob
from pathlib import Path
import time

## Config

global CONFIG_PATH
CONFIG_PATH = "config_lung_immuno.yaml"
configfile: CONFIG_PATH


FIGURES_DIR = config["figures_dir"]
RESULTS_DIR = config["results_dir"]
PANDA_DIR = config["panda_dir"]
DATASETS = config["datasets"]

DATA_DIR = config["data_dir"]
NETWORKS_DIR = config["networks_dir"]
MERGED_NETWORKS_DIR = config["merged_networks_dir"]
BUILD_NETWORK_OUTPUTS = config.get("build_network_outputs", False)
PATHWAY_GMT_FILE = os.path.join("resources", "c2.cp.reactome.v7.1.symbols.gmt")
GMT_FILES = {
    "c2": PATHWAY_GMT_FILE
}
COVARIATES = ["with_covariates", "no_covariates"]
GMT_TYPES = ["c2"]
#SUBTYPES = ["all", "adeno"]
SUBTYPES = ["adeno"]
TREATMENTS = ["pd1"]
# TREATMENTS = ["all", "pd1"]


# Expression and sample files for PANDA
IMMUNE_COHORT_TPM_FILE = os.path.join(
    DATA_DIR,
    "immune_cohort/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct"
)
CLINICAL_FILE_IMMUNO = os.path.join(DATA_DIR, "immune_cohort/clinical_filtered.txt")
CLINICAL_FILE_EXTENDED_IMMUNO = os.path.join(DATA_DIR, "immune_cohort/41588_2023_1355_MOESM3_ESM.xlsx")

DATASET_EXPRESSION_PANDA_FILE = os.path.join(PANDA_DIR,  "{dataset}_expression_for_PANDA.tsv")
DATASET_SAMPLES_PANDA_DATASET = os.path.join(PANDA_DIR, "{dataset}_samples_for_PANDA.tsv")
PCA_PLOTS_EXPRESSION_FILE = os.path.join(
    FIGURES_DIR,
    "{dataset}_pca_plots_expression.pdf"
)

# Network files merging and normalization

MERGED_NETWORKS_FILE_DATASET = os.path.join(MERGED_NETWORKS_DIR, "{dataset}_merged_net_raw.RData")
MERGED_NETWORKS_FILE_DATASET_NORMALIZED = os.path.join(MERGED_NETWORKS_DIR, "{dataset}_merged_net_normalized.RData")
PANDA_NETWORK_FILE_DATASET = os.path.join(NETWORKS_DIR, "{dataset}", "panda_net.txt")
NETWORK_EDGE_FILE_DATASET = os.path.join(PANDA_DIR, "{dataset}_network_edges.txt")
LIONESS_SAMPLE_MAPPING_DATASET = os.path.join(PANDA_DIR, "{dataset}_lioness_sample_mapping.txt")
NETWORK_DIR_DATASET = os.path.join(NETWORKS_DIR, "{dataset}")
MERGED_NETWORKS_FILE_DATASET_NORMALIZED_INPUT = os.path.join(MERGED_NETWORKS_DIR, "{dataset}_merged_net_normalized.RData")

# Basic analysis of networks (indegrees, outdegrees, PCA)
INDEGREES_DATASET = os.path.join(RESULTS_DIR, "{dataset}_indegrees.txt")
OUTDEGREES_DATASET = os.path.join(RESULTS_DIR, "{dataset}_outdegrees.txt")
PCA_INDEGREES_DATASET = os.path.join(FIGURES_DIR, "{dataset}_PCA_indegrees.pdf")
PD1_NETWORK_DATASET = os.path.join(RESULTS_DIR, "{dataset}_network_CD274.RData")

# LIMMA analysis of indegrees and FGSEA
LIMMA_RESULTS_FILE = os.path.join(
    RESULTS_DIR,
    "{dataset}_limma_results_subtype_{subtype}_treatment_{treatment}_{covariates}_{gmt}.txt"
)
FGSEA_RESULTS_FILE = os.path.join(
    RESULTS_DIR,
    "{dataset}_fgsea_results_subtype_{subtype}_treatment_{treatment}_{covariates}_{gmt}.txt"
)
# ALL LIMMA models 

LIMMA_RESULTS_FILE_ALL_MODELS = os.path.join(
    RESULTS_DIR,
    "{dataset}_limma_results_all_models_subtype_{subtype}_treatment_{treatment}_{gmt}.txt"
)
FGSEA_RESULTS_FILE_ALL_MODELS = os.path.join(
    RESULTS_DIR,
    "{dataset}_fgsea_results_all_models_subtype_{subtype}_treatment_{treatment}_{gmt}.txt"
)
ALL_MODELS_SUMMARY_PATHWAY = config.get(
    "all_models_summary_pathway",
    "REACTOME_PD_1_SIGNALING"
)
FGSEA_ALL_MODELS_SUMMARY_FILE = os.path.join(
    RESULTS_DIR,
    "{dataset}_fgsea_summary_all_models_subtype_{subtype}_treatment_{treatment}_{gmt}.txt"
)



BASE_TARGETS = [
    expand(DATASET_EXPRESSION_PANDA_FILE, dataset=DATASETS),
    expand(DATASET_SAMPLES_PANDA_DATASET, dataset=DATASETS),
    expand(PCA_PLOTS_EXPRESSION_FILE, dataset=DATASETS),
    expand(INDEGREES_DATASET, dataset=DATASETS),
    expand(OUTDEGREES_DATASET, dataset=DATASETS),
    expand(PCA_INDEGREES_DATASET, dataset=DATASETS),
    expand(
            LIMMA_RESULTS_FILE,
            dataset=DATASETS,
            subtype=SUBTYPES,
            treatment=TREATMENTS,
            covariates=COVARIATES,
            gmt=GMT_TYPES
        ),
        expand(
            FGSEA_RESULTS_FILE,
            dataset=DATASETS,
            subtype=SUBTYPES,
            treatment=TREATMENTS,
            covariates=COVARIATES,
            gmt=GMT_TYPES
        ),
    expand(
            LIMMA_RESULTS_FILE_ALL_MODELS,
            dataset=DATASETS,
            subtype=["adeno"],
            treatment=["pd1"],
            gmt=GMT_TYPES
        ),
        expand(
            FGSEA_RESULTS_FILE_ALL_MODELS,
            dataset=DATASETS,
            subtype=["adeno"],
            treatment=["pd1"],
            gmt=GMT_TYPES
        ),
        expand(
            FGSEA_ALL_MODELS_SUMMARY_FILE,
            dataset=DATASETS,
            subtype=["adeno"],
            treatment=["pd1"],
            gmt=GMT_TYPES
        )
]

NETWORK_BUILD_TARGETS = [
    expand(LIONESS_SAMPLE_MAPPING_DATASET, dataset=DATASETS),
    expand(MERGED_NETWORKS_FILE_DATASET, dataset=DATASETS),
    expand(MERGED_NETWORKS_FILE_DATASET_NORMALIZED, dataset=DATASETS),
    expand(NETWORK_EDGE_FILE_DATASET, dataset=DATASETS),
]

ALL_TARGETS = BASE_TARGETS + (NETWORK_BUILD_TARGETS if BUILD_NETWORK_OUTPUTS else [])

rule all:
    input: ALL_TARGETS


rule immune_cohort_data_preprocessing:
    input:
        tpm_file = IMMUNE_COHORT_TPM_FILE,
        clinical_file = CLINICAL_FILE_IMMUNO
    output:
        expression_output_file = DATASET_EXPRESSION_PANDA_FILE,
        samples_output_file = DATASET_SAMPLES_PANDA_DATASET,
        pca_plot_file = PCA_PLOTS_EXPRESSION_FILE
    log:
        "logs/immune_cohort_data_preprocessing_{dataset}.log"
    message:
        "Preprocessing immune cohort TPM expression for {wildcards.dataset}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/immune_cohort_data_preprocessing.R \
            --tpm_file {input.tpm_file} \
            --clinical_file {input.clinical_file} \
            --expression_output_file {output.expression_output_file} \
            --samples_output_file {output.samples_output_file} \
            --pca_plot_file {output.pca_plot_file} \
            &> {log}
        """

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
        bin = config["bin"],
    shell:
        """
        Rscript {params.bin}/create_lioness_sample_mapping_file.R \
            --network_dir {input.networks_dir} \
            --samples_panda_file {input.samples_file} \
            --output_file {output.mapping} \
            > {log} 2>&1
        """
rule combine_lioness_networks:
    """
    Save LIONESS sample-specific networks to cancer-specific RData files.
    """
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

# ## Apply quantile normalization on the networks
rule normalize_networks:
    """
    Apply quantile normalization to LIONESS networks.
    """
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

# ## Create a network edge file
rule create_network_edge_file:
    """
    Create a network edge file from the PANDA network file.
    """
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
rule network_analysis:
    """
    Basic network analysis
    
    Output files:
    - Target gene specific network (.RData)
    - Network indegrees (.tsv)
    - PCA plot (.pdf)

    """
    input:
        network = MERGED_NETWORKS_FILE_DATASET_NORMALIZED_INPUT,
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
    log:
        "logs/network_analysis_{dataset}.log"
    message:
        "Performing basic network analysis on {wildcards.dataset} dataset with normalized network"
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

rule lung_immuno_limma:
    input:
        indegree_file = INDEGREES_DATASET,
        clinical_file = CLINICAL_FILE_IMMUNO,
        pathway_gmt = lambda wc: GMT_FILES[wc.gmt]
    output:
        limma_results_file = LIMMA_RESULTS_FILE,
        fgsea_results_file= FGSEA_RESULTS_FILE
    log:
        "logs/limma_analysis_{dataset}_{subtype}_{treatment}_{covariates}_{gmt}.log"
    message:
        "Performing limma analysis for subtype={wildcards.subtype}, treatment={wildcards.treatment}, covariates={wildcards.covariates}, gmt={wildcards.gmt}"
    params:
        bin = config["bin"],
        subtype=lambda wc: wc.subtype,
        treatment=lambda wc: wc.treatment,
        covariates=lambda wc: "TRUE" if wc.covariates == "with_covariates" else "FALSE",

    wildcard_constraints:
        subtype="all|adeno",
        treatment="all|pd1",
        covariates="with_covariates|no_covariates",
        gmt="c2"

    shell:
        """
         Rscript {params.bin}/lung_immuno_limma_indegrees.R \
            --indegree_file {input.indegree_file} \
            --clinical_file {input.clinical_file} \
            --pathway_gmt_file {input.pathway_gmt} \
            --limma_results_file {output.limma_results_file} \
            --fgsea_results_file {output.fgsea_results_file} \
            --subtype_type {params.subtype} \
            --treatment_type {params.treatment} \
            --covariates {params.covariates} \
        &> {log}
        """

rule lung_immuno_indegree_all_models:
    input:
        indegree_file = INDEGREES_DATASET,
        clinical_file = CLINICAL_FILE_IMMUNO,
        clinical_file_extended = CLINICAL_FILE_EXTENDED_IMMUNO,
        expression_file = DATASET_EXPRESSION_PANDA_FILE,
        samples_file = DATASET_SAMPLES_PANDA_DATASET,
        pathway_gmt = lambda wc: GMT_FILES[wc.gmt]
    output:
        limma_results_file = LIMMA_RESULTS_FILE_ALL_MODELS,
        fgsea_results_file = FGSEA_RESULTS_FILE_ALL_MODELS
    log:
        "logs/limma_all_models_{dataset}_{subtype}_{treatment}_{gmt}.log"
    message:
        "Running fixed all-model limma set (M0 + covariate/sensitivity models; no covariates flag) for subtype={wildcards.subtype}, treatment={wildcards.treatment}, gmt={wildcards.gmt}"
    params:
        bin = config["bin"],
        subtype = lambda wc: wc.subtype,
        treatment = lambda wc: wc.treatment
    wildcard_constraints:
        subtype = "adeno",
        treatment = "pd1",
        gmt = "c2"
    shell:
        """
        Rscript {params.bin}/lung_immuno_indegree_all_models.R \\
            --indegree_file {input.indegree_file} \\
            --clinical_file {input.clinical_file} \\
            --clinical_file_extended {input.clinical_file_extended} \\
            --expression_file {input.expression_file} \\
            --samples_file {input.samples_file} \\
            --pathway_gmt_file {input.pathway_gmt} \\
            --limma_results_file {output.limma_results_file} \\
            --fgsea_results_file {output.fgsea_results_file} \\
            --subtype_type {params.subtype} \\
            --treatment_type {params.treatment} \\
        &> {log}
        """
rule summarise_all_models_fgsea:
    input:
        fgsea_results_file = FGSEA_RESULTS_FILE_ALL_MODELS
    output:
        summary_table = FGSEA_ALL_MODELS_SUMMARY_FILE
    log:
        "logs/summarise_all_models_fgsea_{dataset}_{subtype}_{treatment}_{gmt}.log"
    message:
        "Summarising all-model FGSEA for pathway={params.pathway}, subtype={wildcards.subtype}, treatment={wildcards.treatment}, gmt={wildcards.gmt}"
    params:
        bin = config["bin"],
        pathway = ALL_MODELS_SUMMARY_PATHWAY
    wildcard_constraints:
        subtype = "adeno",
        treatment = "pd1",
        gmt = "c2"
    shell:
        """
        Rscript {params.bin}/summarise_all_models_fgsea.R \\
            --fgsea_results_file {input.fgsea_results_file} \\
            --pathway {params.pathway} \\
            --output_table {output.summary_table} \\
            &> {log}
        """

