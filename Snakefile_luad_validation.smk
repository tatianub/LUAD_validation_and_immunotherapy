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
CONFIG_PATH = "config_lung_validation.yaml"
configfile: CONFIG_PATH


FIGURES_DIR = config["figures_dir"]
RESULTS_DIR = config["results_dir"]
PANDA_DIR = config["panda_dir"]
DATASETS = config["datasets"]
MERGED_NETWORKS_DIR = config["merged_networks_dir"]

NETWORKS_DIR = config["networks_dir"]
DATA_DIR = config["data_dir"]
NORM_TYPES =["raw", "norm"]


DATASET_ORIGINAL_EXPRESSION_FILE = os.path.join(DATA_DIR, "{dataset}", "{dataset}_expression.RData")
DATASET_ORIGINAL_PROBE_ANNOTATION_FILE = os.path.join(DATA_DIR, "{dataset}", "{dataset}_probe_annotation.RData")
DATASET_ORIGINAL_CLINICAL_FILE = os.path.join(DATA_DIR, "{dataset}", "{dataset}_clinical.RData")

# Dolgalev preprocessing
DOLGALEV_DIR = os.path.join(DATA_DIR, "LUAD_Dolgalev")
# PROCESSING DOLGALEV DATA
DOLGALEV_ORIGINAL_EXPRESSION_FILE = os.path.join(DOLGALEV_DIR, "rna-counts-normalized-n246.csv")
DOLGALEV_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "LUAD_Dolgalev_PCA_plot_with_outliers.pdf")
DOLGALEV_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "LUAD_Dolgalev_expression_for_PANDA.tsv")
DOLGALEV_SAMPLES_PANDA = os.path.join(PANDA_DIR, "LUAD_Dolgalev_samples_for_PANDA.tsv")

DATASET_EXPRESSION_PANDA_FILE = os.path.join(PANDA_DIR,  "{dataset}_expression_for_PANDA.tsv")
DATASET_SAMPLES_PANDA_DATASET = os.path.join(PANDA_DIR, "{dataset}_samples_for_PANDA.tsv")

# PROCESSING CHEN DATA

CHEN_DIR = os.path.join(DATA_DIR, "LUAD_Chen")
CHEN_ORIGINAL_EXPRESSION_FILE = os.path.join(CHEN_DIR, "Expr_tumor.txt")
CHEN_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "LUAD_Chen_PCA_plot_with_outliers.pdf")
CHEN_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "LUAD_Chen_expression_for_PANDA.tsv")
CHEN_SAMPLES_PANDA = os.path.join(PANDA_DIR, "LUAD_Chen_samples_for_PANDA.tsv")

rule all:
    input:
        expand(DATASET_ORIGINAL_EXPRESSION_FILE, dataset=DATASETS),
        expand(DATASET_ORIGINAL_PROBE_ANNOTATION_FILE, dataset=DATASETS),
        expand(DATASET_ORIGINAL_CLINICAL_FILE, dataset=DATASETS),
        DOLGALEV_OUTPUT_PCA_FILE,
        DOLGALEV_EXPRESSION_PANDA,
        DOLGALEV_SAMPLES_PANDA,
        CHEN_OUTPUT_PCA_FILE,
        CHEN_EXPRESSION_PANDA,
        CHEN_SAMPLES_PANDA

rule download_from_GEO:
    output:
        expression_file = DATASET_ORIGINAL_EXPRESSION_FILE,
        probe_annotation_file = DATASET_ORIGINAL_PROBE_ANNOTATION_FILE,
        clinical_file = DATASET_ORIGINAL_CLINICAL_FILE
        
    log:
        "logs/download_GEO_{dataset}.log"

    params:
        bin = config["bin"],
        gse_id = lambda wildcards: wildcards.dataset
    shell:
        """
        Rscript {params.bin}/download_GEO_data.R \
            --GSE_id {params.gse_id} \
            --expression_file {output.expression_file} \
            --probe_annotation_file {output.probe_annotation_file} \
            --clinical_file {output.clinical_file} \
            > {log} 2>&1
        """

rule process_dolgalev_data:
    input:
        expression_file = DOLGALEV_ORIGINAL_EXPRESSION_FILE
    output:
        pca_file = DOLGALEV_OUTPUT_PCA_FILE,
        expression_panda = DOLGALEV_EXPRESSION_PANDA,
        samples_panda = DOLGALEV_SAMPLES_PANDA
    log:
        "logs/process_expression_data_dolgalev.log"
    message:
        "Processing expression data"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/LUAD_Dolgalev_expression_preprocessing.R \
            --expression_file {input.expression_file} \
            --exp_clean {output.expression_panda} \
            --pca_plot {output.pca_file} \
            --samples_file {output.samples_panda} \
            2> {log}
        """

rule chen_luad_expression_analysis:
    """
    Process LUAD Chen expression data to remove outliers and perform PCA analysis
    """
    input:
        expression_file = CHEN_ORIGINAL_EXPRESSION_FILE,
    output:
        pca_plot = CHEN_OUTPUT_PCA_FILE,
        expression_matrix = CHEN_EXPRESSION_PANDA,
        samples_file = CHEN_SAMPLES_PANDA
    log:
        "logs/process_expression_data_chen.log"
    message:
        "Processing expression data"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/LUAD_Chen_expression_preprocessing.R \
            --expression_file {input.expression_file} \
            --exp_clean {output.expression_matrix} \
            --pca_plot {output.pca_plot} \
            --samples_file {output.samples_file} \
            2> {log}
        """