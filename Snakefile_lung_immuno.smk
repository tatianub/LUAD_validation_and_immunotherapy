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


NETWORKS_DIR = config["networks_dir"]
DATA_DIR = config["data_dir"]


# Expression and sample files for PANDA
IMMUNE_COHORT_TPM_FILE = os.path.join(
    DATA_DIR,
    "immune_cohort/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct"
)
CLINICAL_FILE_IMMUNO = os.path.join(DATA_DIR, "immune_cohort/clinical_filtered.txt")

DATASET_EXPRESSION_PANDA_FILE = os.path.join(PANDA_DIR,  "{dataset}_expression_for_PANDA.tsv")
DATASET_SAMPLES_PANDA_DATASET = os.path.join(PANDA_DIR, "{dataset}_samples_for_PANDA.tsv")
PCA_PLOTS_EXPRESSION_FILE = os.path.join(
    FIGURES_DIR,
    "{dataset}_pca_plots_expression.pdf"
)



rule all:
    input:
        expand(DATASET_EXPRESSION_PANDA_FILE, dataset=DATASETS),
        expand(DATASET_SAMPLES_PANDA_DATASET, dataset=DATASETS),
        expand(PCA_PLOTS_EXPRESSION_FILE, dataset=DATASETS)




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