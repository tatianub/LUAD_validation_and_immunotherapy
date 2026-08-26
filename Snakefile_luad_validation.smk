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


rule all:
    input:
        expand(DATASET_ORIGINAL_EXPRESSION_FILE, dataset=DATASETS),
        expand(DATASET_ORIGINAL_PROBE_ANNOTATION_FILE, dataset=DATASETS),
        expand(DATASET_ORIGINAL_CLINICAL_FILE, dataset=DATASETS)

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
