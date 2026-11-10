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
import re
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

DATASET_COHORTS_FILE = os.path.join(DATA_DIR, "cohorts_info", "{dataset}_cohorts.txt")
# DATASET_FEATURE_FILE = os.path.join(DATA_DIR, "{dataset}", "{dataset}_features.txt")
EXPRESSION_PCA_FIGURE = os.path.join(FIGURES_DIR, "{dataset}_pca_analysis.pdf")
DATASET_EXPRESSION_PANDA_FILE = os.path.join(PANDA_DIR,  "{dataset}_expression_for_PANDA.tsv")
DATASET_SAMPLES_PANDA_DATASET = os.path.join(PANDA_DIR, "{dataset}_samples_for_PANDA.tsv")


# Dolgalev preprocessing
DOLGALEV_DIR = os.path.join(DATA_DIR, "LUAD_Dolgalev")
# PROCESSING DOLGALEV DATA
DOLGALEV_ORIGINAL_EXPRESSION_FILE = os.path.join(DOLGALEV_DIR, "rna-counts-normalized-n246.csv")
DOLGALEV_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "LUAD_Dolgalev_pca_analysis.pdf")
DOLGALEV_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "LUAD_Dolgalev_expression_for_PANDA.tsv")
DOLGALEV_SAMPLES_PANDA = os.path.join(PANDA_DIR, "LUAD_Dolgalev_samples_for_PANDA.tsv")

DATASET_EXPRESSION_PANDA_FILE = os.path.join(PANDA_DIR,  "{dataset}_expression_for_PANDA.tsv")
DATASET_SAMPLES_PANDA_DATASET = os.path.join(PANDA_DIR, "{dataset}_samples_for_PANDA.tsv")

# PROCESSING CHEN DATA

CHEN_DIR = os.path.join(DATA_DIR, "LUAD_Chen")
CHEN_ORIGINAL_EXPRESSION_FILE = os.path.join(CHEN_DIR, "Expr_tumor.txt")
CHEN_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "LUAD_Chen_pca_analysis.pdf")
CHEN_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "LUAD_Chen_expression_for_PANDA.tsv")
CHEN_SAMPLES_PANDA = os.path.join(PANDA_DIR, "LUAD_Chen_samples_for_PANDA.tsv")



# PROCESSING OKAYAMA DATA
OKAYAMA_DIR = os.path.join(DATA_DIR, "LUAD_Okayama")
OKAYAMA_ORIGINAL_EXPRESSION_FILE = os.path.join(OKAYAMA_DIR, "Lung.cancer.ADENO.GSE31210.HGU133Plus2_EntrezCDF.data.RDS")
OKAYAMA_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "LUAD_Okayama_pca_analysis.pdf")
OKAYAMA_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "LUAD_Okayama_expression_for_PANDA.tsv")
OKAYAMA_SAMPLES_PANDA = os.path.join(PANDA_DIR, "LUAD_Okayama_samples_for_PANDA.tsv")

# PROCESSING GSE81089 DATA (custom preprocessing)
GSE81089_DIR = os.path.join(DATA_DIR, "GSE81089")
LUAD_GSE81089_ORIGINAL_EXPRESSION_FILE = os.path.join(GSE81089_DIR, "GSE81089_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
LUAD_GSE81089_ANNOTATION_FILE = os.path.join(DATA_DIR, "cohorts_info", "Human.GRCh38.p13.annot.tsv")
LUAD_GSE81089_SAMPLES_FILE = os.path.join(GSE81089_DIR, "sample_annotation_GSE81089.txt")
LUAD_GSE81089_OUTPUT_PCA_FILE = os.path.join(FIGURES_DIR, "GSE81089_pca_analysis.pdf")
GSE81089_EXPRESSION_PANDA = os.path.join(PANDA_DIR, "GSE81089_expression_for_PANDA.tsv")
GSE81089_SAMPLES_PANDA = os.path.join(PANDA_DIR, "GSE81089_samples_for_PANDA.tsv")

## FEATURE FILES

FEATURES_HGU133PLUS2 = os.path.join(DATA_DIR, "cohorts_info", "HGU133_features.txt")
FEATURES_GSE13213 = os.path.join(DATA_DIR, "cohorts_info", "GSE13213_features.txt")
FEATURES_GSE41271 = os.path.join(DATA_DIR, "cohorts_info", "GSE41271_features.txt")
FEATURES_GSE42127 = os.path.join(DATA_DIR, "cohorts_info", "GSE42127_features.txt")

# Map each dataset to the appropriate feature file (shared files allowed).
DATASET_FEATURE_FILES = {
    "GSE13213": FEATURES_GSE13213,
    "GSE41271": FEATURES_GSE41271,
    "GSE42127": FEATURES_GSE42127,
    "GSE50081": FEATURES_HGU133PLUS2,
    "GSE37745": FEATURES_HGU133PLUS2,
    "GSE30219": FEATURES_HGU133PLUS2
}

PREPROCESS_EXPRESSION_DATASETS = [
    dataset for dataset in DATASETS if dataset in DATASET_FEATURE_FILES
]

PREPROCESS_EXPRESSION_DATASET_PATTERN = "|".join(
    re.escape(dataset) for dataset in PREPROCESS_EXPRESSION_DATASETS
)


def get_dataset_feature_file(wildcards):
    dataset = wildcards.dataset
    if dataset not in DATASET_FEATURE_FILES:
        raise ValueError(
            f"No feature file configured for dataset '{dataset}'. "
            f"Update DATASET_FEATURE_FILES in Snakefile_luad_validation.smk."
        )
    return DATASET_FEATURE_FILES[dataset]

GSE30219_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE30219", "GSE30219_clinical.RData")
GSE30219_PUBLICATION_CLINICAL = os.path.join(DATA_DIR, "GSE30219", "GSE30129_clinical_data_pdf.txt")
GSE13213_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE13213", "GSE13213_clinical.RData")
GSE41271_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE41271", "GSE41271_clinical.RData")
GSE42127_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE42127", "GSE42127_clinical.RData")
GSE50081_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE50081", "GSE50081_clinical.RData")
GSE37745_CLINICAL_FILE = os.path.join(DATA_DIR, "GSE37745", "GSE37745_clinical.RData")
LUAD_GSE81089_SERIES_MATRIX = os.path.join(DATA_DIR, "GSE81089", "GSE81089_series_matrix.txt")
LUAD_Dolgalev_CLINICAL_FILE = os.path.join(DATA_DIR, "LUAD_Dolgalev", "meta_tumor_normal.txt")
LUAD_Okayama_SERIES_MATRIX = os.path.join(DATA_DIR, "LUAD_Okayama", "GSE31210_series_matrix.txt")
LUAD_Chen_CLINICAL_FILE = os.path.join(DATA_DIR, "LUAD_Chen", "meta_tumor.txt")
LUAD_Chen_CLINICAL_FILE_FROM_PUBLICATION = os.path.join(DATA_DIR, "LUAD_Chen", "41588_2019_569_MOESM3_ESM.xlsx")
TCGA_LUAD_CLINICAL_FILE = os.path.join(DATA_DIR, "TCGA_LUAD", "clinical.txt")
PD1_VALIDATION_PLOT = os.path.join(FIGURES_DIR, "gsea_univariate_cox_results_validation.pdf")
CLINICAL_FILE_OUTPUT = os.path.join(
    DATA_DIR,
    "{dataset}",
    "clinical_subset_all.txt"
)

VALIDATION_DATASETS = {
    "GSE30219": {
        "input": {
            "clinical_file": GSE30219_CLINICAL_FILE,
            "publication_clinical": GSE30219_PUBLICATION_CLINICAL,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
            "--clinical_file_from_publication": "publication_clinical",
        },
    },

    "GSE13213": {
        "input": {
            "clinical_file": GSE13213_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },

    "GSE41271": {
        "input": {
            "clinical_file": GSE41271_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },

    "GSE42127": {
        "input": {
            "clinical_file": GSE42127_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },

    "GSE50081": {
        "input": {
            "clinical_file": GSE50081_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },

    "GSE37745": {
        "input": {
            "clinical_file": GSE37745_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },

    "GSE81089": {
        "input": {
            "gse_series_matrix": LUAD_GSE81089_SERIES_MATRIX,
        },
        "arg_map": {
            "--gse_series_matrix_file": "gse_series_matrix",
        },
    },

    "LUAD_Dolgalev": {
        "input": {
            "clinical_file": LUAD_Dolgalev_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },
    "LUAD_Okayama": {
        "input": {
            "gse_series_matrix": LUAD_Okayama_SERIES_MATRIX,
        },
        "arg_map": {
            "--gse_series_matrix_file": "gse_series_matrix",
        },
    },
    "LUAD_Chen": {
        "input": {
            "clinical_file": LUAD_Chen_CLINICAL_FILE,
            "publication_clinical": LUAD_Chen_CLINICAL_FILE_FROM_PUBLICATION,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
            "--clinical_file_from_publication": "publication_clinical",
        },
    },
    "TCGA_LUAD": {
        "input": {
            "clinical_file": TCGA_LUAD_CLINICAL_FILE,
        },
        "arg_map": {
            "--clinical_file": "clinical_file",
        },
    },
}

def make_clinical_args(wildcards, input):
    return " ".join(
        f"{flag} {input[input_name]}"
        for flag, input_name in VALIDATION_DATASETS[wildcards.dataset]["arg_map"].items()
    )




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
        CHEN_SAMPLES_PANDA,
        expand(EXPRESSION_PCA_FIGURE, dataset=DATASETS),
        expand(DATASET_EXPRESSION_PANDA_FILE, dataset=DATASETS),
        expand(DATASET_SAMPLES_PANDA_DATASET, dataset=DATASETS),
        LUAD_GSE81089_OUTPUT_PCA_FILE,
        GSE81089_EXPRESSION_PANDA,
        GSE81089_SAMPLES_PANDA,
        OKAYAMA_OUTPUT_PCA_FILE,
        OKAYAMA_EXPRESSION_PANDA,
        OKAYAMA_SAMPLES_PANDA,
        expand(
            CLINICAL_FILE_OUTPUT,
            dataset=VALIDATION_DATASETS.keys()
        )


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


rule luad_GSE81089_expression_analysis:
    """
    Process LUAD GSE81089 expression data to perform PCA analysis
    """
    input:
        expression_file = LUAD_GSE81089_ORIGINAL_EXPRESSION_FILE,
        annotation_file = LUAD_GSE81089_ANNOTATION_FILE,
        samples_file = LUAD_GSE81089_SAMPLES_FILE
    output:
        pca_plot = LUAD_GSE81089_OUTPUT_PCA_FILE,
        expression_matrix = GSE81089_EXPRESSION_PANDA,
        samples_file = GSE81089_SAMPLES_PANDA
    log:
        "logs/process_expression_LUAD_GSE81089_expression_data.log"
    message:
        "Processing expression data"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/LUAD_GSE81089_expression_preprocessing.R \
            --expression_file {input.expression_file} \
            --annotation_file {input.annotation_file} \
            --samples_file {input.samples_file} \
            --output_pca_file {output.pca_plot} \
            --output_expression_clean_file {output.expression_matrix} \
            --output_samples_file {output.samples_file} \
            2> {log}
        """

rule preprocess_expression:
    wildcard_constraints:
        dataset = PREPROCESS_EXPRESSION_DATASET_PATTERN

    input:
        expression_file = DATASET_ORIGINAL_EXPRESSION_FILE,
        cohorts_file = DATASET_COHORTS_FILE,
        probe_annotation_file = get_dataset_feature_file
    output:
        expression_for_PANDA = DATASET_EXPRESSION_PANDA_FILE,
        samples_for_PANDA = DATASET_SAMPLES_PANDA_DATASET,
        pca_figure = EXPRESSION_PCA_FIGURE
    params:
        bin = config["bin"],
        gse_id = lambda wildcards: wildcards.dataset,
        output_dir = PANDA_DIR
        
    log:
        "logs/preprocess_expression_{dataset}.log"

    shell:
        """
        Rscript {params.bin}/prepare_expression_with_probes.R \
            --GSE_id {params.gse_id} \
            --expression_file {input.expression_file} \
            --cohorts_file {input.cohorts_file} \
            --probe_file {input.probe_annotation_file} \
            --exp_clean {output.expression_for_PANDA} \
            --samples_file {output.samples_for_PANDA} \
            --pca_plot {output.pca_figure} \
            --output_dir {params.output_dir} \
            > {log} 2>&1
        """

rule luad_Okayama_expression_analysis:
    """
    Process LUAD Okayama expression data to perform PCA analysis
    """
    input:
        expression_file = OKAYAMA_ORIGINAL_EXPRESSION_FILE
    output:
        pca_plot = OKAYAMA_OUTPUT_PCA_FILE,
        expression_matrix = OKAYAMA_EXPRESSION_PANDA,
        samples_file = OKAYAMA_SAMPLES_PANDA
    log:
        "logs/process_expression_LUAD_Okayama_expression_data.log"
    message:
        "Processing expression data"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/LUAD_Okayama_expression_preprocessing.R \
            --expression_file {input.expression_file} \
            --output_pca_file {output.pca_plot} \
            --output_expression_clean_file {output.expression_matrix} \
            --output_samples_file {output.samples_file} \
            2> {log}
        """

rule prepare_clinical_data:
    input:
        unpack(lambda wc: VALIDATION_DATASETS[wc.dataset]["input"])
    output:
        clinical_output=CLINICAL_FILE_OUTPUT
    params:
        bin=config["bin"],
        script=lambda wc: f"{wc.dataset}_select_survival_info_clinical.R",
        args=make_clinical_args
    log:
        "logs/prepare_{dataset}_clinical_data.log"
    message:
        "Preparing {wildcards.dataset} clinical data"
    shell:
        """
        Rscript {params.bin}/{params.script} \
            {params.args} \
            --clinical_output_file {output.clinical_output}
        """
