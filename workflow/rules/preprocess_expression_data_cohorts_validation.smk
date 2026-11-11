###############################################################################
# EXPRESSION PREPROCESSING — LUAD VALIDATION COHORTS
###############################################################################
# Preprocess expression data from ten independent LUAD cohorts for PANDA.
# Raw expression files are produced by rules/download_expression_data_cohorts_validation.smk

#
# Outputs for each cohort:
#   - {dataset}_expression_for_PANDA.tsv  (processed gene expression)
#   - {dataset}_samples_for_PANDA.tsv     (samples retained for PANDA)
#   - {dataset}_pca_analysis.pdf         (PCA quality-control figure)
#
# Six standard GEO cohorts use a shared probe-to-gene preprocessing script.
# The remaining four cohorts use dataset-specific R scripts.
#
###############################################################################

###############################################################################
# TARGET: ALL EXPRESSION PREPROCESSING OUTPUTS
###############################################################################

rule all_expression:
    input:
        expand(DATASET_EXPRESSION_PANDA_FILE, dataset=PREPROCESS_EXPRESSION_DATASETS),
        expand(DATASET_SAMPLES_PANDA_DATASET, dataset=PREPROCESS_EXPRESSION_DATASETS),
        expand(EXPRESSION_PCA_FIGURE, dataset=PREPROCESS_EXPRESSION_DATASETS),
        DOLGALEV_OUTPUT_PCA_FILE,
        DOLGALEV_EXPRESSION_PANDA,
        DOLGALEV_SAMPLES_PANDA,
        CHEN_OUTPUT_PCA_FILE,
        CHEN_EXPRESSION_PANDA,
        CHEN_SAMPLES_PANDA,
        LUAD_GSE81089_OUTPUT_PCA_FILE,
        GSE81089_EXPRESSION_PANDA,
        GSE81089_SAMPLES_PANDA,
        OKAYAMA_OUTPUT_PCA_FILE,
        OKAYAMA_EXPRESSION_PANDA,
        OKAYAMA_SAMPLES_PANDA


###############################################################################
# 1. LUAD_DOLGALEV — NORMALIZED RNA-SEQ COUNTS
#    Input: rna-counts-normalized-n246.csv (Figshare)
###############################################################################

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
            > {log} 2>&1
        """


###############################################################################
# 2. LUAD_CHEN — LOSTdb EXPRESSION DATA
#    Input: Expr_tumor.txt
###############################################################################

rule chen_luad_expression_analysis:
    input:
        expression_file = CHEN_ORIGINAL_EXPRESSION_FILE
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
            > {log} 2>&1
        """


###############################################################################
# 3. GSE81089 — RNA-SEQ TPM EXPRESSION
#    Inputs: NCBI TPM matrix, GRCh38 gene annotation, sample metadata
###############################################################################

rule luad_GSE81089_expression_analysis:
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
            > {log} 2>&1
        """


###############################################################################
# 4. STANDARD GEO DATASETS — MICROARRAY EXPRESSION
#    GSE13213, GSE41271, GSE42127, GSE50081, GSE37745, GSE30219
#    Probe annotation is selected by get_dataset_feature_file().
#    GSE50081, GSE37745 and GSE30219 share HGU133_features.txt.
#    Each GEO dataset is processed as one cohort (no cohort-selection file).
###############################################################################

rule preprocess_expression:
    wildcard_constraints:
        dataset = PREPROCESS_EXPRESSION_DATASET_PATTERN
    input:
        expression_file = DATASET_ORIGINAL_EXPRESSION_FILE,
        probe_annotation_file = get_dataset_feature_file
    output:
        expression_for_PANDA = DATASET_EXPRESSION_PANDA_FILE,
        samples_for_PANDA = DATASET_SAMPLES_PANDA_DATASET,
        pca_figure = EXPRESSION_PCA_FIGURE
    params:
        bin = config["bin"],
        gse_id = lambda wildcards: wildcards.dataset
    log:
        "logs/preprocess_expression_{dataset}.log"
    shell:
        """
        Rscript {params.bin}/prepare_expression_with_probes.R \
            --GSE_id {params.gse_id} \
            --expression_file {input.expression_file} \
            --probe_file {input.probe_annotation_file} \
            --exp_clean {output.expression_for_PANDA} \
            --samples_file {output.samples_for_PANDA} \
            --pca_plot {output.pca_figure} \
            > {log} 2>&1
        """


###############################################################################
# 5. LUAD_OKAYAMA (GSE31210) — PRECOG EXPRESSION DATA
#    Input: EntrezCDF-processed expression matrix (.RDS)
###############################################################################

rule luad_Okayama_expression_analysis:
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
            > {log} 2>&1
        """
