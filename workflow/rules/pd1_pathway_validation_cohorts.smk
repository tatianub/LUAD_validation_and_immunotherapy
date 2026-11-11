###############################################################################
# PD-1 PATHWAY VALIDATION — LUAD COHORTS
###############################################################################
# Run univariate Cox regression and pathway GSEA for each validation cohort,
# then summarize the PD-1 pathway results across cohorts in a PDF plot.
###############################################################################

rule all_pd1_validation:
    input:
        PD1_VALIDATION_PLOT


###############################################################################
# 1. UNIVARIATE COX REGRESSION AND PATHWAY GSEA
###############################################################################
rule univariate_cox_and_gsea:
    input:
        indegree = INDEGREES_DATASET,
        clinical = CLINICAL_FILE_OUTPUT,
        pathway_gmt_file = PATHWAY_GMT_FILE
    output:
        results_output_file = GSEA_UNIVARIATE_COX_RESULTS_DATASET
    params:
        bin = config["bin"],
        cores = 10,
        covariates = "age,gender"
    wildcard_constraints:
        dataset = "|".join(re.escape(dataset) for dataset in VALIDATION_DATASETS),
        norm = "norm"
    log:
        "logs/univariate_cox_and_gsea_{dataset}_{norm}.log"
    message:
        "Performing univariate Cox regression and GSEA analysis on "
        "{wildcards.dataset} with {wildcards.norm} network"
    shell:
        """
        Rscript {params.bin}/univariate_cox_gsea_analysis.R \
            --indegree_file {input.indegree} \
            --clinical_file {input.clinical} \
            --pathway_gmt_file {input.pathway_gmt_file} \
            --covariates {params.covariates} \
            --cores {params.cores} \
            --results_output_file {output.results_output_file} \
            > {log} 2>&1
        """


###############################################################################
# 2. PLOT PD-1 PATHWAY VALIDATION RESULTS ACROSS COHORTS
###############################################################################
rule plot_univariate_cox_gsea_validation:
    input:
        result_files = expand(
            GSEA_UNIVARIATE_COX_RESULTS_DATASET,
            dataset=VALIDATION_DATASETS.keys(),
            norm=["norm"]
        )
    output:
        validation_plot_file = PD1_VALIDATION_PLOT
    params:
        bin = config["bin"],
        result_files = lambda wildcards, input: ",".join(input.result_files)
    log:
        "logs/plot_univariate_cox_gsea_validation.log"
    shell:
        """
        Rscript {params.bin}/plot_P1_outcome_validation_results.R \
            --result_files {params.result_files} \
            --validation_plot_file {output.validation_plot_file} \
            > {log} 2>&1
        """
