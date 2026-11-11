###############################################################################
# CLINICAL DATA PREPROCESSING — LUAD VALIDATION COHORTS
###############################################################################
# Prepare cohort-specific survival/clinical data using the corresponding
# R script: bin/{dataset}_select_survival_info_clinical.R
#
# Input files and command-line arguments are selected through VALIDATION_DATASETS
# and make_clinical_args(), both defined in rules/common.smk.
#
# Output: data/{dataset}/clinical_subset_all.txt
#
###############################################################################


rule all_clinical:
    input:
        expand(CLINICAL_FILE_OUTPUT, dataset=VALIDATION_DATASETS.keys())


###############################################################################
# 1. PREPARE COHORT-SPECIFIC CLINICAL AND SURVIVAL DATA
###############################################################################
rule prepare_clinical_data:
    wildcard_constraints:
        dataset = "|".join(re.escape(dataset) for dataset in VALIDATION_DATASETS)
    input:
        unpack(lambda wc: VALIDATION_DATASETS[wc.dataset]["input"])
    output:
        clinical_output = CLINICAL_FILE_OUTPUT
    params:
        bin = config["bin"],
        script = lambda wc: f"{wc.dataset}_select_survival_info_clinical.R",
        args = make_clinical_args
    log:
        "logs/prepare_{dataset}_clinical_data.log"
    message:
        "Preparing {wildcards.dataset} clinical data"
    shell:
        """
        Rscript {params.bin}/{params.script} \
            {params.args} \
            --clinical_output_file {output.clinical_output} \
            > {log} 2>&1
        """


###############################################################################
# COPY TCGA_LUAD CLINICAL DATA FROM RESOURCES
###############################################################################

rule copy_TCGA_LUAD_clinical:
    input:
        clinical_file = TCGA_LUAD_CLINICAL_RESOURCE
    output:
        clinical_file = TCGA_LUAD_CLINICAL_FILE
    message:
        "Copying TCGA_LUAD clinical data from resources"
    shell:
        """
        mkdir -p $(dirname {output.clinical_file:q})
        cp {input.clinical_file:q} {output.clinical_file:q}
        """