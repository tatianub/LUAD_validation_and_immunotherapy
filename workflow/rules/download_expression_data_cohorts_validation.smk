
###############################################################################
# DATA DOWNLOAD — LUAD VALIDATION COHORTS
###############################################################################
#
# This file contains rules for downloading expression and clinical data
# used in the independent LUAD validation cohorts.
#
# DATA SOURCES
#
# 1. GEO datasets (GEOquery)
#    Datasets: GSE13213, GSE41271, GSE42127,
#              GSE50081, GSE37745, GSE30219
#    Source:   NCBI Gene Expression Omnibus (GEO)
#    URL:      https://www.ncbi.nlm.nih.gov/geo/
#    Method:   download_GEO_data.R
#    Files:    Expression, probe annotation and clinical RData
#
# 2. GSE81089
#    Source:   NCBI GEO (RNA-seq)
#    URL:      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089
#    Files:    NCBI-generated TPM expression matrix
#              GEO series matrix (clinical metadata)
#
# 3. LUAD_Okayama (GSE31210)
#    Source:   PRECOG database (Gentles et al., Nature Medicine, 2015)
#    URL:      https://precog.stanford.edu/
#    Files:    EntrezCDF-processed expression matrix (.RDS)
#              GEO series matrix (clinical metadata, from NCBI GEO)
#
# 4. LUAD_Chen
#    Source:   LOSTdb (Chen et al., 2020)
#    URL:      http://lostdbcancer.com:8080/
#    Archive:  Chen-NG_2020_LUAD_RNA.tar.gz
#    Files:    Expr_tumor.txt (expression data)
#              meta_tumor.txt (clinical metadata)
#    Additional clinical data:
#              41588_2019_569_MOESM3_ESM.xlsx
#    Source:   Nature Genetics supplementary material
#    DOI:      https://doi.org/10.1038/s41588-019-0569-6
#
#
# 5. LUAD_Dolgalev
#    Source:   LOSTdb (Dolgalev et al., 2023)
#    URL:      http://lostdbcancer.com:8080/
#    Archive:  Dolgalev_2023_LUAD_RNA.tar.gz
#    File:     meta_tumor_normal.txt (clinical metadata)
#    Expression: Figshare (https://figshare.com/ndownloader/files/37863348)
#    File:     rna-counts-normalized-n246.csv (expression matrix)
#
#
# ADDITIONAL INPUTS
#    Feature annotation files are stored in data/cohorts_info/:
#      - GSE13213_features.txt
#      - GSE41271_features.txt
#      - GSE42127_features.txt
#      - HGU133_features.txt (shared by GSE30219, GSE37745, GSE50081)
#      - Human.GRCh38.p13.annot.tsv (GSE81089)
#
#    Publication-derived clinical data are stored in resources/:
#      - GSE30129_clinical_data_pdf.txt (GSE30219)
#
###############################################################################

###############################################################################
# RULE ALL
###############################################################################

rule all_download:
    input:
        # LUAD_Chen
        CHEN_ORIGINAL_EXPRESSION_FILE,
        LUAD_Chen_CLINICAL_FILE,
        LUAD_Chen_CLINICAL_FILE_FROM_PUBLICATION,
         # GSE81089
        LUAD_GSE81089_ORIGINAL_EXPRESSION_FILE,
        LUAD_GSE81089_SERIES_MATRIX,
        LUAD_GSE81089_SAMPLES_FILE,
        LUAD_GSE81089_ANNOTATION_FILE,
        # LUAD_Okayama
        OKAYAMA_ORIGINAL_EXPRESSION_FILE,
        OKAYAMA_SERIES_MATRIX,
        # GEO datasets
        expand(
            DATASET_ORIGINAL_EXPRESSION_FILE,
            dataset=GEO_DATASETS
        ),
        expand(
            DATASET_ORIGINAL_PROBE_ANNOTATION_FILE,
            dataset=GEO_DATASETS
        ),
        expand(
            DATASET_ORIGINAL_CLINICAL_FILE,
            dataset=GEO_DATASETS
        ),
        expand(
            DATASET_FEATURE_FILE,
            dataset=GEO_FEATURE_DATASETS
        ),
        # LUAD_Dolgalev
        LUAD_Dolgalev_CLINICAL_FILE,
        DOLGALEV_ORIGINAL_EXPRESSION_FILE,
        # HGU133PLUS2 features
        FEATURES_HGU133PLUS2,
        GSE30219_PUBLICATION_CLINICAL


###############################################################################
# 1. LUAD_CHEN — LOSTdb EXPRESSION AND CLINICAL DATA
###############################################################################

rule download_LUAD_Chen_data:
    output:
        expression_file = CHEN_ORIGINAL_EXPRESSION_FILE,
        clinical_file = LUAD_Chen_CLINICAL_FILE
    log:
        "logs/download_LUAD_Chen_data.log"
    params:
        url = "http://lostdbcancer.com:8080/api/search/data/download/Chen-NG_2020_LUAD_RNA.tar.gz",
        outdir = CHEN_DIR
    shell:
        """
        (
            set -e

            mkdir -p {params.outdir:q}

            echo "Start: $(date)"
            echo "Downloading LUAD_Chen expression and clinical data"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                '{params.url}' \
                -o {params.outdir:q}/Chen-NG_2020_LUAD_RNA.tar.gz

            echo "Checking downloaded archive"
            gzip -t {params.outdir:q}/Chen-NG_2020_LUAD_RNA.tar.gz

            echo "Extracting expression and clinical files"
            tar -xzf {params.outdir:q}/Chen-NG_2020_LUAD_RNA.tar.gz \
                -C {params.outdir:q} \
                --strip-components=1 \
                data/Expr_tumor.txt data/meta_tumor.txt

            echo "Validating extracted files"
            test -s {output.expression_file:q}
            test -s {output.clinical_file:q}

            echo "Removing downloaded archive"
            rm {params.outdir:q}/Chen-NG_2020_LUAD_RNA.tar.gz

            echo "Output files:"
            ls -lh {output.expression_file:q} {output.clinical_file:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """

###############################################################################
# 2. LUAD_CHEN — PUBLICATION SUPPLEMENTARY CLINICAL DATA
###############################################################################

rule download_LUAD_Chen_supplementary_clinical:
    output:
        clinical_file = LUAD_Chen_CLINICAL_FILE_FROM_PUBLICATION
    log:
        "logs/download_LUAD_Chen_supplementary_clinical.log"
    params:
        url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0569-6/MediaObjects/41588_2019_569_MOESM3_ESM.xlsx"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.clinical_file:q})

            echo "Start: $(date)"
            echo "Downloading LUAD_Chen supplementary clinical data"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                '{params.url}' \
                -o {output.clinical_file:q}.tmp

            test -s {output.clinical_file:q}.tmp

            mv {output.clinical_file:q}.tmp \
               {output.clinical_file:q}

            echo "Output file:"
            ls -lh {output.clinical_file:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """


###############################################################################
# 3. GSE81089 — TPM EXPRESSION MATRIX
###############################################################################

rule download_GSE81089_TPM:
    output:
        expression_file = LUAD_GSE81089_ORIGINAL_EXPRESSION_FILE
    log:
        "logs/download_GSE81089_TPM.log"
    params:
        url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81089&format=file&type=rnaseq_counts&file=GSE81089_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.expression_file:q})

            echo "Start: $(date)"
            echo "Downloading GSE81089 TPM expression matrix"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                -o {output.expression_file:q}.gz.part \
                '{params.url}'

            echo "Validating gzip archive"
            gzip -t {output.expression_file:q}.gz.part

            echo "Decompressing expression matrix"
            gzip -dc {output.expression_file:q}.gz.part \
                > {output.expression_file:q}.part

            test -s {output.expression_file:q}.part

            mv {output.expression_file:q}.part \
               {output.expression_file:q}

            rm -f {output.expression_file:q}.gz.part

            echo "Output file:"
            ls -lh {output.expression_file:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """


###############################################################################
# 4. GSE81089 — GEO SERIES MATRIX
###############################################################################

rule download_GSE81089_series_matrix:
    output:
        matrix = LUAD_GSE81089_SERIES_MATRIX
    log:
        "logs/download_GSE81089_series_matrix.log"
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/matrix/GSE81089_series_matrix.txt.gz"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.matrix:q})

            echo "Start: $(date)"
            echo "Downloading GSE81089 GEO series matrix"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                -o {output.matrix:q}.gz.part \
                '{params.url}'

            echo "Validating gzip archive"
            gzip -t {output.matrix:q}.gz.part

            echo "Decompressing series matrix"
            gzip -dc {output.matrix:q}.gz.part \
                > {output.matrix:q}.part

            test -s {output.matrix:q}.part

            mv {output.matrix:q}.part {output.matrix:q}

            rm -f {output.matrix:q}.gz.part

            echo "Output file:"
            ls -lh {output.matrix:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """


###############################################################################
# GSE81089 — EXTRACT SAMPLE ANNOTATION FROM GEO SERIES MATRIX
###############################################################################

rule extract_GSE81089_sample_annotation:
    input:
        series_matrix = LUAD_GSE81089_SERIES_MATRIX
    output:
        sample_annotation = LUAD_GSE81089_SAMPLES_FILE
    log:
        "logs/extract_GSE81089_sample_annotation.log"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin:q}/extract_GSE81089_sample_annotation.R \
            --series_matrix {input.series_matrix:q} \
            --output_file {output.sample_annotation:q} \
            > {log:q} 2>&1
        """


rule download_GSE81089_annotation:
    output:
        annotation_file = LUAD_GSE81089_ANNOTATION_FILE
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/geo/rnaseq_counts/annot/Human.GRCh38.p13.annot.tsv.gz"
    log:
        "logs/download_GSE81089_annotation.log"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.annotation_file:q})

            wget -q '{params.url}' \
                -O {output.annotation_file:q}.gz.part

            gzip -t {output.annotation_file:q}.gz.part

            gzip -dc {output.annotation_file:q}.gz.part \
                > {output.annotation_file:q}.part

            mv {output.annotation_file:q}.part \
               {output.annotation_file:q}

            rm -f {output.annotation_file:q}.gz.part

            echo "Annotation downloaded:"
            ls -lh {output.annotation_file:q}
        ) > {log:q} 2>&1
        """

###############################################################################
# 5. LUAD_OKAYAMA — PRECOG EXPRESSION DATA
###############################################################################

rule download_Okayama_PRECOG:
    output:
        expression_file = OKAYAMA_ORIGINAL_EXPRESSION_FILE
    log:
        "logs/download_Okayama_PRECOG.log"
    params:
        url = "https://precog.stanford.edu/precogdata/Data/Lung.cancer.ADENO.GSE31210.HGU133Plus2_EntrezCDF.data.RDS"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.expression_file:q})

            echo "Start: $(date)"
            echo "Downloading LUAD_Okayama expression data from PRECOG"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                -o {output.expression_file:q}.part \
                '{params.url}'

            test -s {output.expression_file:q}.part

            mv {output.expression_file:q}.part \
               {output.expression_file:q}

            echo "Output file:"
            ls -lh {output.expression_file:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """


###############################################################################
# 6. LUAD_OKAYAMA — GEO SERIES MATRIX
###############################################################################

rule download_GSE31210_series_matrix:
    output:
        matrix = OKAYAMA_SERIES_MATRIX
    log:
        "logs/download_GSE31210_series_matrix.log"
    params:
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz"
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.matrix:q})

            echo "Start: $(date)"
            echo "Downloading GSE31210 GEO series matrix"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                -o {output.matrix:q}.gz.part \
                '{params.url}'

            echo "Validating gzip archive"
            gzip -t {output.matrix:q}.gz.part

            echo "Decompressing series matrix"
            gzip -dc {output.matrix:q}.gz.part \
                > {output.matrix:q}.part

            test -s {output.matrix:q}.part

            mv {output.matrix:q}.part {output.matrix:q}

            rm -f {output.matrix:q}.gz.part

            echo "Output file:"
            ls -lh {output.matrix:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """


###############################################################################
# 7. STANDARD GEO DATASETS — GEOquery
###############################################################################

rule download_from_GEO:
    wildcard_constraints:
        dataset = GEO_DATASET_PATTERN
    output:
        expression_file = DATASET_ORIGINAL_EXPRESSION_FILE,
        probe_annotation_file = DATASET_ORIGINAL_PROBE_ANNOTATION_FILE,
        clinical_file = DATASET_ORIGINAL_CLINICAL_FILE
    log:
        "logs/download_GEO_{dataset}.log"
    params:
        bin = config["bin"],
        gse_id = lambda wc: wc.dataset
    shell:
        """
        (
            set -e

            mkdir -p $(dirname {output.expression_file:q})

            echo "Start: $(date)"
            echo "Downloading GEO dataset: {params.gse_id}"
            echo "Using GEOquery"

            Rscript {params.bin:q}/download_GEO_data.R \
                --GSE_id {params.gse_id:q} \
                --expression_file {output.expression_file:q} \
                --probe_annotation_file {output.probe_annotation_file:q} \
                --clinical_file {output.clinical_file:q}

            echo "Validating output files"
            test -s {output.expression_file:q}
            test -s {output.probe_annotation_file:q}
            test -s {output.clinical_file:q}

            echo "Output files:"
            ls -lh \
                {output.expression_file:q} \
                {output.probe_annotation_file:q} \
                {output.clinical_file:q}

            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """



rule prepare_GEO_features:
    input:
        probe_file = DATASET_PROBE_ANNOTATION_FILE
    output:
        features_file = DATASET_FEATURE_FILE
    wildcard_constraints:
        dataset = "|".join(GEO_FEATURE_DATASETS)
    params:
        bin = config["bin"]
    log:
        "logs/prepare_GEO_features_{dataset}.log"
    shell:
        """
        mkdir -p $(dirname {output.features_file:q})

        Rscript {params.bin}/prepare_GEO_features.R \
            --dataset {wildcards.dataset} \
            --probe_file {input.probe_file:q} \
            --output_file {output.features_file:q} \
            > {log:q} 2>&1
        """

rule prepare_HGU133_features:
    output:
        features_file = FEATURES_HGU133PLUS2
    params:
        bin = config["bin"]
    log:
        "logs/prepare_HGU133_features.log"
    shell:
        """
        mkdir -p $(dirname {output.features_file:q})

        Rscript {params.bin}/prepare_HGU133_features.R \
            --output_file {output.features_file:q} \
            > {log:q} 2>&1
        """

###############################################################################
# LUAD_DOLGALEV — LOSTdb CLINICAL METADATA
###############################################################################

rule download_LUAD_Dolgalev_clinical:
    output:
        clinical_file = LUAD_Dolgalev_CLINICAL_FILE
    log:
        "logs/download_LUAD_Dolgalev_clinical.log"
    params:
        url = "http://lostdbcancer.com:8080/api/search/data/download/Dolgalev_2023_LUAD_RNA.tar.gz",
        outdir = DOLGALEV_DIR
    shell:
        """
        (
            set -e

            mkdir -p {params.outdir:q}

            archive={params.outdir:q}/Dolgalev_2023_LUAD_RNA.tar.gz

            echo "Start: $(date)"
            echo "Downloading LUAD_Dolgalev clinical metadata"
            echo "Source: {params.url}"

            curl -fLSs --retry 3 --connect-timeout 30 \
                '{params.url}' \
                -o "$archive"

            echo "Checking archive"
            gzip -t "$archive"

            echo "Locating meta_tumor_normal.txt"

            clinical_member=$(
                tar -tzf "$archive" |
                grep -E '(^|/)meta_tumor_normal[.]txt$' |
                head -n 1
            )

            if [ -z "$clinical_member" ]; then
                echo "ERROR: meta_tumor_normal.txt not found"
                echo "Archive contents:"
                tar -tzf "$archive"
                exit 1
            fi

            echo "Found: $clinical_member"

            echo "Extracting clinical metadata"
            tar -xOzf "$archive" "$clinical_member" \
                > {output.clinical_file:q}.part

            test -s {output.clinical_file:q}.part

            mv {output.clinical_file:q}.part {output.clinical_file:q}

            echo "Removing archive"
            rm "$archive"

            echo "Output:"
            ls -lh {output.clinical_file:q}

            echo "Completed: $(date)"

        ) > {log:q} 2>&1
        """
###############################################################################
# LUAD_DOLGALEV — FIGSHARE EXPRESSION MATRIX
###############################################################################

rule download_LUAD_Dolgalev_expression:
    output:
        expression_file = DOLGALEV_ORIGINAL_EXPRESSION_FILE
    log:
        "logs/download_LUAD_Dolgalev_expression.log"
    params:
        url = "https://figshare.com/ndownloader/files/37863348"
    shell:
        """
        (
            set -e
            mkdir -p $(dirname {output.expression_file:q})

            echo "Start: $(date)"
            echo "Downloading LUAD_Dolgalev expression matrix from Figshare"

            curl -fLSs --retry 3 --connect-timeout 30 \
                '{params.url}' \
                -o {output.expression_file:q}.part

            test -s {output.expression_file:q}.part
            mv {output.expression_file:q}.part {output.expression_file:q}

            ls -lh {output.expression_file:q}
            echo "Completed: $(date)"
        ) > {log:q} 2>&1
        """

# GSE30219 PUBLICATION-DERIVED CLINICAL DATA

rule copy_GSE30219_publication_clinical:
    input:
        clinical_file = GSE30219_PUBLICATION_CLINICAL_RESOURCE
    output:
        clinical_file = GSE30219_PUBLICATION_CLINICAL
    shell:
        """
        mkdir -p $(dirname {output.clinical_file:q})
        cp {input.clinical_file:q} {output.clinical_file:q}
        """