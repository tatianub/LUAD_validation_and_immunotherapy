# LUAD validation and immunotherapy analysis

Validation of PD-1 pathway regulatory heterogeneity in independent LUAD cohorts and analysis of its association with response to immune checkpoint inhibitor therapy.

This repository contains a Snakemake-based computational workflow for integrating transcriptomic and clinical data, inferring regulatory networks, analyzing TF/gene activity, and testing associations between pathway-level dysregulation and immunotherapy response in lung adenocarcinoma.

## Overview

The workflow combines:

- cohort expression data and clinical annotations
- curated immune cell and pathway resources
- TF priors and gene regulatory network inference
- PANDA/LIONESS-style network generation and normalization
- indegree/outdegree network statistics
- LIMMA-based differential analyses
- FGSEA pathway enrichment
- TF and pathway summary plots for immunotherapy-related questions

The current project targets the SU2C-MARK cohort and focuses on LUAD subtypes and PD-1 treatment-related comparisons.

## Repository structure

```text
LUAD_validation_and_immunotherapy/
├── Snakefile_lung_immuno.smk
├── config_lung_immuno.yaml
├── LICENSE
├── README.md
├── workflow/
│   └── bin/
│       ├── combine_networks.R
│       ├── differential_edges_TFs.R
│       ├── immune_cohort_data_preprocessing.R
│       ├── lung_immuno_limma_indegrees.R
│       ├── lung_immuno_indegree_all_models.R
│       ├── make_MYC_four_panel_summary.R
│       ├── network_analysis.R
│       ├── plot_p1_TF_edges_expression.R
│       ├── plot_PD1_signaling_TF_enrichment_results.R
│       ├── quantile_normalize_networks.R
│       └── tf_immune_infiltration_correlations.R
├── data/
│   └── immune_cohort/
│       ├── SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct
│       ├── SU2C-MARK_Harmonized_Curated_Sets_SF_v1.txt
│       ├── clinical_filtered.txt
│       └── ...
├── priors/
├── resources/
│   └── c2.cp.reactome.v7.1.symbols.gmt
├── panda_input/
├── networks/
├── merged_networks/
├── results_lung_immuno/
├── figures_lung_immuno/
├── logs/
├── .gitignore
└── .snakemake/
```

## Pipeline overview

The Snakemake workflow performs the following major stages:

1. Preprocess immune cohort expression data.
2. Create regulatory network inputs for PANDA/LIONESS analysis.
3. Build, merge, and normalize networks.
4. Compute indegrees and outdegrees for TFs/genes.
5. Run LIMMA analyses for subtype and treatment models.
6. Conduct FGSEA pathway enrichment on network-derived signals.
7. Evaluate TF/pathway associations with immune infiltration and clinical features.
8. Produce summary plots and result tables for downstream interpretation.

The pipeline is configured in [config_lung_immuno.yaml](config_lung_immuno.yaml) and initiated from [Snakefile_lung_immuno.smk](Snakefile_lung_immuno.smk).

## Input data

The repository currently uses the SU2C-MARK cohort and expects input resources in:

- `data/immune_cohort/`
- `priors/`
- `resources/`

Key files include:

- TPM expression matrix
- clinical annotation tables
- immune infiltration or curated set files
- TF prior tables
- pathway GMT files

## Requirements

This project depends on:

- Snakemake
- R and Bioconductor-based packages
- Conda environment support for workflow dependencies

The workflow header notes the following environment setup:

```bash
module load snakemake/7.23.1-foss-2022a
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
```

## Running the pipeline

From the project root:

```bash
snakemake --use-conda --conda-frontend conda --cores 1 -s Snakefile_lung_immuno.smk
```

Dry run:

```bash
snakemake --cores 1 -np -s Snakefile_lung_immuno.smk
```

## Outputs

Generated outputs are organized into:

- `results_lung_immuno/`: LIMMA results, FGSEA results, network statistics, and summary tables
- `figures_lung_immuno/`: plots and publication-style figures
- `panda_input/`: PANDA/LIONESS input files
- `networks/`: per-dataset regulatory networks
- `merged_networks/`: merged and normalized network objects
- `logs/`: Snakemake and task logs

## Notes

- Network outputs are already built and stored in the repository, and `build_network_outputs` is set to `false` in the configuration by default.

## License

This project is distributed under the license included in the repository. See [LICENSE](LICENSE).


