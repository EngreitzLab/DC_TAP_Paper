# DC-TAP-seq Validations

> [!WARNING]
> This directory is a work in progress and may not be complete or fully accurate.

This folder contains analysis code and notebook for DC-TAP-seq validation experiments.

## Overview

The analyses are organized in subdirectories named by date (`YYYY-MM-DD`) under the
`analyses/` directory. `README.md` file inside this folder contains a description
of the analysis performed on each date, as well as the code and results.

Results are stored in the `results/` directory, organized by date corresponding to
the analysis that generated them.

## Requirements

### Software Dependencies

| Software         | Version  |
| ---------------- | -------- |
| python           | ≥ 3.12   |
| bedtools         | ≥ 2.31.1 |
| blat             | ≥ 35     |
| bedgraphtobigwig | ≥ 377    |
| bigwigtobedgraph | ≥ 377    |
| fetchchromsizes  | ≥ 377    |

### Environment Setup

For installing the software dependencies listed above,
consider using a conda environment to manage bioinformatic related
software tools. Alternatively, visit
[UCSC's Utilities Downloads page]("https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads")
to download the required software binary files independent of a
package manager.

```bash
conda create -y -n dctap-valid-env \
  python=3.12 bedtools \
  ucsc-bedgraphtobigwig \
  ucsc-bigwigtobedgraph \
  ucsc-fetchchromsizes \
  zlib make cython \
  -c conda-forge -c bioconda
```

It is recommended to use a Python virtual environment.  
To create and activate an environment and install the required libraries:

```bash
cd validations
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Running notebooks

```bash
marimo edit analyses/YYYY-MM-DD/my_notebook.py
```

## Validation screen analyses (2026)

Analyses of the revision validation screen (WTC11, low MOI 1.1/1.5 and high MOI 11.3/13), merged from the EngreitzLab/DC-TAP-paper-validation repository (E. Mattei). R notebooks use `here()`, which resolves to this `validations/` folder (marked by `.here`); package versions are in `renv.lock`. Inputs (`data/`) and outputs (`results/`) are not in the repository; the Cell Ranger outputs are on Google Cloud Storage (set `VALIDATION_CELLRANGER_URI` for `scripts/bash/download_filtered_h5.sh`).


### [2026-05-11_cellranger](analyses/2026-05-11_cellranger/README.md)

Summary: documents the Cell Ranger count outputs used for downstream validation, the helper script that downloads filtered H5 matrices from Google Cloud Storage, and the analysis script that extracts per-sample RNA and CRISPR summary metrics.

Main code:
- `scripts/bash/download_filtered_h5.sh`
- `analyses/2026-05-11_cellranger/extract_metrics_summary.R`

Main outputs:
- `data/filtered_h5/`
- `results/cellranger/cellranger_metrics_summary.tsv`

### [2026-05-11_guide_pool_qc](analyses/2026-05-11_guide_pool_qc/README.md)

Summary: documents guide-pool alignment and QC inputs, the trimming/alignment commands, and the R Markdown used to summarize Bowtie output.

Main code:
- `analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd`

Main outputs:
- `results/guide_pool_qc/`

### [2026-05-11_qc_experiment](analyses/2026-05-11_qc_experiment/README.md)

Summary: documents merged Seurat object creation from filtered H5 inputs, sample and feature metadata integration (including GTF gene types), and QC comparisons against reference scRNA metrics.

Main code:
- `analyses/2026-05-11_qc_experiment/combine_h5_counts.Rmd`
- `analyses/2026-05-11_qc_experiment/qc_experiment.Rmd`
- `scripts/R/utils.R`

Main outputs:
- `results/qc_experiment/`

### [2026-05-11_SCEPTRE](analyses/2026-05-11_SCEPTRE/README.md)

Summary: SCEPTRE analysis framework for single-cell perturbation data across four validation samples split across low- and high-MOI groups. Supports multiple gRNA integration strategies (union, singleton) and replicate merging with sample batch correction.

Main code:
- `analyses/2026-05-11_SCEPTRE/01_create_sceptre_objects_merged_replicates.Rmd` - Parameterized SCEPTRE object creation
- `analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd` - Per-sample discovery pipeline
- `analyses/2026-05-11_SCEPTRE/03_run_sceptre_pipeline_singletons.Rmd` - Per-sample singleton gRNA strategy
- `analyses/2026-05-11_SCEPTRE/04_run_sceptre_pipeline_merge_replicates.Rmd` - Merged replicate analysis

Main outputs:
- `results/2026-05-11_SCEPTRE/objects/` - SCEPTRE objects (per-sample and merged)
- `results/2026-05-11_SCEPTRE/runs/` - Analysis results (discovery, calibration, power)
- `results/2026-05-11_SCEPTRE/tables/` - Summary tables

See [SCEPTRE analysis README](analyses/2026-05-11_SCEPTRE/README.md) for detailed commands to create all required SCEPTRE objects.

## Repository layout

- `analyses/`: analysis-specific notes and code.
- `data/`: local inputs used by the analyses.
- `results/`: generated outputs that are kept locally but ignored by git.
- `scripts/`: helper scripts and command-line entry points.






