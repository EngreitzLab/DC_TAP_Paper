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

## Validation screen and CTCF analyses (2026)

Analyses of the WTC11 validation screen (30 selected element–gene pairs, two replicates each at low and high MOI) and of the CTCF ChIP-seq experiment at the FAM111A and CDC37 CTCF elements.

R notebooks use `here()`, which resolves to this `validations/` folder (marked by `.here`). R package versions are in `renv.lock` (`renv::restore()`); rendering needs pandoc. Inputs (`data/`) and outputs (`results/`) are not tracked: the scripts in `scripts/bash/` download the public inputs from GEO [GSE303901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901), GENCODE and UCSC. The small tables behind each figure are tracked in the `source_data/` folder of each analysis.

### Samples

The code uses the Cell Ranger sample IDs. GEO names the samples by their measured MOI. One row per 10x channel is in `metadata/validation_screen_samples.tsv`.

| Sample ID | Replicate | MOI | 10x channels | GEO sample names |
|---|---|---|---|---|
| `MOI1_1` | 1 | 1.1 (low) | A1–A6 | `validation_moi1_1_rep01_a1` … `a6` |
| `MOI1_5` | 2 | 1.5 (low) | B1–B6 | `validation_moi1_5_rep02_b1` … `b6` |
| `MOI11_3` | 1 | 12.6 (high) | A7–A8 | `validation_moi12_6_rep01_a7`, `a8` |
| `MOI13` | 2 | 14.4 (high) | B7–B8 | `validation_moi14_4_rep02_b7`, `b8` |

The CTCF ChIP-seq libraries are listed in `metadata/ctcf_chip_libraries.tsv`.

### Running the analyses

From `validations/`:

```bash
# Figures and tables from the tracked source data (minutes)
bash scripts/bash/download_references.sh
bash scripts/bash/render_validation_figures.sh

# Full validation-screen pipeline from the GEO count matrices (hours)
bash scripts/bash/download_filtered_h5.sh
bash scripts/bash/run_validation_pipeline.sh 4

# CTCF signal from the GEO bigWig files
python analyses/2026-06-29_CTCF_ChIP_validation/ctcf_signal_from_bigwigs.py --from-geo
```

### [2026-05-11_cellranger](analyses/2026-05-11_cellranger/README.md)

Downloads the filtered Cell Ranger matrices from GEO and plots per-sample RNA and CRISPR metrics.

### [2026-05-11_guide_pool_qc](analyses/2026-05-11_guide_pool_qc/README.md)

Guide representation and alignment mismatches in the two plasmid guide pools.

### [2026-05-11_qc_experiment](analyses/2026-05-11_qc_experiment/README.md)

Merged Seurat object of the 16 channels, cell QC metrics, and comparison of readout-gene UMIs with a reference scRNA-seq dataset.

### [2026-05-11_SCEPTRE](analyses/2026-05-11_SCEPTRE/README.md)

SCEPTRE analysis of the selected pairs: one run per sample, and one per MOI with both replicates (sample as a covariate).

### [2026-05-13_get_statistics](analyses/2026-05-13_get_statistics/README.md)

Heatmap of the 30 selected pairs in the pooled low-MOI analysis, and counts of significant pairs per analysis and element category.

### [2026-05-14_compare_low_high_es](analyses/2026-05-14_compare_low_high_es/README.md)

Effect sizes of the selected hits: screen vs low MOI, screen vs high MOI, and low vs high MOI.

### [2026-06-29_CTCF_ChIP_validation](analyses/2026-06-29_CTCF_ChIP_validation/README.md)

CTCF ChIP-seq signal over the FAM111A and CDC37 CTCF elements in each library.

### [2026-06-30_CTCF_Plots](analyses/2026-06-30_CTCF_Plots/README.md)

CTCF ChIP-seq boxplots and per-guide CRISPRi effects around the CTCF motifs.

### [2026-10-23_ctcf_and_validation_statistics](analyses/2026-10-23_ctcf_and_validation_statistics/README.md)

Statistical tests of the CTCF ChIP-seq signal, and the table of the 30 selected pairs with original and validation estimates.

## Repository layout

- `analyses/`: analysis-specific notes and code.
- `metadata/`: sample sheets, the 30 selected pairs, the guide library and the readout genes.
- `data/`: local inputs used by the analyses.
- `results/`: generated outputs that are kept locally but ignored by git.
- `scripts/`: helper scripts and command-line entry points.
- `workflows/`: WDL workflow that makes the CTCF ChIP-seq bigWig files.






