# DC-TAP Workflow: Analyzing Unbiased TAP-Seq Screens

This repository contains the workflow for analyzing datasets generated from an unbiased DC TAP-Seq screen on K562 and WTC11. The final output file can be found in `results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`. An explanation of each column in that table can be found in `resources/formatting_dc_tap_results/columns_of_Final_DC_TAP_Seq_Results_table.txt`

## Data and code availability

- **Manuscript:** [A systematic survey of distal element-gene regulatory interactions with Direct-Capture Targeted Perturb-seq](https://www.biorxiv.org/content/10.1101/2025.09.16.676677v2)
- **Source code:** https://github.com/EngreitzLab/DC_TAP_Paper (release tag `v1.0.0`).
- **Raw sequencing data:** pilot DC-TAP-seq, validation DC-TAP-seq, and CTCF ChIP-seq experiments are available from the NCBI Gene Expression Omnibus under accession [`GSE303901`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901).
- **Epigenomic annotation data:** ENCODE Portal (https://www.encodeproject.org); file accessions are listed in Supplementary Table 8 of the manuscript.
- **gRNA design files:** IGVF Data Portal — accessions `IGVFFI7272YLVI` (K562) and `IGVFFI0580WJFK` (WTC11).

## Analysis Summary

Setup the CellRanger outputs for SCEPTRE differential expression analysis
- These steps are summarized in `sceptre_setup.smk` and involve processing the CellRanger outputs, modifying the guide_targets.tsv files and setting up SCEPTRE objects for differential expression

Run differential expression analysis and power simulations using SCEPTRE
- These steps are detailed in `sceptre_power_analysis.smk` where most rules are trying to efficiently run the power simulations based on SCEPTRE's package. The results of the power simulations and differential expression analysis are summarized in `rule format_sceptre_output` and then passed for further formatting

Overlap the tested elements with genomic features and liftOver element coordinates
- In `create_encode_output.smk`, the elements are overlapped with genomic features (promoter, gene body) in order to interpret the differential expression results. The results of K562 and WTC11 are combined at the end of this ruleset.

Format all previous steps, and add different flags based on overlap with genomic features and chromatin features
- In `formatting_dc_tap_results.smk`, the results of the previous analyses are summarized, confidence intervals are added based on a dev branch of SCEPTRE, notes about the original design are added, elements are categorized as "distal" or "promoter" and further categorized based on relationship with the tested gene. Thus a category for each element-gene pair is defined. Finally, chromatin categories are assigned to each element.

A complete description of the analysis methodology and the algorithms implemented in this pipeline is provided in the Methods section of the accompanying manuscript (for example, the power-simulation procedure is implemented in the `sceptre_power_analysis` Snakemake ruleset described there).  

## Important Files

Screen Result Files (Created in `sceptre_setup.smk`)
  - Raw Gene Counts Matrix: `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/dge.rds`
  - Raw Guide Counts Matrix: `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/perturb_status.rds`
  - Metadata File: `results/process_validation_datasets/K562_DC_TAP_Seq/metadata.rds`
  - Guide Design File: `results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv`
  - All Pairs Tested for Differential Expression: `results/process_validation_datasets/K562_DC_TAP_Seq/gene_gRNA_group_pairs.rds`
  - Input for SCEPTRE Differential Expression: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/sceptre_diffex_input.rds`

Differential Expression Output Files (Created in `sceptre_power_analysis.smk`)
  - SCEPTRE Differential Expression Results: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_discovery_analysis.rds`
  - SCEPTRE Calibration Check Results (Negative Control Guide Testing): `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_calibration_check.rds`
  - SCEPTRE Object post-differential expression analysis: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/final_sceptre_object.rds`
  - SCEPTRE Differential Expression Summary Statistics: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/analysis_summary.txt`

Singleton Differential Expression Output Files (Created in `sceptre_power_analysis.smk`)
  - All Results - See (2) for descriptions of each file: `results/process_validation_datasets/K562_DC_TAP_Seq/singleton_differential_expression/*`

Power Analysis Output Files (created in `sceptre_power_analysis.smk`)
  - Batched Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/effect_size_*/`
  - Combined Batched Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/combined_power_analysis_results_es_*.tsv`
  - Analyzed Power Sim Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/power_analysis_results_es_*.tsv`
  - Final Power Sim Output: `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv`

Intermediate Labelling of Elements (Created in `create_encode_output.smk`):
  - All Results: `results/create_encode_output/ENCODE/*`
  - The order of execution for this ruleset is `create_encode_dataset` (hg19), `liftover_enhancers` & `liftover_crispr_dataset`, `filter_crispr_dataset` (hg38), `create_ensemble_encode`, `create_ensemble_epbenchmarking`
  - Essentially everything done in this ruleset is represented in the final file

Final Screen Output Files (Created in `formatting_dc_tap_results.smk`):
  - SCEPTRE Differential Expression Results with Confidence Intervals (See point 2 for file descriptions): `results/formatted_dc_tap_results/differential_expression_w_confidence_intervals_K562_DC_TAP_Seq/*`
  - Final Output w/ EG categories and specific pairs modified: `results/formatted_dc_tap_results/results_wo_pos_controls.tsv`
  - Final Output w/ Added columns for 500bp extension and merging: `results/formatted_dc_tap_results/resized_and_merged_input_for_chromatin_categorization_pipeline.tsv`
  - Final Output w/ Chromatin Categories calculated for resized/merged regions: `results/formatted_dc_tap_results/Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`

Supplementary Tables:
  - Summary of element-gene categories: `results/supplementary_tables/summary_of_element_gene_categories.tsv`

## Instructions for Use

### Reproduce the manuscript's figures and tables


- Default pipeline (`rule all` in `workflow/Snakefile`):
  ```
  snakemake --use-conda --cores <N>
  ```
  This builds the final results table (`results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv`), the supplementary tables, and Figures 1–2. The published Table S3 (`Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv`) is this same table with a small set of annotation columns appended in a manual step downstream of `rule all`.
- Figure 3 and the interpretation-analysis figures:
  ```
  snakemake --use-conda --cores <N> results/main_figure_3/craete_main_figure_3.done
  ```
- Figures 4, 5, 6, S3, and S4a are produced by directly running `interpretation_analysis/dc_tap_analysis.run.R` and `interpretation_analysis/dc_tapseq_indirect_effects.Rmd` — see `interpretation_analysis/README.md`.
- CellRanger outputs are available from NCBI GEO under accession `GSE303901` (see **Data and code availability** below); download them into `resources/process_validation_datasets/{K562,WTC11}_DC_TAP_Seq/cell_ranger_output/` to run `sceptre_setup.smk` from scratch.

## System Requirements

**Operating system**

- Developed and run on Linux x86_64 — CentOS Linux 7 (kernel 3.10, `el7`)

**Software dependencies**

- Snakemake (workflow engine; version not pinned anywhere in the repo — recommend Snakemake ≥7 for compatibility with the `--use-conda` per-rule environment pattern used throughout `workflow/rules/*.smk`)
- conda or mamba (to build the 11 per-rule environments in `workflow/envs/`)
- Python 3.8/3.9
- UCSC `liftOver` binary (via `bioconda::ucsc-liftover`, used in `workflow/rules/create_encode_output.smk`)
- `r-sceptre` / `r-sceptre-dev` — installed from the **non-standard `jgalante` Anaconda.org channel**

**Hardware**

- Per-rule memory requests (`resources: mem=` in `workflow/rules/*.smk`) range from 8 GB to 72 GB; several rules request 32–64 GB.
- Per-rule time budgets range up to 12 hours (e.g. `sceptre_differential_expression`, `sceptre_singleton_differential_expression`)

## Installation Guide

1. Clone the repository and check out the tagged release used for the paper:
   ```
   git clone https://github.com/EngreitzLab/DC_TAP_Paper.git
   cd DC_TAP_Paper
   git checkout v1.0.0
   ```
2. Install a conda/mamba distribution, then install Snakemake into a base environment, e.g.:
   ```
   mamba create -n dc_tap_snakemake -c conda-forge -c bioconda "snakemake>=7"
   conda activate dc_tap_snakemake
   ```
3. Each rule's own conda environment (in `workflow/envs/`) is built automatically the first time it runs, via `--use-conda`
4. Typical install time: building the 11 separate conda environments the first time takes tens of minutes.

## Demo

This repository ships the processed count matrices for the K562 screen, so you can re-run the core SCEPTRE differential-expression step directly from the checked-in data:

- `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/dge.rds` — gene-expression count matrix
- `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/perturb_status.rds` — per-cell guide assignments
- `results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv` — guide-to-target design table

The only external file needed is the GENCODE annotation, which the first command below downloads automatically.

**Instructions to run:**

1. Download the GENCODE annotation used for gene–element pairing:
   ```
   snakemake --use-conda --cores 1 results/genome_annotation_files/gencode.v32lift37.annotation.gtf.gz
   ```
2. Build the SCEPTRE input object from the checked-in counts:
   ```
   snakemake --use-conda --cores 1 results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/sceptre_diffex_input.rds
   ```
   This step also regenerates `gene_gRNA_group_pairs.rds`, `gRNA_groups_table.rds`, and `distances.tsv`

3. Run the SCEPTRE differential-expression analysis:
   ```
   snakemake --use-conda --cores 4 results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_discovery_analysis.rds
   ```

**Expected output:** `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_discovery_analysis.rds` — a SCEPTRE results object with one row per tested element–gene pair (test statistic, p-value, and log2 fold-change)

## License

This project is released under the MIT License; see [`LICENSE`](LICENSE).
