# 2026-05-11 Guide Pool QC

## Purpose

This analysis note documents the guide-pool QC workflow used to align MiSeq reads to the guide reference and summarize alignment quality.

## Inputs

FASTQ inputs:

- `gs://<guide-pool-sequencing-bucket>/WALKUP-19964/fastqs/pool1_S1_L001_R1_001.fastq.gz`
- `gs://<guide-pool-sequencing-bucket>/WALKUP-19964/fastqs/pool2_S2_L001_R1_001.fastq.gz`

Reference input:

- `data/guides.fa`

## Code

QC report code:

- `analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd`

Command sequence used for preprocessing and alignment:

```bash
# cutadapt v4.2
cutadapt -a TTTAAX -o ../../../bowtie_pool_control_miseq/pool2_S2_L001_R1_001.trimmed.fastq.gz pool2_S2_L001_R1_001.fastq.gz
cutadapt -a TTTAAX -o ../../../bowtie_pool_control_miseq/pool1_S1_L001_R1_001.trimmed.fastq.gz pool1_S1_L001_R1_001.fastq.gz

# bowtie v2.5.0
bowtie-build guides.fa guides_index
bowtie -v 1 -a --best --strata -x guides_index pool1_S1_L001_R1_001.trimmed.fastq.gz pool1_bowtie_out.txt
bowtie -v 1 -a --best --strata -x guides_index pool2_S2_L001_R1_001.trimmed.fastq.gz pool2_bowtie_out.txt
```

The R Markdown uses `data.table`, `ggplot2`, `stringr`, and `here` to summarize Bowtie output and generate QC plots.

Render the PDF report from the repository root with:

```bash
RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd', output_format = 'pdf_document', output_file = here::here('results', 'guide_pool_qc', 'guide_pool_qc.pdf'))"
```

## Outputs

Primary alignment outputs:

- `results/guide_pool_qc/pool1_bowtie_out.txt`
- `results/guide_pool_qc/pool2_bowtie_out.txt`

Report outputs:

- `results/guide_pool_qc/guide_pool_qc.pdf`
- `results/guide_pool_qc/guide_counts.tsv`

Per-pool QC directories:

- `results/guide_pool_qc/pool1/`
- `results/guide_pool_qc/pool2/`

## Tool versions and settings

- cutadapt: `4.2`
- Bowtie: `2.5.0`
- Alignment flags: `-v 1 -a --best --strata`
- R environment: packages resolved via `renv.lock`