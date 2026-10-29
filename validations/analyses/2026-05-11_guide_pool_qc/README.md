# 2026-05-11 Guide Pool QC

## Purpose

Checks the uniformity of the sgRNA plasmid pools of the validation library. MiSeq reads from the two pools are aligned to the guide sequences, and reads per guide and mismatches per read are summarized.

## Inputs

- MiSeq FASTQs of the two plasmid pools (not public).
- Guide reference: `metadata/validation_library_guides.fa.gz` (582 guides: 482 targeting the 30 elements and 100 non-targeting).
- `source_data/`: per-guide read counts (`guide_counts.tsv`) and per-pool mismatch summaries (`mismatch_counts.tsv`, `mismatch_positions.tsv`). These were computed from the Bowtie alignments and are enough to regenerate the report.

## Code

Trimming and alignment (run outside this repository):

```bash
# cutadapt 4.2
cutadapt -a TTTAAX -o pool1_S1_L001_R1_001.trimmed.fastq.gz pool1_S1_L001_R1_001.fastq.gz
cutadapt -a TTTAAX -o pool2_S2_L001_R1_001.trimmed.fastq.gz pool2_S2_L001_R1_001.fastq.gz

# Bowtie 1
gunzip -c metadata/validation_library_guides.fa.gz > guides.fa
bowtie-build guides.fa guides_index
bowtie -v 1 -a --best --strata -x guides_index pool1_S1_L001_R1_001.trimmed.fastq.gz pool1_bowtie_out.txt
bowtie -v 1 -a --best --strata -x guides_index pool2_S2_L001_R1_001.trimmed.fastq.gz pool2_bowtie_out.txt
```

QC report: `analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd`. If `results/guide_pool_qc/pool{1,2}_bowtie_out.txt` exist, the report is computed from the alignments; otherwise from `source_data/`.

Run from `validations/`:

```bash
Rscript -e "rmarkdown::render('analyses/2026-05-11_guide_pool_qc/pool_QC.Rmd', output_dir = 'results/guide_pool_qc', output_file = 'guide_pool_qc.pdf')"
```

## Outputs

- `results/guide_pool_qc/guide_pool_qc.pdf`
- `results/guide_pool_qc/guide_counts.tsv`
