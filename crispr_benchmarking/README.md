# CRISPR benchmarking resources
This directory contains resources used to benchmark E2G models on the DC-TAP-seq CRISPR data:

- `merge_dc_tap_elements.R`: Resize and merge DC-TAP-seq elements for benchmarking
- `dc_tap_seq_paper.yml`: CRISPR benchmarking pipeline config file specifying input files for
comparisons to run
- `benchmarking_pred_config_with_alpha.tsv`: predictor config file specifying parameters for E2G
predictors

Used E2G prediction files on ENCODE portal:

- ENCODE-rE2G (K562): [ENCFF970QAX](https://www.encodeproject.org/files/ENCFF970QAX/)
- ENCODE-rE2G (WTC11): [ENCFF071ZZU](https://www.encodeproject.org/files/ENCFF071ZZU/)
- ABC (K562): [ENCFF681DDZ](https://www.encodeproject.org/files/ENCFF681DDZ/)
- ABC (WTC11): [ENCFF660VQW](https://www.encodeproject.org/files/ENCFF660VQW/)
