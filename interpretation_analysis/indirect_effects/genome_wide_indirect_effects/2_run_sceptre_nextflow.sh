#!/bin/bash

## Run SCEPTRE Nextflow pipeline to perform all-against-all trans analysis for the Gasperini et al.,
## 2019 dataset. Requires that Nextflow and R are available in the PATH, and the sceptre and ondisc
## R-packages are installed.

# load required modules
ml R/4.2.0
ml java/21.0.4
ml biology
ml nextflow/25.04.7

##########################
# REQUIRED INPUT ARGUMENTS
##########################
data_directory="./sceptre_object_ondisc"
# sceptre object
sceptre_object_fp=$data_directory"/sceptre_object.rds"
# response ODM
response_odm_fp=$data_directory"/response.odm"
# grna ODM
grna_odm_fp=$data_directory"/grna.odm"
# object containing model formula
formula_object="./sceptre_formula.rds"

###################
# OUTPUT DIRECTORY:
##################
output_directory="./sceptre_outputs"

#################
# Invoke pipeline
#################
nextflow run timothy-barry/sceptre-pipeline -r main \
--sceptre_object_fp $sceptre_object_fp \
--response_odm_fp $response_odm_fp \
--grna_odm_fp $grna_odm_fp \
--output_directory $output_directory \
--formula_object $formula_object \
--grna_assignment_method "thresholding" \
--threshold 1 \
--pair_pod_size 250000 \
--run_association_analysis_time_per_pair 0.25s \
--run_association_analysis_memory 32G \
--discovery_pairs trans \
-resume
