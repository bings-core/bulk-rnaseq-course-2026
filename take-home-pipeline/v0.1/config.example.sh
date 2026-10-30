# Copy this file to config.sh and edit it. Every script reads config.sh; nothing else needs editing.
# Plain KEY=VALUE: bash sources it, R reads it with Sys.getenv().
# Every path is relative to the project folder, so run the scripts from there: cd into the folder that holds this file.

# Organism. Choose human or mouse.
SPECIES=human

# References, built by scripts/setup_resources.sh (run it once).
RESOURCE_DIR=resources

# Where results go.
RESULTS_DIR=results

# Your sample sheet and contrast sheet.
SAMPLES_CSV=samples.csv
CONTRASTS_CSV=contrasts.csv

# Threads per tool (STAR and Salmon use this many cores).
THREADS=8

# Cluster account for LSF jobs. Course: acc_BiNGS_bulk. After the course, your own acc_* account.
# Leave empty to run everything interactively (no jobs submitted).
LSF_ACCOUNT=
# Queue for the jobs (course: premium). Change after the course if needed.
LSF_QUEUE=premium
# Optional walltime, e.g. 12:00 (hh:mm). Empty = the queue default.
LSF_WALLTIME=12:00

# Memory per job in GB, by step (human STAR needs about 32 GB; mouse similar).
MEM_FETCH_GB=2
MEM_FASTQC_GB=4
MEM_TRIM_GB=4
MEM_SALMON_GB=8
MEM_STAR_GB=40
MEM_SAMTOOLS_GB=4
MEM_FEATURECOUNTS_GB=8
MEM_QUALIMAP_GB=8
MEM_MULTIQC_GB=4
MEM_DE_GB=16          # DESeq2 and enrichment
MEM_REPORT_GB=8      # HTML report

# Containers, downloaded or built into containers/ (see the README).
R_CONTAINER=containers/rbings_20250925.sif
QUALIMAP_IMAGE=containers/qualimap_2.2.2d.sif

# Optional: run only these sample_ids (space separated). Empty = all samples in samples.csv.
SAMPLE_FILTER=

# MSigDB gene sets for step 3 (comma-separated). H = Hallmark. C2:CP:KEGG = KEGG. C5:GO:BP = GO biological process.
MSIGDB_COLLECTIONS=H

# Project name shown in the final report.
PROJECT_NAME=my_rnaseq_project

# MultiQC 1.35 (Galaxy depot image, downloaded by setup_resources.sh) and the course MultiQC config.
MULTIQC_CONTAINER=containers/multiqc_1.35.sif
MULTIQC_CONFIG=config/multiqc_config.yaml

# Optional QC modules (Session 4 runs both). yes or no.
RUN_FEATURECOUNTS=yes
RUN_QUALIMAP=yes
