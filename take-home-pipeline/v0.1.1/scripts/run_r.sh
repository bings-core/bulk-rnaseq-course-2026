#!/bin/bash
# Run one R step (scripts/02_deg.R, 03_enrichment.R or 04_report.R) inside the course R container.
#   bash scripts/run_r.sh scripts/02_deg.R
source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
module load singularity/3.6.4
export SINGULARITYENV_R_LIBS_USER=/nonexistent     # keep the container on its own R libraries
singularity exec -B /sc/arion -B /hpc "$R_CONTAINER" Rscript "$1"
