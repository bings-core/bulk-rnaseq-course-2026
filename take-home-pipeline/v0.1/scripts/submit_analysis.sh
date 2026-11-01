#!/bin/bash
# Submit steps 2 to 4 as LSF jobs, each waiting for the one before: DE -> enrichment -> HTML report.
#   bash scripts/submit_analysis.sh        (needs step 1 finished for all samples)
# Checks the inputs first (check_inputs.sh) and records the run (run_info.txt). Without LSF_ACCOUNT the three steps run one after another here.
source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
bash scripts/check_inputs.sh analysis || exit 1
write_run_info
de=$(submit "de_deg" "$THREADS" "$MEM_DE_GB" -- bash scripts/run_r.sh scripts/02_deg.R) || exit 1
en=$(submit "de_enrich" "$THREADS" "$MEM_DE_GB" "done($de)" -- bash scripts/run_r.sh scripts/03_enrichment.R) || exit 1
submit "de_report" 1 "$MEM_REPORT_GB" "done($en)" -- bash scripts/run_r.sh scripts/04_report.R
