#!/bin/bash
# Submit step 1 as LSF jobs: one job per sample (all in parallel), then one MultiQC job that waits for them all.
#   bash scripts/submit_preprocess.sh                 # all samples, from the start
#   bash scripts/submit_preprocess.sh --from star     # restart each sample from one step: that step and everything after it
#                                                     # (fastqc, trim, salmon, star, featurecounts, qualimap); earlier outputs must exist
# To restart only some samples, set SAMPLE_FILTER in config.sh (space-separated sample IDs) and use --from.
# Needs LSF_ACCOUNT (and LSF_QUEUE) in config.sh. Without LSF_ACCOUNT it runs the samples one after another here.
# Check progress with bjobs; each job writes logs/tp_<sample>.<jobid>.out and .err.

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
species_refs || exit 1

FROM=fastqc
while [ $# -gt 0 ]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    *) echo "usage: $0 [--from fastqc|trim|salmon|star|featurecounts|qualimap]" >&2; exit 1 ;;
  esac
done
case "$FROM" in
  fastqc|trim|salmon|star|featurecounts|qualimap) ;;
  *) echo "unknown step '$FROM' for --from (use one of: fastqc trim salmon star featurecounts qualimap)" >&2; exit 1 ;;
esac

deps=""
for sid in $(read_samples | awk '{print $1}'); do
  echo "submit sample $sid (from $FROM)"
  jid=$(submit "tp_$sid" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample "$sid" --from "$FROM") || exit 1
  deps="${deps:+$deps && }done($jid)"
done

if [ -n "$LSF_ACCOUNT" ]; then
  submit "tp_multiqc" 1 "$MEM_MULTIQC_GB" "$deps" -- bash scripts/01_preprocess.sh --multiqc
else
  bash scripts/01_preprocess.sh --multiqc
fi
