#!/bin/bash
# Submit step 0 (GEO downloads) as LSF jobs: one job per GSM (in parallel), then one job that writes samples_draft.csv.
#   bash scripts/submit_fetch.sh GSM8020595,GSM8020596      (or a file with one GSM per line)
# Needs LSF_ACCOUNT (and LSF_QUEUE) in config.sh. Without it, the downloads run one after another here.

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
if [ -f "$1" ]; then gsms=$(grep -v '^#' "$1" | tr -s ' \n\t' ' ')
else gsms=$(echo "$1" | tr ',' ' '); fi
[ -n "$gsms" ] || { echo "usage: $0 GSM_LIST_OR_FILE"; exit 1; }

deps=""
for gsm in $gsms; do
  echo "submit download $gsm"
  jid=$(submit "fetch_$gsm" 1 "$MEM_FETCH_GB" -- bash scripts/00_fetch_geo.sh "$gsm") || exit 1
  deps="${deps:+$deps && }done($jid)"
done
if [ -n "$LSF_ACCOUNT" ]; then
  submit "fetch_draft" 1 "$MEM_FETCH_GB" "$deps" -- bash scripts/00_fetch_geo.sh --draft
else
  bash scripts/00_fetch_geo.sh --draft
fi
