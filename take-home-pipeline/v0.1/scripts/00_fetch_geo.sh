#!/bin/bash
# Step 0 (public data only): download the FASTQs of GEO samples.
#   bash scripts/00_fetch_geo.sh GSM8020595,GSM8020596   # comma-separated GSM IDs, one after another (interactive)
#   bash scripts/00_fetch_geo.sh gsm_list.txt            # or one GSM per line in a file
#   bash scripts/00_fetch_geo.sh --draft                 # after downloading: write samples_draft.csv from all records
# For each GSM: find its SRA runs (NCBI E-utilities), get the FASTQ links and MD5s (ENA), download to data/raw/,
# check MD5, and write one record per GSM to data/raw/runs/. Re-running skips files that are already complete.
# submit_fetch.sh runs one LSF job per GSM and then --draft.

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
RAW=data/raw; RUNS_DIR=$RAW/runs
mkdir -p "$RUNS_DIR"

# ---- Section 1: --draft: build samples_draft.csv from the per-GSM records --------------------
if [ "$1" = "--draft" ]; then
  { printf "gsm\tsrr\tlayout\tfile\n"; cat "$RUNS_DIR"/*.tsv 2>/dev/null; } > "$RAW/geo_runs.tsv"
  echo "sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch" > samples_draft.csv
  awk -F'\t' 'NR>1 { g=$1; n[g]++; f[g,n[g]]=$4; lay[g]=$3 }
       END { for (g in n) {
               r1=""; r2=""; paired=(lay[g]=="PAIRED")
               for (i=1;i<=n[g];i++) {
                 if (paired && i%2==0) r2=r2 (r2==""?"":",") f[g,i]
                 else r1=r1 (r1==""?"":",") f[g,i]
               }
               printf "%s,,,%s,%s,auto,\n", g, r1, r2 } }' "$RAW/geo_runs.tsv" | sort >> samples_draft.csv
  echo "Wrote samples_draft.csv. Add condition and replicate for each sample, then save it as samples.csv."
  exit 0
fi

# Lookups. Each prints its answer, or nothing if NCBI/ENA did not answer.
eut_uid()  { curl -s -m 60 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$1&retmode=json" \
             | grep -o '"idlist":\["[0-9]*"' | grep -o '[0-9]*' | head -1; }
eut_srrs() { curl -s -m 60 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=$1" \
             | grep -o 'Run acc="[A-Z]*[0-9]*"' | grep -o 'SRR[0-9]*'; }
ena_row()  { curl -s -m 60 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=library_layout,fastq_ftp,fastq_md5&format=tsv" \
             | tail -1; }

# Run a lookup function with an argument, retrying while it returns nothing (NCBI and ENA rate-limit
# parallel requests). Usage: retry FUNCTION ARG
retry() {
  local out="" i
  for i in 1 2 3 4 5 6 7 8; do
    out=$("$@")
    [ -n "$out" ] && break
    sleep $((i*15))
  done
  printf "%s" "$out"
}

# ---- Section 2: read the GSM list ---------------------------------------------------------
if [ -f "$1" ]; then gsms=$(grep -v '^#' "$1" | tr -s ' \n\t' ' ')
else gsms=$(echo "$1" | tr ',' ' '); fi
[ -n "$gsms" ] || { echo "usage: $0 GSM_LIST_OR_FILE | --draft"; exit 1; }

# ---- Section 3: per GSM, resolve runs and download each FASTQ -----------------------------
for gsm in $gsms; do
  out="$RUNS_DIR/$gsm.tsv"; : > "$out"
  uid=$(retry eut_uid "$gsm")
  if [ -z "$uid" ]; then echo "WARNING: no SRA record for $gsm" >&2; continue; fi
  srrs=$(retry eut_srrs "$uid")
  echo "$gsm -> $(echo $srrs | tr '\n' ' ')"

  for srr in $srrs; do
    row=$(retry ena_row "$srr")
    layout=$(echo "$row" | cut -f2)
    ftps=$(echo "$row" | cut -f3 | tr ';' ' ')
    md5s=$(echo "$row" | cut -f4 | tr ';' ' ')
    if [ -z "$ftps" ]; then echo "WARNING: no ENA FASTQ for $srr ($gsm)" >&2; continue; fi
    set -- $md5s
    i=1
    for ftp in $ftps; do
      f="$RAW/$(basename "$ftp")"; want_md5=$(eval echo \${$i}); i=$((i+1))
      if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want_md5" ]; then
        echo "  have $f"
      else
        # resume until the MD5 matches (curl -C - continues a partial file)
        for attempt in $(seq 1 200); do
          curl -s -m 0 -C - --retry 3 --retry-delay 20 -o "$f" "https://$ftp"
          [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want_md5" ] && break
          echo "  attempt $attempt: $(stat -c %s "$f") bytes so far, retrying" >&2; sleep 15
        done
        [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want_md5" ] || { echo "MD5 mismatch after retries: $f" >&2; continue; }
      fi
      printf "%s\t%s\t%s\t%s\n" "$gsm" "$srr" "$layout" "$f" >> "$out"
    done
  done
done
echo "Done. Next: bash scripts/00_fetch_geo.sh --draft"
