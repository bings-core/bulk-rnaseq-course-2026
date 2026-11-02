#!/bin/bash
# Check the inputs before anything is submitted, so a mistake shows up now and not hours into a run.
#   bash scripts/check_inputs.sh preprocess   # run by submit_preprocess.sh: samples, FASTQs, references, containers
#   bash scripts/check_inputs.sh analysis     # run by submit_analysis.sh: also contrasts, and Salmon outputs for every sample
# Prints every error and warning. Exits 1 on any error, so the submit script stops before submitting anything.
# Warnings (for example fewer than 2 replicates in a condition) do not stop the run.

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
MODE="${1:-preprocess}"
errors=0
err()  { echo "ERROR: $*" >&2; errors=$((errors + 1)); }
warn() { echo "WARNING: $*" >&2; }

# ---- Settings -----------------------------------------------------------------------------------
case "$SPECIES" in human|mouse) ;; *) err "SPECIES must be human or mouse (config.sh has '$SPECIES')" ;; esac
if [ -z "$LSF_ACCOUNT" ]; then
  warn "LSF_ACCOUNT is empty: every step runs here on the login node, one after another. Set it for real runs."
fi
if [ "$MODE" = analysis ]; then
  [ -f "$R_CONTAINER" ] || err "R container not found: $R_CONTAINER (run scripts/setup_resources.sh)"
  TXDB="$RESOURCE_DIR/$SPECIES/txdb/$([ "$SPECIES" = human ] && echo homo_sapiens_grch38_gencode_36.sqlite || echo mus_musculus_grcm38_gencode_M25.sqlite)"
  [ -s "$TXDB" ] || err "TxDb not found: $TXDB (run scripts/setup_resources.sh and wait for setup_txdb)"
else
  [ "$RUN_QUALIMAP" = yes ] && { [ -f "$QUALIMAP_IMAGE" ] || err "Qualimap container not found: $QUALIMAP_IMAGE (run scripts/setup_resources.sh)"; }
  [ -f "$MULTIQC_CONTAINER" ] || err "MultiQC container not found: $MULTIQC_CONTAINER (run scripts/setup_resources.sh)"
  # The index folders appear while they are being built, so check that each build finished
  REF="$RESOURCE_DIR/$SPECIES"
  grep -q "DONE: Genome generation" "$REF/star/2.7.5b/index/Log.out" 2>/dev/null \
    || err "STAR index for $SPECIES is not finished ($REF/star/2.7.5b/index/Log.out has no 'DONE: Genome generation'). Wait for setup_star to finish"
  [ -s "$REF/salmon/1.2.1/index/info.json" ] \
    || err "Salmon index for $SPECIES is missing or unfinished ($REF/salmon/1.2.1/index/info.json). Wait for setup_salmon to finish"
  [ -s "$REF/gtf/annotation.gtf" ] || err "GTF missing: $REF/gtf/annotation.gtf (run scripts/setup_resources.sh)"
  species_refs 2>/dev/null || err "references for $SPECIES are missing or incomplete in $RESOURCE_DIR (run scripts/setup_resources.sh)"
fi

# ---- Sample sheet -------------------------------------------------------------------------------
if [ ! -f "$SAMPLES_CSV" ]; then
  err "samples sheet not found: $SAMPLES_CSV (copy samples.example.csv and edit it)"
  echo "$errors error(s). Nothing was submitted." >&2; exit 1
fi
want_header="sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch"
if [ "$(head -1 "$SAMPLES_CSV" | tr -d '\r')" != "$want_header" ]; then
  err "samples.csv header must be exactly: $want_header"
fi
while IFS= read -r line; do
  case "${line:0:1}" in
    E) err "${line:2}" ;;
    W) warn "${line:2}" ;;
  esac
done < <(awk -F, '
  NR == 1 { next }
  NF == 0 { next }
  {
    n++; r = "row " NR
    if ($1 == "") { print "E|" r ": sample_id is empty"; next }
    if ($1 ~ /[ \t]/) print "E|" r ": sample_id \"" $1 "\" contains spaces"
    if ($1 in seen) print "E|" r ": duplicate sample_id " $1
    seen[$1] = 1
    s = " (" $1 ")"
    if ($2 == "") print "E|" r s ": condition is empty"; else cnt[$2]++
    if ($3 == "") print "E|" r s ": replicate is empty"
    if ($4 == "") { print "E|" r s ": fastq_1 is empty"; next }
    k1 = split($4, f1, ";")
    for (i = 1; i <= k1; i++) if (system("test -s \"" f1[i] "\"") != 0) print "E|" r s ": fastq_1 not found or empty: " f1[i]
    if ($5 != "") {
      k2 = split($5, f2, ";")
      if (k2 != k1) print "E|" r s ": fastq_1 has " k1 " file(s) but fastq_2 has " k2 " (lanes must match)"
      for (i = 1; i <= k2; i++) if (system("test -s \"" f2[i] "\"") != 0) print "E|" r s ": fastq_2 not found or empty: " f2[i]
    }
    if ($6 != "" && $6 !~ /^(auto|unstranded|forward|reverse)$/) print "E|" r s ": strandedness must be auto, unstranded, forward or reverse (got " $6 ")"
    if ($7 != "") nb++; else ne++
  }
  END {
    if (nb > 0 && ne > 0) print "E|batch is filled for some samples only (" nb " filled, " ne " empty). Fill it for all or none"
    for (c in cnt) if (cnt[c] < 2) print "W|condition " c " has " cnt[c] " replicate(s); DESeq2 and the DE tables need at least 2 per condition"
    if (n == 0) print "E|samples.csv has no samples"
  }' "$SAMPLES_CSV")

# SAMPLE_FILTER names must exist
for id in $SAMPLE_FILTER; do
  awk -F, -v id="$id" 'NR > 1 && $1 == id { found = 1 } END { exit !found }' "$SAMPLES_CSV" \
    || err "SAMPLE_FILTER names '$id', which is not in samples.csv"
done

# ---- Contrast sheet (analysis) ------------------------------------------------------------------
if [ "$MODE" = analysis ]; then
  if [ ! -f "$CONTRASTS_CSV" ]; then
    err "contrasts sheet not found: $CONTRASTS_CSV (copy contrasts.example.csv and edit it)"
  else
    [ "$(head -1 "$CONTRASTS_CSV" | tr -d '\r')" = "name,treatment,reference" ] || err "contrasts.csv header must be exactly: name,treatment,reference"
    conds="|$(awk -F, 'NR > 1 && $2 != "" { print $2 }' "$SAMPLES_CSV" | sort -u | tr '\n' '|')"
    while IFS= read -r line; do
      case "${line:0:1}" in E) err "${line:2}" ;; esac
    done < <(awk -F, -v conds="$conds" '
      NR == 1 { next }
      NF == 0 { next }
      {
        r = "contrast row " NR
        if ($1 == "" || $1 ~ /[ \t]/) print "E|" r ": name is empty or contains spaces"
        if ($1 in seen) print "E|" r ": duplicate name " $1
        seen[$1] = 1
        if (index(conds, "|" $2 "|") == 0) print "E|" r " (" $1 "): treatment \"" $2 "\" is not a condition in samples.csv"
        if (index(conds, "|" $3 "|") == 0) print "E|" r " (" $1 "): reference \"" $3 "\" is not a condition in samples.csv"
        if ($2 == $3) print "E|" r " (" $1 "): treatment and reference are the same"
      }' "$CONTRASTS_CSV")
  fi
fi

# ---- Salmon outputs (analysis) --------------------------------------------------------------------
if [ "$MODE" = analysis ] && [ "$errors" -eq 0 ]; then   # only when the sample sheet is valid
  for sid in $(read_samples | awk '{print $1}'); do
    [ -s "$RESULTS_DIR/preprocessed/salmon/$sid/quant.sf" ] || err "no Salmon output for $sid. Finish step 1 (submit_preprocess.sh) first"
  done
fi

if [ "$errors" -gt 0 ]; then
  echo "$errors error(s). Nothing was submitted. Fix them and run again." >&2
  exit 1
fi
echo "Inputs OK for $MODE: $(read_samples | wc -l) sample(s)."
