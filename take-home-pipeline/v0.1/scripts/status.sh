#!/bin/bash
# Show how far each sample got in step 1, the --from step to restart it from, and any failed LSF jobs.
#   bash scripts/status.sh
# A step counts as done when its output file exists. Restart with the step it names:
#   bash scripts/submit_preprocess.sh --from <step>      (set SAMPLE_FILTER in config.sh to restart only some samples)

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
P="$RESULTS_DIR/preprocessed"
STEPS=(fastqc trim salmon star featurecounts qualimap)

# Output file that shows a step finished, for one sample and layout (SE or PE)
step_file() {
  local step="$1" sid="$2" layout="$3"
  case "$step" in
    fastqc)        [ "$layout" = SE ] && echo "$P/fastqc_raw/${sid}_fastqc.zip" || echo "$P/fastqc_raw/${sid}_R1_fastqc.zip" ;;
    trim)          [ "$layout" = SE ] && echo "$P/trim_galore_trimmed/${sid}_trimmed.fq.gz" || echo "$P/trim_galore_trimmed/${sid}_R1_val_1.fq.gz" ;;
    salmon)        echo "$P/salmon/$sid/quant.sf" ;;
    star)          echo "$P/star/${sid}_Aligned.sortedByCoord.out.bam" ;;
    featurecounts) echo "$P/featurecounts/${sid}.featureCounts.txt" ;;
    qualimap)      echo "$P/qualimap/$sid/rnaseq_qc_results.txt" ;;
  esac
}

printf "%-22s %-6s" "sample" "layout"
for s in "${STEPS[@]}"; do printf " %-13s" "$s"; done
printf " %s\n" "restart from"
incomplete=0
while read -r sid r1 r2 strand; do
  if [ "$r2" = "-" ]; then layout=SE; else layout=PE; fi
  printf "%-22s %-6s" "$sid" "$layout"
  restart="complete"
  for s in "${STEPS[@]}"; do
    f=$(step_file "$s" "$sid" "$layout")
    if [ -s "$f" ] || { [ "$s" = featurecounts ] && [ "$RUN_FEATURECOUNTS" != yes ]; } || { [ "$s" = qualimap ] && [ "$RUN_QUALIMAP" != yes ]; }; then
      printf " %-13s" "done"
    else
      printf " %-13s" "missing"
      [ "$restart" = complete ] && restart="$s"
    fi
  done
  printf " %s\n" "$restart"
  [ "$restart" = complete ] || incomplete=$((incomplete + 1))
done < <(read_samples)
echo
echo "$incomplete sample(s) not complete."

# LSF jobs from submit_preprocess.sh
if [ -n "$LSF_ACCOUNT" ] && command -v bjobs >/dev/null 2>&1; then
  jobs=$(bjobs -a -w 2>/dev/null | awk '$7 ~ /^tp_/ { print $1, $3, $7 }')
  failed=$(echo "$jobs" | awk '$2 == "EXIT" { print }')
  waiting=$(echo "$jobs" | awk '$2 == "PEND" && $3 == "tp_multiqc" { print $1 }')
  if [ -n "$failed" ]; then
    echo "Failed LSF jobs (job id, status, name):"
    echo "$failed" | sed 's/^/  /'
    echo "Logs: logs/<name>.<jobid>.err"
  fi
  if [ -n "$failed" ] && [ -n "$waiting" ]; then
    echo "The MultiQC job ($waiting) waits for the failed jobs and will never start. Remove it with: bkill $waiting"
  fi
  running=$(echo "$jobs" | awk '$2 == "RUN" || $2 == "PEND"' | wc -l)
  echo "Sample and MultiQC jobs still pending or running: $running"
fi
