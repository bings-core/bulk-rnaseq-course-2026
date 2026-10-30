#!/bin/bash
# Step 1: preprocessing, as in Session 4. For one sample (or all, when run interactively):
#   FastQC (raw) -> TrimGalore -> Salmon (library type detected) -> STAR -> samtools -> featureCounts -> Qualimap
# and, with --multiqc, one MultiQC report over all samples (the same patched 1.35 container and config as the course pipeline).
#   bash scripts/01_preprocess.sh                 # every sample, one after another (interactive)
#   bash scripts/01_preprocess.sh --sample KO_1   # one sample (used by submit_preprocess.sh, one LSF job each)
#   bash scripts/01_preprocess.sh --multiqc       # the MultiQC report (run after all samples are done)
#   bash scripts/01_preprocess.sh --sample KO_1 --qualimap   # only the Qualimap step for one sample (after STAR)
#   bash scripts/01_preprocess.sh --from star     # restart from a step: that step and everything after it (fastqc, trim, salmon, star, featurecounts, qualimap)
# Works for single-end (fastq_2 empty) and paired-end samples.
# With --from, the earlier steps are not run; their outputs in results/ must already exist for the samples.

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1
MODE=samples
STEPS=(fastqc trim salmon star featurecounts qualimap)
FROM=fastqc
while [ $# -gt 0 ]; do
  case "$1" in
    --sample)  SAMPLE_FILTER="$2"; shift 2 ;;
    --from)    FROM="$2"; shift 2 ;;
    --multiqc) MODE=multiqc; shift ;;
    --qualimap) MODE=qualimap; shift ;;
    *) echo "usage: $0 [--sample SAMPLE_ID | --multiqc | --from STEP]" >&2; exit 1 ;;
  esac
done
START=-1
for i in "${!STEPS[@]}"; do [ "${STEPS[$i]}" = "$FROM" ] && START=$i; done
if [ "$START" -lt 0 ]; then echo "unknown step '$FROM' for --from (use one of: ${STEPS[*]})" >&2; exit 1; fi
# A step runs when it is the --from step or comes after it.
run_from() {
  local i
  for i in "${!STEPS[@]}"; do [ "${STEPS[$i]}" = "$1" ] && [ "$i" -ge "$START" ] && return 0; done
  return 1
}

# ---- Section 1: tools and reference ----------------------------------------------------
module load fastqc/0.11.8 trim_galore/0.6.6 star/2.7.5b salmon/1.2.1 samtools/1.11 subread/2.0.1 singularity/3.6.4
species_refs || exit 1
P="$RESULTS_DIR/preprocessed"                  # folder names match the course MultiQC config
mkdir -p "$P"/{fastqc_raw,trim_galore_trimmed,salmon,star,samtools,featurecounts,qualimap,strandedness}

# Qualimap for one sample (BAM from STAR), in the Qualimap container (see setup_resources.sh).
run_qualimap() {
  local sid="$1" layout="$2" bam="$3" pe=""
  [ "$layout" = PE ] && pe="-pe"
  run singularity exec -B /sc/arion -B /hpc "$QUALIMAP_IMAGE" qualimap rnaseq \
      --java-mem-size="${MEM_QUALIMAP_GB:-8}G" -bam "$bam" -gtf "$GTF" -outdir "$P/qualimap/$sid" $pe
}

# ---- Section 2: per-sample work ----------------------------------------------------------
if [ "$MODE" = samples ]; then
[ -n "$SAMPLE_FILTER" ] && set -e      # a sample job stops at the first failed step
while read -r sid r1 r2 strand; do
  if [ "$r2" = "-" ]; then layout=SE; else layout=PE; fi
  echo "== $sid ($layout), from $FROM"

  # 2.0 Give the input reads the sample_id as their name, so FastQC/TrimGalore files (and MultiQC) show it.
  # Several lanes are concatenated into one file per read (gzip concatenation is valid).
  mkdir -p "$RESULTS_DIR/input"
  link_reads() {   # link_reads SRC_LIST DEST  (comma-separated list = lanes)
    local dest="$2"
    if [[ "$1" == *,* ]]; then cat ${1//,/ } > "$dest"; else ln -sf "$(readlink -f "$1")" "$dest"; fi
  }
  if [ "$layout" = SE ]; then
    run_from trim && link_reads "$r1" "$RESULTS_DIR/input/${sid}.fastq.gz"; r1="$RESULTS_DIR/input/${sid}.fastq.gz"
  else
    run_from trim && link_reads "$r1" "$RESULTS_DIR/input/${sid}_R1.fastq.gz"; r1="$RESULTS_DIR/input/${sid}_R1.fastq.gz"
    run_from trim && link_reads "$r2" "$RESULTS_DIR/input/${sid}_R2.fastq.gz"; r2="$RESULTS_DIR/input/${sid}_R2.fastq.gz"
  fi

  # 2a. FastQC on the raw reads (comma-separated lanes become a list)
  if run_from fastqc; then
    run fastqc -t "$THREADS" -o "$P/fastqc_raw" "$r1" ${r2:+"$r2"} 2>/dev/null
  fi

  # 2b. TrimGalore (PE uses --paired); FastQC of the trimmed reads goes next to them
  if run_from trim; then
    if [ "$layout" = SE ]; then
      run trim_galore --fastqc --cores "$THREADS" --gzip -o "$P/trim_galore_trimmed" "$r1"
    else
      run trim_galore --paired --fastqc --cores "$THREADS" --gzip -o "$P/trim_galore_trimmed" "$r1" "$r2"
    fi
  fi
  if [ "$layout" = SE ]; then
    t1=$(ls "$P/trim_galore_trimmed/${sid}"*trimmed.fq.gz | head -1); t2=""
  else
    t1=$(ls "$P/trim_galore_trimmed/${sid}_R1"*val_1.fq.gz | head -1)
    t2=$(ls "$P/trim_galore_trimmed/${sid}_R2"*val_2.fq.gz | head -1)
  fi

  # 2c. Salmon with automatic library-type detection (-l A)
  if run_from salmon; then
    if [ "$layout" = SE ]; then reads=(-r "$t1"); else reads=(-1 "$t1" -2 "$t2"); fi
    run salmon quant -i "$SALMON_INDEX" -l A "${reads[@]}" -p "$THREADS" \
        --validateMappings --seqBias --gcBias -o "$P/salmon/$sid"
  fi

  # 2d. Strandedness: Salmon's inferred library type -> unstranded / forward / reverse, recorded per sample
  if [ ! -f "$P/salmon/$sid/lib_format_counts.json" ]; then
    echo "ERROR: no Salmon output for $sid. Run from salmon or earlier." >&2; exit 1
  fi
  inferred=$(grep -o '"expected_format": *"[^"]*"' "$P/salmon/$sid/lib_format_counts.json" | sed 's/.*"\([^"]*\)"$/\1/')
  case "$inferred" in
    IU|U)   st=unstranded; fc=0 ;;
    ISR|SR) st=reverse;    fc=2 ;;
    ISF|SF) st=forward;    fc=1 ;;
    *)      st=undetermined; fc=2 ;;
  esac
  if [ "$strand" != auto ]; then st="$strand"; case "$strand" in unstranded) fc=0;; forward) fc=1;; reverse) fc=2;; esac; fi
  if run_from salmon; then   # one line per sample: replace an earlier line from a previous run
    sed -i "/^${sid}\t/d" "$P/strandedness/strandedness_report.tsv" 2>/dev/null
    printf "%s\t%s\t%s\n" "$sid" "$inferred" "$st" >> "$P/strandedness/strandedness_report.tsv"
  fi

  # 2e. STAR alignment (sorted BAM and gene counts), then 2f. samtools
  if [ "$layout" = SE ]; then fq="$t1"; else fq="$t1 $t2"; fi
  bam="$P/star/${sid}_Aligned.sortedByCoord.out.bam"
  if run_from star; then
    run STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX" --readFilesIn $fq \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --quantMode GeneCounts --outFileNamePrefix "$P/star/${sid}_"
    # samtools index, flagstat, stats, idxstats (MultiQC reads these)
    run samtools index "$bam"
    run samtools flagstat "$bam" > "$P/samtools/${sid}.flagstat"
    run samtools stats "$bam" > "$P/samtools/${sid}.stats"
    run samtools idxstats "$bam" > "$P/samtools/${sid}.idxstats"
  fi

  # 2g. featureCounts on exons (strand from 2d; -p for paired-end)
  if run_from featurecounts && [ "$RUN_FEATURECOUNTS" = yes ]; then
    pe=""; [ "$layout" = PE ] && pe="-p"
    run featureCounts -T "$THREADS" $pe -s "$fc" -t exon -g gene_id -a "$GTF" \
        -o "$P/featurecounts/${sid}.featureCounts.txt" "$bam"
  fi

  # 2h. Qualimap rnaseq (optional)
  if run_from qualimap && [ "$RUN_QUALIMAP" = yes ]; then run_qualimap "$sid" "$layout" "$bam"; fi
done < <(read_samples)

fi

# Qualimap only: reuse the BAMs already made by STAR
if [ "$MODE" = qualimap ]; then
set -e
while read -r sid r1 r2 strand; do
  if [ "$r2" = "-" ]; then layout=SE; else layout=PE; fi
  [ "$RUN_QUALIMAP" = yes ] && run_qualimap "$sid" "$layout" "$P/star/${sid}_Aligned.sortedByCoord.out.bam"
done < <(read_samples)
fi

# ---- Section 3: MultiQC, the same container and config as the course pipeline --------------
if [ "$MODE" = multiqc ]; then
set -e
module load singularity/3.6.4
OUT="$P/multiqc_report"; mkdir -p "$OUT/.home"          # folder name starts with multiqc_ so MultiQC skips its own output
run singularity exec --cleanenv --bind /sc/arion,/hpc \
    --env HOME="$OUT/.home" --env TMPDIR="$OUT/.home" --env MPLCONFIGDIR="$OUT/.home" --env XDG_CACHE_HOME="$OUT/.home" \
    "$MULTIQC_CONTAINER" multiqc --no-version-check -f --config "$MULTIQC_CONFIG" \
    --filename multiqc_report.html -o "$OUT" "$P"
echo "Done. Open $OUT/multiqc_report.html"
fi
