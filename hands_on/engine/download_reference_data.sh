#!/usr/bin/env bash
# =============================================================================
# download_reference_data.sh — restore the large binary blobs the .gitignore
# holds back from GitHub. On Minerva (which is where the course runs) this is
# just a set of symlinks to the shared reference tree; nothing is actually
# downloaded.
#
# USAGE
#   cd /sc/arion/projects/BiNGS_bulk/$USER/hands_on   # your own clone
#   bash engine/download_reference_data.sh
# =============================================================================
set -euo pipefail

SHARED=/sc/arion/projects/BiNGS/bings_omics/data/bings/2026/bulk_course_data/hands_on
STUDENT="$(pwd)"

if [[ ! -d "$SHARED" ]]; then
  echo "ERROR: shared reference not found at $SHARED"
  echo "This script assumes Minerva. Contact the instructor if you are elsewhere."
  exit 1
fi

if [[ "$STUDENT" -ef "$SHARED" ]]; then
  echo "You are inside the shared reference tree already — nothing to link."
  exit 0
fi

echo "▶ Linking large reference blobs from  $SHARED"
echo "                                into  $STUDENT"
echo

# ─── 1. Singularity container (~7 GB) ────────────────────────────────────────
mkdir -p "$STUDENT/engine/singularity_containers"
for sif in "$SHARED"/engine/singularity_containers/*.sif; do
  [[ -e "$sif" ]] || continue
  ln -sf "$sif" "$STUDENT/engine/singularity_containers/$(basename "$sif")"
  echo "  link  engine/singularity_containers/$(basename "$sif")"
done

# ─── 2. STAR + Salmon indices (~2.8 GB combined) ─────────────────────────────
mkdir -p "$STUDENT/engine/annotation"
for idx in star salmon; do
  if [[ -d "$SHARED/engine/annotation/$idx" ]]; then
    ln -sfn "$SHARED/engine/annotation/$idx" "$STUDENT/engine/annotation/$idx"
    echo "  link  engine/annotation/$idx/"
  fi
done

# ─── 3. chr5 genome FASTA (177 MB) ───────────────────────────────────────────
mkdir -p "$STUDENT/engine/annotation/data"
if [[ -f "$SHARED/engine/annotation/data/GRCh38.p13.chr5.genome.fa" ]]; then
  ln -sf "$SHARED/engine/annotation/data/GRCh38.p13.chr5.genome.fa" \
         "$STUDENT/engine/annotation/data/GRCh38.p13.chr5.genome.fa"
  echo "  link  engine/annotation/data/GRCh38.p13.chr5.genome.fa"
fi

# ─── 4. Raw chr5 fastqs (1.1 GB total) ───────────────────────────────────────
mkdir -p "$STUDENT/data_rna/reduced_chr5/raw/fastq"
for fq in "$SHARED"/data_rna/reduced_chr5/raw/fastq/*.fastq*; do
  [[ -e "$fq" ]] || continue
  ln -sf "$fq" "$STUDENT/data_rna/reduced_chr5/raw/fastq/$(basename "$fq")"
done
echo "  link  data_rna/reduced_chr5/raw/fastq/*.fastq  ($(ls "$SHARED"/data_rna/reduced_chr5/raw/fastq/*.fastq 2>/dev/null | wc -l) files)"

echo
echo "✓ All large-file symlinks in place — you can now run session 4."
