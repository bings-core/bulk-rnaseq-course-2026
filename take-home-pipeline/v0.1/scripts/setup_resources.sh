#!/bin/bash
# Set up the references and the Qualimap container for this project. Run once, from the project root:
#   bash scripts/setup_resources.sh
# It downloads GENCODE (each file checked against GENCODE's MD5SUMS), builds the STAR and Salmon indexes and the TxDb,
# and downloads the Qualimap and MultiQC containers. With LSF_ACCOUNT set, each part is its own job and the builds wait for the download.
# Output: resources/<species>/ and containers/qualimap_2.2.2d.sif, containers/multiqc_1.35.sif
#   bash scripts/setup_resources.sh --step download|star|salmon|txdb|qualimap|multiqc   # one part (the jobs call this)

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1

# GENCODE release per species, matching the engine references the course uses.
case "$SPECIES" in
  human) GENCODE_DIR=Gencode_human; RELEASE=release_36; PREFIX=gencode.v36
         GENOME=GRCh38.primary_assembly.genome.fa.gz; TXDB_FILE=homo_sapiens_grch38_gencode_36.sqlite ;;
  mouse) GENCODE_DIR=Gencode_mouse; RELEASE=release_M25; PREFIX=gencode.vM25
         GENOME=GRCm38.primary_assembly.genome.fa.gz; TXDB_FILE=mus_musculus_grcm38_gencode_M25.sqlite ;;
  *) echo "ERROR: SPECIES must be human or mouse (config.sh has '$SPECIES')" >&2; exit 1 ;;
esac
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/$GENCODE_DIR/$RELEASE"
REF="$RESOURCE_DIR/$SPECIES"
QUALIMAP_URL="https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--hdfd78af_2"
MULTIQC_URL="https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2"
R_URL="https://bings-aws-public-data.s3.us-east-1.amazonaws.com/ulukag01/rbings_20250925.sif"

# Download the three GENCODE files into $REF/download, checking each against MD5SUMS (retrying partial downloads).
download() {
  local dl="$REF/download"
  mkdir -p "$dl"
  run curl -sfL -o "$dl/MD5SUMS" "$BASE/MD5SUMS" || return 1
  local f tries
  for f in "$PREFIX.annotation.gtf.gz" "$PREFIX.transcripts.fa.gz" "$GENOME"; do
    tries=0
    until (cd "$dl" && grep " $f\$" MD5SUMS | md5sum -c --status -); do
      tries=$((tries + 1))
      [ "$tries" -gt 5 ] && { echo "ERROR: $f failed its MD5 check after 5 downloads" >&2; return 1; }
      curl -sfL -C - -o "$dl/$f" "$BASE/$f" || true
    done
  done
}

# Unpack the downloads into the layout the pipeline reads (plain FASTA and GTF, no compression).
prepare() {
  local dl="$REF/download" out
  mkdir -p "$REF/gtf" "$REF/genome" "$REF/transcripts"
  out="$REF/gtf/annotation.gtf"
  [ -s "$out" ] || { gunzip -c "$dl/$PREFIX.annotation.gtf.gz" > "$out.part" && mv "$out.part" "$out"; }
  out="$REF/genome/genome.fa"
  [ -s "$out" ] || { gunzip -c "$dl/$GENOME" > "$out.part" && mv "$out.part" "$out"; }
  out="$REF/transcripts/transcripts.fa"
  [ -s "$out" ] || { gunzip -c "$dl/$PREFIX.transcripts.fa.gz" > "$out.part" && mv "$out.part" "$out"; }
}

# STAR index (2.7.5b, sjdbOverhang 100 as in the engine references).
build_star() {
  module load star/2.7.5b
  mkdir -p "$REF/star/2.7.5b/index"
  run STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$REF/star/2.7.5b/index" \
      --genomeFastaFiles "$REF/genome/genome.fa" --sjdbGTFfile "$REF/gtf/annotation.gtf" --sjdbOverhang 100
}

# Salmon index (1.2.1, transcripts only, as in the engine references).
build_salmon() {
  module load salmon/1.2.1
  run salmon index -t "$REF/transcripts/transcripts.fa" -i "$REF/salmon/1.2.1/index" -p "$THREADS"
}

# TxDb used by 02_deg.R (gene and transcript models), built from the GTF in the R container.
build_txdb() {
  mkdir -p "$REF/txdb"
  run bash scripts/run_r.sh scripts/lib/make_txdb.R
}

# Qualimap container from the Galaxy depot (Qualimap 2.2.1 is not on the depot; 2.2.2d is the nearest release).
get_qualimap() {
  mkdir -p containers
  local img="$QUALIMAP_IMAGE"
  [ -s "$img" ] && return 0
  run curl -sfL -o "$img.part" "$QUALIMAP_URL" && mv "$img.part" "$img"
}

# MultiQC 1.35 container from the Galaxy depot (the stock release; the course image adds only front-end changes).
get_multiqc() {
  mkdir -p containers
  local img="$MULTIQC_CONTAINER"
  [ -s "$img" ] && return 0
  run curl -sfL -o "$img.part" "$MULTIQC_URL" && mv "$img.part" "$img"
}

# R container (R, DESeq2, clusterProfiler, msigdbr), the same image the course uses, from its public S3 location.
get_rcontainer() {
  mkdir -p containers
  local img="$R_CONTAINER"
  [ -s "$img" ] && return 0
  run curl -sfL -C - -o "$img.part" "$R_URL" && mv "$img.part" "$img"
}

if [ "$1" = "--step" ]; then
  case "$2" in
    download) download && prepare ;;
    star)     build_star ;;
    salmon)   build_salmon ;;
    txdb)     build_txdb ;;
    qualimap) get_qualimap ;;
    multiqc)  get_multiqc ;;
    rcontainer) get_rcontainer ;;
    *) echo "usage: $0 [--step download|star|salmon|txdb|qualimap|multiqc|rcontainer]" >&2; exit 1 ;;
  esac
  exit
fi

# Whole setup. With LSF: download first, then the four jobs that need it. Without LSF: one part after another.
mkdir -p "$REF"
if [ -z "$LSF_ACCOUNT" ]; then
  for step in rcontainer download star salmon txdb qualimap multiqc; do bash "$0" --step "$step" || exit 1; done
  exit
fi
dl=$(submit "setup_download" 1 "$MEM_FETCH_GB" -- bash "$0" --step download) || exit 1
rc=$(submit "setup_rcontainer" 1 "$MEM_FETCH_GB" -- bash "$0" --step rcontainer) || exit 1
submit "setup_star"     "$THREADS" "$MEM_STAR_GB"    "done($dl)" -- bash "$0" --step star
submit "setup_salmon"   "$THREADS" "$MEM_SALMON_GB"  "done($dl)" -- bash "$0" --step salmon
submit "setup_txdb"     1          "$MEM_DE_GB"      "done($dl) && done($rc)" -- bash "$0" --step txdb
submit "setup_qualimap" 1          "$MEM_FETCH_GB"   -- bash "$0" --step qualimap
submit "setup_multiqc"  1          "$MEM_FETCH_GB"   -- bash "$0" --step multiqc
echo "Submitted. Check with bjobs; the references are ready when the four setup_ jobs are DONE."
