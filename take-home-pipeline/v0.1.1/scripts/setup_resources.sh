#!/bin/bash
# Set up the references and the Qualimap container for this project. Run once, from the project root:
#   bash scripts/setup_resources.sh
# It downloads GENCODE (each file checked against GENCODE's MD5SUMS), builds the STAR and Salmon indexes and the TxDb,
# and downloads the Qualimap and MultiQC containers. With LSF_ACCOUNT set, each part is its own job and the builds wait for the download.
# Output: resources/<species>/ and containers/qualimap_2.2.2d.sif, containers/multiqc_1.35.sif
#   bash scripts/setup_resources.sh --step download|star|salmon|txdb|qualimap|multiqc   # one part (the jobs call this)

source "$(dirname "$0")/lib/helpers.sh"
load_config || exit 1

# Reference files per species, matching the engine references the course uses.
# Rat has no GENCODE release. Its references are NCBI RefSeq mRatBN7.2, the same assembly as UCSC rn7 (the engine's rat
# assembly). NCBI gives the genome, the GTF and the transcript FASTA; all three are downloaded (see download_rat).
RAT_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2
case "$SPECIES" in
  human) GENCODE_DIR=Gencode_human; RELEASE=release_36; PREFIX=gencode.v36
         GENOME=GRCh38.primary_assembly.genome.fa.gz; TXDB_FILE=homo_sapiens_grch38_gencode_36.sqlite
         GTF_GZ=$PREFIX.annotation.gtf.gz; TX_GZ=$PREFIX.transcripts.fa.gz ;;
  mouse) GENCODE_DIR=Gencode_mouse; RELEASE=release_M25; PREFIX=gencode.vM25
         GENOME=GRCm38.primary_assembly.genome.fa.gz; TXDB_FILE=mus_musculus_grcm38_gencode_M25.sqlite
         GTF_GZ=$PREFIX.annotation.gtf.gz; TX_GZ=$PREFIX.transcripts.fa.gz ;;
  rat)   PREFIX=GCF_015227675.2_mRatBN7.2
         GENOME=${PREFIX}_genomic.fna.gz; GTF_GZ=${PREFIX}_genomic.gtf.gz; TX_GZ=${PREFIX}_rna.fna.gz
         TXDB_FILE=rattus_norvegicus_mratbn7_2_refseq.sqlite ;;
  *) echo "ERROR: SPECIES must be human, mouse or rat (config.sh has '$SPECIES')" >&2; exit 1 ;;
esac
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/$GENCODE_DIR/$RELEASE"
REF="$RESOURCE_DIR/$SPECIES"
QUALIMAP_URL="https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--hdfd78af_2"
MULTIQC_URL="https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2"
R_URL="https://bings-aws-public-data.s3.us-east-1.amazonaws.com/ulukag01/rbings_20250925.sif"
# SHA-256 of the three containers, recorded from the test downloads (the R container is byte-identical to the course copy).
QUALIMAP_SHA256=b42a6e6b718b6428cb7018b8c3ac26df70d8df443bf0901e24072aec23d4040e
MULTIQC_SHA256=05e8af004f5f05a238bf66556c04ca137272ad578f858221bbfaf17d36732261
R_SHA256=009f0201527c4da421d6ff04699e1e998f2a8507e414a93d31c88f36dffac75e

# Rat: download the genome, GTF and transcript FASTA from NCBI into $REF/download, each checked against the MD5 that
# NCBI publishes in md5checksums.txt (retrying partial downloads, as for GENCODE).
download_rat() {
  local dl="$REF/download" f want tries
  mkdir -p "$dl"
  run curl -sfL -o "$dl/md5checksums.txt" "$RAT_BASE/md5checksums.txt" || return 1
  for f in "$GENOME" "$GTF_GZ" "$TX_GZ"; do
    want=$(awk -v f="./$f" '$2 == f {print $1}' "$dl/md5checksums.txt")
    [ -n "$want" ] || { echo "ERROR: $f is not listed in md5checksums.txt" >&2; return 1; }
    tries=0
    until (cd "$dl" && echo "$want  $f" | md5sum -c --status - 2>/dev/null); do
      tries=$((tries + 1))
      [ "$tries" -gt 5 ] && { echo "ERROR: $f failed its MD5 check after 5 downloads" >&2; return 1; }
      curl -sfL -C - -o "$dl/$f" "$RAT_BASE/$f" || true
    done
  done
}

# Download the three GENCODE files into $REF/download, checking each against MD5SUMS (retrying partial downloads).
download() {
  if [ "$SPECIES" = rat ]; then download_rat; return; fi
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
  # NCBI names the gene symbol with the attribute gene "..."; 02_deg.R reads gene_name, so rat rows get gene_name too
  gtf_fix() { if [ "$SPECIES" = rat ]; then sed 's/ gene "\([^"]*\)";/ gene "\1"; gene_name "\1";/'; else cat; fi; }
  [ -s "$out" ] || { gunzip -c "$dl/$GTF_GZ" | gtf_fix > "$out.part" && mv "$out.part" "$out"; }
  out="$REF/genome/genome.fa"
  [ -s "$out" ] || { gunzip -c "$dl/$GENOME" > "$out.part" && mv "$out.part" "$out"; }
  out="$REF/transcripts/transcripts.fa"
  [ -s "$out" ] || { gunzip -c "$dl/$TX_GZ" > "$out.part" && mv "$out.part" "$out"; }
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

# Download one container and check its SHA-256 (recorded from the test downloads). A file that does not match is
# downloaded again; a finished download that does not match is kept as .part and the step fails.
fetch_image() {
  local url="$1" img="$2" want="$3"
  mkdir -p containers
  if [ -s "$img" ] && echo "$want  $img" | sha256sum -c --status -; then echo "ok: $img (checksum matches)"; return 0; fi
  [ -e "$img" ] && { echo "checksum does not match for $img, downloading it again" >&2; rm -f "$img"; }
  curl -sfL -C - -o "$img.part" "$url" || echo "curl did not finish cleanly for $img; checking the file" >&2
  if echo "$want  $img.part" | sha256sum -c --status -; then
    mv "$img.part" "$img"; echo "ok: $img"; return 0
  fi
  echo "ERROR: checksum mismatch for $img.part (expected $want). Delete $img.part and run the step again." >&2
  return 1
}

# Qualimap container from the Galaxy depot (Qualimap 2.2.1 is not on the depot; 2.2.2d is the nearest release).
get_qualimap() { fetch_image "$QUALIMAP_URL" "$QUALIMAP_IMAGE" "$QUALIMAP_SHA256"; }

# MultiQC 1.35 container from the Galaxy depot (the stock release; the course image adds only front-end changes).
get_multiqc() { fetch_image "$MULTIQC_URL" "$MULTIQC_CONTAINER" "$MULTIQC_SHA256"; }

# R container (R, DESeq2, clusterProfiler, msigdbr), the same image the course uses, from its public S3 location.
get_rcontainer() { fetch_image "$R_URL" "$R_CONTAINER" "$R_SHA256"; }

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
