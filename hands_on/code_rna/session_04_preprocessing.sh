#!/usr/bin/env bash
# =============================================================================
#  BiNGS Bulk RNA-seq Course — Session 4
#  Preprocessing & quantification, run interactively on a reserved Minerva node
# =============================================================================
#
#  HOW TO USE THIS FILE
#  --------------------
#   Each section has two parts:
#     (a) an intro — what the tool does and why
#     (b) the command(s) — every command is filled in and runnable
#     (c) an EXPLORE OUTPUT block — commands you run to see what the tool
#         produced, what each file means, and how to read it
#   The "Explore output" blocks are the point of this exercise. Anyone can
#   run `fastqc *.fastq`; the pipeline skill is knowing what came out and
#   how to read it.
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# Section 1 — Environment
# ─────────────────────────────────────────────────────────────────────────────
# Define all the paths we'll use throughout the script. Each variable is a
# named location that would otherwise be repeated dozens of times.

# Convention: run this script from the participant's OWN hands_on/ dir, i.e.
#   cd /sc/arion/projects/BiNGS_bulk/$USER/hands_on
#   bash code_rna/session_04_preprocessing.sh
# ref_dir is the current working directory; all outputs land inside
# ${ref_dir}/data_rna (i.e., alongside the raw inputs). Each participant's own
# hands_on/ clone gives them their own private write area — no need for a
# separate BiNGS_bulk/${USER}/data_rna/ location.
ref_dir="$(pwd)"
annotation_dir="${ref_dir}/engine/annotation"
fastq_dir="${ref_dir}/data_rna/reduced_chr5/raw/fastq"   # chr5 fastqs
code_dir="${ref_dir}/code_rna"
data_dir="${ref_dir}/data_rna"

# Convention (mirrors production pipeline):
#   data_rna/reduced_chr5/raw/             — chr5 fastqs + metadata (input for this session)
#   data_rna/reduced_chr5/preprocessed/    — session-4 outputs (populated below)
#   data_rna/{reduced_chr5,all_chr}/processed/  — session-5 and -6 outputs (all_chr also
#                                            holds fridls01 full-genome inputs for session 6)
preprocessed_dir="${data_dir}/reduced_chr5/preprocessed"

# One dir per tool under preprocessed/ — keeps outputs tidy and makes
# MultiQC's job trivial later.
fastqc_dir="${preprocessed_dir}/fastqc"
trimgalore_dir="${preprocessed_dir}/trimgalore"
star_dir="${preprocessed_dir}/star"
samtools_dir="${preprocessed_dir}/samtools"
qualimap_dir="${preprocessed_dir}/qualimap"
subread_dir="${preprocessed_dir}/subread"
deeptools_dir="${preprocessed_dir}/deeptools"
salmon_dir="${preprocessed_dir}/salmon"
multiqc_dir="${preprocessed_dir}/multiqc"

mkdir -p \
  "$preprocessed_dir" \
  "$fastqc_dir" "$trimgalore_dir" \
  "$star_dir" "$samtools_dir" "$qualimap_dir" \
  "$subread_dir" "$deeptools_dir" \
  "$salmon_dir" "$multiqc_dir"



# ─────────────────────────────────────────────────────────────────────────────
# Section 2 — Concept refresher + a peek at the raw FASTQ
# ─────────────────────────────────────────────────────────────────────────────
# Pipeline overview:
#     raw fastq ─┬─ FastQC ────────────────────┐
#                └─ TrimGalore ─┬─ FastQC ─────┤
#                               ├─ STAR ─ samtools QC ─┤
#                               └─ Salmon ─────────────┤
#                                                      └─ MultiQC (one HTML)
#
# Before we run anything, look at what a raw FASTQ actually IS.
# FASTQ format: 4 lines per read.
#   Line 1: @<read_id>
#   Line 2: <bases: A/C/G/T/N>
#   Line 3: +
#   Line 4: <per-base quality, one ASCII char per base (Phred+33)>

ls -lh "$fastq_dir"
#   → four .fastq files, ~200–320 MB each

head -8 "${fastq_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1.fastq"
#   → the first TWO reads, all 8 lines

# How many reads per file?
for f in "${fastq_dir}"/*.fastq; do
  echo "$(basename "$f"): $(( $(wc -l < "$f") / 4 )) reads"
done
#   → ~1.1–1.6 M reads each (chr5-only subset of a full experiment)


# ─────────────────────────────────────────────────────────────────────────────
# Section 3 — FastQC on RAW reads          (wall-clock ~1 min/sample)
# ─────────────────────────────────────────────────────────────────────────────
# FastQC scans each FASTQ and produces an HTML report with per-base quality,
# adapter content, duplication rate, and GC bias. Run it on raw reads first.

module load fastqc/0.11.8

for f in "${fastq_dir}"/*_chr5_R1.fastq; do
  name=$(basename "$f" _RNA_chr5_R1.fastq)
  echo "▶ FastQC on ${name}"
  start=$(date +%s)

  fastqc --threads 4 "$f" --outdir "$fastqc_dir"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore FastQC output ────────────────────────────────────────────────────
# Each sample produced two things:
#   *_fastqc.html   — human-readable interactive report (open in browser)
#   *_fastqc.zip    — machine-readable version: raw text metrics + PNGs
ls -lh "$fastqc_dir"

# Unpack ONE zip to see the internals
cd "$fastqc_dir"
unzip -o SKMel147_ARID2WT_R1_RNA_chr5_R1_fastqc.zip >/dev/null
ls SKMel147_ARID2WT_R1_RNA_chr5_R1_fastqc/
#   → summary.txt         one PASS/WARN/FAIL per module (this is what MultiQC parses)
#     fastqc_data.txt     every metric as text — the source of truth
#     fastqc_report.html  the HTML (same as the top-level .html)
#     Images/             PNGs shown in the report
#     Icons/              status icons

# PASS/WARN/FAIL summary — always the first thing to look at
cat SKMel147_ARID2WT_R1_RNA_chr5_R1_fastqc/summary.txt

# The "Basic Statistics" block from the raw data (total reads, %GC, seq length)
awk '/>>Basic Statistics/,/>>END_MODULE/' \
    SKMel147_ARID2WT_R1_RNA_chr5_R1_fastqc/fastqc_data.txt

# Cross-sample summary — status of every module for every sample, on one screen
echo
echo "=== FastQC status matrix (raw reads) ==="
for z in *_fastqc.zip; do
  base=$(basename "$z" .zip)
  echo "--- $base ---"
  unzip -p "$z" "$base/summary.txt"
done
cd - >/dev/null

# What does each of the 12 FastQC modules test?
#   Basic Statistics                — reads, length, GC%, encoding
#   Per base sequence quality       — Phred quality by cycle
#   Per tile sequence quality       — flowcell-position artifacts (Illumina only)
#   Per sequence quality scores     — average-quality distribution per read
#   Per base sequence content       — %A/C/G/T by cycle (should be flat after cycle ~5)
#   Per sequence GC content         — should be roughly normal around genome mean
#   Per base N content              — %N by cycle (should be ~0)
#   Sequence Length Distribution    — variable = trimmed data, constant = raw
#   Sequence Duplication Levels     — technical vs. biological duplication
#   Overrepresented sequences       — flags adapters and rRNA contamination
#   Adapter Content                 — adapter dimer / read-through
#   Kmer Content                    — over-represented short motifs


# ─────────────────────────────────────────────────────────────────────────────
# Section 4 — TrimGalore                    (wall-clock ~2 min/sample)
# ─────────────────────────────────────────────────────────────────────────────
# TrimGalore auto-detects the adapter and trims it (plus low-quality bases at
# read ends) using cutadapt under the hood. `--fastqc` re-runs FastQC on the
# trimmed reads so we can compare before/after.

module load trim_galore/0.6.6

for f in "${fastq_dir}"/*_chr5_R1.fastq; do
  name=$(basename "$f" _RNA_chr5_R1.fastq)
  echo "▶ TrimGalore on ${name}"
  start=$(date +%s)

  trim_galore \
    --fastqc \
    --cores 4 \
    --gzip \
    --output_dir "$trimgalore_dir" \
    "$f"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore TrimGalore output ────────────────────────────────────────────────
# Each sample produced FOUR things:
#   *_trimmed.fq.gz                             gzipped trimmed reads
#   *_RNA_chr5_R1.fastq_trimming_report.txt     cutadapt log
#   *_trimmed_fastqc.html + .zip                FastQC on the trimmed reads
ls "$trimgalore_dir"

# Peek at the trimmed FASTQ (it's gzipped — use zcat, not cat)
zcat "${trimgalore_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_trimmed.fq.gz" | head -12
#   → same 4-lines-per-read format as raw, but reads are variable-length now

# How many reads survived trimming?
raw=$(( $(wc -l < "${fastq_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1.fastq") / 4 ))
trm=$(( $(zcat "${trimgalore_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_trimmed.fq.gz" | wc -l) / 4 ))
echo "Raw:     $raw reads"
echo "Trimmed: $trm reads"
awk -v r=$raw -v t=$trm 'BEGIN{printf "Dropped: %.2f%%\n", 100*(r-t)/r}'

# Read the full trimming report — key sections
cat "${trimgalore_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1.fastq_trimming_report.txt"
# Things to notice:
#   - "Adapter type" — should say "Illumina TruSeq, auto-detected"
#   - "Total reads processed" vs "Reads written (passing filters)"
#   - "Reads with adapters" — 20-70% is normal, high just means many short fragments
#   - "Removed sequences" section — how many dropped for too-short vs too-many-N

# One-line summary per sample
echo
echo "=== TrimGalore summary (all samples) ==="
for r in "${trimgalore_dir}"/*_trimming_report.txt; do
  base=$(basename "$r" _RNA_chr5_R1.fastq_trimming_report.txt)
  pass=$(grep -oP 'Reads written \(passing filters\):\s+\K[\d,]+' "$r" | head -1)
  adap=$(grep -oP 'Reads with adapters:\s+[\d,]+ \(\K[\d.]+' "$r" | head -1)
  printf "  %-30s  passed=%-12s  adapter%%=%s%%\n" "$base" "$pass" "$adap"
done


# ─────────────────────────────────────────────────────────────────────────────
# Section 5 — STAR index (OPTIONAL — already built for you)
# ─────────────────────────────────────────────────────────────────────────────
# The chr5 STAR index is already at:
#   ${annotation_dir}/star/2.7.5b/chr5_index
# For reference, this is how it was generated:
#
#   module load star/2.7.5b
#   STAR --runMode genomeGenerate \
#        --genomeDir  "${annotation_dir}/star/2.7.5b/chr5_index" \
#        --genomeFastaFiles "${annotation_dir}/data/GRCh38.p13.chr5.genome.fa" \
#        --sjdbGTFfile "${annotation_dir}/data/gencode.v36.annotation_chr5.gtf" \
#        --genomeSAindexNbases 12 \
#        --runThreadN 16
#   module purge
#
# What lives inside a STAR index?
ls "${annotation_dir}/star/2.7.5b/chr5_index"
#   → Genome, SA, SAindex          the suffix-array based genome index (big binary)
#     chrNameLength.txt            chr → length lookup
#     exonInfo.tab, geneInfo.tab   parsed GTF metadata
#     sjdbList.fromGTF.out.tab     splice junctions extracted from the GTF
#     transcriptInfo.tab           transcript metadata
#     Log.out, genomeParameters.txt  how the index was built


# ─────────────────────────────────────────────────────────────────────────────
# Section 6 — STAR alignment                (wall-clock ~3 min/sample)
# ─────────────────────────────────────────────────────────────────────────────
# STAR is a splice-aware aligner: it can align a read across an intron by
# consulting the GTF. We output BAM sorted by coordinate (needed for indexing
# and for IGV) and turn on --quantMode GeneCounts (a free gene-count table
# that we won't actually use, since Salmon is our count source, but it's a
# nice sanity check).

module load star/2.7.5b

star_index="${annotation_dir}/star/2.7.5b/chr5_index"

for f in "${trimgalore_dir}"/*_trimmed.fq.gz; do
  name=$(basename "$f" _RNA_chr5_R1_trimmed.fq.gz)
  echo "▶ STAR align ${name}"
  start=$(date +%s)

  STAR \
    --runThreadN 4 \
    --genomeDir "$star_index" \
    --readFilesIn "$f" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${star_dir}/${name}_RNA_chr5_R1_" \
    --outSAMunmapped Within \
    --quantMode GeneCounts

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore STAR output ──────────────────────────────────────────────────────
# For each sample STAR wrote (prefix = <name>_RNA_chr5_R1_):
#   *_Aligned.sortedByCoord.out.bam    the alignments, sorted by position
#   *_Log.final.out                    QC-oriented summary — READ THIS
#   *_Log.out                          full verbose run log
#   *_Log.progress.out                 per-second progress during the run
#   *_SJ.out.tab                       high-confidence splice junctions found
#   *_ReadsPerGene.out.tab             per-gene read counts (from --quantMode)
ls "$star_dir" | grep ARID2WT_R1

# The summary log — the one number-heavy file per alignment
cat "${star_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_Log.final.out"
# Line-by-line:
#   Number of input reads                        reads in the trimmed FASTQ
#   Uniquely mapped reads number / %             THE gold-standard mapping metric
#   Number of splices: Total                     splice events called
#   Mismatch rate per base                       <1% for good data
#   % of reads mapped to multiple loci           moderate OK; >30% = paralogs/rRNA
#   % of reads unmapped: too short               most common failure mode

# One-line mapping-rate summary per sample
echo
echo "=== STAR uniquely-mapped % across samples ==="
for lg in "${star_dir}"/*_Log.final.out; do
  base=$(basename "$lg" _RNA_chr5_R1_Log.final.out)
  pct=$(grep 'Uniquely mapped reads %' "$lg" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
  printf "  %-30s  %s\n" "$base" "$pct"
done

# Splice junctions found (annotated + novel)
head -5 "${star_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_SJ.out.tab"
#   Columns: chr  intron_start  intron_end  strand  motif  annotated  unique_reads  multi_reads  max_overhang
#   annotated = 1 → junction was in the GTF; 0 → novel splice event STAR discovered
wc -l "${star_dir}"/*_SJ.out.tab

# STAR's own gene counts (we use Salmon downstream, but this is a free QC)
head -8 "${star_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_ReadsPerGene.out.tab"
#   The FIRST 4 rows are library-wide stats:
#     N_unmapped     N_multimapping     N_noFeature     N_ambiguous
#   Then one row per gene:
#     gene_id  unstranded_count  forward_stranded  reverse_stranded

# Look inside the BAM (requires samtools — we load it here temporarily)
module load samtools/1.11
echo
echo "=== BAM header (first 15 lines) ==="
samtools view -H "${star_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_Aligned.sortedByCoord.out.bam" | head -15
#   @HD  file-format version + sort order
#   @SQ  one line per chromosome (name + length)
#   @PG  the STAR command that produced this BAM

echo
echo "=== First 3 alignments ==="
samtools view "${star_dir}/SKMel147_ARID2WT_R1_RNA_chr5_R1_Aligned.sortedByCoord.out.bam" | head -3
# Columns: read_id  flag  chr  pos  MAPQ  CIGAR  mate_chr  mate_pos  tlen  seq  qual  <tags>
# CIGAR examples: 76M = 76 matches; 40M100N36M = 40M, 100bp intron (N), 36M (spliced!)
module purge


# ─────────────────────────────────────────────────────────────────────────────
# Section 6.5 — subread featureCounts        (wall-clock ~10 s/sample)
# ─────────────────────────────────────────────────────────────────────────────
# featureCounts is the alignment-based counter from the subread suite. It
# takes a STAR BAM + the GTF and returns per-gene read counts. This is a
# DIFFERENT way of getting counts than Salmon (§10): featureCounts assigns
# each read to at most one gene based on where it aligned in the genome;
# Salmon assigns each read to transcripts probabilistically using a
# selective-alignment model. Both are legitimate DESeq2 inputs.
#
# Note: this is the same tool the BiNGS A3 pipeline runs (see
#   /sc/arion/projects/BiNGS/bings_analysis/code/R/modules/subread.R).
# A3 uses `-g gene_type` to build BIOTYPE counts for MultiQC. We use
# `-g gene_id` here because the pedagogical goal is a gene-level count
# matrix comparable to Salmon.

module load subread/2.0.1

genome_gtf="${annotation_dir}/data/gencode.v36.annotation_chr5.gtf"

for bam in "${star_dir}"/*_Aligned.sortedByCoord.out.bam; do
  name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  echo "▶ featureCounts on ${name}"
  start=$(date +%s)

  featureCounts \
    -a "$genome_gtf" \
    -t exon \
    -g gene_id \
    -s 0 \
    -T 4 \
    -o "${subread_dir}/${name}.featureCounts.txt" \
    "$bam"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore featureCounts output ─────────────────────────────────────────────
# Two output files per sample:
#   *.featureCounts.txt          the count matrix (7 columns)
#   *.featureCounts.txt.summary  a stats table used by MultiQC
ls "$subread_dir"

# The main output — 7 columns; row 1 is a comment, row 2 is a header, then data
head -3 "${subread_dir}/SKMel147_ARID2WT_R1.featureCounts.txt"
#   Column layout:
#     Geneid  Chr  Start  End  Strand  Length  <path-to-bam>
#     ENSG…   chr5 1234   5678 +       4444    123
# The last column name is the BAM path — awkward but that's what makes
# multi-sample matrices merge-able.

# How many genes got a non-zero count?
awk 'NR>2 && $NF>0' "${subread_dir}/SKMel147_ARID2WT_R1.featureCounts.txt" | wc -l

# The summary — where reads went and why
cat "${subread_dir}/SKMel147_ARID2WT_R1.featureCounts.txt.summary"
# Rows:
#   Assigned                — reads counted toward a gene (what DESeq2 sees)
#   Unassigned_Ambiguity    — read overlaps ≥2 gene features (dropped)
#   Unassigned_NoFeatures   — read aligned to intergenic space
#   Unassigned_MultiMapping — read aligned to >1 place in the genome
#   Unassigned_Unmapped     — read was unmapped in the BAM (0 here — filtered upstream)

# Per-sample assignment rate
echo
echo "=== featureCounts assignment rate across samples ==="
for s in "${subread_dir}"/*.featureCounts.txt.summary; do
  base=$(basename "$s" .featureCounts.txt.summary)
  assigned=$(awk '$1=="Assigned"{print $2}' "$s")
  total=$(awk 'NR>1 {sum+=$2} END{print sum}' "$s")
  awk -v n="$base" -v a="$assigned" -v t="$total" \
      'BEGIN{printf "  %-30s  assigned=%-10d  rate=%.2f%%\n", n, a, 100*a/t}'
done

# ── Compare featureCounts vs Salmon for a couple of genes ───────────────────
# (Deferred until after §10 Salmon runs — the check makes more sense once
#  both count sources exist. If you re-run this section after §10, you can
#  pick any high-count gene from featureCounts and grep the transcripts
#  in Salmon's quant.sf via the tx2gene mapping to compare orders of
#  magnitude.)


# ─────────────────────────────────────────────────────────────────────────────
# Section 7 — samtools QC: index / flagstat / idxstats / stats   (~10 s/sample)
# ─────────────────────────────────────────────────────────────────────────────
# Every BAM used with IGV or with per-region tools needs a companion .bai
# index. `flagstat` counts alignment flags. `idxstats` shows reads-per-chr.
# `stats` gives insert-size distribution and other library-level metrics.

module load samtools/1.11

for bam in "${star_dir}"/*_Aligned.sortedByCoord.out.bam; do
  name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  echo "▶ samtools QC ${name}"
  start=$(date +%s)

  samtools index    "$bam"
  samtools flagstat "$bam" > "${samtools_dir}/${name}_Aligned.sortedByCoord.out.flagstat"
  samtools idxstats "$bam" > "${samtools_dir}/${name}_Aligned.sortedByCoord.out.idxstats"
  samtools stats    "$bam" > "${samtools_dir}/${name}_Aligned.sortedByCoord.out.stats"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore samtools QC output ───────────────────────────────────────────────
# `samtools index` produced a .bai next to each .bam (NOT in samtools_dir):
ls "${star_dir}"/*.bai
#   → the .bai is a binary index that lets any tool (IGV, deepTools, samtools
#     view chr5:1000-2000, ...) jump straight to a region without reading the
#     whole BAM

ls "$samtools_dir"
#   → three files per sample: .flagstat  .idxstats  .stats

# ── flagstat: a 13-line tally of BAM flags ───────────────────────────────────
cat "${samtools_dir}/SKMel147_ARID2WT_R1_Aligned.sortedByCoord.out.flagstat"
# Key rows:
#   in total       — should ≈ trimmed read count + multi-mapping alignments
#   mapped (X%)    — ≥99% on our chr5-reduced data (bad reads never got in)
#   secondary      — extra alignments for multi-mapping reads
#   supplementary  — chimeric / split alignments (very rare in RNA-seq)
#   duplicates     — 0 here because we didn't run MarkDuplicates

# ── idxstats: reads per chromosome ───────────────────────────────────────────
head "${samtools_dir}/SKMel147_ARID2WT_R1_Aligned.sortedByCoord.out.idxstats"
#   Columns: chr  chr_length  mapped_reads  unmapped_reads
#   → almost all reads on chr5; every other chr should be ~0 (sanity check that
#     the chr5-only reduction worked)

# ── stats: 30+ metrics in "SN" (summary numbers) lines ───────────────────────
echo
echo "=== samtools stats — SN section ==="
grep '^SN' "${samtools_dir}/SKMel147_ARID2WT_R1_Aligned.sortedByCoord.out.stats" | head -25
# Notable lines:
#   raw total sequences, reads mapped, reads unmapped
#   error rate, average length, average quality
#   insert size average/std dev (PE only — 0 for our SE data)


# ─────────────────────────────────────────────────────────────────────────────
# Section 7.5 — deeptools bamCoverage → bigwig    (wall-clock ~30 s/sample)
# ─────────────────────────────────────────────────────────────────────────────
# bamCoverage converts a BAM into a bigWig (.bw) — a compact, indexed,
# per-bin coverage track. This is the exact format that IGV, UCSC, and
# the BiNGS A3 UCSC track hub all consume for visualization.
#
# CPM normalization = counts per million mapped reads. It makes the y-axis
# of the coverage track roughly comparable across samples with different
# library sizes.
#
# For stranded libraries A3 produces TWO bigwigs per sample
# (--filterRNAstrand forward + reverse). Our data was auto-detected as
# unstranded by Salmon in §10 (or you can override; see
#   /sc/arion/projects/BiNGS/bings_analysis/code/R/modules/deeptools2.R),
# so a single unstranded bigwig is fine.
#
# Requires: samtools index has already run (from §7) — bamCoverage seeks
# into the BAM via the .bai.

module load macs/2.1.0     # bundles deeptools binaries on Minerva

for bam in "${star_dir}"/*_Aligned.sortedByCoord.out.bam; do
  name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  echo "▶ bamCoverage on ${name}"
  start=$(date +%s)

  bamCoverage \
    --bam "$bam" \
    --outFileName "${deeptools_dir}/${name}_unstranded.bw" \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing CPM \
    --numberOfProcessors 4

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore bamCoverage / bigwig output ─────────────────────────────────────
# One bigwig per sample:
ls -lh "$deeptools_dir"
#   → ~1–5 MB each — much smaller than the BAM they came from,
#     because .bw stores per-bin summaries instead of per-read alignments

# What does a bigwig actually contain? Ask UCSC utilities:
module load ucsc_utils/2020_03_17 2>/dev/null || module load ucsc-tools 2>/dev/null || true
bigWigInfo "${deeptools_dir}/SKMel147_ARID2WT_R1_unstranded.bw" 2>/dev/null | head -15
#   → version, isCompressed, isSwapped, primaryDataSize/primaryIndexSize,
#     zoomLevels (indexed at multiple resolutions for fast whole-chromosome views),
#     chromCount, basesCovered, mean/min/max/std of the coverage signal

# Show coverage of a specific region (chr5:1,000,000–1,001,000) as text
bigWigToBedGraph -chrom=chr5 -start=1000000 -end=1001000 \
    "${deeptools_dir}/SKMel147_ARID2WT_R1_unstranded.bw" \
    /dev/stdout 2>/dev/null | head -10 || \
  echo "  (bigWigToBedGraph not available in loaded ucsc_utils version — skip)"
#   → four columns: chr  start  end  coverage_value

module purge

# What can you DO with these bigwigs?
#   - Load them into IGV: File → Load from File → *.bw (chr5 view)
#   - Serve them from an httpd + drop a UCSC "custom track" line pointing at
#     the URL (this is what A3's ucscgenomebrowser_prepare_hub_rna does,
#     with AWS S3 as the httpd — we don't run that in this course)
#   - Feed them to deeptools plotProfile / plotHeatmap for average-signal
#     plots over gene bodies or peak sets
#   - Feed them to pyGenomeTracks for figure-quality browser panels


# ─────────────────────────────────────────────────────────────────────────────
# Section 8 — Qualimap RNA-seq QC (OPTIONAL, ~1 min/sample)
# ─────────────────────────────────────────────────────────────────────────────
# Qualimap adds RNA-specific QC: 5'→3' coverage bias, junction saturation,
# fraction of reads in exons/introns/intergenic. Skip if short on time.

module load qualimap/2.2.1
module load java/1.8.0_211    # override auto-loaded java/21; qualimap 2.2.1 needs Java 8
genome_gtf="${annotation_dir}/data/gencode.v36.annotation_chr5.gtf"
unset DISPLAY   # avoids X11 error on compute nodes

for bam in "${star_dir}"/*_Aligned.sortedByCoord.out.bam; do
  name=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  echo "▶ qualimap rnaseq ${name}"
  start=$(date +%s)

  qualimap rnaseq \
    -bam "$bam" \
    -gtf "$genome_gtf" \
    -outdir "${qualimap_dir}/${name}" \
    --java-mem-size=4G || echo "  (qualimap failed on ${name} — skipping)"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore Qualimap output ──────────────────────────────────────────────────
ls "${qualimap_dir}/SKMel147_ARID2WT_R1"
#   → qualimapReport.html         the interactive report (open in browser)
#     rnaseq_qc_results.txt       plain-text version (grep-friendly)
#     css/, images_qualimapReport/, raw_data_qualimapReport/

# The plain-text summary — everything the HTML shows, but as text
head -60 "${qualimap_dir}/SKMel147_ARID2WT_R1/rnaseq_qc_results.txt"
# Key sections:
#   Reads alignment           total/mapped/unique/multi-mapped
#   Reads genomic origin      exonic vs intronic vs intergenic (should be exon-heavy)
#   Transcript coverage       5'→3' bias metric (mean coverage at each transcript position)
#   Junction analysis         reads on annotated vs novel junctions


# ─────────────────────────────────────────────────────────────────────────────
# Section 9 — Salmon index (OPTIONAL — already built for you)
# ─────────────────────────────────────────────────────────────────────────────
# The chr5 Salmon index is already at:
#   ${annotation_dir}/salmon/1.2.1/chr5_index
# It was built with the genome as a decoy sequence — cleaner than a plain
# transcriptome index because it lets Salmon "reject" reads coming from
# unannotated intergenic regions.
#
# What lives inside a Salmon index?
ls "${annotation_dir}/salmon/1.2.1/chr5_index"
#   → seq.bin, pos.bin, ctable.bin, ctg_offsets.bin  binary index structures
#     info.json                                       parameters + version info
#     duplicate_clusters.tsv                          identical-sequence transcripts
#     pre_indexing.log, ref_indexing.log              build logs


# ─────────────────────────────────────────────────────────────────────────────
# Section 10 — Salmon quant                 (wall-clock ~2 min/sample)
# ─────────────────────────────────────────────────────────────────────────────
# Salmon quantifies transcript abundances. `-l A` auto-detects library type;
# `--validateMappings` upgrades pure pseudo-alignment to selective alignment
# (much more robust). `--seqBias --gcBias` correct for coverage biases.

module load salmon/1.2.1
salmon_index="${annotation_dir}/salmon/1.2.1/chr5_index"

for f in "${trimgalore_dir}"/*_trimmed.fq.gz; do
  name=$(basename "$f" _RNA_chr5_R1_trimmed.fq.gz)
  echo "▶ Salmon quant ${name}"
  start=$(date +%s)

  salmon quant \
    -i "$salmon_index" \
    -l A \
    -r "$f" \
    -p 4 \
    --validateMappings \
    --seqBias \
    --gcBias \
    -o "${salmon_dir}/${name}"

  echo "  ⏱  $(( $(date +%s) - start ))s"
done

module purge

# ── Explore Salmon output ────────────────────────────────────────────────────
# One directory per sample; each contains:
#   quant.sf                   THE per-transcript quant table (tximport reads this)
#   cmd_info.json              the salmon command that produced this dir
#   lib_format_counts.json     inferred library type + counts per orientation
#   aux_info/meta_info.json    mapping rate, num processed/mapped, EM iterations
#   logs/salmon_quant.log      the run log
#   libParams/                 parameter details for reproducibility
ls "${salmon_dir}/SKMel147_ARID2WT_R1"
ls "${salmon_dir}/SKMel147_ARID2WT_R1/aux_info"

# ── quant.sf — the main output ──────────────────────────────────────────────
head -5 "${salmon_dir}/SKMel147_ARID2WT_R1/quant.sf"
#   Name             transcript ID (Ensembl ENST…)
#   Length           full transcript length in bp
#   EffectiveLength  transcript length adjusted for fragment-length distribution
#   TPM              Transcripts Per Million — length-normalized; comparable across genes
#   NumReads         estimated reads assigned to this transcript — GOES INTO DESEQ2
wc -l "${salmon_dir}/SKMel147_ARID2WT_R1/quant.sf"
#   → number of transcripts in the chr5 salmon index

# How many transcripts got a non-zero count?
awk 'NR>1 && $5>0' "${salmon_dir}/SKMel147_ARID2WT_R1/quant.sf" | wc -l

# Top-10 most highly expressed transcripts in this sample (by NumReads)
echo
echo "=== Top 10 transcripts by NumReads — SKMel147_ARID2WT_R1 ==="
(head -1 "${salmon_dir}/SKMel147_ARID2WT_R1/quant.sf" ; \
 sort -k5,5nr "${salmon_dir}/SKMel147_ARID2WT_R1/quant.sf" | head -10) | column -t

# ── lib_format_counts.json — the inferred library type ──────────────────────
cat "${salmon_dir}/SKMel147_ARID2WT_R1/lib_format_counts.json"
#   expected_format:  what -l A guessed (SR = single-end reverse-stranded, SF = single-end forward, U = unstranded)
#   compatible_fragments: reads matching the guess
#   assigned_fragments:   reads actually used for quantification

# ── meta_info.json — the "mapping rate" you'd quote in a report ─────────────
head -20 "${salmon_dir}/SKMel147_ARID2WT_R1/aux_info/meta_info.json"
# Notable keys:
#   num_processed     total reads read from the trimmed FASTQ
#   num_mapped        reads Salmon could confidently assign to a transcript
#   percent_mapped    num_mapped / num_processed × 100
#   library_types     the format(s) matched — matches lib_format_counts.json
#   num_decoy_fragments  reads whose best mapping was to a decoy (dropped)

echo
echo "=== Mapping-rate summary across samples ==="
for s in "${salmon_dir}"/SKMel147_*; do
  name=$(basename "$s")
  proc=$(grep -oP '"num_processed":\s*\K\d+' "$s/aux_info/meta_info.json")
  mapd=$(grep -oP '"num_mapped":\s*\K\d+' "$s/aux_info/meta_info.json")
  rate=$(grep -oP '"percent_mapped":\s*\K[\d.]+' "$s/aux_info/meta_info.json")
  printf "  %-30s  processed=%-10s mapped=%-10s rate=%s%%\n" "$name" "$proc" "$mapd" "$rate"
done


# ─────────────────────────────────────────────────────────────────────────────
# Section 11 — MultiQC                       (wall-clock ~30 s)
# ─────────────────────────────────────────────────────────────────────────────
# MultiQC crawls the data dir, recognizes each tool's outputs, and produces
# ONE interactive HTML report. This is what you actually show a wet-lab
# collaborator when handing off preprocessing results.

# NOTE: Minerva's `module load python/3.7.3` ships a MultiQC that crashes on
# a cryptography/pyOpenSSL bug. Use the user-installed multiqc 1.35 from
# ~/.local/bin instead (installed once via `pip install --user multiqc`
# under python/3.10.17 — see the launch instructions).
export PATH="$HOME/.local/bin:$PATH"

cd "$preprocessed_dir"

multiqc \
  --config  "${ref_dir}/code_rna/multiqc_config_rnaseq.yaml" \
  --title   "BiNGS Bulk RNA-seq Course — ${USER}" \
  --filename "multiqc_report_${USER}.html" \
  --outdir  "$multiqc_dir" \
  fastqc/ trimgalore/ star/ samtools/ subread/ deeptools/ qualimap/ salmon/

# ── Explore MultiQC output ───────────────────────────────────────────────────
ls "$multiqc_dir"
#   → multiqc_report_${USER}.html          the one thing you actually open
#     multiqc_report_${USER}_data/          all parsed values as TSV (great for scripts)

ls "${multiqc_dir}/multiqc_report_${USER}_data"
#   → multiqc_general_stats.txt   one-line-per-sample summary (numbers from every tool)
#     multiqc_fastqc.txt          every FastQC metric across all samples
#     multiqc_star.txt            every STAR Log.final.out across all samples
#     multiqc_salmon.txt          every Salmon meta_info metric
#     multiqc_samtools_*.txt      flagstat/stats/idxstats aggregated
#     multiqc_sources.txt         WHICH FILE each metric came from
#     multiqc.log, multiqc_data.json

# The most useful table — one row per sample, columns from every tool
echo
echo "=== MultiQC general stats (all samples, all tools) ==="
column -t -s $'\t' "${multiqc_dir}/multiqc_report_${USER}_data/multiqc_general_stats.txt"

# Where did MultiQC pull each metric from? (debugging aid)
head "${multiqc_dir}/multiqc_report_${USER}_data/multiqc_sources.txt"

# Now open the HTML in OnDemand Files. As a session we'll walk through:
#   - General Statistics section (compare to the TSV above)
#   - STAR bar chart: uniquely-mapped % per sample
#   - Salmon bar chart: assigned reads + inferred fragment length
#   - FastQC (raw) vs FastQC (trimmed): did trimming actually improve quality?
#   - Software Versions section at the bottom (reproducibility)


# ─────────────────────────────────────────────────────────────────────────────
# End of session 4
# ─────────────────────────────────────────────────────────────────────────────
# Your Salmon quant.sf files under ${salmon_dir} are what session 5 will read
# into R with tximport. 