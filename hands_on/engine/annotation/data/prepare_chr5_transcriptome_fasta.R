# =============================================================================
#  How the chr5 transcriptome FASTA was built
# =============================================================================
#
#  Extract every transcript_id present in the chr5 GTF, then subset the
#  full-genome GENCODE transcripts FASTA to just those transcripts. The
#  resulting `regions.txt` is a plain list of FASTA headers (one per line)
#  that a downstream `samtools faidx` / `seqkit grep` invocation can use to
#  cut a chr5-only transcriptome FASTA out of the full one.
#
#  This script is documentation for how we generated
#  `gencode.v36.chr5.transcripts.fa` — the Salmon index in
#  `engine/annotation/salmon/1.2.1/chr5_index/` is built on top of it. You
#  do NOT need to run this again unless the reference GENCODE version bumps
#  or the chr5 subset changes.
#
#  CONVENTION: launch from the hands_on/ dir, i.e.
#     cd /sc/arion/projects/BiNGS/bings_omics/data/bings/2026/bulk_course_data/hands_on
#     R -f engine/annotation/data/prepare_chr5_transcriptome_fasta.R
#
#  Inputs  (from engine/annotation/data/):
#     gencode.v36.annotation_chr5.gtf   — chr5-subset GENCODE GTF (already staged)
#     gencode.v36.transcripts.fa        — FULL GENCODE v36 transcripts FASTA
#                                         (NOT staged in the hands_on tree — grab
#                                          the full-genome fasta from wherever the
#                                          reference lives; the download link is
#                                          in the GENCODE v36 release notes)
#
#  Output:
#     regions.txt                       — chr5-only FASTA headers, one per line
# =============================================================================

# Resolve paths off the current working directory (assumed = hands_on/)
data_dir <- if (dir.exists("engine/annotation/data")) {
  "engine/annotation/data"
} else {
  # Fallback: script may have been sourced from its own dir, or from an
  # unrelated cwd — reach for the absolute path in the shared reference tree
  "/sc/arion/projects/BiNGS/bings_omics/data/bings/2026/bulk_course_data/hands_on/engine/annotation/data"
}

gtf_path <- file.path(data_dir, "gencode.v36.annotation_chr5.gtf")
fa_path  <- file.path(data_dir, "gencode.v36.transcripts.fa")   # full-genome tx FASTA
out_path <- file.path(data_dir, "regions.txt")

message("GTF        : ", gtf_path)
message("Transcripts: ", fa_path)
message("Output     : ", out_path)

# Read chr5 GTF and pull out transcript rows
gtf_chr5    <- rtracklayer::import(gtf_path)
gtf_chr5_df <- as.data.frame(gtf_chr5)
gtf_chr5_df_tx <- gtf_chr5_df[gtf_chr5_df$type == "transcript", ]
message("chr5 transcripts in GTF: ", nrow(gtf_chr5_df_tx))

# Read all FASTA headers from the full-genome transcript FASTA
fa_all    <- readLines(fa_path)
fa_all    <- grep("^>", fa_all, value = TRUE)
fa_all    <- gsub("^>", "", fa_all)
fa_all_tx <- sapply(strsplit(fa_all, split = "|", fixed = TRUE), "[[", 1)

# Keep only the headers whose transcript_id is in the chr5 GTF
keep_idx <- fa_all_tx %in% gtf_chr5_df_tx$transcript_id
message("Matched to full-genome FASTA: ", sum(keep_idx), " / ", length(fa_all))

writeLines(fa_all[keep_idx], out_path)
message("Wrote ", out_path)
