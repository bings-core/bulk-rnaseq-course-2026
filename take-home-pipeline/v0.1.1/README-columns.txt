samples.csv
  sample_id      unique name, no spaces
  condition      group label used in contrasts
  replicate      replicate number within condition
  fastq_1        path to R1 (or the only file for single-end). Several lanes: comma-separated, concatenated in order
  fastq_2        path to R2. Leave EMPTY for single-end data
  strandedness   auto (default), unstranded, forward or reverse. auto is detected from Salmon and written to the report
  batch          optional. If filled for any sample, batch is added to the DESeq2 design

contrasts.csv
  name           label for the comparison (used in output file names)
  treatment      condition in the numerator
  reference      condition in the denominator (log2FC = treatment vs reference)
