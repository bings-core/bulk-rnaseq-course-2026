# =============================================================================
#  BiNGS Bulk RNA-seq Course — Session 5
#  Quantification (chr5) + Differential expression (full genome)
# =============================================================================


# ─── Section 0 — Session overview ────────────────────────────────────────────
#
#  Two-part session:
#
#  PART A — CHR5 QUANTIFICATION (Sections 3-7)
#    Take the four `quant.sf` files from session 4 (chr5-only, quick), import
#    with tximport, build a DESeqDataSet, and save the count tables under
#    reduced_chr5/processed/gene_expression/tables/ (raw + normalized medianratios
#    + VST + sample_metadata). Educational: shows the quantification pipeline
#    end-to-end on a small, fast dataset.
#
#  PART B — FULL-GENOME DEG (Sections 8-20)   ← real biology lives here
#    Switch to completed full-genome pipeline results (already
#    computed upstream). We reuse raw + normalized counts and DEG
#    table (all under all_chr/processed/) to run the interpretation half of
#    the session — dispersion, PCA, sample distances, volcano, MA plot,
#    ComplexHeatmap top-50 DEGs, genome-wide log2FC circos, per-gene boxplots
#    (with p.adj + LFC overlaid), patchwork combined volcano/MA figure, and
#    gene-list export. Chr5-only DEGs are too noisy (n=2,600 genes) for
#    meaningful biology; the full-genome dataset (22,717 tested genes,
#    8 samples) is what you'd interpret in practice.
#
#  Rationale: preprocessing and quantification are the pedagogically-important
#  mechanics; the ARID2 biology is what participants should walk away with. Doing
#  quantification on chr5 keeps session fast; doing downstream on all_chr keeps
#  the results real.


# ─── Section 1 — Environment ────────────────────────────────────────────────

# On a headless compute node R has no X11 display, so ANY base-R plot call
# would open the fallback pdf() device and dump into "Rplots.pdf" in cwd.
# We always route that fallback to /dev/null — every real save goes through
# explicit png()/pdf() opens inside save_figure_both / save_base_plot_both,
# so nothing "user-visible" is suppressed. Guarantees no stray Rplots.pdf
# whether the script is sourced interactively OR run via `R -f`.
pdf(NULL)

# Layout convention:
#   data_rna/reduced_chr5/raw/                 — chr5 fastqs + sample metadata
#   data_rna/reduced_chr5/preprocessed/        — session-4 outputs
#   data_rna/reduced_chr5/processed/           — session-5 outputs, split into:
#     ├── gene_expression/{tables,figures}/    (matrices, PCA, distances, dispersion)
#     └── differential_expression/{tables,figures}/  (DEG CSVs, volcano, MA, heatmap, boxplots)
#   data_rna/all_chr/processed/                — full-genome inputs
#                                                AND session-6 outputs (co-located)

user <- Sys.getenv("USER")
message("Running as user: ", user)
# data_dir lives inside the participant's own hands_on/ clone — the SAME tree
# ref_dir points to. Every session writes into data_rna/ next to its inputs.
data_dir <- file.path(getwd(), "data_rna")

# Session-4 salmon quants live under reduced_chr5/preprocessed/salmon/
salmon_input_dir <- file.path(data_dir, "reduced_chr5", "preprocessed", "salmon")

# Convention: run this script from the participant's OWN hands_on/ dir, i.e.
#   cd /sc/arion/projects/BiNGS_bulk/$USER/hands_on
#   ml singularity/3.6.4
#   singularity exec engine/singularity_containers/rbings_20250925.sif R
#   source("code_rna/session_05_deg_analysis.R") 
# ref_dir is the current working directory (fall back to the shared path if
# the script was launched from elsewhere).
.hands_on_default <- "/sc/arion/projects/BiNGS/bings_omics/data/bings/2026/bulk_course_data/hands_on"
ref_dir <- if (dir.exists(file.path(getwd(), "engine", "annotation"))) getwd() else .hands_on_default
if (!dir.exists(salmon_input_dir)) {
  warning("No salmon dir under ", salmon_input_dir,
          " — falling back to shared reference at ", ref_dir,
          "/data_rna/reduced_chr5/preprocessed/salmon")
  salmon_input_dir <- file.path(ref_dir, "data_rna", "reduced_chr5", "preprocessed", "salmon")
}
message("Reading Salmon quants from: ", salmon_input_dir)

# Chr5 sample metadata (read from the shared reference CSV)
sample_metadata_path <- file.path(ref_dir, "data_rna", "reduced_chr5", "raw",
                                  "sample_metadata", "sample_metadata_rna.csv")
stopifnot(file.exists(sample_metadata_path))

# Where the chr5 reference annotation lives:
annotation_dir <- file.path(ref_dir, "engine", "annotation")

# Two output areas under reduced_chr5/processed/, each with tables/ + figures/
ge_dir  <- file.path(data_dir, "reduced_chr5", "processed", "gene_expression")
deg_dir <- file.path(data_dir, "reduced_chr5", "processed", "differential_expression")
ge_tables  <- file.path(ge_dir,  "tables")
ge_figures <- file.path(ge_dir,  "figures")
deg_tables <- file.path(deg_dir, "tables")
deg_figures <- file.path(deg_dir, "figures")
for (d in c(ge_tables, ge_figures, deg_tables, deg_figures)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
message("Chr5 gene-expression outputs → ", ge_dir,  " (tables/)")

# All_chr areas — used starting Section 8 for the full-genome downstream.
# Inputs (raw / normalized counts, DEG csv, sample_metadata) live
# in the STUDENT dir under all_chr/processed/ — the LSF runner pre-populates
# them from the shared REF at the start of the run so a participant's local
# workspace mirrors what would be handed to them at session time.
# Outputs (VST table, PCA / dispersion / distance / volcano / MA / heatmap /
# boxplot figures, gene lists) also go under the STUDENT dir so the final
# rsync mirrors them cleanly back to REF without wiping the inputs.
allchr_root       <- file.path(data_dir, "all_chr", "processed")
allchr_ge_tables  <- file.path(allchr_root, "gene_expression",         "tables")
allchr_ge_figures <- file.path(allchr_root, "gene_expression",         "figures")
allchr_de_tables  <- file.path(allchr_root, "differential_expression", "tables")
allchr_de_figures <- file.path(allchr_root, "differential_expression", "figures")

# Fall back to the shared reference if this participant's all_chr inputs are
# missing (parallels the salmon_input_dir fallback above).
allchr_ref_ge_tables <- file.path(ref_dir, "data_rna", "all_chr", "processed",
                                   "gene_expression", "tables")
allchr_ref_de_tables <- file.path(ref_dir, "data_rna", "all_chr", "processed",
                                   "differential_expression", "tables")
input_probe <- file.path(allchr_de_tables, "deseq_comparison_KO_vs_CTRL_all_genes.csv")
if (!file.exists(input_probe)) {
  warning("all_chr inputs missing under ", allchr_root,
          " — falling back to shared reference at ", ref_dir,
          "/data_rna/all_chr/processed/")
  allchr_ge_tables <- allchr_ref_ge_tables
  allchr_de_tables <- allchr_ref_de_tables
}

for (d in c(allchr_ge_figures, allchr_de_figures)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
message("All_chr downstream → ", allchr_root, "/{gene_expression,differential_expression}/")


# ─── Section 2 — Libraries ──────────────────────────────────────────────────

library(GenomicFeatures)   # TxDb objects for genome annotation
library(AnnotationDbi)     # SQLite-backed annotation lookup
library(tximport)          # import Salmon (and other) quantification into R
library(DESeq2)            # differential expression
library(apeglm)            # log2FC shrinkage estimator
library(ggplot2)           # plotting
library(ggrepel)           # non-overlapping text labels
library(pheatmap)          # pretty heatmaps
library(RColorBrewer)      # color palettes
library(ComplexHeatmap)    # DEG heatmap split by direction (Section 15)
library(circlize)          # colorRamp2 + genome-wide log2FC circos (Section 16)
library(patchwork)         # combined "Figure 1" panel (Section 18)
suppressPackageStartupMessages({
  library(org.Hs.eg.db)                       # ENSG → ENTREZ for circos (Section 16)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # gene coordinates for circos
})

# Shared enlarged-text theme so titles/axes/legends are legible on a slide.
bump_theme <- function(base_size = 14) {
  theme_linedraw(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      axis.title    = element_text(size = base_size + 1),
      axis.text     = element_text(size = base_size),
      legend.title  = element_text(size = base_size),
      legend.text   = element_text(size = base_size - 1),
      plot.title.position = "plot"
    )
}


# ─── Section 2.5 — save_figure_panel helper ─────────────────────────────────
#
#  Every plot section below saves via this helper so figures land on disk under
#  `<figures_dir>/<output_type>/<figure_name>.<ext>` and the file layout matches
#  what the rest of the pipeline produces.
#
#  Usage (ggplot / pheatmap):
#      save_figure_panel(p, output_directory=figures_dir,
#                        figure_name="pca_arid2", output_type="png",
#                        width=800, height=600)
#  Usage (base R plot — plotDispEsts, plotMA, etc.):
#      plotDispEsts(dds)                  # draws to current device
#      p_disp <- recordPlot()             # capture the drawing
#      save_figure_panel(p_disp, output_directory=figures_dir,
#                        figure_name="dispersion", output_type="png",
#                        width=800, height=600)

# Convenience helper: save the same figure as BOTH png (raster) and pdf
# (vector) at IDENTICAL physical dimensions so text/elements look the same.
#
# Unit convention:
#   * `width` / `height` are PIXELS for the PNG at `dpi` dots-per-inch.
#   * Physical size = width/dpi × height/dpi INCHES.
#   * We pass `res = dpi` to png() (default is 72!) so the raster's physical
#     size actually equals width/dpi inches, not width/72. Without this the
#     PDF (width/dpi in) would look smaller than the PNG (width/72 in),
#     causing point-sized text to appear proportionally bigger in the PDF.
#   * pdf() width/height are already in inches, so we pass width/dpi directly.
save_figure_both <- function(figure_panel, output_directory, figure_name,
                             width, height, dpi = 100, ...) {
  # PNG — width/height in pixels, res=dpi so the physical size matches the pdf
  save_figure_panel(figure_panel, output_directory, figure_name,
                    output_type = "png",
                    width = width, height = height, res = dpi, ...)
  # PDF — inches
  width_in  <- width  / dpi
  height_in <- height / dpi
  save_figure_panel(figure_panel, output_directory, figure_name,
                    output_type = "pdf",
                    width = width_in, height = height_in, ...)
  message("  ↳ saved  ", file.path(output_directory, "png", paste0(figure_name, ".png")))
}

# Base-R plot companion: takes an unevaluated expression that draws to the
# current graphics device, and saves both png+pdf under
# <output_directory>/<png|pdf>/<figure_name>.<ext>.
save_base_plot_both <- function(expr, output_directory, figure_name,
                                width, height, dpi = 100, pointsize = 12) {
  for (ext in c("png", "pdf")) {
    dir.create(file.path(output_directory, ext),
               showWarnings = FALSE, recursive = TRUE)
  }
  png_path <- file.path(output_directory, "png", paste0(figure_name, ".png"))
  pdf_path <- file.path(output_directory, "pdf", paste0(figure_name, ".pdf"))
  # PNG (pixels)
  png(png_path, width = width, height = height, pointsize = pointsize)
  eval(expr, envir = parent.frame())
  invisible(dev.off())
  # PDF (inches)
  pdf(pdf_path, width = width / dpi, height = height / dpi)
  eval(expr, envir = parent.frame())
  invisible(dev.off())
  message("  ↳ saved  ", png_path)
}

save_figure_panel = function(figure_panel, output_directory, figure_name, output_type = "html", use_orca = FALSE, self_contained = FALSE, background = "white", title = NULL, knitr_options = list(), ...) {
    figure_plot_folder = file.path(output_directory, output_type)
    if (file.exists(figure_plot_folder) == FALSE) {
        dir.create(figure_plot_folder, recursive = TRUE)
    }

    if (output_type == "html") {
        figure_plot_folder = file.path(figure_plot_folder, ifelse(self_contained, "self_contained", "dependent"))
        figure_data_folder = file.path(figure_plot_folder, "lib")
        if (file.exists(figure_data_folder) == FALSE) {
            dir.create(figure_data_folder, recursive = TRUE)
        }
        figure_path = file.path(figure_plot_folder, paste0(figure_name, ".", output_type))
        htmlwidgets::saveWidget(widget = figure_panel,
                            file = figure_path,
                            selfcontained = self_contained,
                            libdir = figure_data_folder,
                            background = background,
                            title = ifelse(is.null(title), session(figure_panel)[[1]], title),
                            knitrOptions = knitr_options)
    } else if (output_type %in% c("png", "jpg", "jpeg", "eps", "svg", "pdf")) {
        # Save an html panel (e.g. plotly) as one of these data formats using orca
        if (use_orca == TRUE) {
            current_working_directory = getwd()
            setwd(figure_plot_folder)
            figure_path = paste0(figure_name, ".", output_type)
            res = tryCatch({
                plotly::orca(figure_panel, figure_path, more_args = c('--disable-gpu'), ...)
            }, error = function(err) {
                print(paste0("ERROR: Orca cannot export the figure panel to: ", figure_path, " !!!"))
            }, finally = {
                setwd(current_working_directory)
            })
            figure_path = file.path(figure_plot_folder, paste0(figure_name, ".", output_type))
        } else {
            figure_path = file.path(figure_plot_folder, paste0(figure_name, ".", output_type))
            if (output_type == "png") {
                png(filename = figure_path, bg = background, ...)
                print(figure_panel)
                dev.off()
            } else if (output_type %in% c("jpg", "jpeg")) {
                jpeg(filename = figure_path, bg = background, ...)
                print(figure_panel)
                dev.off()
            } else if (output_type == "eps") {
                cairo_ps(filename = figure_path, bg = background, ...)
                print(figure_panel)
                dev.off()
            } else if (output_type == "svg") {
                svg(filename = figure_path, bg = background, ...)
                print(figure_panel)
                dev.off()
            } else if (output_type == "pdf") {
                pdf(file = figure_path, title = title, bg = background, ...)
                print(figure_panel)
                dev.off()
            } else {
                figure_path = NULL
                stop(paste0("Invalid output type: ", output_type, "!!!"))
            }
        }
    } else {
        figure_path = NULL
        stop(paste0("Figure output type not supported: ", output_type))
    }
    return(figure_path)
}


# ─── Section 3 — The experiment ─────────────────────────────────────────────
#
#  4 RNA-seq libraries from SKmel147 (melanoma):
#    - ARID2 WT: control, two biological replicates
#    - ARID2 KO: CRISPR knockout of the PBAF/SWI-SNF subunit ARID2, two reps
#  Only reads that mapped to chromosome 5 (~181 Mb, ~2,600 genes) were kept.

# Read the chr5 sample metadata from the shared reference CSV. Schema:
#   sample_id, file_name, replicate, genotype, cellline, species
sample_metadata_df <- read.csv(sample_metadata_path, stringsAsFactors = FALSE)

# Give the display frame the nicer column names the rest of the script
# expects — File Names / Sample Names / Genotype / Cell Line
colnames(sample_metadata_df) <- c("Sample Names", "File Names", "Replicate",
                                  "Genotype", "Cell Line", "Species")
rownames(sample_metadata_df) <- sample_metadata_df[["Sample Names"]]

# ── Explore ──
print(sample_metadata_df)


# ─── Section 4 — Annotation data ────────────────────────────────────────────
#
#  Genome annotation (GTF → TxDb) tells us which transcripts belong to which
#  genes and where each exon/intron/UTR is on the genome. The TxDb SQLite is
#  pre-built so nothing needs to be re-parsed at run time.

gtf_path  <- file.path(annotation_dir, "data", "gencode.v36.annotation_chr5.gtf")
txdb_path <- file.path(annotation_dir, "data", "gencode.v36.annotation_chr5.sqlite")
message("GTF:  ", gtf_path)
message("TxDb: ", txdb_path)

# For reference only — do NOT run. This is how the TxDb was built:
#
#   txdb <- GenomicFeatures::makeTxDbFromGFF(
#     file       = gtf_path,
#     dataSource = "gencode - grch38 - v36 - homo_sapiens",
#     organism   = "Homo sapiens"
#   )
#   AnnotationDbi::saveDb(txdb, txdb_path)

# Load the pre-built TxDb
txdb <- AnnotationDbi::loadDb(txdb_path)

# Pull a compact gene-level metadata table from the GTF (gene symbols,
# biotypes, HGNC IDs — useful when interpreting a DEG list).
gtf    <- rtracklayer::import(gtf_path)
gtf_df <- as.data.frame(gtf)
gene_metadata <- gtf_df[gtf_df$type == "gene",
                        c("gene_id", "gene_type", "gene_name", "havana_gene", "hgnc_id")]
rownames(gene_metadata) <- gene_metadata$gene_id

# ── Explore ──
head(gene_metadata, 20)
table(gene_metadata$gene_type)      # count by biotype (protein_coding, lncRNA, …)


# ─── Section 5 — Import Salmon quantifications ──────────────────────────────
#
#  Salmon writes per-TRANSCRIPT abundance. DESeq2 wants per-GENE counts, so:
#    (a) build a transcript→gene map from the TxDb
#    (b) let tximport aggregate transcript counts up to gene level.

sample_metadata_df$file_path <- file.path(
  salmon_input_dir,
  gsub("_RNA_chr5_R1.fastq$", "", sample_metadata_df[["File Names"]]),
  "quant.sf"
)

# Sanity: every quant.sf exists
stopifnot(all(file.exists(sample_metadata_df$file_path)))

# Save the metadata alongside the DEG results for provenance
write.csv(sample_metadata_df, file.path(ge_tables, "sample_metadata.csv"),
          row.names = FALSE)

files <- setNames(sample_metadata_df$file_path,
                  sample_metadata_df[["Sample Names"]])

# Transcript → gene map from the TxDb
k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Transcript-level import first (mostly for interest — we won't use it)
txi_tx   <- tximport(files, type = "salmon", txOut = TRUE)

# Aggregate to gene-level counts (this is what DESeq2 wants)
txi_gene <- summarizeToGene(txi_tx, tx2gene)

# ── Explore ──
str(txi_tx,   max.level = 1)
str(txi_gene, max.level = 1)

raw_gene_counts <- as.matrix(round(txi_gene$counts))
mode(raw_gene_counts) <- "integer"
write.csv(raw_gene_counts,
          file.path(ge_tables, "raw_gene_expression.csv"),
          row.names = TRUE)

cat("✓ Count matrix:", nrow(raw_gene_counts), "genes ×",
                       ncol(raw_gene_counts), "samples\n")


# ─── Section 6 — Build the DESeq2 object ────────────────────────────────────
#
#  DESeqDataSetFromTximport() bundles counts + metadata + a "design formula"
#  that says "model expression as a function of Genotype". The design formula
#  is what tells DESeq2 which contrasts are possible.

# Guardrail: metadata rows and count-matrix columns must line up
stopifnot(all(colnames(txi_gene$counts) == rownames(sample_metadata_df)))

dds <- DESeqDataSetFromTximport(
  txi_gene,
  colData = sample_metadata_df,
  design  = ~ Genotype
)
# Attach the gene-level metadata (biotype, symbol, etc.) to the dds
mcols(dds) <- gene_metadata[rownames(dds), ]

# Pre-filter very lowly expressed genes: ≥10 reads in ≥2 samples (= smallest group)
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]
message("After pre-filter: ", nrow(dds), " genes retained")

# Force the reference level so the direction of change is intuitive
# (positive log2FC → higher in KO than WT)
dds$Genotype <- relevel(dds$Genotype, ref = "ARID2_WT")

# Fit the model (estimates size factors, dispersions, and Wald tests all at once)
dds <- DESeq(dds)

# 🎓 The design formula in one sentence:
#   ~ Genotype  tells DESeq2 to fit
#     log(mean_expression) = intercept + β_Genotype × I(sample is KO)
#   β_Genotype is what we test for significance. If we had batch effects,
#   we'd add + Batch to model them out.


# ─── Section 7 — Exploratory analysis: normalization ────────────────────────
#
#  Two ways to normalize:
#    - Median-of-ratios: DESeq2's default *count-scale* normalization.
#      Good for looking at individual genes.
#    - VST (variance-stabilizing transform): a *log-scale* transform whose
#      variance is roughly constant across the mean. Good for PCA, clustering,
#      heatmaps.

norm_medratios <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_medratios),
          file.path(ge_tables, "normalized_gene_expression_medianratios.csv"),
          row.names = TRUE)

norm_vst <- vst(dds, blind = TRUE, nsub = 500)   # nsub only needed on tiny (chr5) datasets
write.csv(as.data.frame(assay(norm_vst)),
          file.path(ge_tables, "normalized_gene_expression_vst.csv"),
          row.names = TRUE)

# ── Explore ──
head(norm_medratios, 5)
head(assay(norm_vst), 5)


# ═══════════════════════════════════════════════════════════════════════════
#  PART B — FULL-GENOME DEG (upstream pipeline results)
# ═══════════════════════════════════════════════════════════════════════════
#
#  Chr5-only DEG (Part A) produces very few significant hits — after fitting
#  dispersion on ~2,600 chr5 genes with n=2 replicates per genotype, most
#  candidates lose significance to FDR. That's a fact about the reduced
#  dataset, not about the biology.
#
#  For the interpretation half of the session we switch to the full-genome
#  pipeline results (already computed upstream, 22,717 tested
#  genes across 8 samples: 4 KO + 4 WT). Everything from here on — dispersion,
#  PCA, sample distances, volcano, MA, DEG heatmap, per-gene boxplots — reads
#  from all_chr/processed/ and writes back under the same tree.


# ─── Section 8 — Load all_chr data ──────────────────────────────────────────
#
#  Load precomputed count matrices + DEG table, build a fresh
#  DESeqDataSet for the all_chr samples so we can call plotDispEsts / vst /
#  plotPCA on it, and (once) compute + save the VST matrix that Part B's
#  PCA / sample-distance / heatmap sections reuse.

# Read DEG table (all genes, with status column already assigned)
deg_all_path <- file.path(allchr_de_tables,
                          "deseq_comparison_KO_vs_CTRL_all_genes.csv")
stopifnot(file.exists(deg_all_path))
deg_all <- read.csv(deg_all_path, stringsAsFactors = FALSE)

# Recode "up_regulated"/"down_regulated"/"not_deg" strings into
# the display factor we'll use everywhere in Part B.
deg_all$deg_status <- factor(
  ifelse(deg_all$status == "up_regulated",   "Up Regulated",
  ifelse(deg_all$status == "down_regulated", "Down Regulated",
                                             "Not DEG")),
  levels = c("Down Regulated", "Not DEG", "Up Regulated")
)
message("all_chr DEG counts (KO vs CTRL):")
print(table(DEG_Status = deg_all$deg_status))

# Read raw counts, median-of-ratios normalized counts, sample metadata
counts_all_raw <- read.csv(file.path(allchr_ge_tables, "raw_gene_expression.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
counts_all_norm <- read.csv(file.path(allchr_ge_tables,
                                      "normalized_gene_expression_medianratios.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
meta_all <- read.csv(file.path(allchr_ge_tables, "sample_metadata.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

# Align columns: keep only the sample columns that appear in metadata, in
# metadata's row order
sample_cols_all <- intersect(colnames(counts_all_raw), meta_all$sample_id)
meta_all <- meta_all[match(sample_cols_all, meta_all$sample_id), , drop = FALSE]
rownames(meta_all) <- meta_all$sample_id

count_mat_all <- as.matrix(counts_all_raw[, sample_cols_all])
rownames(count_mat_all) <- counts_all_raw$gene_id
mode(count_mat_all) <- "integer"

norm_mat_all <- as.matrix(counts_all_norm[, sample_cols_all])
rownames(norm_mat_all) <- counts_all_norm$gene_id
mode(norm_mat_all) <- "numeric"

# Look-up: gene_id → gene_name (for row labels on heatmap / boxplots)
gene_name_map <- setNames(counts_all_raw$gene_name, counts_all_raw$gene_id)

# Set factor levels so CTRL is the reference; palette matches Part A
# (green control, orange KO).
meta_all$condition <- factor(meta_all$condition, levels = c("CTRL", "KO"))
allchr_cond_cols   <- c("CTRL" = "#708238", "KO" = "#F46D43")

# Build DESeqDataSet + prefilter low-count genes (matches cutoff:
# keep genes with ≥10 counts in ≥ smallest-group-size samples).
dds_all <- DESeqDataSetFromMatrix(count_mat_all, colData = meta_all,
                                  design = ~ condition)
smallest_group <- min(table(meta_all$condition))
keep_all <- rowSums(counts(dds_all) >= 10) >= smallest_group
dds_all  <- dds_all[keep_all, ]
message("Genes after prefilter (≥10 counts in ≥",
        smallest_group, " samples): ", nrow(dds_all))

# Fit DESeq2 (needed for plotDispEsts). We don't use its DEG table here —
# is authoritative — but the dispersion fit itself is the
# interesting object for Section 9.
dds_all <- DESeq(dds_all)

# VST for PCA / sample distances (blind = TRUE for exploratory viz)
norm_vst_all <- vst(dds_all, blind = TRUE, nsub = min(1000, nrow(dds_all)))

# Save the VST matrix — the one all_chr GE table missing vs reduced_chr5
vst_mat_all <- as.data.frame(assay(norm_vst_all))
vst_mat_all <- cbind(gene_id = rownames(vst_mat_all), vst_mat_all)
write.csv(vst_mat_all,
          file.path(allchr_ge_tables, "normalized_gene_expression_vst.csv"),
          row.names = FALSE)
message("Wrote all_chr VST table to ", allchr_ge_tables,
        "/normalized_gene_expression_vst.csv (", nrow(vst_mat_all), " genes)")


# ─── Section 9 — Dispersion plot (all_chr) ──────────────────────────────────
#
#  Each dot is one gene. Red curve = mean-dispersion trend DESeq2 fits;
#  blue dots = shrunk final dispersion estimates that get used in testing.
#  Textbook: most blue dots hug the red curve; some outliers up top don't
#  get shrunk (real outliers). With 22k+ genes and n=4 per group, the fit is
#  visibly cleaner than what you'd see on chr5-only Part A.

plotDispEsts(dds_all)

# Base-R plots (plotDispEsts, plotMA) can't be captured cleanly with
# recordPlot() → save_figure_panel in a non-interactive R session, so we open
# a png device directly and re-draw.
save_base_plot_both(
  quote({
    par(cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, mar = c(5, 5, 4, 2))
    plotDispEsts(dds_all,
                 main = "Dispersion estimates — ARID2 KO vs CTRL (all_chr)")
  }),
  output_directory = allchr_ge_figures, figure_name = "dispersion_estimates",
  width = 1000, height = 750, pointsize = 16)


# ─── Section 10 — PCA (all_chr) ─────────────────────────────────────────────
#
#  PCA on VST-transformed counts, top 500 most variable genes. With 4 reps
#  per group we expect PC1 to separate KO from WT and replicates to cluster
#  by genotype.

pcaData    <- plotPCA(norm_vst_all, intgroup = "condition", ntop = 500,
                      returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = allchr_cond_cols) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_label_repel(aes(label = name),
                   size = 5, box.padding = 0.5, point.padding = 0.5,
                   segment.color = "grey50", show.legend = FALSE) +
  coord_fixed() +
  labs(title = "PCA — ARID2 KO vs CTRL (all_chr)",
       subtitle = "top 500 most-variable genes on VST-normalized counts") +
  bump_theme()
print(p_pca)
save_figure_both(p_pca, output_directory = allchr_ge_figures,
                 figure_name = "pca_arid2ko_vs_ctrl",
                 width = 950, height = 800)


# ─── Section 11 — Sample-to-sample distances (all_chr) ──────────────────────
#
#  Euclidean distance on VST counts (top 500 variable genes). Same-genotype
#  pairs should be darker (closer) than cross-genotype pairs.

tmp    <- assay(norm_vst_all)
rv     <- rowVars(tmp)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

sampleDists      <- dist(t(assay(norm_vst_all[select, ])))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(norm_vst_all)
colnames(sampleDistMatrix) <- colnames(norm_vst_all)

p_sample_dist <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         border_color = "grey80",
         fontsize = 14, fontsize_row = 15, fontsize_col = 15,
         main = "Sample-to-sample Euclidean distances (all_chr)\nVST, top 500 variable genes")
save_figure_both(p_sample_dist, output_directory = allchr_ge_figures,
                 figure_name = "sample_to_sample_distances",
                 width = 900, height = 750)


# ─── Section 12 — Shared "top labeled" DEG list (all_chr) ───────────────────
#
#  Volcano, MA, DEG-heatmap annotation, and per-gene boxplots ALL use the
#  same 10 genes: top 5 up + top 5 down by adjusted p-value from
#  DEG table.

deg_cols <- c("Down Regulated" = "#1052bd",
              "Not DEG"        = "#c0c0c0",
              "Up Regulated"   = "#cc212f")

.up  <- deg_all[deg_all$deg_status == "Up Regulated",   , drop = FALSE]
.dn  <- deg_all[deg_all$deg_status == "Down Regulated", , drop = FALSE]
top_up_labeled   <- head(.up[order(.up$padj), ], n = 5)
top_down_labeled <- head(.dn[order(.dn$padj), ], n = 5)
top_labeled      <- rbind(top_down_labeled, top_up_labeled)
rm(.up, .dn)
cat("Genes labeled in volcano / MA / boxplots (all_chr):",
    paste(top_labeled$gene_name, collapse = ", "), "\n")


# ─── Section 13 — Volcano plot (all_chr) ────────────────────────────────────

pvalue_threshold <- 0.05
log2fc_threshold <- 1

volcano_df <- deg_all
volcano_df$neg_log10_pval <- -log10(volcano_df$padj)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_pval,
                                    color = deg_status)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = deg_cols) +
  xlab(expression(log[2]~"Fold Change")) +
  ylab(expression(-log[10]~"adjusted p-value")) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "longdash", color = "grey40") +
  geom_hline(yintercept = -log10(pvalue_threshold),
             linetype = "longdash", color = "grey40") +
  # Force a segment from every labeled point to its label.
  # Tuning notes: bumped `force` (label↔label + label↔point repulsion) and
  # padding way up so gene labels don't overlap with each other or the dense
  # dot cloud; `seed` fixed for reproducible layouts across re-runs.
  geom_text_repel(data = transform(top_labeled,
                                   neg_log10_pval = -log10(padj)),
                  aes(label = gene_name),
                  size = 5, fontface = "bold",
                  min.segment.length = 0,
                  segment.color = "grey40", segment.size = 0.4,
                  box.padding   = 1.2,
                  point.padding = 0.6,
                  force         = 8,
                  force_pull    = 0.5,
                  max.overlaps  = Inf,
                  seed          = 42,
                  show.legend   = FALSE) +
  labs(colour = "DEG status",
       title = "Volcano — ARID2 KO vs CTRL (all_chr)") +
  bump_theme()
print(p_volcano)
save_figure_both(p_volcano, output_directory = allchr_de_figures,
                 figure_name = "volcano_arid2ko_vs_ctrl",
                 width = 1100, height = 800)


# ─── Section 14 — MA plot (all_chr) ─────────────────────────────────────────
#
#  x-axis = baseMean (log₁₀), y-axis = log₂FC. DEGs colored red/blue, non-DEGs
#  grey. Same 10 top-hit genes as volcano are labeled with segments.

ma_df <- deg_all

p_ma <- ggplot(ma_df, aes(x = baseMean, y = log2FoldChange, color = deg_status)) +
  geom_point(alpha = 0.75, size = 3) +
  scale_color_manual(values = deg_cols) +
  scale_x_log10() +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "longdash", color = "dodgerblue") +
  # Stronger repulsion than the volcano — MA panel is thin vertically (LFC
  # clipped to ±4) so labels cluster tightly at ±3 and easily collide. Nudge
  # labels vertically toward the outer edge (top for up-regulated, bottom for
  # down-regulated) so they clear both the dot cloud and each other.
  geom_text_repel(data = top_labeled, aes(label = gene_name),
                  size = 5, fontface = "bold",
                  min.segment.length = 0,
                  segment.color = "grey40", segment.size = 0.4,
                  box.padding   = 1.5,
                  point.padding = 0.7,
                  force         = 15,
                  force_pull    = 0.3,
                  nudge_y       = ifelse(top_labeled$deg_status == "Up Regulated",
                                          0.9, -0.9),
                  direction     = "both",
                  max.overlaps  = Inf,
                  max.iter      = 20000,
                  seed          = 42,
                  show.legend   = FALSE) +
  labs(x = "baseMean (log₁₀)", y = expression(log[2]~"Fold Change"),
       colour = "DEG status",
       title = "MA plot — ARID2 KO vs CTRL (all_chr)") +
  coord_cartesian(ylim = c(-4, 4)) +
  bump_theme()
print(p_ma)
save_figure_both(p_ma, output_directory = allchr_de_figures,
                 figure_name = "ma_arid2ko_vs_ctrl",
                 width = 1000, height = 750)


# ─── Section 15 — DEG heatmap (ComplexHeatmap, split by direction) ─────────
#
#  Row-scaled log2 expression, capped at the top 50 DEGs by padj 
#  (full-genome DEG table has 269 down + 290 up — 559 rows would overflow a
#  slide). ComplexHeatmap adds a row_split by DEG direction so up- vs down-
#  regulated blocks are visually separated, and column annotation shows the
#  KO / CTRL condition of each sample.

sig <- deg_all[deg_all$deg_status != "Not DEG", , drop = FALSE]
sig <- sig[order(sig$padj), ]
top_n_heatmap <- min(50, nrow(sig))
sig_top <- sig[seq_len(top_n_heatmap), ]

# Row-scale log2(counts + 1) so each row's colors show *relative* change
row_mat_ch <- t(scale(t(log2(norm_mat_all[sig_top$gene_id, sample_cols_all,
                                          drop = FALSE] + 1))))
rownames(row_mat_ch) <- sig_top$gene_name

col_ann_ch <- ComplexHeatmap::HeatmapAnnotation(
  Condition = meta_all$condition,
  col = list(Condition = allchr_cond_cols),
  annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"))
row_ann_ch <- ComplexHeatmap::rowAnnotation(
  Direction = sig_top$status,
  col = list(Direction = c(up_regulated   = "#cc212f",
                           down_regulated = "#1052bd")))

ht_top50 <- ComplexHeatmap::Heatmap(
  row_mat_ch,
  name             = "z-score",
  top_annotation   = col_ann_ch,
  right_annotation = row_ann_ch,
  row_split        = sig_top$status,
  cluster_columns  = TRUE,
  cluster_row_slices = FALSE,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#1052bd", "white", "#cc212f")),
  row_names_gp    = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 12),
  column_title    = paste0("Top ", top_n_heatmap,
                           " differentially-expressed genes (padj, all_chr)"))

save_base_plot_both(
  quote(ComplexHeatmap::draw(ht_top50)),
  output_directory = allchr_de_figures,
  figure_name = "complexheatmap_top50_degs",
  width = round(900 * 0.65), height = round(1100 * 0.65))   # 585 × 715


# ─── Section 16 — Circos: genome-wide DEG log2FC (circlize) ─────────────────
#
#  Whole-genome view of where DE regulation lives, plotted on an ideogram.
#  Three tracks: (1) ALL tested genes as dark-gray dots at their log2FC;
#  (2) significant up-regulated (padj<0.05, red); (3) significant
#  down-regulated (padj<0.05, blue).
#
#  Requires mapping ENSG → chromosome coordinates via org.Hs.eg.db + TxDb.
#  Some ENSG IDs have version suffixes (ENSG00000123.14) — strip for lookup.

# ENSG → ENTREZ
eids <- sub("\\..*$", "", deg_all$gene_id)
suppressMessages({
  ens2entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = eids, keytype = "ENSEMBL",
    column = "ENTREZID", multiVals = "first")
})

# ENTREZ → GRanges (single-strand-only avoids multi-mapping duplicates)
txdb    <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
gene_gr <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
mask    <- !is.na(ens2entrez) & ens2entrez %in% names(gene_gr)
gr      <- gene_gr[ens2entrez[mask]]

circos_all <- data.frame(chr   = as.character(GenomicRanges::seqnames(gr)),
                         start = GenomicRanges::start(gr),
                         end   = GenomicRanges::end(gr),
                         value = deg_all$log2FoldChange[mask],
                         sig   = !is.na(deg_all$padj[mask]) &
                                 deg_all$padj[mask] < 0.05)
# Autosomes + X + Y only
autosomes  <- paste0("chr", c(1:22, "X", "Y"))
circos_all <- circos_all[circos_all$chr %in% autosomes, ]
circos_sig <- circos_all[circos_all$sig, ]

# Y-range for the "all genes" track — clip cleanly around ±6
y_all <- c(-6, 6)
draw_circos <- quote({
  circos.clear()
  # Give the top margin extra room so the 2-line title has clearance
  # (default mar=c(1,1,1,1) clips the first row of a "\n"-separated title).
  par(mar = c(1, 1, 5, 1))
  circos.par("track.height" = 0.15, "gap.after" = 2)
  circos.initializeWithIdeogram(species = "hg38",
                                chromosome.index = autosomes)
  # Track 1: all tested genes, dark gray, larger dots
  circos.genomicTrack(circos_all, ylim = y_all,
    panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, pch = 16, cex = 0.45,
                           col = "grey30")
    })
  # Track 2: significant up-regulated (red)
  circos.genomicTrack(circos_sig[circos_sig$value > 0, ],
    ylim = c(0, max(circos_sig$value, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, pch = 16, cex = 0.7,
                           col = "#cc212f")
    })
  # Track 3: significant down-regulated (blue)
  circos.genomicTrack(circos_sig[circos_sig$value < 0, ],
    ylim = c(min(circos_sig$value, na.rm = TRUE), 0),
    panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, pch = 16, cex = 0.7,
                           col = "#1052bd")
    })
  title("Chromosome-level DEG landscape\n(grey = all tested, red = up-in-KO, blue = down-in-KO)")
})
save_base_plot_both(draw_circos,
                    output_directory = allchr_de_figures,
                    figure_name = "circlize_deg_lfc_genome",
                    width = 900, height = 900)


# ─── Section 17 — Per-gene boxplots (with p.adj + LFC overlay) ──────────────
#
#  Same 10 top-hit genes labeled in the volcano and MA plots (top 5 up +
#  top 5 down by padj). Y-axis is log₂ normalized expression. Each panel
#  overlays the DE p.adjust and log₂FoldChange from DESeq2 result
#  so the tested statistics are anchored directly to the figure.

plot_gene_box <- function(g) {
  df <- data.frame(log2_norm_exp = log2(norm_mat_all[g$gene_id,
                                                     sample_cols_all] + 1),
                   condition     = meta_all$condition)
  # Position the annotation just below the top of the y range
  y_top <- max(df$log2_norm_exp, na.rm = TRUE)
  y_rng <- diff(range(df$log2_norm_exp, na.rm = TRUE))
  ann_y <- y_top + y_rng * 0.10   # a bit above the highest point
  arrow <- if (g$deg_status == "Up Regulated") "↑" else "↓"

  ggplot(df, aes(condition, log2_norm_exp, color = condition)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 4) +
    scale_color_manual(values = allchr_cond_cols) +
    # Anchor the DE stats to the panel with an annotation label
    annotate("label",
             x = 1.5, y = ann_y,
             label = paste0("padj = ",
                            format(g$padj, digits = 2, scientific = TRUE),
                            "\nlog2FC = ",
                            format(round(g$log2FoldChange, 2), nsmall = 2)),
             fill = "grey95", color = "black",
             fontface = "bold", size = 4.5,
             label.padding = unit(0.35, "lines")) +
    labs(y     = expression(log[2] ~ "(normalized expression + 1)"),
         title = paste0(arrow, " ", g$gene_name, "  (", g$gene_id, ")")) +
    coord_cartesian(clip = "off") +
    bump_theme() +
    theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))
}

# Filename-safe symbol (replace anything that isn't a letter/digit/_/- with _)
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

for (i in seq_len(nrow(top_labeled))) {
  g <- top_labeled[i, ]
  dir_  <- if (g$deg_status == "Up Regulated") "up" else "down"
  p_box <- plot_gene_box(g)
  print(p_box)
  save_figure_both(p_box, output_directory = allchr_de_figures,
                   figure_name = paste0("boxplot_", dir_, "_", safe_name(g$gene_name)),
                   width = 700, height = 650)
}


# ─── Section 18 — Combined Figure 1 (patchwork: volcano + MA) ───────────────
#
#  Publication-style two-panel figure that stitches the volcano and MA plots
#  side-by-side. Reuses `p_volcano` (Section 13) and `p_ma` (Section 14)
#  built above.

combined_figure <- (p_volcano | p_ma) +
  patchwork::plot_annotation(
    title    = "ARID2 KO vs CTRL — differential expression (all_chr)",
    subtitle = paste0(sum(deg_all$deg_status == "Up Regulated"),
                      " up · ",
                      sum(deg_all$deg_status == "Down Regulated"),
                      " down (padj < 0.05, |log2FC| ≥ 1)"),
    theme    = theme(plot.title    = element_text(size = 18, face = "bold"),
                     plot.subtitle = element_text(size = 13)))
print(combined_figure)
save_figure_both(combined_figure, output_directory = allchr_de_figures,
                 figure_name = "combined_figure_volcano_ma",
                 width = 1800, height = 800)


# ─── Section 19 — Export gene lists (all_chr) ──────────────────────────────
#
#  Write the up- and down-regulated gene names to plain-text files (one
#  symbol per line). These files are the standard input format for downstream
#  functional-analysis tools — session 6 picks them up.

upregulated_genes   <- deg_all$gene_name[deg_all$deg_status == "Up Regulated"]
downregulated_genes <- deg_all$gene_name[deg_all$deg_status == "Down Regulated"]

writeLines(upregulated_genes,   file.path(allchr_de_tables, "upregulated_degs.txt"))
writeLines(downregulated_genes, file.path(allchr_de_tables, "downregulated_degs.txt"))

cat("Up-regulated (",   length(upregulated_genes),   "):\n", sep = "")
cat(paste(head(upregulated_genes,   30), collapse = "\n"), "\n", sep = "")
if (length(upregulated_genes) > 30) cat("... [", length(upregulated_genes) - 30,
                                        " more]\n", sep = "")
cat("\n")
cat("Down-regulated (", length(downregulated_genes), "):\n", sep = "")
cat(paste(head(downregulated_genes, 30), collapse = "\n"), "\n", sep = "")
if (length(downregulated_genes) > 30) cat("... [", length(downregulated_genes) - 30,
                                          " more]\n", sep = "")


# ─── Section 20 — Session info ──────────────────────────────────────────────

sessionInfo()

# End of session 5.
