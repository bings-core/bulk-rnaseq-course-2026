# =============================================================================
#  BiNGS Bulk RNA-seq Course — Session 6
#  Functional analysis: from a DEG list to biological meaning
# =============================================================================
#
#
# ─── Section 0 — Session overview ────────────────────────────────────────────
#
#  Session 5 produced a table of differentially expressed genes (DEGs).
#  Session 6 asks: "so what?" — do these genes have a shared BIOLOGICAL theme?
#  This is called functional analysis / pathway enrichment.
#
#  We use two complementary methods that are the standard everywhere in the
#  field:
#
#    (1) ORA — Over-Representation Analysis
#        Take a fixed list (your up- or down-regulated genes). Ask: which
#        pre-defined gene sets (pathways) contain MORE of these genes than
#        you'd expect by chance?
#        - Statistical test: hypergeometric (Fisher's exact).
#        - Input: a THRESHOLDED gene list (e.g. padj < 0.05 & |log2FC| ≥ 1).
#
#    (2) GSEA — Gene Set Enrichment Analysis
#        Rank ALL genes by their DE score (not just DEGs). Ask: does a given
#        pathway's genes tend to be near the top (or bottom) of the ranked
#        list?
#        - Statistical test: Kolmogorov–Smirnov-style running sum.
#        - Input: a RANKED gene list (all genes, no threshold).
#
#  Same DEG table, different questions, complementary answers.
#
#  GENE SETS: MSigDB collections, fetched at run time via the msigdbr package.
#    H          Hallmark  (50 curated pathways — start here)
#    C2 KEGG    KEGG canonical pathways
#    C2 REACTOME  Reactome canonical pathways
#    C2 WIKIPATHWAYS  Community-curated
#    C5 GO_BP   Gene Ontology — Biological Process
#    C5 GO_CC   Gene Ontology — Cellular Component
#    C5 GO_MF   Gene Ontology — Molecular Function


# ─── Section 1 — Environment ────────────────────────────────────────────────

# On a headless compute node R has no X11 display, so ANY base-R plot call
# would open the fallback pdf() device and dump into "Rplots.pdf" in cwd.
# We always route that fallback to /dev/null — every real save goes through
# explicit png()/pdf() opens inside save_figure_both / save_base_plot_both,
# so nothing "user-visible" is suppressed. Guarantees no stray Rplots.pdf
# whether the script is sourced interactively OR run via `R -f`.
pdf(NULL)

# Layout convention:
#   data_rna/reduced_chr5/raw/                  — chr5 fastqs + metadata (session-4 input)
#   data_rna/reduced_chr5/preprocessed/         — session-4 outputs
#   data_rna/reduced_chr5/processed/            — session-5 outputs
#   data_rna/all_chr/processed/                 — full-genome inputs (from fridls01
#                                                  upstream pipeline) AND session-6 outputs,
#                                                  co-located under gene_expression/,
#                                                  differential_expression/, functional_analysis/

user <- Sys.getenv("USER")
message("Running as user: ", user)
# data_dir lives inside the participant's own hands_on/ clone — same tree as
# ref_dir. Session 6 both READS fridls01 inputs from all_chr/processed/ AND
# WRITES its own functional-analysis outputs there.
data_dir <- file.path(getwd(), "data_rna")

# The shared reference tree — has the FULL-genome ARID2 KO vs CTRL DEG,
# gene-expression counts, and sample metadata from fridls01's real experiment.
# Lives under data_rna/all_chr/processed/ per the new convention (these files
# were already processed upstream by fridls01's pipeline; from the course's
# perspective they are inputs, but stored alongside our own processed outputs).
# Convention: run this script from the participant's OWN hands_on/ dir, i.e.
#   cd /sc/arion/projects/BiNGS_bulk/$USER/hands_on
#   ml singularity/3.6.4
#   singularity exec engine/singularity_containers/rbings_20250925.sif R
#   source("code_rna/session_06_functional_analysis.R") 
# ref_dir is the current working directory (fall back to the shared path if
# the script was launched from elsewhere).
.hands_on_default <- "/sc/arion/projects/BiNGS/bings_omics/data/bings/2026/bulk_course_data/hands_on"
ref_dir     <- if (dir.exists(file.path(getwd(), "engine", "annotation"))) getwd() else .hands_on_default
input_root  <- file.path(ref_dir, "data_rna", "all_chr", "processed")

deg_path         <- file.path(input_root, "differential_expression", "tables",
                              "deseq_comparison_KO_vs_CTRL_all_genes.csv")
counts_raw_path  <- file.path(input_root, "gene_expression", "tables",
                              "raw_gene_expression.csv")
counts_norm_path <- file.path(input_root, "gene_expression", "tables",
                              "normalized_gene_expression_medianratios.csv")
sample_meta_path <- file.path(input_root, "gene_expression", "tables",
                              "sample_metadata.csv")
stopifnot(all(file.exists(c(deg_path, counts_raw_path, counts_norm_path, sample_meta_path))))
message("Full-genome inputs read from: ", input_root)

# Two output areas under all_chr/processed/, each with tables/ + figures/
ge_dir  <- file.path(data_dir, "all_chr", "processed", "gene_expression")
fa_dir  <- file.path(data_dir, "all_chr", "processed", "functional_analysis")
ge_tables  <- file.path(ge_dir, "tables")
ge_figures <- file.path(ge_dir, "figures")
fa_tables  <- file.path(fa_dir, "tables")
fa_figures <- file.path(fa_dir, "figures")
for (d in c(ge_tables, ge_figures, fa_tables, fa_figures)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
message("Gene-expression outputs   → ", ge_dir, " (tables/, figures/)")
message("Functional-analysis output → ", fa_dir, " (tables/, figures/)")


# ─── Section 2 — Libraries ──────────────────────────────────────────────────

library(clusterProfiler)   # GSEA + ORA wrappers
library(msigdbr)           # MSigDB gene sets, downloaded once by the package
library(enrichplot)        # dotplot / gseaplot2 / cnetplot for enrichment results
library(org.Hs.eg.db)      # human gene ID annotation (used implicitly by enricher)
library(DOSE)              # helpers used by enrichplot

library(ggplot2)           # plotting
library(ggrepel)           # non-overlapping labels
library(dplyr)             # tidy data wrangling
library(tidyr)             # pivot_wider / pivot_longer
library(forcats)           # reorder factors for barplots
library(stringr)           # string manipulation for tidy plot labels
library(scales)            # comma(), pretty axis breaks
library(pheatmap)          # heatmaps
library(RColorBrewer)      # color palettes


# ─── Section 2.5 — save_figure_panel helper ─────────────────────────────────
#  Copied verbatim from /sc/arion/projects/BiNGS/bings_analysis/code/R/visualization/utils.R:94
#  Every ggplot / pheatmap below saves via this into
#      <figures_dir>/<output_type>/<figure_name>.<ext>

# Convenience helper: save the same figure as BOTH png (raster) and pdf
# (vector) at IDENTICAL physical dimensions. See session_05 for the full
# rationale — the crucial bit is res=dpi on the png call, so a 500-px-wide
# png at dpi=100 is physically 5 in, matching the 5-in-wide pdf. Without
# res=dpi the png defaults to 72 dpi and ends up physically bigger, which
# makes point-sized text look proportionally smaller in the png / bigger
# in the pdf.
save_figure_both <- function(figure_panel, output_directory, figure_name,
                             width, height, dpi = 100, ...) {
  save_figure_panel(figure_panel, output_directory, figure_name,
                    output_type = "png",
                    width = width, height = height, res = dpi, ...)
  save_figure_panel(figure_panel, output_directory, figure_name,
                    output_type = "pdf",
                    width = width / dpi, height = height / dpi, ...)
  message("  ↳ saved  ", file.path(output_directory, "png", paste0(figure_name, ".png")))
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
                png(filename = figure_path, bg = background, ...); print(figure_panel); dev.off()
            } else if (output_type %in% c("jpg", "jpeg")) {
                jpeg(filename = figure_path, bg = background, ...); print(figure_panel); dev.off()
            } else if (output_type == "eps") {
                cairo_ps(filename = figure_path, bg = background, ...); print(figure_panel); dev.off()
            } else if (output_type == "svg") {
                svg(filename = figure_path, bg = background, ...); print(figure_panel); dev.off()
            } else if (output_type == "pdf") {
                pdf(file = figure_path, title = title, bg = background, ...); print(figure_panel); dev.off()
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


# ─── Section 3 — Load DEG + gene expression + sample metadata ───────────────
#
#  This is the SAME biology as session 5 (SKmel147 ARID2 KO vs WT/CTRL) but
#  on the FULL genome (22,717 tested genes / 60,240 total gene models) — from
#  fridls01's completed pipeline run. Session 5 used chr5-only reads for
#  speed; here we use the real thing so the enrichment numbers are meaningful.

deg <- read.csv(deg_path, stringsAsFactors = FALSE)

# Normalized counts (median-of-ratios) — used later for heatmap of leading-edge genes
counts_norm <- read.csv(counts_norm_path, stringsAsFactors = FALSE)

# Sample metadata: sample_id, replicate, condition (KO / CTRL), etc.
sample_meta <- read.csv(sample_meta_path, stringsAsFactors = FALSE)

# ── Explore ──
cat("DEG rows:                ", nrow(deg), "\n")
cat("Gene-expression rows:    ", nrow(counts_norm), "\n")
cat("Samples:                 ", nrow(sample_meta), "\n\n")
print(table(deg$status))
cat("\nSample metadata (subset):\n")
print(sample_meta[, c("sample_id","replicate","condition")])

cat("\nTop 5 up-regulated by padj:\n")
print(deg |> filter(status == "up_regulated") |> arrange(padj) |> head(5) |>
              select(gene_name, log2FoldChange, padj))
cat("\nTop 5 down-regulated by padj:\n")
print(deg |> filter(status == "down_regulated") |> arrange(padj) |> head(5) |>
              select(gene_name, log2FoldChange, padj))


# ─── Section 4 — Prepare enrichment inputs ──────────────────────────────────
#
#  ORA needs:  a HIT LIST (up-regulated symbols; down-regulated symbols)
#              a UNIVERSE (all tested gene symbols) — the background for the
#              hypergeometric test.
#
#  GSEA needs: a RANKED gene list = named numeric vector, decreasing.
#              We use "combined score" = sign(log2FC) × -log10(padj) — this
#              weights BOTH direction and confidence and matches what the
#              BiNGS pipeline (B3) uses.

# Deduplicate by gene name (keep the row with the smallest padj if any duplicates)
deg_dedup <- deg |>
  filter(!is.na(gene_name), gene_name != "", !is.na(padj), !is.na(log2FoldChange)) |>
  arrange(padj) |>
  distinct(gene_name, .keep_all = TRUE)

# --- Ranked list for GSEA ---
ranked <- deg_dedup |>
  mutate(combined_score = sign(log2FoldChange) *
                          pmin(-log10(pmax(padj, 1e-300)), 300)) |>
  arrange(desc(combined_score))
gene_list <- setNames(ranked$combined_score, ranked$gene_name)
cat("Ranked list length:", length(gene_list), "\n")
cat("Head:\n"); print(head(gene_list, 5))
cat("Tail:\n"); print(tail(gene_list, 5))

# --- Hit lists + universe for ORA ---
universe_genes <- deg_dedup$gene_name
up_genes       <- deg_dedup$gene_name[deg_dedup$status == "up_regulated"]
down_genes     <- deg_dedup$gene_name[deg_dedup$status == "down_regulated"]
cat("\nUniverse:", length(universe_genes),
    "| up:", length(up_genes),
    "| down:", length(down_genes), "\n")


# ─── Section 5 — MSigDB collections via msigdbr ─────────────────────────────
#
#  msigdbr::msigdbr(species, category, subcategory) returns a tidy data.frame
#  of gene sets (long format: one row per gene per set). We reduce to a
#  TERM2GENE map (columns: gs_name, gene_symbol) which is the input format
#  clusterProfiler wants.

collections <- list(
  "Hallmark (H)"            = list(category = "H",  subcategory = NULL),
  "KEGG (C2 CP)"            = list(category = "C2", subcategory = "CP:KEGG"),
  "Reactome (C2 CP)"        = list(category = "C2", subcategory = "CP:REACTOME"),
  "WikiPathways (C2 CP)"    = list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
  "GO — BP (C5)"            = list(category = "C5", subcategory = "GO:BP"),
  "GO — CC (C5)"            = list(category = "C5", subcategory = "GO:CC"),
  "GO — MF (C5)"            = list(category = "C5", subcategory = "GO:MF")
)

fetch_t2g <- function(x) {
  args <- c(list(species = "Homo sapiens", category = x$category),
            if (!is.null(x$subcategory)) list(subcategory = x$subcategory))
  m <- do.call(msigdbr::msigdbr, args)
  m |> dplyr::select(gs_name, gene_symbol) |> distinct()
}

t2g_all <- lapply(collections, fetch_t2g)

# ── Explore ──
data.frame(
  Collection = names(t2g_all),
  Pathways   = sapply(t2g_all, \(x) length(unique(x$gs_name))),
  Total_gene_memberships = sapply(t2g_all, nrow)
) |> print()


# ─── Section 6 — GSEA — Hallmark first, then loop over all collections ──────

set.seed(42)
run_gsea <- function(t2g, gene_list) {
  clusterProfiler::GSEA(
    geneList     = gene_list,
    TERM2GENE    = t2g,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 1,          # keep all pathways; filter later by is_significant
    eps          = 0,          # exact p-values (slower but accurate)
    verbose      = FALSE,
    seed         = TRUE
  )
}

gsea_results <- lapply(t2g_all, run_gsea, gene_list = gene_list)

# ── Explore — top hits per collection ──
for (cname in names(gsea_results)) {
  cat("\n=== GSEA:", cname, "===\n")
  d <- as.data.frame(gsea_results[[cname]])
  if (nrow(d) == 0) { cat("  (none)\n"); next }
  print(d |> arrange(p.adjust) |> head(6) |>
              select(ID, NES, p.adjust, setSize))
}

# Save every collection's GSEA table
for (cname in names(gsea_results)) {
  safe_c <- gsub("[^A-Za-z0-9]+", "_", cname)
  write.csv(as.data.frame(gsea_results[[cname]]),
            file.path(fa_tables, paste0("gsea_", safe_c, ".csv")),
            row.names = FALSE)
}


# ─── Section 7 — ORA — up-regulated and down-regulated across all collections

run_ora <- function(hits, universe, t2g) {
  clusterProfiler::enricher(
    gene         = hits,
    universe     = universe,
    TERM2GENE    = t2g,
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

ora_up   <- lapply(t2g_all, run_ora, hits = up_genes,   universe = universe_genes)
ora_down <- lapply(t2g_all, run_ora, hits = down_genes, universe = universe_genes)

# ── Explore — top hits per collection ──
for (cname in names(ora_up)) {
  cat("\n=== ORA up:", cname, "===\n")
  d <- if (is.null(ora_up[[cname]])) NULL else as.data.frame(ora_up[[cname]])
  if (is.null(d) || nrow(d) == 0) { cat("  (none)\n"); next }
  print(d |> arrange(p.adjust) |> head(5) |>
              select(ID, Count, GeneRatio, p.adjust))
}

# Save every collection's ORA tables
for (cname in names(ora_up)) {
  safe_c <- gsub("[^A-Za-z0-9]+", "_", cname)
  for (dir_ in c("up", "down")) {
    obj <- if (dir_ == "up") ora_up[[cname]] else ora_down[[cname]]
    if (is.null(obj)) next
    write.csv(as.data.frame(obj),
              file.path(fa_tables, paste0("ora_", dir_, "_", safe_c, ".csv")),
              row.names = FALSE)
  }
}


# ─── Section 8 — GSEA visualization: dotplot + barplot per collection ───────
#
#  Two complementary views of the top 10 pathways per direction (activated =
#  NES > 0 = up in KO; suppressed = NES < 0 = up in CTRL).
#
#  For GSEA the right axes are:
#     x     = NES         (Normalized Enrichment Score, signed)
#     size  = setSize     (dotplot only — pathway gene count)
#     color = p.adjust    (BH-adjusted p-value, clamped to [0, 0.2] so
#                          significance differences visible near the cutoff)

top_n_gsea <- 10

# Shared p.adjust color scale for all functional-analysis plots (0 – 0.2)
# Shared p.adjust palette + LEGEND BREAKS so GSEA / ORA / any downstream plot
# using either the color or fill aesthetic gets a visually identical legend
# (same tick positions, same endpoint colors, same clamping).
.pval_breaks <- seq(0, 0.2, by = 0.05)
.pval_guide  <- guide_colorbar(barheight = unit(55, "mm"),
                               barwidth  = unit(5,  "mm"),
                               ticks.colour = "black",
                               frame.colour = "black")

pval_color_scale <- function() {
  scale_color_gradient(low = "#cc212f", high = "#3182bd",
                       limits = c(0, 0.2), oob = scales::squish,
                       breaks = .pval_breaks, labels = .pval_breaks,
                       guide = .pval_guide,
                       name = "p.adjust")
}
pval_fill_scale <- function() {
  scale_fill_gradient(low = "#cc212f", high = "#3182bd",
                      limits = c(0, 0.2), oob = scales::squish,
                      breaks = .pval_breaks, labels = .pval_breaks,
                      guide = .pval_guide,
                      name = "p.adjust")
}

# Shared enlarged-text theme for all functional-analysis bar/dot plots.
# theme_linedraw(base_size=14) scales axis text ~14 pt, then we bump the title
# and pathway-label sizes further for slide-legible output.
fa_theme <- function(base_size = 14) {
  theme_linedraw(base_size = base_size) +
    theme(
      plot.title        = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle     = element_text(size = base_size),
      axis.title        = element_text(size = base_size + 1),
      axis.text.y       = element_text(size = base_size + 1),
      axis.text.x       = element_text(size = base_size),
      legend.title      = element_text(size = base_size),
      legend.text       = element_text(size = base_size - 1),
      # Facet-strip labels ("suppressed" / "activated"): force BLACK BOLD text
      # on a light grey background so they are legible on a slide.
      strip.text        = element_text(size = base_size + 2, face = "bold", color = "black"),
      strip.background  = element_rect(fill = "grey85", color = "grey40"),
      plot.title.position = "plot"
    )
}

# GSEA barplot only — x = NES, sorted by NES, filled by p.adjust.
# (GSEA dotplot removed per user request; barplot conveys the same info more
# compactly and stays consistent with the ORA / cnetplot styling below.)
plot_gsea_barplot <- function(gsea_obj, cname, top_n = 10) {
  d <- as.data.frame(gsea_obj)
  if (nrow(d) == 0) return(NULL)
  d_up   <- d |> filter(NES > 0) |> arrange(p.adjust) |> head(top_n)
  d_down <- d |> filter(NES < 0) |> arrange(p.adjust) |> head(top_n)
  d2 <- bind_rows(d_up, d_down)
  if (nrow(d2) == 0) return(NULL)
  d2 <- d2 |>
    mutate(label = str_trunc(gsub("_", " ",
                                  sub("^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOCC_|^GOMF_", "", ID)),
                             55)) |>
    arrange(NES) |>
    mutate(label = factor(label, levels = unique(label)))
  ggplot(d2, aes(x = NES, y = label, fill = p.adjust)) +
    geom_col(color = "black", size = 0.25) +
    pval_fill_scale() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    labs(x = "Normalized Enrichment Score (NES)  →  positive = up in KO",
         y = NULL,
         title = paste0("GSEA — ", cname),
         subtitle = paste0("top ", top_n, " pathways per direction, sorted by NES; ",
                           "fill = p.adjust (0–0.2)")) +
    fa_theme()
}

for (cname in names(gsea_results)) {
  safe_c <- gsub("[^A-Za-z0-9]+", "_", cname)
  p_bar <- plot_gsea_barplot(gsea_results[[cname]], cname, top_n = top_n_gsea)
  if (is.null(p_bar)) next
  print(p_bar)
  save_figure_both(p_bar, output_directory = fa_figures,
                   figure_name = paste0("gsea_barplot_", safe_c),
                   width = 1950, height = 750)     # 650 × 1.5
}


# ─── Section 9 — GSEA visualization: running enrichment (gseaplot2) ─────────
#
#  For the SINGLE top pathway from Hallmark, plot the classic GSEA running
#  enrichment curve — the visual signature of the algorithm.

hallmark_gsea <- gsea_results[["Hallmark (H)"]]
top_hallmark  <- as.data.frame(hallmark_gsea) |> arrange(p.adjust) |> head(1)
cat("Top Hallmark pathway:", top_hallmark$ID, "\n")

# Plot a GSEA running-enrichment curve for a specific Hallmark pathway ID.
# pvalue_table = TRUE would show 3 cols (pvalue, p.adjust, NES). We only want
# NES + p.adjust — so turn that table off and put those two numbers in the
# title instead.
plot_gsea_running <- function(gsea_obj, pathway_id) {
  row <- as.data.frame(gsea_obj)
  row <- row[row$ID == pathway_id, , drop = FALSE]
  if (nrow(row) == 0) {
    message("  gseaplot2: pathway '", pathway_id, "' not found in results — skipping.")
    return(invisible(NULL))
  }
  title <- paste0(
    row$ID,
    "   |   NES = ",     formatC(row$NES,      format = "f", digits = 2),
    "   p.adjust = ",    formatC(row$p.adjust, format = "e", digits = 2)
  )
  p <- enrichplot::gseaplot2(gsea_obj, geneSetID = pathway_id,
                             title = title,
                             pvalue_table = FALSE,
                             base_size = 15)
  print(p)
  save_figure_both(p, output_directory = fa_figures,
                   figure_name = paste0("gsea_running_enrichment_",
                                        gsub("[^A-Za-z0-9]+","_", pathway_id)),
                   width = 1650, height = 850)      # 550 × 1.5
  invisible(p)
}

# Two Hallmark pathways of interest, both plotted by explicit name so the
# figure files are named after the biology (not "the current top pathway"):
#   - MYC_TARGETS_V1           — proliferation program (top by p.adjust on this dataset)
#   - INTERFERON_GAMMA_RESPONSE — canonical ARID2-KO signal via HLA session II derepression
plot_gsea_running(hallmark_gsea, "HALLMARK_MYC_TARGETS_V1")
plot_gsea_running(hallmark_gsea, "HALLMARK_INTERFERON_GAMMA_RESPONSE")


# ─── Section 10 — ORA visualization: dotplot + barplot ──────────────────────
#
#  ORA dotplot: dot size = Count (overlap), color = p.adjust, y-axis = pathway.
#  ORA barplot: bar length = Count, split by adjusted p-value bucket.

# Build the ORA dotplot from scratch (not enrichplot::dotplot) so we can
# fully control the p.adjust color scale AND its legend — enrichplot's
# internal scale replaces guide= silently, giving inconsistent legends
# across GSEA barplot vs ORA dotplot.
plot_ora_dot <- function(obj, cname, direction, top_n = 12) {
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df <- df |> arrange(p.adjust) |> head(top_n)
  # Parse "16/223" GeneRatio → 0.0717 numeric
  gr <- strsplit(df$GeneRatio, "/")
  df$GeneRatioNum <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df$label <- str_trunc(gsub("_", " ", sub("^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOCC_|^GOMF_", "", df$ID)), 55)
  # Order y-axis by GeneRatio (largest at top). factor() needs unique levels
  # — two pathways can truncate to the same label after str_trunc, so unique().
  df$label <- factor(df$label,
                     levels = unique(df$label[order(df$GeneRatioNum)]))

  ggplot(df, aes(x = GeneRatioNum, y = label, size = Count, color = p.adjust)) +
    geom_point() +
    pval_color_scale() +
    scale_size_continuous(range = c(3, 10)) +
    labs(x = "GeneRatio", y = NULL,
         title = paste0("ORA — ", cname, "  (", direction, "-regulated in KO)"),
         subtitle = "dot size = overlap (Count), color = p.adjust (0–0.2)") +
    fa_theme()
}

for (cname in names(ora_up)) {
  safe_c <- gsub("[^A-Za-z0-9]+", "_", cname)
  for (dir_ in c("up", "down")) {
    obj <- if (dir_ == "up") ora_up[[cname]] else ora_down[[cname]]
    if (is.null(obj) || nrow(as.data.frame(obj)) == 0) next
    p <- plot_ora_dot(obj, cname, dir_)
    if (is.null(p)) next
    print(p)
    save_figure_both(p, output_directory = fa_figures,
                     figure_name = paste0("ora_dotplot_", dir_, "_", safe_c),
                     width = 1876, height = 750)   # 625 × 1.5
  }
}


# ─── Section 11 — Gene-concept network for the top Hallmark ORA hit ─────────
#
#  cnetplot shows the WHICH GENES drive WHICH pathways. Great for spotting
#  hub genes (a single gene appearing in many enriched pathways).

# One log2FC vector shared by both directions — colors gene nodes in the network.
fc_vec <- setNames(deg_dedup$log2FoldChange, deg_dedup$gene_name)

make_cnet <- function(ora_obj, direction_label) {
  p <- tryCatch(
    enrichplot::cnetplot(ora_obj, showCategory = 5, foldChange = fc_vec,
                         node_label = "all",
                         cex_category       = 3.0,   # pathway (yellow) node dots — much bigger
                         cex_gene           = 2.5,   # gene node dots            — much bigger
                         cex_label_category = 2.6,   # pathway labels
                         cex_label_gene     = 2.0) + # gene names
      # Diverging blue-white-red palette centered at 0, clamped to [-2, 2]
      scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                            midpoint = 0, limits = c(-2, 2),
                            oob = scales::squish,
                            name = expression(log[2]~FC)) +
      labs(title = paste0("Hallmark ORA (", direction_label,
                          " in KO) — gene-concept network"),
           subtitle = "gene node color = log2FC (KO vs CTRL, clamped to ±2); pathway nodes in gold") +
      theme(plot.title    = element_text(size = 34, face = "bold"),
            plot.subtitle = element_text(size = 22),
            legend.title  = element_text(size = 22, face = "bold"),
            legend.text   = element_text(size = 20),
            legend.key.size = unit(1.4, "cm")),
    error = function(e) { message("cnetplot failed (", direction_label, "): ",
                                  conditionMessage(e)); NULL }
  )
  if (is.null(p)) return(NULL)
  # enrichplot::cnetplot has no `fontface` arg for pathway-node labels — post-
  # hoc bold every text layer (gene labels get bold too, which reads fine).
  for (i in seq_along(p$layers)) {
    l <- p$layers[[i]]
    if (inherits(l$geom, "GeomText") ||
        inherits(l$geom, "GeomTextRepel") ||
        inherits(l$geom, "GeomLabel")) {
      p$layers[[i]]$aes_params$fontface <- "bold"
    }
  }
  p
}

# Both directions — Hallmark up-regulated AND down-regulated
for (dir_ in c("up", "down")) {
  ora_obj <- if (dir_ == "up") ora_up[["Hallmark (H)"]]
             else               ora_down[["Hallmark (H)"]]
  if (is.null(ora_obj) || nrow(as.data.frame(ora_obj)) == 0) next
  p_cnet <- make_cnet(ora_obj, dir_)
  if (is.null(p_cnet)) next
  print(p_cnet)
  save_figure_both(p_cnet, output_directory = fa_figures,
                   figure_name = paste0("ora_hallmark_", dir_, "_cnetplot"),
                   width = 2476, height = 1350)   # 825 × 1.5, 450 × 1.5
}


# ─── Section 12 — Expression heatmap of leading-edge genes ──────────────────
#
#  Bridge from PATHWAYS back to raw EXPRESSION. Take the top Hallmark GSEA
#  pathway, extract its leading-edge (core-enrichment) genes, and plot their
#  z-scored normalized expression across the 8 samples. Rows should cluster
#  cleanly into two groups (KO vs CTRL) if the pathway really drives the
#  observed transcriptional program.

# Helper — plot LE heatmap for a specific pathway ID from any GSEA result
sample_cols <- setdiff(colnames(counts_norm), c("gene_id", "gene_name"))
col_ann <- data.frame(Condition = sample_meta$condition[match(sample_cols,
                                                              sample_meta$sample_id)])
rownames(col_ann) <- sample_cols

plot_leading_edge_heatmap <- function(gsea_obj, pathway_id) {
  row <- as.data.frame(gsea_obj)
  row <- row[row$ID == pathway_id, , drop = FALSE]
  if (nrow(row) == 0) {
    message("  leading-edge heatmap: pathway '", pathway_id,
            "' not found — skipping.")
    return(invisible(NULL))
  }
  le_genes <- strsplit(row$core_enrichment, "/")[[1]]
  cat("Leading-edge genes for", pathway_id, ":", length(le_genes), "\n")
  cat("  First 10:", paste(head(le_genes, 10), collapse = ", "), "\n")

  le_mat <- counts_norm |>
    filter(gene_name %in% le_genes) |>
    distinct(gene_name, .keep_all = TRUE)
  rownames(le_mat) <- le_mat$gene_name
  le_mat <- as.matrix(le_mat[, sample_cols])
  mode(le_mat) <- "numeric"
  le_mat <- log2(le_mat + 1)     # log2(x+1) then row-scale in pheatmap

  p_le_hm <- pheatmap(le_mat, scale = "row",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      annotation_col = col_ann,
                      annotation_colors = list(Condition = c(KO   = "#F46D43",
                                                             CTRL = "#708238")),
                      show_rownames = nrow(le_mat) < 80,
                      fontsize = 14, fontsize_row = 9, fontsize_col = 14,
                      border_color = NA,
                      main = paste0(gsub("_", " ", pathway_id), "\n",
                                    "Leading edge genes (", nrow(le_mat), ")\n",
                                    "log2 normalized counts"))
  save_figure_both(p_le_hm, output_directory = fa_figures,
                   figure_name = paste0("leading_edge_heatmap_",
                                        gsub("[^A-Za-z0-9]+","_", pathway_id)),
                   width = 700, height = 900)
  invisible(p_le_hm)
}

# Both pathways — MYC targets (top overall) + INTERFERON_GAMMA_RESPONSE
plot_leading_edge_heatmap(hallmark_gsea, "HALLMARK_MYC_TARGETS_V1")
plot_leading_edge_heatmap(hallmark_gsea, "HALLMARK_INTERFERON_GAMMA_RESPONSE")



# ─── Section 13 — Session info ──────────────────────────────────────────────

sessionInfo()

# End of session 6.
