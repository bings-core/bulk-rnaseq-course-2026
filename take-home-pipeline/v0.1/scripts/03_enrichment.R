# Step 3: functional enrichment, as in Session 6.
#   For each contrast: over-representation (up and down gene lists) and GSEA on the ranked genes, using MSigDB gene sets via msigdbr.
# Run from the project root after step 2, inside the container:
#   source config.sh; singularity exec -B /sc/arion -B /hpc "$R_CONTAINER" Rscript scripts/03_enrichment.R
# Gene sets come from MSIGDB_COLLECTIONS in config.sh: H = Hallmark (default); C2:CP:KEGG = KEGG; C5:GO:BP = GO biological process.

suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr); library(enrichplot); library(ggplot2); library(readr); library(dplyr) })

env <- function(k) Sys.getenv(k)
species   <- env("SPECIES")
results_dir <- env("RESULTS_DIR")
contrasts <- read.csv(env("CONTRASTS_CSV"), stringsAsFactors = FALSE)
de_dir    <- file.path(results_dir, "processed", "differential_expression", "tables")
fa_dir    <- file.path(results_dir, "processed", "functional_analysis")
fig_dir   <- file.path(fa_dir, "figures"); dir.create(file.path(fig_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "pdf"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(fa_dir, "tables"), showWarnings = FALSE)
msig_species <- if (species == "human") "Homo sapiens" else "Mus musculus"

# ---- Section 1: gene sets (msigdbr, offline) --------------------------------------------------
collections <- strsplit(env("MSIGDB_COLLECTIONS"), ",")[[1]]
if (!length(collections) || !nzchar(collections[1])) collections <- "H"
t2g <- bind_rows(lapply(collections, function(cc) {
  parts <- strsplit(trimws(cc), ":")[[1]]
  if (length(parts) == 1) x <- msigdbr(species = msig_species, category = parts[1])
  else x <- msigdbr(species = msig_species, category = parts[1], subcategory = paste(parts[-1], collapse = ":"))
  x[, c("gs_name", "gene_symbol")]
}))
t2g <- distinct(t2g)
cat("gene sets:", length(unique(t2g$gs_name)), "from", paste(collections, collapse = ", "), "\n")
save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, "png", paste0(name, ".png")), p, width = w, height = h, dpi = 200)
  ggsave(file.path(fig_dir, "pdf", paste0(name, ".pdf")), p, width = w, height = h)
}

# ---- Section 2: per contrast, ORA and GSEA ----------------------------------------------------
set.seed(42)
for (i in seq_len(nrow(contrasts))) {
  cn  <- contrasts$name[i]
  tab <- read_csv(file.path(de_dir, paste0(cn, "_all_genes.csv")), show_col_types = FALSE)
  tab <- tab[!is.na(tab$gene_name) & nzchar(tab$gene_name), ]
  universe <- unique(tab$gene_name[!is.na(tab$padj)])
  up   <- unique(tab$gene_name[tab$status == "up_regulated"])
  down <- unique(tab$gene_name[tab$status == "down_regulated"])

  # 2a. over-representation, up and down lists separately
  for (dir in c("up", "down")) {
    genes <- if (dir == "up") up else down
    if (length(genes) < 5) { cat("skip ORA", dir, cn, ": fewer than 5 genes\n"); next }
    ora <- enricher(genes, universe = universe, TERM2GENE = t2g, pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500)
    res <- as.data.frame(ora)
    write_csv(res, file.path(fa_dir, "tables", paste0("ora_", dir, "_", cn, ".csv")))
    if (nrow(res)) save_plot(dotplot(ora, showCategory = 15) + ggtitle(paste("ORA", dir, cn)), paste0("ora_", dir, "_", cn))
  }

  # 2b. GSEA on the ranked list: log2FC * -log10(p), ties broken by a small jitter
  rank_stat <- tab$log2FoldChange * -log10(tab$pvalue)
  rank_stat[!is.finite(rank_stat)] <- 0
  rank_stat <- rank_stat + runif(length(rank_stat), -1e-6, 1e-6)
  gl <- sort(setNames(rank_stat, tab$gene_name), decreasing = TRUE); gl <- gl[!duplicated(names(gl))]   # GSEA needs a decreasing ranking
  gsea <- GSEA(gl, TERM2GENE = t2g, pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, verbose = FALSE, seed = TRUE)
  gres <- as.data.frame(gsea)
  write_csv(gres, file.path(fa_dir, "tables", paste0("gsea_", cn, ".csv")))
  if (nrow(gres)) save_plot(dotplot(gsea, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) + ggtitle(paste("GSEA", cn)), paste0("gsea_", cn))
}
cat("Done. Tables in", file.path(fa_dir, "tables"), "and figures in", file.path(fa_dir, "figures"), "\n")
