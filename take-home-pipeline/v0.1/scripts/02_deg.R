# Step 2: differential expression, as in Session 5.
#   Salmon quant.sf -> tximport -> DESeq2 -> PCA, correlation heatmap, volcano, MA, top-50 heatmaps, DE tables.
# Run from the project root after step 1, inside the container (submit_analysis.sh does this):
#   bash scripts/run_r.sh scripts/02_deg.R
# Settings come from config.sh (through the environment) and from samples.csv / contrasts.csv.

suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(GenomicFeatures); library(AnnotationDbi)
  library(rtracklayer); library(ggplot2); library(ggrepel); library(pheatmap); library(readr); library(dplyr)
})

# ---- Section 1: settings and sample sheets ---------------------------------------------
env <- function(k) Sys.getenv(k)
species      <- env("SPECIES")
resource_dir <- file.path(env("RESOURCE_DIR"), species)
results_dir  <- env("RESULTS_DIR")
samples      <- read.csv(env("SAMPLES_CSV"), stringsAsFactors = FALSE)
contrasts    <- read.csv(env("CONTRASTS_CSV"), stringsAsFactors = FALSE)
filt         <- strsplit(env("SAMPLE_FILTER"), " +")[[1]]
if (length(filt) && any(nzchar(filt))) samples <- samples[samples$sample_id %in% filt, ]
samples$condition <- factor(samples$condition, levels = unique(c(contrasts$reference, samples$condition)))
use_batch <- any(!is.na(samples$batch) & nzchar(as.character(samples$batch)))   # an empty batch column reads as NA
design_formula <- if (use_batch) ~ batch + condition else ~ condition
de_dir  <- file.path(results_dir, "processed", "differential_expression")
exp_dir <- file.path(results_dir, "processed", "gene_expression")
fig_dir <- file.path(de_dir, "figures")
for (d in c(file.path(de_dir, "tables"), exp_dir, file.path(fig_dir, c("png", "pdf")))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 7, h = 6) {          # ggplot objects
  ggsave(file.path(fig_dir, "png", paste0(name, ".png")), p, width = w, height = h, dpi = 200)
  ggsave(file.path(fig_dir, "pdf", paste0(name, ".pdf")), p, width = w, height = h)
}
save_base <- function(draw, name, w = 7, h = 7) {        # pheatmap drawing functions
  png(file.path(fig_dir, "png", paste0(name, ".png")), width = w, height = h, units = "in", res = 200); draw(); dev.off()
  pdf(file.path(fig_dir, "pdf", paste0(name, ".pdf")), width = w, height = h); draw(); dev.off()
}

# One colour per condition, the same palette in the PCA, the correlation heatmap and the contrast heatmaps.
cond_levels <- levels(samples$condition)
pal <- setNames(scales::hue_pal()(length(cond_levels)), cond_levels)
write.csv(data.frame(condition = cond_levels, colour = unname(pal)),
          file.path(de_dir, "tables", "condition_colors.csv"), row.names = FALSE)
TOP_N <- 500

# ---- Section 2: gene names and counts (tximport from Salmon) -------------------------------
db <- loadDb(file.path(resource_dir, "txdb", paste0(
        if (species == "human") "homo_sapiens_grch38_gencode_36" else "mus_musculus_grcm38_gencode_M25", ".sqlite")))
tx2gene <- suppressMessages(AnnotationDbi::select(db, AnnotationDbi::keys(db, "TXNAME"), "GENEID", "TXNAME"))
tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)   # Salmon transcript IDs have no version number
gene_names <- as.data.frame(import(file.path(resource_dir, "gtf", "annotation.gtf"), feature.type = "gene"))
gene_names <- unique(gene_names[, c("gene_id", "gene_name")])
quant_files <- file.path(results_dir, "preprocessed", "salmon", samples$sample_id, "quant.sf")
names(quant_files) <- samples$sample_id
missing <- quant_files[!file.exists(quant_files)]
if (length(missing)) stop("Missing Salmon output: ", paste(missing, collapse = ", "), ". Run step 1 first.")
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene[, c("TXNAME", "GENEID")], ignoreTxVersion = TRUE)
counts <- as.data.frame(txi$counts); counts$gene_id <- rownames(counts)
counts <- left_join(counts, gene_names, by = "gene_id")
write_csv(counts, file.path(exp_dir, "raw_gene_counts.csv"))

# ---- Section 3: DESeq2 model and normalised expression ----------------------------------------
coldata <- data.frame(row.names = samples$sample_id, condition = samples$condition, batch = factor(samples$batch))
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = design_formula)
dds <- dds[rowSums(counts(dds)) >= 5, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
write_csv(as.data.frame(assay(vsd)) |> tibble::rownames_to_column("gene_id") |> left_join(gene_names, by = "gene_id"),
          file.path(exp_dir, "normalized_gene_expression_vst.csv"))

# ---- Section 4: clustering figures (PCA and correlation heatmap, top variable genes) -------------
n_top <- min(TOP_N, nrow(vsd))
pca <- plotPCA(vsd, intgroup = "condition", ntop = n_top, returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
save_plot(ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
            geom_point(size = 3) +
            geom_text_repel(size = 3, box.padding = 0.6, point.padding = 0.5, force = 5,
                            min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
            scale_colour_manual(values = pal) +
            labs(x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"),
                 title = paste0("PCA (vst, top ", n_top, " most variable genes)")) + theme_bw(),
          "pca_vst")

top_var <- order(apply(assay(vsd), 1, var), decreasing = TRUE)[seq_len(n_top)]
cor_mat <- cor(assay(vsd)[top_var, ], method = "pearson")
ann <- data.frame(condition = samples$condition, row.names = samples$sample_id)
ann_colors <- list(condition = pal)
save_base(function() print(pheatmap(cor_mat, annotation_col = ann, annotation_row = ann,
                                    annotation_colors = ann_colors,
                                    main = paste0("Sample correlation (Pearson, vst, top ", n_top, " genes)"))),
          "correlation_heatmap", 7.5, 7.5)

# ---- Section 5: one result table, volcano, MA and top-50 heatmap per contrast ------------------
for (i in seq_len(nrow(contrasts))) {
  cn    <- contrasts$name[i]
  treat <- contrasts$treatment[i]; ref <- contrasts$reference[i]
  res <- results(dds, contrast = c("condition", treat, ref), alpha = 0.05)
  tab <- as.data.frame(res); tab$gene_id <- rownames(tab)
  tab <- left_join(tab, gene_names, by = "gene_id")
  tab$status <- with(tab, ifelse(!is.na(padj) & padj < 0.05 & log2FoldChange >= 1, "up_regulated",
                          ifelse(!is.na(padj) & padj < 0.05 & log2FoldChange <= -1, "down_regulated", "not_deg")))
  write_csv(tab, file.path(de_dir, "tables", paste0(cn, "_all_genes.csv")), na = "NA")
  writeLines(na.omit(tab$gene_name[tab$status == "up_regulated"]), file.path(de_dir, "tables", paste0(cn, "_upregulated.txt")))
  writeLines(na.omit(tab$gene_name[tab$status == "down_regulated"]), file.path(de_dir, "tables", paste0(cn, "_downregulated.txt")))

  # Plot the grey (not DE) points first, so the DE points sit on top.
  plot_tab <- tab[!is.na(tab$padj) & !is.na(tab$log2FoldChange), ]
  plot_tab <- plot_tab[order(plot_tab$status != "not_deg"), ]
  cols <- c(up_regulated = "#d7191c", down_regulated = "#2c7bb6", not_deg = "grey70")
  save_plot(ggplot(plot_tab, aes(log2FoldChange, -log10(padj), colour = status)) + geom_point(size = 0.8) +
            scale_colour_manual(values = cols) +
            labs(title = paste("Volcano:", cn), y = "-log10 adjusted p-value (padj)") + theme_bw(), paste0("volcano_", cn))
  save_plot(ggplot(plot_tab, aes(log10(baseMean + 1), log2FoldChange, colour = status)) + geom_point(size = 0.8) +
            scale_colour_manual(values = cols) + labs(title = paste("MA:", cn)) + theme_bw(), paste0("ma_", cn))

  # Heatmap: only the two conditions, reference samples first, then treatment samples.
  cols_ref <- samples$sample_id[samples$condition == ref]
  cols_trt <- samples$sample_id[samples$condition == treat]
  sample_order <- c(cols_ref, cols_trt)
  top <- head(tab[order(tab$padj), ], 50)
  top <- top[top$gene_id %in% rownames(vsd), ]
  m <- assay(vsd)[top$gene_id, sample_order, drop = FALSE]
  rownames(m) <- top$gene_name[match(rownames(m), top$gene_id)]
  keep <- !duplicated(rownames(m)); m <- m[keep, , drop = FALSE]
  m <- m - rowMeans(m)
  col_ann <- data.frame(condition = factor(samples$condition[match(sample_order, samples$sample_id)], levels = c(ref, treat)),
                        row.names = sample_order)
  save_base(function() print(pheatmap(m, cluster_cols = FALSE, annotation_col = col_ann,
                                      annotation_colors = ann_colors,
                                      main = paste("Top 50 DEGs:", treat, "vs", ref))),
            paste0("heatmap_top50_", cn), 7, 9)
}
cat("Done. Tables in", de_dir, "and figures in", fig_dir, "\n")
