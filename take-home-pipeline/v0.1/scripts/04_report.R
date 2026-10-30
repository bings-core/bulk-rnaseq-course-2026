# Step 4: one self-contained HTML report (figures and the MultiQC report embedded, no external files) with the main results.
#   Run through submit_analysis.sh, or: bash scripts/run_r.sh scripts/04_report.R
# Output: results/report/rnaseq_report.html   (open it in any browser; copy the single file to share it)

suppressPackageStartupMessages({ library(base64enc); library(readr); library(dplyr); library(htmltools); library(ggplot2) })

env <- function(k) Sys.getenv(k)
results_dir <- env("RESULTS_DIR"); project <- env("PROJECT_NAME"); species <- env("SPECIES")
samples     <- read.csv(env("SAMPLES_CSV"), stringsAsFactors = FALSE)
contrasts   <- read.csv(env("CONTRASTS_CSV"), stringsAsFactors = FALSE)
pre  <- file.path(results_dir, "preprocessed")
proc <- file.path(results_dir, "processed")
de   <- file.path(proc, "differential_expression", "figures", "png")
fa   <- file.path(proc, "functional_analysis")
rep_dir <- file.path(results_dir, "report"); dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Section 1: helpers -------------------------------------------------------------------
img <- function(path, caption) {
  if (!file.exists(path)) return(tags$p(tags$em(paste("Missing figure:", basename(path)))))
  tags$figure(tags$img(src = paste0("data:image/png;base64,", base64encode(path)), style = "max-width:100%;"),
              tags$figcaption(caption))
}
tbl <- function(df, digits = 3) {
  if (!nrow(df)) return(tags$p(tags$em("No rows.")))
  df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, digits) else x)
  tags$table(class = "t", tags$thead(tags$tr(lapply(names(df), tags$th))),
             tags$tbody(lapply(seq_len(nrow(df)), function(i) tags$tr(lapply(df[i, ], function(v) tags$td(as.character(v)))))))
}
read_or_empty <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE) else data.frame()
# Top 10 enrichment terms. ORA: dot plot (x = gene ratio, dot size = gene set size, colour = adjusted p).
# GSEA: NES bars for the 10 terms with the largest |NES|, highest NES at the top (red = enriched in treatment, blue = in reference).
enrich_plot <- function(df, kind, title) {
  if (!nrow(df)) return(tags$p(tags$em("No terms.")))
  if (kind == "GSEA") df <- df[order(-abs(df$NES)), ][seq_len(min(10, nrow(df))), ]
  if (kind == "GSEA") df <- df[order(-df$NES), ]
  df <- head(df, 10)
  df$term <- factor(df$ID, levels = rev(df$ID))
  p <- if (kind == "ORA") {
    df$gene_ratio <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    df$set_size <- as.numeric(sub("/.*", "", df$BgRatio))
    ggplot(df, aes(gene_ratio, term, size = set_size, colour = p.adjust)) + geom_point() +
      scale_colour_gradient(low = "#c0392b", high = "#2166ac", name = "adjusted p") +
      scale_size_continuous(name = "set size") + labs(x = "Gene ratio")
  } else {
    ggplot(df, aes(NES, term, fill = NES > 0)) + geom_col() +
      scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2166ac"), guide = "none") + labs(x = "NES")
  }
  p <- p + labs(y = NULL, title = title) + theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 7))
  f <- tempfile(fileext = ".png")
  ggsave(f, p, width = 6, height = 3.6, dpi = 150)
  tags$figure(style = "margin:0;", tags$img(src = paste0("data:image/png;base64,", base64encode(f)), style = "max-width:100%;"))
}
# A table with its plot beside it (wraps under the table on narrow screens).
side_by_side <- function(table_el, plot_el) {
  tags$div(style = "display:flex;gap:16px;align-items:flex-start;flex-wrap:wrap;margin-bottom:12px;",
           tags$div(style = "flex:1 1 380px;min-width:0;overflow-x:auto;", table_el),
           tags$div(style = "flex:1 1 380px;min-width:0;", plot_el))
}
# First number after `sep` on the first line matching `pattern` (NA if the file or line is missing).
count_from <- function(f, pattern, sep) {
  if (!file.exists(f)) return(NA_real_)
  line <- grep(pattern, readLines(f, warn = FALSE), value = TRUE)[1]
  if (is.na(line)) return(NA_real_)
  as.numeric(gsub(",", "", sub(paste0(".*", sep, "\\s*([0-9,]+).*"), "\\1", line)))
}

# ---- Section 2: sample summary (layout, strandedness) and read counts ----------------------
strand <- if (file.exists(file.path(pre, "strandedness", "strandedness_report.tsv")))
  read.delim(file.path(pre, "strandedness", "strandedness_report.tsv"), header = FALSE,
             col.names = c("sample_id", "salmon_inferred", "strandedness_used")) else data.frame()
sample_tab <- data.frame(sample_id = samples$sample_id, condition = samples$condition, replicate = samples$replicate,
                         layout = ifelse(!is.na(samples$fastq_2) & nzchar(samples$fastq_2), "paired", "single"),
                         strandedness = if (nrow(strand)) strand$strandedness_used[match(samples$sample_id, strand$sample_id)] else NA)

# Reads left after each step: raw input, after trimming, STAR unique alignments, reads assigned to genes.
read_rows <- lapply(seq_len(nrow(samples)), function(i) {
  sid <- samples$sample_id[i]
  raw      <- count_from(file.path(pre, "trim_galore_trimmed", paste0(sid, ".fastq.gz_trimming_report.txt")), "Total reads processed", ":")
  trimmed  <- count_from(file.path(pre, "trim_galore_trimmed", paste0(sid, ".fastq.gz_trimming_report.txt")), "Reads written", ":")
  star_u   <- count_from(file.path(pre, "star", paste0(sid, "_Log.final.out")), "Uniquely mapped reads number", "\\|")
  assigned <- count_from(file.path(pre, "featurecounts", paste0(sid, ".featureCounts.txt.summary")), "^Assigned\t", "\t")
  data.frame(sample = sid, raw_reads = raw, after_trimming = trimmed, uniquely_mapped = star_u,
             assigned_to_genes = assigned, pct_assigned = round(100 * assigned / raw, 1))
})
read_tab <- bind_rows(read_rows)
read_fmt <- read_tab
for (col in c("raw_reads", "after_trimming", "uniquely_mapped", "assigned_to_genes")) read_fmt[[col]] <- format(read_tab[[col]], big.mark = ",", trim = TRUE)
read_fmt$pct_assigned <- paste0(read_tab$pct_assigned, "%")
names(read_fmt) <- c("Sample", "Raw reads", "After trimming", "Uniquely mapped (STAR)", "Assigned to genes (featureCounts)", "Assigned, % of raw")

# DE gene counts per contrast (padj < 0.05 and |log2FC| >= 1)
de_counts <- bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
  tab <- read_or_empty(file.path(proc, "differential_expression", "tables", paste0(contrasts$name[i], "_all_genes.csv")))
  up <- sum(tab$status == "up_regulated", na.rm = TRUE); dn <- sum(tab$status == "down_regulated", na.rm = TRUE)
  data.frame(Contrast = paste0(contrasts$treatment[i], " vs ", contrasts$reference[i]), Up = up, Down = dn, Total = up + dn)
}))

# ---- Section 3: page sections ---------------------------------------------------------------
css <- "body{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#1c2027;line-height:1.45}
h1{color:#1a3a5c;margin-bottom:4px}.subtitle{font-size:16px;color:#1a3a5c;margin:0 0 6px 0}
h2{border-bottom:2px solid #c09434;padding-bottom:4px;margin-top:36px}h3{margin-top:24px;color:#1a3a5c}
table.t{border-collapse:collapse;font-size:13px;margin:8px 0}table.t th,table.t td{border:1px solid #d0d4da;padding:4px 8px;text-align:left}
table.t th{background:#eef1f5}figure{margin:16px 0}figcaption{font-size:13px;color:#4b5460}.note{font-size:13px;color:#4b5460}
select{font-size:15px;padding:4px 8px;margin:8px 0 16px 0}
"

multiqc_html <- file.path(pre, "multiqc_report", "multiqc_report.html")
multiqc_link <- if (file.exists(multiqc_html))
  tags$p(tags$a(href = paste0("data:text/html;base64,", base64encode(multiqc_html)), download = "multiqc_report.html",
                "Download the MultiQC report (HTML)"), " Open the downloaded file in any browser.") else
  tags$p(tags$em("Missing MultiQC report: run step 1 with --multiqc first."))

sections <- list(
  tags$h2("1. Samples and quality control"),
  tags$p("Layout and library strandedness for each sample. Strandedness is detected by Salmon, unless set in samples.csv."),
  tbl(sample_tab),
  tags$p(class = "note", "Strandedness used for the report: ", paste(unique(sample_tab$strandedness), collapse = ", "), "."),
  tags$h3("Read counts per sample"),
  tags$p("Reads left after each step, from the raw FASTQ to the reads assigned to genes that enter the differential expression analysis. Paired-end samples are counted per read pair."),
  tbl(read_fmt),
  tags$p(class = "note", "Raw reads: input FASTQ (TrimGalore). After trimming: reads passing the quality and adapter filters. Uniquely mapped: STAR unique alignments. Assigned to genes: reads counted for genes by featureCounts, the final usable number of reads."),
  tags$h3("MultiQC report (raw reads, trimming, alignment)"),
  tags$p("The full MultiQC report covers FastQC on raw reads, TrimGalore, STAR and Salmon mapping, and featureCounts."),
  multiqc_link,
  tags$h2("2. Clustering"),
  tags$p("Sample relationships are computed on variance-stabilised expression of the 500 most variable genes."),
  img(file.path(de, "pca_vst.png"), "PCA of variance-stabilised expression (top 500 most variable genes). Samples of the same condition should cluster together."),
  img(file.path(de, "correlation_heatmap.png"), "Pearson correlation between samples (top 500 most variable genes). Row and column colours show the condition, using the same colours as the PCA."),
  tags$h2("3. Differential expression"),
  tags$p("A gene is called DE when padj < 0.05 and |log2 fold change| >= 1. Positive log2FC means higher in the treatment group."),
  tbl(de_counts, digits = 7),
  tags$p(class = "note", "Choose a contrast from the menu below for its volcano, MA, heatmap and enrichment results.")
)

de_blocks <- list()
for (i in seq_len(nrow(contrasts))) {
  cn    <- contrasts$name[i]
  lab   <- paste0(contrasts$treatment[i], " vs ", contrasts$reference[i])
  tab   <- read_or_empty(file.path(proc, "differential_expression", "tables", paste0(cn, "_all_genes.csv")))
  n_up  <- sum(tab$status == "up_regulated", na.rm = TRUE); n_dn <- sum(tab$status == "down_regulated", na.rm = TRUE)
  ora_up <- read_or_empty(file.path(fa, "tables", paste0("ora_up_", cn, ".csv")))
  ora_dn <- read_or_empty(file.path(fa, "tables", paste0("ora_down_", cn, ".csv")))
  gsea   <- read_or_empty(file.path(fa, "tables", paste0("gsea_", cn, ".csv")))
  keep <- function(d) if (nrow(d)) d[, intersect(c("ID", "GeneRatio", "p.adjust", "NES"), names(d)), drop = FALSE] else d
  de_blocks[[i]] <- tags$div(class = "contrast", `data-contrast` = lab,
    tags$h3(paste0("Contrast: ", lab)),
    tags$p(paste0(n_up, " up-regulated and ", n_dn, " down-regulated genes.")),
    img(file.path(de, paste0("volcano_", cn, ".png")), "Volcano plot: log2 fold change against -log10 adjusted p-value. Red = up, blue = down, grey = not DE. DE genes are drawn on top of the grey points."),
    img(file.path(de, paste0("ma_", cn, ".png")), "MA plot: mean expression against log2 fold change. DE genes are drawn on top of the grey points."),
    img(file.path(de, paste0("heatmap_top50_", cn, ".png")), "Top 50 genes by padj, row-centred vst expression. Only the two conditions in this contrast: reference samples first, then treatment samples. Column colours show the condition."),
    tags$p(tags$strong("Enrichment (MSigDB, top 10 each)")),
    tags$p(tags$em("Over-representation, up-regulated genes")),
    side_by_side(tbl(head(keep(ora_up), 10)), enrich_plot(ora_up, "ORA", "Up-regulated genes")),
    tags$p(tags$em("Over-representation, down-regulated genes")),
    side_by_side(tbl(head(keep(ora_dn), 10)), enrich_plot(ora_dn, "ORA", "Down-regulated genes")),
    tags$p(tags$em("GSEA on the ranked gene list (NES > 0: enriched in the treatment group)")),
    side_by_side(tbl(head(keep(gsea), 10)), enrich_plot(gsea, "GSEA", "GSEA (NES)"))
  )
}
menu <- tags$p(tags$label(`for` = "contrast", "Contrast: "),
               tags$select(id = "contrast", onchange = "showContrast(this.value)",
                           lapply(seq_len(nrow(contrasts)), function(i) {
                             lab <- paste0(contrasts$treatment[i], " vs ", contrasts$reference[i])
                             tags$option(value = lab, lab)
                           })))
js <- "function showContrast(lab){document.querySelectorAll('.contrast').forEach(function(d){d.hidden = d.dataset.contrast !== lab;});}
showContrast(document.getElementById('contrast').value);"
sections <- c(sections, list(menu), de_blocks, list(tags$script(HTML(js))))

sections <- c(sections, list(
  tags$h2("Software"),
  tags$p(class = "note", "FastQC, TrimGalore, STAR, Salmon, samtools, featureCounts, MultiQC, tximport, DESeq2, clusterProfiler, msigdbr. Versions are listed in the README."),
  tags$h2("Acknowledgements"),
  tags$p(class = "note", "The development of the Bioinformatics for Next Generation Sequencing (BiNGS) shared resource facility is partially supported by the NCI P30 (P30CA196521) Cancer Center support grant, the ISMMS Skin Biology and Disease Resource-based Center NIAMS P30 support grant (AR079200), and the Black Family Stem Cell Institute. This work was also supported in part through the computational resources and staff expertise provided by Scientific Computing at the Icahn School of Medicine at Mount Sinai and supported by the Clinical and Translational Science Awards (CTSA) grant UL1TR004419 from the National Center for Advancing Translational Sciences. Research reported in this paper was supported by the Office of Research Infrastructure of the National Institutes of Health under award number S10OD026880. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health."),
  tags$p(class = "note", "We asked that our work be acknowledged in publications and presentations supported by BiNGS. Please also consider including our bioinformaticians in the authors list in cases where they contributed significantly.")
))

# ---- Section 4: build the page --------------------------------------------------------------
page <- tags$html(
  tags$head(tags$meta(charset = "utf-8"), tags$title(paste("RNA-seq report:", project)), tags$style(css)),
  tags$body(
    tags$h1(paste("RNA-seq report:", project)),
    tags$p(class = "subtitle", "BiNGS Bulk RNA-seq Course 2026 take-home pipeline v0.1, developed by Gulay B. Ulukaya from BiNGS."),
    tags$p(class = "note", "If outputs of this pipeline are used in presentations or publications, please give credit to BiNGS."),
    tags$p(class = "note", paste0("Species: ", species, ". Samples: ", nrow(samples), ". Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), ".")),
    sections
  )
)
out <- file.path(rep_dir, "rnaseq_report.html")
save_html(page, out)
cat("Report written to", out, "\n")
