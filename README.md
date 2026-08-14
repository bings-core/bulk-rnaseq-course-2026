# Bulk RNA-seq Analysis Course

Welcome! Everything for the ten-week BiNGS bulk RNA-seq course lives here — setup guides,
practice data, the analysis pipeline, and reading. Ask questions any time in the
[Issue](../../issues) tab.

> 🚨 **Urgent day-of issues** (can't access the building, can't find the room, etc.): call **Dan Hasson at 646-824-6449**.

- **Meets:** Wednesdays, 4:00–5:30 PM · Sep 9 – Nov 11, 2026 (10 sessions)
- **Platform:** Minerva HPC (project `BiNGS_bulk`, `/sc/arion/projects/BiNGS_bulk/`)
- **Organized by:** BiNGS Core · Tisch Cancer Center, Mount Sinai

## Asking questions

Please open an [Issue](../../issues). 
TAs and instructors monitor them between sessions — office hours can be set up upon request.

## Schedule

| #  | Date         | Time         | Location   | Topic                                          | 
|----|--------------|--------------|------------|------------------------------------------------|
| 01 | Wed · Sep 9  | 4:00-5:30 PM | Hess 8-101 | Introduction to bulk RNA-seq & BiNGS           | 
| 02 | Wed · Sep 16 | 4:00-5:30 PM | Hess 5-101 | Introduction to HPC & R in biostatistics       | 
| 03 | Wed · Sep 23 | 4:00-5:30 PM | Hess 5-101 | QC metrics & tools                             | 
| 04 | Wed · Sep 30 | 4:00-5:30 PM | Hess 5-101 | Preprocessing                                  |
| 05 | Wed · Oct 7  | 4:00-5:30 PM | Hess 5-101 | Differential expression analysis               |
| 06 | Wed · Oct 14 | 4:00-5:30 PM | James 12   | Functional enrichment analysis                 | 
| 07 | Wed · Oct 21 | 4:00-5:30 PM | James 12   | Interpreting results · public datasets · single-cell intro  |
| 08 | Wed · Oct 28 | 4:00-5:30 PM | Hess 5-101 | Integration: ATAC-seq / ChIP-seq & epigenetics |
| 09 | Wed · Nov 4  | 4:00-5:30 PM | Hess 5-101 | Bring your own data · full pipeline            | 
| 10 | Wed · Nov 11 | 4:00-5:30 PM | Hess 8-101 | Bring your own data · full pipeline            |

## Syllabus

### Session 1 — Introduction to bulk RNA-seq & BiNGS
- **Style:** Lecture
- **Content:** Overview of bulk RNA-seq: project design, library preparation, sequencing, alignment, and analysis. Participants learn what questions bulk RNA-seq is best suited to answer and when to reach for it, plus a brief tour of adjacent transcriptomic modalities (splicing analysis, long-read RNA-seq). Covers the most important outputs of the analysis and its key takeaways. The session also introduces the course project and the BiNGS Core.

### Session 2 — Introduction to HPC & R in biostatistics
- **Style:** Lecture + Hands-on
- **Content:** Combined tour of Mount Sinai's HPC (Minerva) and R for RNA-seq. On the HPC side: navigating the Minerva file system, basic Linux commands, submitting batch jobs (LSF), and what Singularity containers are for and why the pipeline uses them. On the R side: loading packages, running commands interactively vs from scripts, and preparing the packages the DEG and functional-analysis sessions will need. Participants finish the session with a working environment for the rest of the course.
- **Tools:** Linux, LSF, Singularity containers, OnDemand Code Server, OnDemand RStudio Server, R

### Session 3 — QC metrics & tools
- **Style:** Lecture
- **Content:** What "quality" means for a bulk RNA-seq experiment and how we measure it. Walkthrough of the QC metrics that matter at each stage — raw-read quality, adapter/base bias, alignment rate, duplication, gene-body coverage, strandedness, rRNA contamination, sample correlation, PCA, variance patterns — and the tools that emit them. Participants learn how to read a MultiQC report and what should trigger dropping or reprocessing a sample before differential-expression analysis.
- **Tools:** FastQC, MultiQC, Qualimap, samtools, deepTools, PCA / correlation heatmaps

### Session 4 — Preprocessing
- **Style:** Lecture + Hands-on
- **Content:** Raw FASTQs → gene- and transcript-level counts. Participants run the full preprocessing chain on chr5-subset reads: read QC, adapter/quality trimming, alignment to the reference genome, alignment-level QC, alignment-based counting, and lightweight-mapping quantification. The session ends with a MultiQC report aggregating every tool run and a comparison of gene-level (featureCounts) vs transcript-level (Salmon) count matrices, both valid DESeq2 inputs.
- **Tools:** FastQC, TrimGalore, STAR, subread/featureCounts, samtools, deepTools bamCoverage, Qualimap, Salmon, MultiQC

### Session 5 — Differential expression analysis
- **Style:** Lecture + Hands-on
- **Content:** From count matrix to a ranked list of differentially expressed genes. Participants define the sample-condition design, import Salmon quantifications, normalize with DESeq2, fit the model, and extract DEGs with shrunken log2 fold changes. Outputs include the DEG table, a top-DEG expression heatmap, a volcano plot, an MA plot, per-gene expression boxplots, and a genome-wide karyogram of DEG log2FC — the standard figure panel for a DEG story.
- **Tools:** tximport, DESeq2, apeglm, ggplot2, ggrepel, ComplexHeatmap, circlize, patchwork

### Session 6 — Functional enrichment analysis
- **Style:** Lecture + Hands-on
- **Content:** From a DEG list to biology. Participants run both over-representation analysis (ORA) on up- and down-regulated DEGs and Gene Set Enrichment Analysis (GSEA) on the ranked full gene list against MSigDB Hallmark gene sets, then interpret the results as dotplots, cnetplots, running-enrichment plots, and leading-edge heatmaps for the top pathways (MYC targets, interferon-gamma response). The session closes by connecting the enrichment story back to the biological question introduced in Session 1.
- **Tools:** clusterProfiler, msigdbr, enrichplot, ggplot2, Enrichr

### Session 7 — Presenting & interpreting · single-cell intro
- **Style:** Lecture
- **Content:** How to tell the story once the analysis is done: which QC metrics belong in a figure, how to justify sample removal from PCA, framing the biological question and background, presenting DEGs and highlighting the most relevant functional-analysis output, and writing the methods section for the whole pipeline. Also covers where to find public bulk RNA-seq datasets (GEO, SRA, TCGA/cBioPortal, UCSC Genome Browser) to complement or benchmark your own, and closes with a brief look at how the same experimental question would be tackled with single-cell RNA-seq.
- **Tools:** GEO, SRA, TCGA / cBioPortal

### Session 8 — Integration: ATAC-seq / ChIP-seq & epigenetics
- **Style:** Lecture + Hands-on
- **Content:** How chromatin-accessibility (ATAC-seq) and TF/histone-mark (ChIP-seq, CUT&RUN) data can be layered onto RNA-seq results to explain *why* genes are differentially expressed. Participants work with matched ATAC and ChIP peak sets on the course dataset, cross-reference peaks with DEGs, visualize coverage at DEG loci, and annotate peak-gene associations. The session covers when to reach for each assay and how to design an experiment that supports integrated analysis from the start.
- **Tools:** deepTools, ChIPseeker, bedtools, IGV, UCSC Genome Browser

### Sessions 9 & 10 — Bring your own data · full pipeline
- **Style:** Hands-on
- **Content:** Two open sessions where participants run the full production BiNGS bulk RNA-seq pipeline on data they bring themselves. TAs and instructors work through it alongside each participant end-to-end: sample-metadata setup, preprocessing, DEG analysis, functional enrichment, and (where relevant) epigenetic integration. Goal: leave the course with a reproducible analysis of your own data, not just the course dataset.
- **Tools:** the full BiNGS bulk RNA-seq pipeline (all tools from Sessions 4–8)

## Getting started

Before sessions begin, complete these three steps from [`minerva/README.md`](minerva/README.md):

1. [Request a Minerva account](minerva/README.md#request-a-minerva-account) — request an account and once you get a confirmation from HPC that your minerva account is activated, send us your username so we can add you to the course project (`BiNGS_bulk`).
2. [Set up Microsoft Authenticator](minerva/README.md#set-up-microsoft-authenticator) — install and register the MFA app if you haven't before. Mount Sinai requires this for every Minerva login and VPN tunnel.
3. [One-time access check](minerva/README.md#one-time-access-check) — one-line command that creates your own folder under the shared project and confirms your access is live.

Then at the start of each Wednesday session, launch [OnDemand Code Server](minerva/README.md#ondemand-code-server) (bash sessions) or [OnDemand RStudio Server](minerva/README.md#ondemand-rstudio-server) (R sessions) with reservation `BINGS_1`. Outside of session times, leave reservation empty.

## Repository layout

| Path | Contents |
|------|----------|
| `resources/` | Tutorials, docs & papers |
| `minerva/` | Minerva account, VPN, OnDemand Code Server / RStudio Server |
| `organizers/` | Course planning & admin |
| `instructors-tas/` | Docs for the course team |

