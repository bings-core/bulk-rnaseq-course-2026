# Resources

Tutorials, documentation, and papers, organized by pipeline step. _Fill in links as
materials are finalized; keep this file promotable (no internal notes)._

## Papers — bulk RNA-seq technology & analysis

Four foundational review papers, in this folder as PDFs. Read them in the
order below to go from *what RNA-seq is* → *what today's best-practice
pipeline looks like* → *the analyst's toolkit*.

- **[Wang, Gerstein & Snyder (2009). *RNA-Seq: a revolutionary tool for transcriptomics*. Nat Rev Genet 10:57–63.](nrg2484[49].pdf)**
  The seminal RNA-Seq review — introduces the technology, contrasts it with
  microarrays and EST/SAGE, and lays out the read-mapping + expression-quantification
  workflow that every modern pipeline still descends from. Useful for the "why
  are we doing this instead of microarrays" framing in Session 1.

- **[Conesa *et al.* (2016). *A survey of best practices for RNA-seq data analysis*. Genome Biology 17:13.](s13059-016-0881-8[6].pdf)**
  Practical, step-by-step review of RNA-seq analysis: experimental design,
  library prep, QC, alignment, quantification, DE testing, functional
  interpretation, and integration with other omics. Pairs directly with the
  Session 3-7 hands-on flow. The reference every senior lab uses when
  arguing about defaults.

- **[Van den Berge *et al.* (2019). *RNA Sequencing Data: Hitchhiker's Guide to Expression Analysis*. Annu Rev Biomed Data Sci 2:139–173.](annurev-biodatasci-072018-021255[31].pdf)**
  Deep dive into the statistical machinery underneath count-based expression
  analysis — normalization models, dispersion estimation, DE testing,
  differential transcript usage, single-cell adjustments. This is the paper
  to read alongside the Session 5 (DESeq2) walkthrough.

- **[Stark, Grzelak & Hadfield (2019). *RNA sequencing: the teenage years*. Nat Rev Genet 20:631–656.](s41576-019-0150-2[5].pdf)**
  A 10-years-on update of the 2009 Wang review. Covers where the field has
  gone since: long-read tech, spatial transcriptomics, single-cell, direct
  RNA sequencing, and the analytical challenges each of those brings.
  Good scoping read for Session 7 (single-cell teaser) and Session 8
  (multi-omic integration).

## Minerva training

Scientific Computing runs free Minerva training, and two of this autumn's
live sessions are close to ideal preparation for this course. **We
recommend attending both.** Thursdays, 1:00–2:00 PM, in person with a Zoom
option — you register per session to get the link.

| Date | Session | Covers |
|------|---------|--------|
| **Thu Sep 17** | Introduction to Minerva | Minerva resources, account and logging in, the user software environment, other services |
| **Thu Sep 24** | Load Sharing Facility (LSF) Job Scheduler | Job submission and monitoring, parallel jobs and GPUs, job arrays and the self-scheduler |

The timing lines up: each falls the day after our own Session 2 (Sep 16,
HPC and R) and Session 3 (Sep 23). If anything in ours moves too fast,
theirs covers the same ground from scratch.

Their first session (Sep 10) is about the Minerva *Restricted* Cluster,
which this course doesn't use — you can skip it.

Registration and the full schedule are on the landing page:

- 📎 **[Minerva User Group & Training Sessions](https://labs.icahn.mssm.edu/minervalab/resources/the-minerva-user-group-and-training-classes/)**

If you'd rather not wait, or you miss one, the archive has recordings of
past runs of the same material:

- **February 2, 2026 — Introduction to Minerva.** File systems, module system,
  batch vs interactive jobs, quotas. Foundation for everything the course does
  on Minerva — start here.
- **October 3, 2025 — Access Minerva via web browser Open OnDemand.** How to
  open an RStudio / code-server / interactive shell from your browser without
  touching the terminal. This is the entry point most course participants use.
- **February 26, 2026 — Load Sharing Facility (LSF) Job Scheduler.** How
  `bsub` / `bjobs` / interactive-vs-batch job submission works on Minerva.
  Needed for anything the course dispatches to compute nodes.

## Course project background — why this dataset matters

The dataset you will preprocess, quantify, and analyze across Sessions 4-6
is **not** a synthetic teaching set. It is the actual **SKMel147 ARID2-KO vs
ARID2-WT** melanoma RNA-seq that underpins the Bernstein lab's 2022 Cell
Reports study on ARID2 loss and melanoma metastasis — the paper below. In
session, you re-run the transcriptomic half of that study end-to-end.

- 📰 **News — [*Researchers Discover a Mutated Gene That Triggers More Aggressive Melanoma*](https://reports.mountsinai.org/article/tisch2023-05-arid2)** (Tisch Cancer Institute Report, 2023).
  Plain-language write-up of the finding: mutations in **ARID2** turn out to
  make melanoma more aggressive by rewiring the tumor cell's chromatin
  landscape, opening a route to better therapy strategies.

- 📄 **Paper — Carcamo S, Nguyen CB, Grossi E, Filipescu D, Alpsoy A, Dhiman A, Sun D, Narang S, Imig J, Martin TC, Parsons R, Aifantis I, Tsirigos A, Aguirre-Ghiso JA, Dykhuizen EC, Hasson D, Bernstein E. [*Altered BAF occupancy and transcription factor dynamics in PBAF-deficient melanoma*](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)00389-8). Cell Reports 39(1):110637 (5 April 2022). doi:10.1016/j.celrep.2022.110637 · PMID 35385731 · PMCID PMC9013128.**

  ARID2 is the most recurrently mutated SWI/SNF chromatin-remodeler subunit
  in melanoma. The authors model ARID2 loss in melanoma cells (**SKMel147** —
  the same line you'll analyze) and show it destabilizes the PBAF complex,
  redistributes BAF genome-wide, and shifts chromatin accessibility at
  enhancers driving invasion-related gene expression. ARID2-deficient cells
  gain metastatic capacity in multiple animal models. The RNA-seq arm of
  that study is exactly the assay you preprocess in Session 4, run
  differential expression on in Session 5, and interpret functionally in
  Session 6.

  **Why this matters for you.** Every hands-on artifact you make in session —
  QC report, count matrix, DEG table, pathway enrichment — is one panel of
  a published, biologically consequential Mount Sinai paper. The first
  author, Saul Carcamo, sits inside the same BiNGS core that built this
  course. You are not learning RNA-seq analysis on a toy example; but
  reproducing (and exploring) real published findings.
