# BiNGS Bulk RNA-seq Course — Hands-on Materials (Fall 2026)

Reference materials for the **hands-on interactive sessions** of the BiNGS bulk RNA-seq course (Mount Sinai, Fall 2026). Sessions 4-6 run live in session against this tree; every preprocessing / quantification / DE / functional step is sized to finish in **seconds-to-minutes** so participants can step through the full analysis interactively rather than wait on long LSF jobs.

Course structure:
- **[Session 4 — Preprocessing & Quantification](SESSION4.md)** (bash, ~5 min): FastQC → TrimGalore → STAR → featureCounts → samtools QC → bamCoverage → Qualimap → Salmon → MultiQC
- **[Session 5 — Quantification + Differential Expression](SESSION5.md)** (R, ~50 s): tximport → DESeq2 chr5 quant + full-genome DEG (dispersion, PCA, sample distances, volcano, MA, ComplexHeatmap, circos, boxplots, patchwork combined)
- **[Session 6 — Functional / Pathway Analysis](SESSION6.md)** (R, ~70 s): clusterProfiler + msigdbr GSEA + ORA across 7 MSigDB collections + enrichplot visualizations + leading-edge heatmaps

---

## Biological context

Small subset of a real melanoma experiment from the Bernstein lab (SKMel147 ARID2 KO vs WT):

| Sample                | Cell line | Genotype  | Purpose         |
|-----------------------|-----------|-----------|-----------------|
| `SKMel147_ARID2WT_R1` | SKMel147  | ARID2 WT  | Control rep 1   |
| `SKMel147_ARID2WT_R2` | SKMel147  | ARID2 WT  | Control rep 2   |
| `SKMel147_ARID2KO_R1` | SKMel147  | ARID2 KO  | Perturbed rep 1 |
| `SKMel147_ARID2KO_R2` | SKMel147  | ARID2 KO  | Perturbed rep 2 |

**ARID2** is a subunit of the PBAF SWI/SNF chromatin-remodeling complex; it is frequently mutated in melanoma. Knocking it out in SKMel147 alters chromatin accessibility and, downstream, gene expression. The teaching goal: identify which genes change expression between ARID2 KO and WT.

**Reduction to chromosome 5.** Raw sequencing yielded ~30-50 M reads per sample genome-wide. To make the classroom exercises fast we extracted **only the reads that map to chromosome 5** (~1.1-1.6 M reads per sample) and pre-built chr5 references (genome fasta, GTF, STAR index, Salmon index with decoys). Chr5 is ~181 Mb / ~2,600 genes — small enough to align + quantify in seconds per sample, large enough to give real biology.

> ⚠️ **Reduced-chromosome caveat.** Because chr5 discards the rest of the genome, library-size normalization and dispersion estimation are less stable than on a full dataset. Session 5 handles this by using chr5 only for the **quantification** half, then switching to a pre-computed **full-genome** DEG result from an upstream pipeline run for the interpretation half (dispersion, PCA, volcano, MA, heatmap, boxplots, functional analysis).

---

## Quickstart on Minerva

```bash
cd /sc/arion/projects/BiNGS_bulk/$USER/hands_on

# Run one session at a time, following its README
bash code_rna/session_04_preprocessing.sh          # session 4: preprocessing
module load singularity/3.6.4
singularity exec engine/singularity_containers/rbings_20250925.sif R
# then inside R:
> source("code_rna/session_05_deg_analysis.R")     # session 5: DEG
> source("code_rna/session_06_functional_analysis.R")   # session 6: functional
```

Full step-by-step walkthroughs in [SESSION4.md](SESSION4.md), [SESSION5.md](SESSION5.md), [SESSION6.md](SESSION6.md).

---

## Directory layout

```
hands_on/                                    ← repo root (this repo)
├── README.md                                ← you are here
├── SESSION4.md                              ← preprocessing walkthrough
├── SESSION5.md                              ← DEG walkthrough
├── SESSION6.md                              ← functional analysis walkthrough
│
├── code_rna/
│   ├── session_04_preprocessing.sh          ← bash preprocessing pipeline
│   ├── session_05_deg_analysis.R            ← 20-section DEG script
│   ├── session_06_functional_analysis.R     ← 13-section functional-analysis script
│   └── multiqc_config_rnaseq.yaml
│
├── engine/                                  ← REFERENCE MATERIAL (participant-read-only)
│   ├── download_reference_data.sh           ← symlink the big blobs back on Minerva
│   ├── annotation/data/
│   │   ├── gencode.v36.annotation_chr5.gtf         ← chr5 gene models
│   │   ├── gencode.v36.annotation_chr5.sqlite      ← pre-built TxDb
│   │   ├── gencode.v36.chr5.transcripts.fa         ← chr5 transcripts
│   │   ├── prepare_chr5_transcriptome_fasta.R      ← helper (how the tx fasta was made)
│   │   └── GRCh38.p13.chr5.genome.fa               ← 177 MB — .gitignored, symlinked
│   ├── annotation/star/                     ← 2.0 GB STAR index — .gitignored, symlinked
│   ├── annotation/salmon/                   ← 800 MB Salmon index — .gitignored, symlinked
│   └── singularity_containers/
│       └── rbings_20250925.sif              ← 7.1 GB — .gitignored, symlinked
│                                                 R 4.3.1, DESeq2, clusterProfiler, msigdbr,
│                                                 enrichplot, ComplexHeatmap, circlize, patchwork
│
└── data_rna/                                ← RNA-seq pipeline data (scope-first)
    ├── reduced_chr5/                        ← chr5-only teaching subset
    │   ├── raw/                             ← chr5 fastqs (.gitignored) + sample_metadata_rna.csv
    │   ├── preprocessed/                    ← session-4 outputs (fastqc, trimgalore, star,
    │   │                                       subread, samtools, deeptools, qualimap, salmon,
    │   │                                       multiqc)
    │   └── processed/                       ← session-5 Part A outputs (chr5 quant tables)
    │       └── gene_expression/{tables,figures}/
    │
    └── all_chr/                             ← full-genome (upstream pipeline inputs + session-5/6 outputs)
        └── processed/
            ├── gene_expression/{tables,figures}/
            ├── differential_expression/{tables,figures}/
            └── functional_analysis/{tables,figures}/
```

---

## Compute environment

The session runs on Mount Sinai's Minerva HPC using the reserved node/queue via OnDemand:

| Resource         | Value                        |
|------------------|------------------------------|
| Reservation ID   | `BINGS_1`                    |
| Project account  | `acc_BiNGS_bulk`             |
| Queue            | `premium`                    |

The reservation `BINGS_1` is live **continuously from Sep 2 to Nov 12**, not
only during Wednesday class hours — so use it whenever you come back to redo a
session or work through the material on your own time. See
[Using the reservation](../minerva/README.md#using-the-reservation).

All R sessions run using the `rbings_20250925.sif` singularity container (R 4.3.1 with DESeq2 1.42, clusterProfiler 4.10, msigdbr 7.5, enrichplot 1.22, org.Hs.eg.db 3.18, plus ComplexHeatmap, circlize, patchwork, and every other package the three sessions need).

---

## Data provenance

- **chr5 fastqs** — chr5-subset extracts of Bernstein-lab SKMel147 ARID2-KO / -WT RNA-seq.
- **Full-genome DEG + gene expression** used in Session 5 Part B + Session 6. Comparison: ARID2 KO (4 reps: ARID2KO1-4) vs CTRL (4 reps: ARID2NT1-4). Filtered to 22,717 tested genes; 269 down + 290 up at padj < 0.05 & |log2FC| ≥ 1.
- **Container** — `rbings_20250925.sif` is a BiNGS R build; the version referenced here is pinned for reproducibility.
