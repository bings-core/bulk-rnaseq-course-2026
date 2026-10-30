# Take-home bulk RNA-seq pipeline (v0.1)

From raw FASTQ to differential expression, enrichment and one HTML report, on Minerva (LSF).
Developed by Gulay B. Ulukaya from BiNGS. **If outputs of this pipeline are used in presentations or publications, please give credit to BiNGS.**

The pipeline follows Sessions 4 to 6 of the BiNGS Bulk RNA-seq course:

| Step | Script | What it does |
|---|---|---|
| setup | `scripts/setup_resources.sh` | Downloads GENCODE, builds the STAR and Salmon indexes and the TxDb, and downloads the containers (run once) |
| 0 (optional) | `scripts/00_fetch_geo.sh`, `scripts/submit_fetch.sh` | Downloads FASTQs for GEO samples from a list of GSM IDs |
| 1 | `scripts/01_preprocess.sh`, `scripts/submit_preprocess.sh` | FastQC, TrimGalore, Salmon (library type detected), STAR, samtools, featureCounts, Qualimap, MultiQC |
| 2 | `scripts/02_deg.R` | DESeq2: PCA, sample correlation, volcano, MA, top-50 heatmaps, DE tables |
| 3 | `scripts/03_enrichment.R` | Over-representation (up and down genes) and GSEA with MSigDB gene sets |
| 4 | `scripts/04_report.R` | One self-contained HTML report with the figures, tables and MultiQC report |

Steps 2 to 4 are submitted together by `scripts/submit_analysis.sh`.

## Requirements

- Minerva with LSF. Course users: account `acc_BiNGS_bulk`. After the course, use your own `acc_*` account.
- Disk: about 60 GB for the references and containers, plus the results.
- Memory: STAR needs about 32 to 40 GB per job. Plan 8 cores per sample.
- Internet from the compute nodes (GENCODE, the Galaxy depot, the public S3 bucket, and NCBI/ENA for step 0).
- Site modules (available to all Minerva users): `fastqc/0.11.8`, `trim_galore/0.6.6`, `star/2.7.5b`, `salmon/1.2.1`, `samtools/1.11`, `subread/2.0.1`, `singularity/3.6.4`.

## 1. Get the pipeline

Clone only this folder:

```bash
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

Or download the repository and copy `take-home-pipeline/v0.1` to where you want to work, for example:

```bash
cp -r bulk-rnaseq-course-2026/take-home-pipeline/v0.1 /sc/arion/projects/<your_group>/<you>/rnaseq_project
cd /sc/arion/projects/<your_group>/<you>/rnaseq_project
```

Every path in the pipeline is relative to this folder. Run all commands from the folder that holds `config.sh`.

## 2. Configure

```bash
cp config.example.sh config.sh
```

Edit `config.sh`. The settings you are most likely to change:

| Setting | What it does |
|---|---|
| `SPECIES` | `human` or `mouse` |
| `LSF_ACCOUNT` | Your LSF account. Leave empty to run everything here, one step after another, with no jobs submitted |
| `LSF_QUEUE`, `LSF_WALLTIME` | Queue and walltime for the jobs |
| `THREADS` | Cores per tool (STAR and Salmon) |
| `MEM_*_GB` | Memory per job, by step |
| `MSIGDB_COLLECTIONS` | Gene sets for step 3, comma-separated: `H` (Hallmark), `C2:CP:KEGG`, `C5:GO:BP` |
| `RUN_FEATURECOUNTS`, `RUN_QUALIMAP` | `yes` or `no` for the optional QC steps |
| `SAMPLE_FILTER` | Run only these sample IDs (space-separated). Empty = all samples |
| `PROJECT_NAME` | Name shown at the top of the report |

Paths (`RESOURCE_DIR`, `RESULTS_DIR`, `SAMPLES_CSV`, `CONTRASTS_CSV`, the containers) are relative to this folder. You normally do not need to change them.

### Customizing for your data

| You have | What to do |
|---|---|
| Single-end reads | Leave `fastq_2` empty for that sample |
| Paired-end reads | Fill `fastq_2` with the R2 file. Single and paired samples can be mixed in one project |
| Several lanes per sample | Put the files in `fastq_1` (and `fastq_2`) as a comma-separated list, in lane order. They are concatenated |
| Known strandedness | Set `strandedness` to `unstranded`, `forward` or `reverse`. The default `auto` detects it with Salmon |
| A batch effect (for example, sequencing run) | Fill the `batch` column. DESeq2 then uses `~ batch + condition` |
| Human or mouse | `SPECIES=human` or `SPECIES=mouse`. Other species are not supported in v0.1 |
| More or other gene sets | Change `MSIGDB_COLLECTIONS` |
| Only some samples | Set `SAMPLE_FILTER` |

## 3. Set up the references and containers (once)

```bash
bash scripts/setup_resources.sh
```

With `LSF_ACCOUNT` set, this submits one job per part: the GENCODE download (each file is checked against GENCODE's MD5SUMS), the R container, the STAR index, the Salmon index, the TxDb, the Qualimap container and the MultiQC container. The builds start after the download. Check the four `setup_` jobs with `bjobs`; the references are ready when they are all `DONE`. Without `LSF_ACCOUNT` the same parts run one after another in your shell.

The setup writes `resources/<species>/` and `containers/`. These are not part of the repository. In our test the STAR index took under an hour.

## 4. Prepare your samples and contrasts

Start from the two templates:

```bash
cp samples.example.csv samples.csv
cp contrasts.example.csv contrasts.csv
```

`samples.csv`:

| Column | Meaning |
|---|---|
| `sample_id` | Unique name, no spaces |
| `condition` | Group label used in the contrasts |
| `replicate` | Replicate number within the condition |
| `fastq_1` | R1 file, or the only file for single-end data. Comma-separated for several lanes |
| `fastq_2` | R2 file. Leave empty for single-end data |
| `strandedness` | `auto` (default), `unstranded`, `forward` or `reverse` |
| `batch` | Optional. If any sample has a value, `batch` goes into the DESeq2 design |

`contrasts.csv`:

| Column | Meaning |
|---|---|
| `name` | Label for the comparison, used in output file names |
| `treatment` | Condition in the numerator |
| `reference` | Condition in the denominator. log2FC = treatment vs reference |

Your FASTQ paths in `samples.csv` can be relative to this folder (for example `data/raw/...`) or absolute.

**Public GEO data (optional).** `scripts/submit_fetch.sh GSM8020595,GSM8020596,...` downloads the FASTQs of GEO samples (one LSF job per GSM) and writes `samples_draft.csv`. The draft fills in only what GEO provides, so check the condition and replicate columns yourself before you use it as `samples.csv`.

## 5. Run

Step 1, one job per sample and then MultiQC:

```bash
bash scripts/submit_preprocess.sh
```

Check progress with `bjobs`. When every sample job is `DONE` and MultiQC is `DONE`, run steps 2 to 4:

```bash
bash scripts/submit_analysis.sh
```

The report is ready when the `de_report` job is `DONE`.

Without `LSF_ACCOUNT`, run the same steps one at a time:

```bash
bash scripts/01_preprocess.sh --sample Ras_Ctrl_rep1    # one sample
bash scripts/01_preprocess.sh --multiqc                  # after all samples
bash scripts/run_r.sh scripts/02_deg.R
bash scripts/run_r.sh scripts/03_enrichment.R
bash scripts/run_r.sh scripts/04_report.R
```

## 6. Outputs

```
results/
  preprocessed/
    fastqc_raw/  trim_galore_trimmed/  salmon/  star/  samtools/  featurecounts/  qualimap/  strandedness/
    multiqc_report/multiqc_report.html
  processed/
    gene_expression/                       raw counts and VST-normalised expression
    differential_expression/tables/        <contrast>_all_genes.csv, _upregulated.txt, _downregulated.txt
    differential_expression/figures/       PCA, correlation, volcano, MA, top-50 heatmaps (png and pdf)
    functional_analysis/tables/            ora_up_, ora_down_ and gsea_ tables per contrast
  report/rnaseq_report.html                the report: open it in any browser and share it as one file
```

The report has three sections: samples and quality control (with the MultiQC report as a download), clustering, and differential expression (one contrast at a time, from a menu).

## Known differences and limits (v0.1)

- **Salmon is not reproducible run to run.** Two identical Salmon runs differ slightly, so DEG counts can change by a few genes between reruns. STAR and featureCounts results are reproducible.
- **Qualimap** is version 2.2.2d from the Galaxy depot. The course uses 2.2.1, which is not on the depot.
- **MultiQC** is the stock 1.35 image from the Galaxy depot. The course uses a build with the same analysis code and different front-end files (jQuery, and AI features disabled), so the tables and plots should match.
- **Human and mouse only.** One species per project.
- Results depend on your data, so check the QC section of the report before interpreting the DE results.
