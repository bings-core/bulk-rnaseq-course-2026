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

Worked, step-by-step examples with real commands and outputs, from simplest to most complicated:

1. [Mouse from GEO, single-end, three contrasts](../examples/01_crc_mouse_geo/RUN.md)
2. [Human patient RNA-seq, single-end, two conditions](../examples/02_nup98_nsd1_human/RUN.md)
3. [Mouse with a batch term, and restarts from each step](../examples/03_crc_mouse_geo_variants/RUN.md)
4. [Human single-end options: multi-lane, strandedness, cutoffs, optional steps, sample subset](../examples/04_nup98_nsd1_human_variants/RUN.md)
5. [Human paired-end with a batch term (airway)](../examples/05_airway_pe_batch/RUN.md)
6. [Mixed single-end and paired-end, and a killed job](../examples/06_mixed_se_pe/RUN.md)

## Requirements

- Minerva with LSF. Course users: account `acc_BiNGS_bulk`. After the course, use your own `acc_*` account.
- Disk: about 45 GB per species for the references and containers. Each sample needs about 7 GB for a 90-million-read sample (trimmed FASTQs and the BAM), so plan for that plus the results.
- Memory: each sample job asks for 40 GB in total on 8 cores (`MEM_STAR_GB`). Its peak use was about 30 GB. The STAR index build peaked at about 35 GB.
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
| `MEM_*_GB` | Total memory per job, by step. LSF on Minerva counts memory per core, so the script divides this by `THREADS` |
| `PADJ_CUTOFF`, `LFC_CUTOFF` | A gene is DE when padj < `PADJ_CUTOFF` and \|log2FC\| >= `LFC_CUTOFF` (defaults 0.05 and 1). `PADJ_CUTOFF` is also the enrichment cutoff in step 3 |
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
| Several lanes per sample | Put the files in `fastq_1` (and `fastq_2`) as a semicolon-separated list (`;`, not a comma, which would split the CSV field), in lane order. They are concatenated |
| Known strandedness | Set `strandedness` to `unstranded`, `forward` or `reverse`. The default `auto` detects it with Salmon |
| A batch effect (for example, sequencing run) | Fill the `batch` column. DESeq2 then uses `~ batch + condition` |
| Human or mouse | `SPECIES=human` or `SPECIES=mouse`. Other species are not supported in v0.1 |
| More or other gene sets | Change `MSIGDB_COLLECTIONS` |
| Only some samples | Set `SAMPLE_FILTER` |
| Other DE or enrichment cutoffs | Set `PADJ_CUTOFF` and `LFC_CUTOFF` |

## 3. Set up the references and containers (once)

### Before December 2026: copy them from the shared folder

The references and containers were built with this pipeline and are in the shared test folders. Copy them into your
project folder instead of building them. Choose your species.

Human (`SPECIES=human`):

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/airway_pe_batch/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
cp -r $SHARED/resources .
mkdir -p containers
cp $SHARED/containers/*.sif containers/
ls resources/human containers
```

Mouse (`SPECIES=mouse`):

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/mouse_geo/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
cp -r $SHARED/resources .
mkdir -p containers
cp $SHARED/containers/*.sif containers/
ls resources/mouse containers
```

The copy takes a while because the references are large. Before you start step 1, `scripts/check_inputs.sh` checks
that the STAR and Salmon indexes are complete.

### From December 2026: set up the references and containers

```bash
bash scripts/setup_resources.sh
```


```bash
bash scripts/setup_resources.sh
```

With `LSF_ACCOUNT` set, this submits one job per part: the GENCODE download (each file is checked against GENCODE's MD5SUMS), the R container, the STAR index, the Salmon index, the TxDb, the Qualimap container and the MultiQC container. The builds start after the download. That is seven jobs. Check them with `bjobs`. The references are ready when `setup_download`, `setup_star`, `setup_salmon` and `setup_txdb` are all `DONE`; the three container jobs finish earlier. Without `LSF_ACCOUNT` the same parts run one after another in your shell.

The input check (`check_inputs.sh`) refuses to submit until the STAR and Salmon indexes are finished, so wait for those jobs before you submit.

The setup writes `resources/<species>/` and `containers/`. These are not part of the repository. In our test the STAR index took under an hour. Each container is checked against its SHA-256 after download; a file that does not match is downloaded again.

## 4. Prepare your samples and contrasts

Start from the two templates:

```bash
cp samples.example.csv samples.csv
cp contrasts.example.csv contrasts.csv
```

`samples.csv` has one row per sample. The template looks like this (paired-end example):

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
CTRL_1,CTRL,1,data/raw/CTRL_1_R1.fastq.gz,data/raw/CTRL_1_R2.fastq.gz,auto,
CTRL_2,CTRL,2,data/raw/CTRL_2_R1.fastq.gz,data/raw/CTRL_2_R2.fastq.gz,auto,
KO_1,KO,1,data/raw/KO_1_R1.fastq.gz,data/raw/KO_1_R2.fastq.gz,auto,
KO_2,KO,2,data/raw/KO_2_R1.fastq.gz,data/raw/KO_2_R2.fastq.gz,auto,
```

Columns:

- `sample_id`: unique name, no spaces.
- `condition`: group label used in the contrasts.
- `replicate`: replicate number within the condition.
- `fastq_1`: R1 file, or the only file for single-end data. Comma-separated for several lanes.
- `fastq_2`: R2 file. Leave empty for single-end data.
- `strandedness`: `auto` (default), `unstranded`, `forward` or `reverse`.
- `batch`: optional. If any sample has a value, `batch` goes into the DESeq2 design.

`contrasts.csv` has one row per comparison. The template looks like this:

```csv
name,treatment,reference
KO_vs_CTRL,KO,CTRL
```

Columns:

- `name`: label for the comparison, used in output file names.
- `treatment`: condition in the numerator.
- `reference`: condition in the denominator. log2FC = treatment vs reference.

Your FASTQ paths in `samples.csv` can be relative to this folder (for example `data/raw/...`) or absolute.

**Public GEO data (optional).** `scripts/submit_fetch.sh GSM8020595,GSM8020596,...` downloads the FASTQs of GEO samples (one LSF job per GSM) and writes `samples_draft.csv`. The draft fills in only what GEO provides, so check the condition and replicate columns yourself before you use it as `samples.csv`.

## 5. Run

Step 1, one job per sample and then MultiQC:

```bash
bash scripts/submit_preprocess.sh
```

Before anything is submitted, `scripts/check_inputs.sh` checks the sample sheet (IDs, FASTQ files, strandedness values, batch column), the references and the containers. It lists every problem and submits nothing if there is an error. It also warns about conditions with fewer than 2 replicates.

Check progress with `bash scripts/status.sh`. It shows which steps each sample has finished, the `--from` step to restart from, and any failed jobs. Or use `bjobs`. When every sample job is `DONE` and MultiQC is `DONE`, run steps 2 to 4:

```bash
bash scripts/submit_analysis.sh
```

The report is ready when the `de_report` job is `DONE`. `submit_analysis.sh` checks the contrasts and the Salmon outputs first.

### If a job fails

1. Run `bash scripts/status.sh`. It names the failed job, and the `--from` step to restart each sample from.
2. Read the job's error log, `logs/<name>.<jobid>.err`.
3. Fix the cause. Restart the failed samples with `--from` (see below).
4. If the MultiQC job is waiting for a failed sample, `status.sh` tells you to remove it with `bkill`. It cannot start until every sample is `DONE`.

Without `LSF_ACCOUNT`, every step runs on the login node, one after another. Use that only for one sample or a quick test. Real runs need `LSF_ACCOUNT`.

Without `LSF_ACCOUNT`, run the same steps one at a time:

```bash
bash scripts/01_preprocess.sh --sample Ras_Ctrl_rep1    # one sample
bash scripts/01_preprocess.sh --multiqc                  # after all samples
bash scripts/run_r.sh scripts/02_deg.R
bash scripts/run_r.sh scripts/03_enrichment.R
bash scripts/run_r.sh scripts/04_report.R
```

### Restart from a step

If a step fails for some samples, restart those samples from that step instead of from the beginning. `--from` runs the named step and every step after it; earlier outputs must already exist.

```bash
SAMPLE_FILTER="KO_1 KO_2"  # set in config.sh: only these samples
bash scripts/submit_preprocess.sh --from star
```

The steps are `fastqc`, `trim`, `salmon`, `star`, `featurecounts` and `qualimap`. MultiQC always runs at the end.

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
  run_info.txt                             the run record: date, host, pipeline commit, containers, and a copy of config.sh
  tool_versions.tsv                        tool versions for the run (also shown in the report's Software section)
```

The report has three sections: samples and quality control (with the MultiQC report as a download), clustering, and differential expression (one contrast at a time, from a menu).
