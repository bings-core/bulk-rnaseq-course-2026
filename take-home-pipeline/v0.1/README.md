# Take-home bulk RNA-seq pipeline (v0.1)

From raw FASTQ to differential expression, enrichment and one HTML report, on Minerva (LSF).
Developed by Gulay B. Ulukaya from BiNGS. **If outputs of this pipeline are used in presentations or publications, please give credit to BiNGS.**

The pipeline follows Sessions 4 to 6 of the BiNGS Bulk RNA-seq course:

| Step | Script | What it does |
|---|---|---|
| setup | `scripts/setup_resources.sh` | Downloads GENCODE, builds the STAR and Salmon indexes and the TxDb, and downloads the containers (run once per species) |
| 0 (optional) | `scripts/00_fetch_geo.sh`, `scripts/submit_fetch.sh` | Downloads FASTQs for GEO samples from a list of GSM IDs |
| 1 | `scripts/01_preprocess.sh`, `scripts/submit_preprocess.sh` | FastQC, TrimGalore, Salmon (library type detected), STAR, samtools, featureCounts, Qualimap, MultiQC |
| 2 | `scripts/02_deg.R` | DESeq2: PCA, sample correlation, volcano, MA, top-50 heatmaps, DE tables |
| 3 | `scripts/03_enrichment.R` | Over-representation (up and down genes) and GSEA with MSigDB gene sets |
| 4 | `scripts/04_report.R` | One self-contained HTML report with the figures, tables and MultiQC report |

Steps 2 to 4 are submitted together by `scripts/submit_analysis.sh`.

The pipeline is one folder. To analyse a new dataset you do not edit the scripts. You copy the data and the references in, fill in two sample sheets and set a few lines of `config.sh`. Section 5 shows this for three datasets.

Worked, step-by-step examples with real commands and outputs, from simplest to most complicated:

1. [Mouse from GEO, single-end, three contrasts](../examples/01_crc_mouse_geo/RUN.md)
2. [Human patient RNA-seq, single-end, two conditions](../examples/02_nup98_nsd1_human/RUN.md)
3. [Mouse with a batch term, and restarts from each step](../examples/03_crc_mouse_geo_variants/RUN.md)
4. [Human single-end options: multi-lane, strandedness, cutoffs, optional steps, sample subset](../examples/04_nup98_nsd1_human_variants/RUN.md)
5. [Human paired-end with a batch term (airway)](../examples/05_airway_pe_batch/RUN.md)
6. [Mixed single-end and paired-end, and a killed job](../examples/06_mixed_se_pe/RUN.md)

## Requirements

- Minerva with LSF. Course users: account `acc_BiNGS_bulk`. After the course, use your own `acc_*` account.
- Disk: about 45 GB per species for the references and containers. Each sample needs about 7 GB for a 90-million-read sample (trimmed FASTQs and the BAM), so plan for that plus the results. The FASTQs you copy in count too.
- Memory: each sample job asks for 40 GB in total on 8 cores (`MEM_STAR_GB`). Its peak use was about 30 GB. The STAR index build peaked at about 35 GB.
- Internet from the compute nodes (GENCODE, the Galaxy depot, the public S3 bucket, and NCBI/ENA for step 0).
- Site modules (available to all Minerva users): `fastqc/0.11.8`, `trim_galore/0.6.6`, `star/2.7.5b`, `salmon/1.2.1`, `samtools/1.11`, `subread/2.0.1`, `singularity/3.6.4`.

## 1. Get the pipeline

Clone only this folder into your project area:

```bash
mkdir -p /sc/arion/projects/<your_group>/<you>/rnaseq_project && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

Every path in the pipeline is relative to this folder. Run all commands from the folder that holds `config.sh`.

For a second dataset, copy this folder to a new project folder and work there, so each project keeps its own `config.sh`, sheets, data, references and results:

```bash
cp -r /sc/arion/projects/<your_group>/<you>/rnaseq_project/bulk-rnaseq-course-2026/take-home-pipeline/v0.1 /sc/arion/projects/<your_group>/<you>/rnaseq_project_2
cd /sc/arion/projects/<your_group>/<you>/rnaseq_project_2
```

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
| `SAMPLE_FILTER` | Run only these sample IDs (space-separated, in quotes). Empty = all samples |
| `PROJECT_NAME` | Name shown at the top of the report |

Paths (`RESOURCE_DIR`, `RESULTS_DIR`, `SAMPLES_CSV`, `CONTRASTS_CSV`, the containers) are relative to this folder. You normally do not need to change them.

Check the file after editing. A value with a space must be quoted, for example `SAMPLE_FILTER="KO_1 KO_2"`. An unquoted second ID is run as a command.

### Customizing for your data

| You have | What to do |
|---|---|
| Single-end reads | Leave `fastq_2` empty for that sample |
| Paired-end reads | Fill `fastq_2` with the R2 file. Single and paired samples can be mixed in one project |
| Several lanes per sample | Put the files in `fastq_1` (and `fastq_2`) as a semicolon-separated list (`;`, not a comma, which would split the CSV field), in lane order. They are concatenated |
| Known strandedness | Set `strandedness` to `unstranded`, `forward` or `reverse`. The default `auto` detects it with Salmon |
| A batch effect (for example, sequencing run or cell line) | Fill the `batch` column. DESeq2 then uses `~ batch + condition` |
| Human or mouse | `SPECIES=human` or `SPECIES=mouse`. Other species are not supported in v0.1 |
| More or other gene sets | Change `MSIGDB_COLLECTIONS` |
| Only some samples | Set `SAMPLE_FILTER` |
| Other DE or enrichment cutoffs | Set `PADJ_CUTOFF` and `LFC_CUTOFF` |

## 3. Get the data

The FASTQs go into `data/raw/` of your project folder. Each sample's `fastq_1` and `fastq_2` in `samples.csv` point there (`data/raw/<file>`), or to any absolute path.

```bash
mkdir -p data/raw
```

Choose one of the four sources.

### a. Copy the shared test datasets

Two datasets are already on the shared project folder. Copy them rather than downloading them again.

Mouse colorectal organoids, GEO GSE239754, nine single-end samples:

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/mouse_geo/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
cp $SHARED/data/raw/*.fastq.gz data/raw/
cp -r $SHARED/data/raw/runs data/raw/
ls data/raw/*.fastq.gz | wc -l    # expected: 9
```

Human airway smooth muscle cells, GEO GSE52778, eight paired-end samples (16 files):

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/airway_pe_batch/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
cp $SHARED/data/raw/*_[12].fastq.gz data/raw/
cp -r $SHARED/data/raw/runs data/raw/
ls data/raw/*_[12].fastq.gz | wc -l    # expected: 16
```

The `*_[12].fastq.gz` pattern copies only the R1 and R2 files. The shared folder also has an unpaired `SRR<id>.fastq.gz` next to some pairs, which this pattern skips.

### b. Copy your own FASTQs

Copy them in, keeping their names:

```bash
cp /sc/arion/projects/<your_group>/<you>/fastq/*.fastq.gz data/raw/
```

Check the files before you start. Each one must be gzip-compressed (`gzip -t data/raw/<file>.fastq.gz`). Paired-end files need matching R1 and R2 names in the sample sheet.

### c. Download public GEO FASTQs

`scripts/submit_fetch.sh` downloads the FASTQs of GEO samples (one LSF job per GSM) into `data/raw/` and writes `samples_draft.csv`. Each file is checked against the MD5 that ENA publishes:

```bash
bash scripts/submit_fetch.sh GSM7670963,GSM7670964,GSM7670965
```

Wait for the `fetch_` jobs to finish (`bjobs`). The draft has the file names but no condition or replicate. Write `samples.csv` yourself (section 5), using the GEO sample titles for the groups.

### d. Lanes and mixed layouts

If one sample was sequenced on several lanes, copy every lane's file in and list them in lane order. Section 5 shows the sheet.

## 4. Set up the references and containers (once per species)

The references (GENCODE annotation, the STAR and Salmon indexes and the TxDb) and the containers are the same for every project of a species. Build them once, or copy them from the shared folder.

### Before December 2026: copy them from the shared folder

The references and containers were built with this pipeline and are in the shared test folders. Copy them into your project folder instead of building them. Choose your species.

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

The copy takes a while because the references are large. Instead of copying, you can link the folder into several projects with `ln -s $SHARED/resources resources`. Before you start step 1, `scripts/check_inputs.sh` checks that the STAR and Salmon indexes are complete.

### From December 2026: set up the references and containers

```bash
bash scripts/setup_resources.sh
```

With `LSF_ACCOUNT` set, this submits one job per part: the GENCODE download (each file is checked against GENCODE's MD5SUMS), the R container, the STAR index, the Salmon index, the TxDb, the Qualimap container and the MultiQC container. The builds start after the download. That is seven jobs. Check them with `bjobs`. The references are ready when `setup_download`, `setup_star`, `setup_salmon` and `setup_txdb` are all `DONE`; the three container jobs finish earlier. Without `LSF_ACCOUNT` the same parts run one after another in your shell.

The input check (`check_inputs.sh`) refuses to submit until the STAR and Salmon indexes are finished, so wait for those jobs before you submit.

The setup writes `resources/<species>/` and `containers/`. These are not part of the repository. In our test the STAR index took under an hour. Each container is checked against its SHA-256 after download; a file that does not match is downloaded again.

## 5. Describe your samples and contrasts

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
- `fastq_1`: R1 file, or the only file for single-end data. Semicolon-separated for several lanes.
- `fastq_2`: R2 file. Leave empty for single-end data.
- `strandedness`: `auto` (default), `unstranded`, `forward` or `reverse`.
- `batch`: optional. If any sample has a value, `batch` goes into the DESeq2 design, and every sample needs one.

`contrasts.csv` has one row per comparison. The template looks like this:

```csv
name,treatment,reference
KO_vs_CTRL,KO,CTRL
```

Columns:

- `name`: label for the comparison, used in output file names.
- `treatment`: condition in the numerator.
- `reference`: condition in the denominator. log2FC = treatment vs reference, so a positive log2FC is higher in the treatment.

### Three datasets, set up the same way

Each dataset below needs only the `config.sh` lines shown, the two sheets, and the data from section 3. The references from section 4 are the same for every dataset of one species.

**Dataset A: mouse, single-end, nine samples (GEO GSE239754).** Three groups, three contrasts. Copy the data from section 3a (mouse) and the references from section 4 (mouse).

`config.sh`:

```bash
SPECIES=mouse
LSF_ACCOUNT=acc_BiNGS_bulk
PROJECT_NAME=crc_mouse_geo
```

`samples.csv`:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
WT_rep1,WT,1,data/raw/SRR25474517.fastq.gz,,auto,
WT_rep2,WT,2,data/raw/SRR25474516.fastq.gz,,auto,
WT_rep3,WT,3,data/raw/SRR25474515.fastq.gz,,auto,
APC_rep1,APC,1,data/raw/SRR25474514.fastq.gz,,auto,
APC_rep2,APC,2,data/raw/SRR25474513.fastq.gz,,auto,
APC_rep3,APC,3,data/raw/SRR25474512.fastq.gz,,auto,
AKSP_rep1,AKSP,1,data/raw/SRR25474511.fastq.gz,,auto,
AKSP_rep2,AKSP,2,data/raw/SRR25474510.fastq.gz,,auto,
AKSP_rep3,AKSP,3,data/raw/SRR25474509.fastq.gz,,auto,
```

`contrasts.csv`:

```csv
name,treatment,reference
APC_vs_WT,APC,WT
AKSP_vs_WT,AKSP,WT
AKSP_vs_APC,AKSP,APC
```

**Dataset B: human, paired-end, eight samples with a batch (GEO GSE52778).** Untreated and dexamethasone-treated cells from four donors. The donor goes into `batch`, so DESeq2 corrects for it. Copy the data from section 3a (airway) and the references from section 4 (human).

`config.sh`:

```bash
SPECIES=human
LSF_ACCOUNT=acc_BiNGS_bulk
PROJECT_NAME=airway_dex
```

`samples.csv`:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
N61311_untreated,untreated,1,data/raw/SRR1039508_1.fastq.gz,data/raw/SRR1039508_2.fastq.gz,auto,N61311
N61311_Dex,Dex,1,data/raw/SRR1039509_1.fastq.gz,data/raw/SRR1039509_2.fastq.gz,auto,N61311
N052611_untreated,untreated,2,data/raw/SRR1039512_1.fastq.gz,data/raw/SRR1039512_2.fastq.gz,auto,N052611
N052611_Dex,Dex,2,data/raw/SRR1039513_1.fastq.gz,data/raw/SRR1039513_2.fastq.gz,auto,N052611
N080611_untreated,untreated,3,data/raw/SRR1039516_1.fastq.gz,data/raw/SRR1039516_2.fastq.gz,auto,N080611
N080611_Dex,Dex,3,data/raw/SRR1039517_1.fastq.gz,data/raw/SRR1039517_2.fastq.gz,auto,N080611
N061011_untreated,untreated,4,data/raw/SRR1039520_1.fastq.gz,data/raw/SRR1039520_2.fastq.gz,auto,N061011
N061011_Dex,Dex,4,data/raw/SRR1039521_1.fastq.gz,data/raw/SRR1039521_2.fastq.gz,auto,N061011
```

`contrasts.csv`:

```csv
name,treatment,reference
Dex_vs_untreated,Dex,untreated
```

**Dataset C: your own data.** A starting point for a project with paired-end samples, one sample on two lanes, and one sample with known strandedness. Replace the names with yours and keep the same columns.

`config.sh`:

```bash
SPECIES=human                 # or mouse
LSF_ACCOUNT=acc_<your_account>
PROJECT_NAME=my_project
PADJ_CUTOFF=0.05
LFC_CUTOFF=1
MSIGDB_COLLECTIONS=H,C2:CP:KEGG
```

`samples.csv`:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
Ctrl_1,Control,1,data/raw/Ctrl_1_R1.fastq.gz,data/raw/Ctrl_1_R2.fastq.gz,auto,run1
Ctrl_2,Control,2,data/raw/Ctrl_2_R1.fastq.gz,data/raw/Ctrl_2_R2.fastq.gz,auto,run2
Tx_1,Treated,1,data/raw/Tx_1_L001_R1.fastq.gz;data/raw/Tx_1_L002_R1.fastq.gz,data/raw/Tx_1_L001_R2.fastq.gz;data/raw/Tx_1_L002_R2.fastq.gz,auto,run1
Tx_2,Treated,2,data/raw/Tx_2_R1.fastq.gz,data/raw/Tx_2_R2.fastq.gz,reverse,run2
```

`contrasts.csv`:

```csv
name,treatment,reference
Treated_vs_Control,Treated,Control
```

Things to check in this example: `Tx_1` has two lanes, so both R1 files and both R2 files are listed in lane order, separated by `;`. `Tx_2` has its strandedness set to `reverse`, so the pipeline does not detect it. Every sample has a `batch` value, so DESeq2 uses `~ batch + condition`.

To look at a subset first, set `SAMPLE_FILTER="Ctrl_1 Tx_1"` in `config.sh`. Remove the filter when you are ready for all samples.

### Reusing a project for a new dataset

1. Copy the pipeline folder to a new project folder (section 1).
2. Set up its `resources/` and `containers/` (section 4), or link them from an existing project.
3. Copy the FASTQs into `data/raw/` (section 3).
4. Write `samples.csv` and `contrasts.csv` (section 5).
5. Set `SPECIES`, `LSF_ACCOUNT` and `PROJECT_NAME` in `config.sh`, and any cutoffs you want.

Results go into `results/` of that project folder only, so projects do not overwrite each other.

## 6. Run

Step 0, before anything is submitted, check the sample sheet (IDs, FASTQ files, strandedness values, batch column), the references and the containers. It lists every problem and submits nothing if there is an error. It also warns about conditions with fewer than 2 replicates.

```bash
bash scripts/check_inputs.sh preprocess
```

Fix anything it reports before you run submit_preprocess.

Step 1, one job per sample and then MultiQC:

```bash
bash scripts/submit_preprocess.sh
```

Check progress with `bash scripts/status.sh`. You can run it over and over again until you see everything is done. It shows which steps each sample has finished, the `--from` step to restart from, and any failed jobs. Or use `bjobs`. 

```bash
bash scripts/status.sh
```

When every sample job is `DONE` and MultiQC is `DONE`, run steps 2 to 4:

```bash
bash scripts/submit_analysis.sh
```

The report is ready when the `de_report` job is `DONE`. `submit_analysis.sh` checks the contrasts and the Salmon outputs first.

### If a job fails

1. Run `bash scripts/status.sh`. It names the failed job, and the `--from` step to restart each sample from.
2. Read the job's error log, `logs/<name>.<jobid>.err`.
3. Fix the cause. Restart the failed samples with `--from` (see below).
4. If the MultiQC job is waiting for a failed sample, `status.sh` tells you to remove it with `bkill`. It cannot start until every sample is `DONE`.

### Without LSF

Leave `LSF_ACCOUNT` empty to run the same steps in your shell, one after another. Use that only for one sample or a quick test on the login node. Real runs need `LSF_ACCOUNT`.

```bash
bash scripts/01_preprocess.sh --sample WT_rep1    # one sample
bash scripts/01_preprocess.sh --multiqc           # after all samples
bash scripts/run_r.sh scripts/02_deg.R
bash scripts/run_r.sh scripts/03_enrichment.R
bash scripts/run_r.sh scripts/04_report.R
```

### Restart from a step

If a step fails for some samples, restart those samples from that step instead of from the beginning. `--from` runs the named step and every step after it; earlier outputs must already exist.

```bash
# in config.sh:  SAMPLE_FILTER="KO_1 KO_2"
bash scripts/submit_preprocess.sh --from star
```

The steps are `fastqc`, `trim`, `salmon`, `star`, `featurecounts` and `qualimap`. MultiQC always runs at the end. Remove `SAMPLE_FILTER` from `config.sh` before you submit the whole project again.

## 7. Outputs

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
