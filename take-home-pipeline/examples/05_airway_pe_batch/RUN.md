# Example run: paired-end human RNA-seq with a batch term (airway, untreated vs Dex)

Dexamethasone (Dex) treatment of airway smooth muscle cells from four donors, untreated vs Dex, paired-end reads from
ENA (project PRJNA229998). Each donor is a cell line, which is used as the `batch` column. This run was done on Minerva on
2026-11-01 with the take-home pipeline (v0.1) from the latest `main`.

## What you will run

| Step | Command | Time on this run |
|---|---|---|
| Clone the pipeline | `git clone --sparse ...` | 1 min |
| Download the FASTQs from ENA | `submit_fetch.sh` | about 35 min for 8 samples in parallel (19 GB) |
| Set up references and containers | `setup_resources.sh` | about 1 hour (STAR index about 45 min) |
| Step 1: preprocessing | `submit_preprocess.sh` | about 30 to 40 min per sample, all in parallel |
| Steps 2 to 4: DE, enrichment, report | `submit_analysis.sh` | about 10 min |

Requirements: Minerva with LSF, your `acc_*` account (the course account is `acc_BiNGS_bulk`), about 80 GB of free space,
and write access to `/sc/arion/projects/BiNGS_bulk/$USER`.

## About the dataset

- **Study:** Human airway smooth muscle transcriptome changes in response to asthma medications. GEO series GSE52778, SRA
  project PRJNA229998, from the Department of Biostatistics, Epidemiology and Informatics, University of Pennsylvania.
- **Biology:** primary human airway smooth muscle (HASM) cells from four white male donors, each under four conditions:
  untreated, albuterol (a beta-2 agonist, 1 uM for 18 h), dexamethasone (Dex, a glucocorticoid, 1 uM for 18 h), and
  albuterol plus Dex. This example compares untreated with Dex only. The donors are the four cell lines, used as the batch.
- **Sequencing:** Illumina HiSeq 2000 with the TruSeq assay, 75 bp paired-end reads.
- **Reference:** the original analysis used hg19 and the Tuxedo tools. This pipeline uses GRCh38 (GENCODE 36), so the
  counts and DE results will not match the published analysis exactly.

## Samples

Eight paired-end samples: four donors, each untreated and Dex. The GEO samples and runs are:

| Donor (batch) | Untreated: GEO, run | Dex: GEO, run |
|---|---|---|
| N61311 | GSM1275862, SRR1039508 | GSM1275863, SRR1039509 |
| N052611 | GSM1275866, SRR1039512 | GSM1275867, SRR1039513 |
| N080611 | GSM1275870, SRR1039516 | GSM1275871, SRR1039517 |
| N061011 | GSM1275874, SRR1039520 | GSM1275875, SRR1039521 |

Some runs (SRR1039513, SRR1039516, SRR1039520 and SRR1039521) also have an unpaired file on ENA. The pipeline uses only
the `_1` and `_2` files of each run.

## 1. Clone the pipeline from GitHub

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/airway_pe_batch && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

If `mkdir` fails, you do not have write access to `BiNGS_bulk`. Ask the BiNGS team before going on.

## 2. Configure

Copy the template to `config.sh`:

```bash
cp config.example.sh config.sh
```

Open `config.sh` in a text editor and set these three values. Leave everything else as it is.

| Setting | Value |
|---|---|
| `SPECIES` | `human` |
| `LSF_ACCOUNT` | `acc_BiNGS_bulk` |
| `PROJECT_NAME` | `airway_pe_batch` |

## 3. Get the FASTQs

### Before December 2026: copy them from the shared folder

Copy the sixteen paired-end FASTQ files (`_1` and `_2` for each of the eight runs) and the run records from the shared test
folder into `data/raw/`:

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/airway_pe_batch/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
mkdir -p data/raw
cp $SHARED/data/raw/*_[12].fastq.gz data/raw/
cp -r $SHARED/data/raw/runs data/raw/
ls data/raw/*_[12].fastq.gz | wc -l    # expected: 16
```

The four unpaired files that ENA lists for some runs (for example `SRR1039513.fastq.gz`) are not copied, because no sample
uses them. The copy is about 24 GB and takes a few minutes. Then go on to step 4.

### From December 2026: download from ENA

```bash
bash scripts/submit_fetch.sh GSM1275862,GSM1275863,GSM1275866,GSM1275867,GSM1275870,GSM1275871,GSM1275874,GSM1275875
```

This submits one LSF job per GSM (`fetch_<GSM>`) and then `fetch_draft`. The FASTQ files are written to `data/raw/`. The download
needs internet access from the compute nodes (NCBI and ENA). Wait until all nine jobs are `DONE` (`bjobs -w | grep fetch_`).

The draft in `samples_draft.csv` lists the `_1` and `_2` file of each run as `fastq_1` and `fastq_2`. The sample sheet
below uses these names.

## 4. Write the sample and contrast sheets

Copy the template, then replace its contents with the CSV below:

```bash
cp samples.example.csv samples.csv
```

The content of `samples.csv` is:

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

- Both `fastq_1` and `fastq_2` are filled: the samples are paired-end.
- `batch` is the donor (cell line), so DESeq2 uses `~ batch + condition`.
- `replicate` is the donor number within each condition.

Copy the template, then replace its contents with the CSV below:

```bash
cp contrasts.example.csv contrasts.csv
```

The content of `contrasts.csv` is:

```csv
name,treatment,reference
Dex_vs_untreated,Dex,untreated
```

Check the sheets. Before the references exist, this reports them as missing; that is expected:

```bash
bash scripts/check_inputs.sh preprocess
```

## 5. Set up the references and containers

### Before December 2026: copy them from the shared folder

Copy the human references and the three containers from the shared test folder. The references were built with this
pipeline, so you do not need to build them:

```bash
SHARED=/sc/arion/projects/BiNGS_bulk/shared/pipeline_tests/airway_pe_batch/bulk-rnaseq-course-2026/take-home-pipeline/v0.1
cp -r $SHARED/resources .
mkdir -p containers
cp $SHARED/containers/*.sif containers/
ls resources/human containers
bash scripts/check_inputs.sh preprocess
# Inputs OK for preprocess: 8 sample(s).
```

The copy is about 43 GB and takes a while. `check_inputs.sh` confirms that the STAR and Salmon indexes are complete.

### From December 2026: set up the references and containers

```bash
bash scripts/setup_resources.sh
```

Wait for `setup_star` and `setup_salmon` to finish (see the README), then check again:

```bash
bash scripts/check_inputs.sh preprocess
# Inputs OK for preprocess: 8 sample(s).
```

## 6. Step 1: preprocessing

```bash
bash scripts/submit_preprocess.sh
bash scripts/status.sh
```

The end of the status output should look like this, with every sample `complete` and `PE` in the layout column:

```
N061011_Dex        PE  done done done done done done complete
...
0 sample(s) not complete.
```

On this run the strandedness detection gave `unstranded` for all eight samples. That matches the library type of this
dataset. Check `results/preprocessed/strandedness/strandedness_report.tsv`.

## 7. Steps 2 to 4: DE, enrichment and the report

```bash
bash scripts/submit_analysis.sh
```

The design is `~ batch + condition`, because the sheet has a batch column.

Results on this run (padj < 0.05 and |log2FC| >= 1):

| Contrast | Up | Down |
|---|---|---|
| Dex_vs_untreated | 477 | 409 |

Known Dex response genes are up in Dex vs untreated (values from the `log2FoldChange` column of `Dex_vs_untreated_all_genes.csv`;
padj < 1e-18 for all four):

| Gene | log2FC |
|---|---|
| FKBP5 | 4.04 |
| DUSP1 | 2.96 |
| PER1 | 3.13 |
| TSC22D3 | 3.19 |

The report's sample table lists the layout as `paired` for all eight samples. Outputs to look at first:

- `results/report/rnaseq_report.html`
- `results/processed/differential_expression/tables/Dex_vs_untreated_all_genes.csv`
- `results/processed/functional_analysis/tables/`: the ORA and GSEA tables.
- `results/run_info.txt` and `results/tool_versions.tsv`.
