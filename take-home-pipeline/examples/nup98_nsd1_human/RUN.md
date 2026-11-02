# Example run: human RNA-seq, Differentiated vs Immature (NUP98-NSD1 patient samples)

This is a worked example of the take-home pipeline (v0.1) on a real dataset, from a fresh clone to the HTML report. Every
command below was run on Minerva on 2026-11-01. The data is stored under the `elvinwagenblast` folder and owned by `poddan01`. Please ask the owner before reusing it.

Pipeline: `take-home-pipeline/v0.1` from the latest `main`.

## What you will run

| Step | Command | Time on this run |
|---|---|---|
| Clone the pipeline | `git clone --sparse ...` | 1 min |
| Copy the FASTQs | `cp` from the source folder | a few min (5 GB) |
| Configure and check the sheets | `check_inputs.sh` | seconds |
| Set up references and containers | `setup_resources.sh` | about 1 hour (STAR index about 45 min) |
| Step 1: preprocessing | `submit_preprocess.sh` | about 7 min per sample, all samples in parallel |
| Steps 2 to 4: DE, enrichment, report | `submit_analysis.sh` | about 5 min |

Requirements: Minerva with LSF, your `acc_*` account (the course account is `acc_BiNGS_bulk`), about 60 GB of free space
(references, containers, FASTQ copies and results), read access to the source folder, and write access to
`/sc/arion/projects/BiNGS_bulk/$USER`.

## 1. Clone the pipeline from GitHub

Run these from your own project folder. Every path below is relative to `take-home-pipeline/v0.1`.

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/human_inhouse && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

If `mkdir` fails, you do not have write access to `BiNGS_bulk`. Ask the BiNGS team before going on.

## 2. Copy the FASTQ files into `data/raw/`

Check that you can read the source folder, then copy the six files:

```bash
ls /sc/arion/projects/BiNGS/bings_omics/data/bings/2026/elvinwagenblast/poddan01/nup98_nsd1_patient_rnaseq/raw/rna/homo_sapiens
mkdir -p data/raw
cp /sc/arion/projects/BiNGS/bings_omics/data/bings/2026/elvinwagenblast/poddan01/nup98_nsd1_patient_rnaseq/raw/rna/homo_sapiens/*/*.fastq.gz data/raw/
ls -l data/raw
```

`ls -l data/raw` should list six `.fastq.gz` files. The files are single-end. `R1T` and `R2T` in the names are not read
pairs; they are just labels. Step 4 shows which file belongs to which sample.

## 3. Configure

Copy the template to `config.sh`:

```bash
cp config.example.sh config.sh
```

Open `config.sh` in a text editor and set these three values. Leave everything else as it is.

| Setting | Value |
|---|---|
| `SPECIES` | `human` |
| `LSF_ACCOUNT` | `acc_BiNGS_bulk` |
| `PROJECT_NAME` | `nup98_nsd1_patient_rnaseq` |

## 4. Write the sample and contrast sheets

Copy the template, then replace its contents with the CSV below:

```bash
cp samples.example.csv samples.csv
```

The content of `samples.csv` is:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
Differentiated_rep1,Differentiated,1,data/raw/KQLGC3_7_R1T_1.fastq.gz,,auto,
Differentiated_rep2,Differentiated,2,data/raw/KQLGC3_8_R1T_2.fastq.gz,,auto,
Differentiated_rep3,Differentiated,3,data/raw/KQLGC3_9_R1T_3.fastq.gz,,auto,
Immature_rep1,Immature,1,data/raw/KQLGC3_10_R2T_1.fastq.gz,,auto,
Immature_rep2,Immature,2,data/raw/KQLGC3_11_R2T_2.fastq.gz,,auto,
Immature_rep3,Immature,3,data/raw/KQLGC3_12_R2T_3.fastq.gz,,auto,
```

- Empty `fastq_2` and `batch` fields (the two commas with nothing between them): the samples are single-end and there is no batch.
- `strandedness` is `auto`, so Salmon detects it.
- Column names must match the header exactly.

Copy the template, then replace its contents with the CSV below:

```bash
cp contrasts.example.csv contrasts.csv
```

The content of `contrasts.csv` is:

```csv
name,treatment,reference
Differentiated_vs_Immature,Differentiated,Immature
```

Log2FC is positive when a gene is higher in Differentiated.

Check the sheets. This also checks the references and containers, which do not exist yet, so the check reports them as
missing. That is expected before setup:

```bash
bash scripts/check_inputs.sh preprocess
```

## 5. Set up the references and containers (once, about 1 hour)

```bash
bash scripts/setup_resources.sh
```

This submits seven LSF jobs. The references are ready when `setup_download`, `setup_star`, `setup_salmon` and `setup_txdb`
are all `DONE`. Check with `bjobs`. The three containers (`setup_rcontainer`, `setup_qualimap`, `setup_multiqc`) finish first.

The check also refuses to go on while the STAR or Salmon index is still being built, and says so. Wait for `setup_star`
and `setup_salmon` to finish, then check again:

```bash
bash scripts/check_inputs.sh preprocess
# Inputs OK for preprocess: 6 sample(s).
```

## 6. Step 1: preprocessing

```bash
bash scripts/submit_preprocess.sh
```

This checks the inputs, writes `results/run_info.txt`, and submits one job per sample and then the MultiQC job.

Check progress:

```bash
bash scripts/status.sh
```

Expected at the end:

```
Immature_rep1  SE  done done done done done done complete
...
0 sample(s) not complete.
```

Strandedness on this run: all six samples were detected as `SF` (forward), with about 99.99% of fragments consistent with that
library type. Check `results/preprocessed/strandedness/strandedness_report.tsv`.

## 7. Steps 2 to 4: DE, enrichment and the report

```bash
bash scripts/submit_analysis.sh
```

The script checks the contrasts and the Salmon outputs, then submits three jobs in order: `de_deg`, `de_enrich`,
`de_report`. The report is ready when `de_report` is `DONE`.

Results on this run:

| Contrast | Up | Down |
|---|---|---|
| Differentiated vs Immature | 1784 | 1143 |

(padj < 0.05 and |log2FC| >= 1, from the DE tables.)

Outputs to look at first:

- `results/report/rnaseq_report.html`: open it in a browser. It has the samples and QC (with the MultiQC report as a download), clustering (PCA and correlation), and the contrast with a menu.
- `results/processed/differential_expression/tables/Differentiated_vs_Immature_all_genes.csv`
- `results/processed/functional_analysis/tables/`: the ORA and GSEA tables.
- `results/run_info.txt` and `results/tool_versions.tsv`: what was run, with which config and tool versions.

## 8. If something fails

1. `bash scripts/status.sh` shows which samples are incomplete, the `--from` step to restart each one from, and the failed jobs that still matter.
2. Read the error log, `logs/<job name>.<job id>.err`.
3. Fix the cause, then restart with `--from` (only the affected samples, using `SAMPLE_FILTER` in `config.sh` if needed).
