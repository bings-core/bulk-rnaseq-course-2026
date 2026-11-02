# Example run: mouse RNA-seq from GEO (WT, APC and AKSP colorectal samples)

This is a worked example of the take-home pipeline (v0.1) on a public dataset, from a fresh clone to the HTML report.
Every command below was run on Minerva on 2026-11-01. The FASTQ files come from GEO, so you can repeat the run as written.

Pipeline: `take-home-pipeline/v0.1` from the latest `main` of this repository.

## About the dataset

- **Study:** RNA-sequencing of genetically engineered intestinal organoid lines modeling colorectal cancer (CRC) progression.
  GEO series GSE239754, SRA project PRJNA1000918, submitted by the Icahn School of Medicine at Mount Sinai.
- **Biology:** mouse intestinal organoids (tumoroids) along the CRC progression series: wild type (WT), `ApcKO` (A), and
  `ApcKO::KrasG12D::Smad4KO::Trp53KO` (AKSP). Here the `ApcKO` group is called APC. Organoids were collected at confluence and
  sequenced in bulk.
- **Design:** three replicates per line, nine samples in total.
- **Sequencing:** Illumina NovaSeq 6000, cDNA libraries, single-end reads (about 46 million reads for WT Rep1, for example).
- **Question:** how transcription changes across the malignancy continuum, from normal organoids to the APC-loss and
  AKSP tumoroids.

## Samples

Nine single-end mouse samples from GEO, three per group:

| Sample ID | GEO sample | SRA run | Group |
|---|---|---|---|
| WT_rep1 | GSM7670963 | SRR25474517 | WT |
| WT_rep2 | GSM7670964 | SRR25474516 | WT |
| WT_rep3 | GSM7670965 | SRR25474515 | WT |
| APC_rep1 | GSM7670966 | SRR25474514 | APC |
| APC_rep2 | GSM7670967 | SRR25474513 | APC |
| APC_rep3 | GSM7670968 | SRR25474512 | APC |
| AKSP_rep1 | GSM7670969 | SRR25474511 | AKSP |
| AKSP_rep2 | GSM7670970 | SRR25474510 | AKSP |
| AKSP_rep3 | GSM7670971 | SRR25474509 | AKSP |

Contrasts (treatment vs reference, so log2FC > 0 means higher in the first group named):

- `APC_vs_WT`
- `AKSP_vs_WT`
- `AKSP_vs_APC`

## What you will run

| Step | Command | Time on this run |
|---|---|---|
| Clone the pipeline | `git clone --sparse ...` | 1 min |
| Download the FASTQs from GEO | `submit_fetch.sh` | about 15 min (9 samples in parallel) |
| Set up references and containers | `setup_resources.sh` | about 1 hour (STAR index about 20 min) |
| Step 1: preprocessing | `submit_preprocess.sh` | about 15 to 20 min per sample (all in parallel) |
| Steps 2 to 4: DE, enrichment, report | `submit_analysis.sh` | about 5 min |

Requirements: Minerva with LSF, your `acc_*` account (the course account is `acc_BiNGS_bulk`), about 60 GB of free space,
write access to `/sc/arion/projects/BiNGS_bulk/$USER`, and internet access from the compute nodes (NCBI, ENA and GENCODE).

## 1. Clone the pipeline from GitHub

Run these from your own project folder. Every path below is relative to `take-home-pipeline/v0.1`.

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/mouse_geo && cd $_
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
| `SPECIES` | `mouse` |
| `LSF_ACCOUNT` | `acc_BiNGS_bulk` |
| `PROJECT_NAME` | `crc_mouse_geo` |

## 3. Download the FASTQs from GEO

```bash
bash scripts/submit_fetch.sh GSM7670963,GSM7670964,GSM7670965,GSM7670966,GSM7670967,GSM7670968,GSM7670969,GSM7670970,GSM7670971
```

This submits one LSF job per GSM (`fetch_<GSM>`), then one job that writes `samples_draft.csv`. The files go to
`data/raw/`, named after their SRA run (for example `data/raw/SRR25474517.fastq.gz`). Each file is checked against the
MD5 that ENA publishes.

Wait until the downloads are finished. This takes about 15 minutes. Check with `bjobs`: every `fetch_` job and
`fetch_draft` should show `DONE`. Then check that all nine files arrived:

```bash
bjobs -w | grep fetch_
ls data/raw/*.fastq.gz | wc -l    # expected: 9
```

The draft has the file names but no condition or replicate. Write `samples.csv` yourself (step 4). The draft does not
know the groups; the GEO sample titles do (WT, APC, AKSP, replicates 1 to 3).

## 4. Write the sample and contrast sheets

Copy the template, then replace its contents with the CSV below:

```bash
cp samples.example.csv samples.csv
```

The content of `samples.csv` is:

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
APC_vs_WT,APC,WT
AKSP_vs_WT,AKSP,WT
AKSP_vs_APC,AKSP,APC
```

Log2FC is positive when a gene is higher in the treatment group.

Check the sheets. Before the references exist, this reports them as missing; that is expected:

```bash
bash scripts/check_inputs.sh preprocess
```

## 5. Set up the references and containers (once, about 1 hour)

```bash
bash scripts/setup_resources.sh
```

This submits seven LSF jobs: `setup_download`, `setup_star`, `setup_salmon`, `setup_txdb`, `setup_rcontainer`,
`setup_qualimap` and `setup_multiqc`. The references are ready when `setup_download`, `setup_star`, `setup_salmon` and
`setup_txdb` are `DONE`. Check with `bjobs`.

Wait for `setup_star` before going on. The check refuses to submit while the STAR or Salmon index is unfinished, and it
says so. Then check again:

```bash
bash scripts/check_inputs.sh preprocess
# Inputs OK for preprocess: 9 sample(s).
```

## 6. Step 1: preprocessing

```bash
bash scripts/submit_preprocess.sh
```

This checks the inputs, writes `results/run_info.txt`, and submits one job per sample and then the MultiQC job. On this
run all nine samples and MultiQC finished on the first attempt.

Check progress:

```bash
bash scripts/status.sh
```

The end of the output should look like this:

```
AKSP_rep3   SE   done done done done done done complete
...
0 sample(s) not complete.
```

If a sample fails, `status.sh` names the failed job and the `--from` step to restart it from, and tells you about the
MultiQC job if it is waiting. See the README section "If a job fails". Restart only the affected samples with
`SAMPLE_FILTER` in `config.sh` and `submit_preprocess.sh --from <step>`.

Strandedness on this run: all nine samples were detected as `SR` (reverse), so use `strandedness_report.tsv` to confirm it
before you interpret the DE results.

## 7. Steps 2 to 4: DE, enrichment and the report

```bash
bash scripts/submit_analysis.sh
```

The script checks the contrasts and the Salmon outputs, then submits `de_deg`, `de_enrich` and `de_report` in order. The
report is ready when `de_report` is `DONE`.

Results on this run (padj < 0.05 and |log2FC| >= 1):

| Contrast | Up | Down |
|---|---|---|
| APC_vs_WT | 2128 | 2411 |
| AKSP_vs_WT | 2898 | 3189 |
| AKSP_vs_APC | 660 | 793 |

Outputs to look at first:

- `results/report/rnaseq_report.html`: the samples and QC (MultiQC as a download), the PCA and correlation, and the three contrasts from the menu.
- `results/processed/differential_expression/tables/<contrast>_all_genes.csv`
- `results/processed/functional_analysis/tables/`: ORA and GSEA per contrast.
- `results/run_info.txt` and `results/tool_versions.tsv`.
- To look at the report on your own computer, copy `results/report/rnaseq_report.html` (for example with `scp` or your file browser) and open it in a browser.

## 8. If something fails

1. `bash scripts/status.sh`: incomplete samples, the `--from` step, and the failed jobs that still matter.
2. The error log: `logs/<job name>.<job id>.err`.
3. Fix the cause and restart with `--from`.
