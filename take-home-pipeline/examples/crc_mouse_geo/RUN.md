# Example run: mouse RNA-seq from GEO (WT, APC and AKSP colorectal samples)

This is a worked example of the take-home pipeline (v0.1) on a public dataset, from a fresh clone to the HTML report.
Every command below was run on Minerva on 2026-11-01. The FASTQ files come from GEO, so you can repeat the run as written.

Pipeline: `take-home-pipeline/v0.1` from the latest `main` of this repository.

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

Requirements: Minerva with LSF, your `acc_*` account (the course account is `acc_BiNGS_bulk`), and about 60 GB of free space.

## 1. Clone the pipeline from GitHub

```bash
mkdir -p /sc/arion/scratch/$USER/pipeline_tests/mouse_geo && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

Run all commands from this folder.

## 2. Configure

```bash
cp config.example.sh config.sh
```

Changed in `config.sh` (everything else stays at the defaults):

```bash
SPECIES=mouse
LSF_ACCOUNT=acc_BiNGS_bulk
PROJECT_NAME=crc_mouse_geo
```

## 3. Download the FASTQs from GEO

```bash
bash scripts/submit_fetch.sh GSM7670963,GSM7670964,GSM7670965,GSM7670966,GSM7670967,GSM7670968,GSM7670969,GSM7670970,GSM7670971
```

This submits one LSF job per GSM (`fetch_<GSM>`), then one job that writes `samples_draft.csv`. The files go to
`data/raw/`, named after their SRA run (for example `data/raw/SRR25474517.fastq.gz`). Each file is checked against the
MD5 that ENA publishes. Check with `bjobs`: all `fetch_` jobs and `fetch_draft` should be `DONE`.

The draft has the file names but no condition or replicate. Write `samples.csv` yourself (step 4). The draft does not
know the groups; the GEO sample titles do (WT, APC, AKSP, replicates 1 to 3).

## 4. Write the sample and contrast sheets

`samples.csv` (single-end, strandedness detected by Salmon):

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

## 8. If something fails

1. `bash scripts/status.sh`: incomplete samples, the `--from` step, and the failed jobs that still matter.
2. The error log: `logs/<job name>.<job id>.err`.
3. Fix the cause and restart with `--from`.

## Notes

- Mouse samples here are single-end and reverse-stranded. Check the strandedness report before you interpret the results.
- The contrasts were chosen by the person running the example, with WT as the reference for the first two. Pick your own.
