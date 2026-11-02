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
(references, containers, FASTQ copies and results), and read access to the source folder.

## 1. Clone the pipeline from GitHub

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/human_inhouse && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
```

Every path in the pipeline is relative to this folder, so run all commands from here.

## 2. Copy the FASTQ files into `data/raw/`

The source folder is read-only for you. Copy the six files, then check each copy against its source:

```bash
SRC=/sc/arion/projects/BiNGS/bings_omics/data/bings/2026/elvinwagenblast/poddan01/nup98_nsd1_patient_rnaseq/raw/rna/homo_sapiens
mkdir -p data/raw
for g in Differentiated_1 Differentiated_2 Differentiated_3 Immature_1 Immature_2 Immature_3; do
  cp $SRC/$g/*.fastq.gz data/raw/
done
for f in data/raw/*.fastq.gz; do
  [ "$(stat -c %s $f)" = "$(stat -c %s $(find $SRC -name $(basename $f)))" ] && gzip -t $f && echo "ok $(basename $f)"
done
```

Expected: six `ok` lines. The files are single-end. `R1T` and `R2T` in the names are not read pairs; they are just labels.

## 3. Configure

```bash
cp config.example.sh config.sh
```

Changed in `config.sh` (everything else stays at the defaults):

```bash
SPECIES=human
LSF_ACCOUNT=acc_BiNGS_bulk
PROJECT_NAME=nup98_nsd1_patient_rnaseq
```

## 4. Write the sample and contrast sheets

`samples.csv` (single-end, strandedness detected by Salmon):

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
Differentiated_rep1,Differentiated,1,data/raw/KQLGC3_7_R1T_1.fastq.gz,,auto,
Differentiated_rep2,Differentiated,2,data/raw/KQLGC3_8_R1T_2.fastq.gz,,auto,
Differentiated_rep3,Differentiated,3,data/raw/KQLGC3_9_R1T_3.fastq.gz,,auto,
Immature_rep1,Immature,1,data/raw/KQLGC3_10_R2T_1.fastq.gz,,auto,
Immature_rep2,Immature,2,data/raw/KQLGC3_11_R2T_2.fastq.gz,,auto,
Immature_rep3,Immature,3,data/raw/KQLGC3_12_R2T_3.fastq.gz,,auto,
```

`contrasts.csv` (treatment is Differentiated, reference is Immature, so log2FC > 0 means higher in Differentiated):

```csv
name,treatment,reference
Differentiated_vs_Immature,Differentiated,Immature
```

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

## Notes for next time

- The first `check_inputs.sh` run is expected to fail until `setup_resources.sh` has finished.
- With `LSF_ACCOUNT` empty, every step runs on the login node. Do not do that for a real dataset.
- Use a fresh clone for each new project, so the scripts and the run record match the commit you ran.
