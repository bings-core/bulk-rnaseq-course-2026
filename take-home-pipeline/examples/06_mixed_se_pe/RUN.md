# Example run: single-end and paired-end samples in one project, and a killed sample job

This example mixes two single-end and two paired-end human samples in one sample sheet. It checks that each sample is processed
with its own layout, and that a sample job killed while it runs can be restarted. It was run on Minerva on 2026-11-01.
The biology is not meaningful here (the groups are chosen only to test the layouts): use it as a check of the mechanics.

## Setup

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/mixed_se_pe && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
cp config.example.sh config.sh
```

Open `config.sh` in a text editor and set these three values:

| Setting | Value |
|---|---|
| `SPECIES` | `human` |
| `LSF_ACCOUNT` | `acc_BiNGS_bulk` |
| `PROJECT_NAME` | `mixed_se_pe` |

Put the FASTQs in `data/raw/` with these names:

- `Differentiated_rep1.fastq.gz` and `Immature_rep1.fastq.gz`: single-end, from the human run example (`KQLGC3_7_R1T_1` and `KQLGC3_10_R2T_1`), renamed.
- `SRR1039508_1.fastq.gz`, `SRR1039508_2.fastq.gz`, `SRR1039509_1.fastq.gz` and `SRR1039509_2.fastq.gz`: paired-end, from the airway run example.

Set up the references as described in the README (or reuse them, as the other examples did).

## 1. Sample and contrast sheets

Copy the template, then replace its contents with the CSV below:

```bash
cp samples.example.csv samples.csv
```

The content of `samples.csv` is:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
Differentiated_rep1,A,1,data/raw/Differentiated_rep1.fastq.gz,,auto,
N61311_Dex,A,2,data/raw/SRR1039509_1.fastq.gz,data/raw/SRR1039509_2.fastq.gz,auto,
Immature_rep1,B,1,data/raw/Immature_rep1.fastq.gz,,auto,
N61311_untreated,B,2,data/raw/SRR1039508_1.fastq.gz,data/raw/SRR1039508_2.fastq.gz,auto,
```

Copy the template, then replace its contents with the CSV below:

```bash
cp contrasts.example.csv contrasts.csv
```

The content of `contrasts.csv` is:

```csv
name,treatment,reference
A_vs_B,A,B
```

Check the sheets, then run step 1:

```bash
bash scripts/check_inputs.sh preprocess
bash scripts/submit_preprocess.sh
```

## 2. Check the layouts

```bash
bash scripts/status.sh
```

The layout column should show `SE` for `Differentiated_rep1` and `Immature_rep1`, and `PE` for `N61311_Dex` and
`N61311_untreated`. On this run it did.

## 3. Killed sample job, then restart

While the sample jobs ran, the job for `N61311_Dex` was killed with `bkill`. Then:

```bash
bash scripts/status.sh
```

`status.sh` reported `N61311_Dex` as incomplete, with the step to restart from (`salmon`, the step that was running), and it
listed the exited job as failed. Because the MultiQC job was waiting on the failed job, `status.sh` also named it, so that it
could be removed:

```bash
bkill <the MultiQC job id that status.sh names>
```

Restart only that sample from the step that `status.sh` named. Set `SAMPLE_FILTER` in `config.sh` as a quoted string, restart,
then set it back to empty:

```bash
# in config.sh:  SAMPLE_FILTER="N61311_Dex"
bash scripts/submit_preprocess.sh --from salmon
# then set SAMPLE_FILTER back to empty
bash scripts/status.sh          # wait until all four samples say complete
```

Result: `N61311_Dex` completed after the restart from `salmon`, and the MultiQC job ran after it.

## 4. Steps 2 to 4

```bash
bash scripts/submit_analysis.sh
```

Result: all three jobs finished without errors. The report's sample table lists each sample with its layout (`single` or
`paired`). The contrast `A_vs_B` has 0 DE genes: with two samples per group and a mixed design, nothing passes the cutoffs.
This is expected for a mechanics test.
