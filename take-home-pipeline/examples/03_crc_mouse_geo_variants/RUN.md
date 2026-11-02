# Example run: mouse RNA-seq with a batch term, and restarting from each step

This example runs the take-home pipeline (v0.1) on the nine single-end mouse samples from the mouse run example (WT, APC and
AKSP), with the replicate number as a batch term. It then checks the restart from each step of step 1. It was run on Minerva on
2026-11-01. The starting point is the mouse run example, so download the nine GEO samples and set up the references as
described there.

## Setup

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/mouse_se_variants && cd $_
git clone --depth 1 --filter=blob:none --sparse https://github.com/bings-core/bulk-rnaseq-course-2026.git
cd bulk-rnaseq-course-2026
git sparse-checkout set take-home-pipeline/v0.1
cd take-home-pipeline/v0.1
cp config.example.sh config.sh
```

Open `config.sh` in a text editor and set these three values:

| Setting | Value |
|---|---|
| `SPECIES` | `mouse` |
| `LSF_ACCOUNT` | `acc_BiNGS_bulk` |
| `PROJECT_NAME` | `mouse_se_variants` |

Do steps 3 and 5 of the mouse run example: the GEO download of the nine samples into `data/raw/` (`submit_fetch.sh`), then
`setup_resources.sh`. Wait for the references as described there. Section 1 below replaces the sample sheet from the mouse example.

## 1. Sample and contrast sheets with a batch column

Copy the template, then replace its contents with the CSV below:

```bash
cp samples.example.csv samples.csv
```

The content of `samples.csv` is:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
WT_rep1,WT,1,data/raw/SRR25474517.fastq.gz,,auto,1
WT_rep2,WT,2,data/raw/SRR25474516.fastq.gz,,auto,2
WT_rep3,WT,3,data/raw/SRR25474515.fastq.gz,,auto,3
APC_rep1,APC,1,data/raw/SRR25474514.fastq.gz,,auto,1
APC_rep2,APC,2,data/raw/SRR25474513.fastq.gz,,auto,2
APC_rep3,APC,3,data/raw/SRR25474512.fastq.gz,,auto,3
AKSP_rep1,AKSP,1,data/raw/SRR25474511.fastq.gz,,auto,1
AKSP_rep2,AKSP,2,data/raw/SRR25474510.fastq.gz,,auto,2
AKSP_rep3,AKSP,3,data/raw/SRR25474509.fastq.gz,,auto,3
```

- `batch` is the replicate number, so each batch holds one sample of each group.
- DESeq2 uses `~ batch + condition` because the batch column is filled.

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

Check the sheets, then run step 1 and steps 2 to 4:

```bash
bash scripts/check_inputs.sh preprocess
bash scripts/submit_preprocess.sh
bash scripts/status.sh          # wait until every sample says complete
bash scripts/submit_analysis.sh
```

Results on this run (padj < 0.05 and |log2FC| >= 1), with the batch term:

| Contrast | Up | Down |
|---|---|---|
| APC_vs_WT | 2307 | 2384 |
| AKSP_vs_WT | 3201 | 3194 |
| AKSP_vs_APC | 724 | 883 |

The counts differ from the mouse run example, which has no batch term.

## 2. Restart from each step

Each step of step 1 was tested on its own. For each test, one output was deleted, `status.sh` was checked, and the sample was
restarted from the step it names:

| Sample | Output deleted | Step named by `status.sh` |
|---|---|---|
| WT_rep1 | FastQC report (`fastqc_raw/WT_rep1_fastqc.zip`) | `fastqc` |
| WT_rep2 | trimmed reads (`trim_galore_trimmed/WT_rep2_trimmed.fq.gz`) | `trim` |
| WT_rep3 | Salmon output folder (`salmon/WT_rep3`) | `salmon` |
| APC_rep1 | STAR BAM (`star/APC_rep1_Aligned.sortedByCoord.out.bam`) | `star` |
| APC_rep2 | featureCounts file (`featurecounts/APC_rep2.featureCounts.txt`) | `featurecounts` |
| APC_rep3 | Qualimap folder (`qualimap/APC_rep3`) | `qualimap` |

Each sample was restarted as its own LSF job, so that the work runs on the cluster and not on the login node:

```bash
source scripts/lib/helpers.sh && load_config
submit "from_WT_rep1" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample WT_rep1 --from fastqc
submit "from_WT_rep2" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample WT_rep2 --from trim
submit "from_WT_rep3" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample WT_rep3 --from salmon
submit "from_APC_rep1" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample APC_rep1 --from star
submit "from_APC_rep2" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample APC_rep2 --from featurecounts
submit "from_APC_rep3" "$THREADS" "$MEM_STAR_GB" -- bash scripts/01_preprocess.sh --sample APC_rep3 --from qualimap
```

To restart a set of samples from one step in a single command, set `SAMPLE_FILTER` in `config.sh` and use
`submit_preprocess.sh --from <step>` (see the README). The outputs were deleted only to test the restarts: you do not need to
delete anything to follow this example.

Result: each of the six samples reached `complete`, and each restart ran its step and everything after it. The
`status.sh` output after the restarts shows `0 sample(s) not complete.`
