# Example run: human single-end options (multi-lane, strandedness, cutoffs, optional steps, sample subset)

This example runs the take-home pipeline (v0.1) on the six single-end human samples used in the human run example
(Differentiated vs Immature), and changes one option at a time: multi-lane input, strandedness overrides, cutoffs and gene
sets, optional QC steps off, and a subset of samples. It was run on Minerva on 2026-11-01. The starting point is the human
run example, so set up the references and the sample copies as described there.

## Setup

```bash
mkdir -p /sc/arion/projects/BiNGS_bulk/$USER/pipeline_tests/human_se_variants && cd $_
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
| `PROJECT_NAME` | `human_se_variants` |

Copy the six single-end FASTQs into `data/raw/`, as in the human run example (`cp` from the same source folder, then `ls -l data/raw`).
Then run `bash scripts/setup_resources.sh` and wait for the references (see the README).

## Variant 1: multi-lane input

One sample is split into two lane files, so that one sample is given as two files. The two parts are created at a read
boundary, so that concatenating them gives the original file. Run these inside `data/raw/`:

```bash
cd data/raw
lines=$(zcat R1T_1.fastq.gz | wc -l); half=$(( (lines/8)*4 ))
zcat R1T_1.fastq.gz | head -n $half | gzip -1 > R1T_1_lane1.fastq.gz
zcat R1T_1.fastq.gz | tail -n +$((half+1)) | gzip -1 > R1T_1_lane2.fastq.gz
cd ../..
```

In `samples.csv`, the two lanes are separated by a semicolon in `fastq_1`. A comma would split the CSV field:

```csv
sample_id,condition,replicate,fastq_1,fastq_2,strandedness,batch
Differentiated_rep1,Differentiated,1,data/raw/R1T_1_lane1.fastq.gz;data/raw/R1T_1_lane2.fastq.gz,,auto,
Differentiated_rep2,Differentiated,2,data/raw/R1T_2.fastq.gz,,auto,
Differentiated_rep3,Differentiated,3,data/raw/R1T_3.fastq.gz,,auto,
Immature_rep1,Immature,1,data/raw/R2T_1.fastq.gz,,auto,
Immature_rep2,Immature,2,data/raw/R2T_2.fastq.gz,,auto,
Immature_rep3,Immature,3,data/raw/R2T_3.fastq.gz,,auto,
```

`contrasts.csv`:

```csv
name,treatment,reference
Differentiated_vs_Immature,Differentiated,Immature
```

Run step 1 and steps 2 to 4 as in the README (`submit_preprocess.sh`, then `submit_analysis.sh`).

Result: the lanes are read as one sample. `Differentiated_rep1` has 14,090,067 assigned reads, the same as the unsplit file.
The DE result on this run was 1786 up and 1146 down (padj < 0.05, |log2FC| >= 1).

## Variant 2: strandedness overrides

In `samples.csv`, set `Immature_rep1` to `unstranded` and `Immature_rep2` to `reverse` (the other samples stay `auto`):

```csv
Immature_rep1,Immature,1,data/raw/R2T_1.fastq.gz,,unstranded,
Immature_rep2,Immature,2,data/raw/R2T_2.fastq.gz,,reverse,
```

Restart only those two samples from `featurecounts`, because their strandedness is used in that step. Set the sample
filter in `config.sh` as a quoted string, then restart:

```bash
# in config.sh:  SAMPLE_FILTER="Immature_rep1 Immature_rep2"
bash scripts/submit_preprocess.sh --from featurecounts
# then set SAMPLE_FILTER back to empty
```

Check `results/preprocessed/strandedness/strandedness_report.tsv`. The third column shows the strandedness used:
`Immature_rep1` is `unstranded`, `Immature_rep2` is `reverse`, and the others are `forward`. The second column still shows
what Salmon inferred (`SF`).

Assigned reads from featureCounts on this run:

| Sample | Strandedness used | Assigned reads |
|---|---|---|
| Immature_rep1 | unstranded | 13,774,290 |
| Immature_rep2 | reverse | 1,284,497 |
| Immature_rep3 | auto (forward) | 14,354,192 |

Reversing a forward library drops the assigned reads by about 90 percent, as expected.

## Variant 3: cutoffs and gene sets

In `config.sh` set `PADJ_CUTOFF=0.01`, `LFC_CUTOFF=1.5` and `MSIGDB_COLLECTIONS=H,C2:CP:KEGG`, then run steps 2 to 4:

```bash
bash scripts/submit_analysis.sh
```

Wait until the jobs have started before you change `config.sh` back, because each job reads the config when it starts.

Results on this run:

- The report text says `padj < 0.01` and `|log2 fold change| >= 1.5`.
- DE: 877 up and 483 down.
- The enrichment ran on 236 gene sets (50 Hallmark and 186 KEGG). No KEGG term reached significance at padj < 0.01 in this comparison.

## Variant 4: optional steps off and a subset of samples

In `config.sh` set `RUN_FEATURECOUNTS=no`, `RUN_QUALIMAP=no` and `SAMPLE_FILTER="Differentiated_rep1 Differentiated_rep2 Immature_rep1 Immature_rep2"`. Keep the defaults for the cutoffs and gene sets. Then run steps 2 to 4:

```bash
bash scripts/submit_analysis.sh
bash scripts/status.sh
```

`status.sh` counts featureCounts and Qualimap as done, because they are not run. The analysis uses the four samples
(`Inputs OK for analysis: 4 sample(s)`), and the report's sample and read-count tables list only those four. The
"Assigned to genes" column is empty, because featureCounts was not run for this analysis.

Results on this run: DE 1664 up and 1007 down.

Restore `config.sh` to its defaults afterwards.
