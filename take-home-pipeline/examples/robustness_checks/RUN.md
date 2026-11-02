# Example checks: running without LSF, and a container that fails its checksum

Two checks of the pipeline's behaviour outside the normal LSF run. Both were run on Minerva on 2026-11-01 in the human
single-end test folder (see the human variants example). The killed-job check is in the mixed single-end and paired-end
example.

## 1. Run steps 2 to 4 without LSF

Without `LSF_ACCOUNT`, the scripts run every step in the shell they are started from, and submit nothing. A shell on the
login node should not run real data, so this check ran the steps inside one LSF job.

In `config.sh`, empty `LSF_ACCOUNT` (the line reads `LSF_ACCOUNT=`), then submit the analysis as one job:

```bash
bsub -P acc_BiNGS_bulk -q premium -n 8 -R "span[hosts=1]" -R "rusage[mem=5120]" -M 5120 -W 4:00 \
     -J nolsf_test -o logs/nolsf.%J.out -e logs/nolsf.%J.err \
     "cd $PWD && bash scripts/submit_analysis.sh"
```

Keep `LSF_ACCOUNT` empty until the job has started, because the job reads `config.sh` when it starts. Then set it back to
`acc_BiNGS_bulk`.

Result: the job ran the input check (`Inputs OK for analysis: 6 sample(s).`), then steps 2 to 4 one after another, wrote the
report and the run record, and submitted no further jobs. The job's output file has the report and the run-record messages.

Without `LSF_ACCOUNT` set, the check prints a warning (`LSF_ACCOUNT is empty: every step runs here on the login node`).

## 2. A container that fails its checksum

`setup_resources.sh` checks the SHA-256 of each container after the download. The check was tested with a deliberately wrong
expected checksum, in a scratch copy of the scripts (the real script is not changed):

```bash
sed -i 's/^MULTIQC_SHA256=.*/MULTIQC_SHA256=0000000000000000000000000000000000000000000000000000000000000000/' scripts/setup_resources.sh
bash scripts/setup_resources.sh --step multiqc
```

Result: the download finished, the checksum did not match, and the script stopped with:

```
ERROR: checksum mismatch for containers/multiqc_1.35.sif.part (expected 0000...0000). Delete containers/multiqc_1.35.sif.part and run the step again.
```

The script exited with status 1. The downloaded file was kept only as `multiqc_1.35.sif.part`. It was not saved as
`multiqc_1.35.sif`, so no step can use it. Delete the `.part` file, and run the step again with the correct checksum.
