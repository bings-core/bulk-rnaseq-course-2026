# Bulk RNA-seq Analysis Course

Welcome! Everything for the ten-week BiNGS bulk RNA-seq course lives here — setup guides,
practice data, the analysis pipeline, and reading. Ask questions any time in the
[Issue](../../issues) tab.

> 🚨 **Urgent day-of issues** (can't access the building, can't find the room, etc.): call **Dan Hasson at 646-824-6449**.

- **Meets:** Wednesdays, 4:00–5:30 PM · Sep 9 – Nov 11, 2026 (10 sessions)
- **Platform:** Minerva HPC (project `BiNGS_bulk`, `/sc/arion/projects/BiNGS_bulk/`)
- **Organized by:** BiNGS Core · Tisch Cancer Center, Mount Sinai

## Before Session 1

Four things, in this order. Budget about half an hour in total. Please finish by
end of day **Friday September 5** so we have enough time to fix anything broken.

### 1. Set up Microsoft Authenticator

Mount Sinai requires multi-factor authentication on **every** Minerva and VPN
login, so nothing else on this list works without it. If you already approve a
phone prompt when you log into Sinai systems, you're done — skip to step 2.

1. Install **Microsoft Authenticator** on your phone, from the
   [App Store](https://apps.apple.com/app/microsoft-authenticator/id983156458) or
   [Google Play](https://play.google.com/store/apps/details?id=com.azure.authenticator).
2. Register your Mount Sinai account with the app by following the official
   Sinai walkthrough:
   **[itsecurity.mssm.edu/ms-authenticator](https://itsecurity.mssm.edu/ms-authenticator/)**

From then on, every Minerva and VPN login sends a prompt to your phone that you
approve, or shows a number you type into the app.

### 2. Get onto the Sinai network

Minerva is not reachable from the open internet. Your laptop has to be on the
Sinai network first — **on campus that means the wifi, anywhere else it means
the VPN**.

**On campus — connect to `MSMC-Green`**

In your wireless settings, pick `MSMC-Green` and sign in. The username format
depends on which side of Mount Sinai you're affiliated with, and this is the
part that trips people up:

| Affiliation | Username format | Example |
|-------------|-----------------|---------|
| **School** (faculty, staff, students) | `mssmcampus\<your-username>` | `mssmcampus\doej01` |
| **Hospital** (employees) | `msnyuhealth\<your-network-ID>` | `msnyuhealth\doej01` |

Note the backslash, and that this is your Sinai username, not your email
address. Full guide with per-device screenshots:
[Mount Sinai Wifi Instructions (PDF)](https://icahn.mssm.edu/files/ISMMS/Assets/About%20the%20School/Computer%20Services/Mount-Sinai-Wifi-Instructions.pdf).

**Off campus — install and connect the F5 VPN**

*One time,* install the F5 client for your operating system:
[Mac](https://itsecurity.mssm.edu/vpn-tunnel-for-mac/) ·
[Windows](https://itsecurity.mssm.edu/vpn-tunnel-for-windows/)

*Every time you work off campus,* connect through it:

1. Go to the VPN portal for your side of Mount Sinai —
   **school:** [dcsmsvpn.mssm.edu](https://dcsmsvpn.mssm.edu) ·
   **hospital:** [msvpn.mountsinai.org](https://msvpn.mountsinai.org)
2. Log in and approve the Microsoft Authenticator prompt.
3. Under **Network Access**, click **Tunnel**. The F5 client opens in its own
   window; when the indicator turns green and reads **Connected**, you're on.
4. Leave that window open while you work — closing it drops the tunnel.

Screenshots of both steps are in
[`minerva/README.md`](minerva/README.md#off-campus-vpn).

### 3. Run the one-time access check

**This is the important one.** It takes about two minutes and it's how you find
out whether your access actually works while there's still time to fix it.

**a. Open a terminal on your own laptop.**
On **Mac** that's the built-in **Terminal** app. On **Windows** it's
**Windows Terminal** (already installed on Win 11; free from the Microsoft
Store on Win 10) — it ships with `ssh`, so nothing else to install.

**b. Log into Minerva.** Replace `<your-minerva-username>` with your own:

```bash
ssh <your-minerva-username>@minerva.hpc.mssm.edu
```

You'll be asked for your **Sinai password**, then your **Microsoft
Authenticator** code or push approval. The password prompt shows nothing at
all as you type — no dots, no asterisks. That's normal; type it and press
Enter.

*(If your Sinai email ends in `@mountsinai.org` rather than `@mssm.edu`, use
`minerva-org.hpc.mssm.edu` instead.)*

When you see a prompt like `[jdoe01@li03c04 ~]$`, you're on Minerva.

**c. Run the check.** Copy this one line, paste it, press Enter:

```bash
bash /sc/arion/projects/BiNGS_bulk/shared/access_check.sh
```

It creates your personal course folder and submits a short test job, so give it
a few seconds. It's safe to run more than once.

**d. Check the output.** You should see all five checks pass:

```
============================================================
  BiNGS Bulk RNA-seq Course - Minerva access check
============================================================
  user: jdoe01    date: Wed Sep  2 17:39:34 EDT 2026
  login node: lc07e55
============================================================

[1/5] Your Minerva account
  OK       You are logged in as: jdoe01

[2/5] Course group membership
  OK       You are in the BiNGS_bulk group.

[3/5] Your personal course folder
  OK       Ready: /sc/arion/projects/BiNGS_bulk/jdoe01

[4/5] Submitting a test job (usually ~10 seconds, please wait)
  OK       Job ran under project account acc_BiNGS_bulk.

[5/5] Classroom reservation
  OK       Job landed on reserved node lc07e04 - reservation BINGS_1 works.
           R on the compute node: R version 4.2.0 (2022-04-22)

------------------------------------------------------------
  ALL CHECKS PASSED

  Your course folder:  /sc/arion/projects/BiNGS_bulk/jdoe01
  You are ready for Session 1. Nothing else to do.
------------------------------------------------------------
RESULT: PASS user=jdoe01 node=lc07e04 date=2026-09-02 17:39:38

  A copy of this output was saved to:
    /sc/arion/projects/BiNGS_bulk/jdoe01/access_check_20260902_173934.log

  You can now type 'exit' to log out of Minerva.
```

If you see **ALL CHECKS PASSED**, you're completely set up — account, folder,
job scheduler, classroom reservation and R.

The script also saves a copy of everything it printed into your course folder,
as `access_check_<date>_<time>.log`. Leave it there — it's how we confirm on our
side that everyone is set up, so you don't need to email us to say it worked.

If anything says **FAILED**, please
[open an Issue](../../issues/new?template=02-question.md) and paste the output
in. Don't wait and don't try to work around it — most of these are ours to fix
and need a day or two, which is exactly why we're asking you to do this now.

**e. Log out.** Type `exit` (or press **Ctrl-D**).

<sub>Curious what it did? The script is
[`minerva/access_check.sh`](minerva/access_check.sh) — it checks your group
membership, creates `/sc/arion/projects/BiNGS_bulk/<your-username>/`, and
submits a five-minute job to the classroom reservation that loads R and reports
which node it ran on.</sub>

### 4. Introduce yourself

**[Open an introduction issue →](../../issues/new?template=01-introduce-yourself.md)**

Tell us who you are, what you work on, and what you want to get out of the
course. This isn't busywork: your answers about prior programming experience and
whether you have your own data shape how we pitch the hands-on sessions and how
we plan Sessions 9 and 10. It also means that when you post a question in week
four, everyone already knows who you are. You'll need a free
[GitHub account](https://github.com/signup) if you don't have one.

## Asking questions

**All course questions go in the [Issues](../../issues) tab.** Not email — Issues.

This is deliberate. When you hit a problem, three other people usually have the
same one, and an issue answers all four of you at once and stays searchable for
the rest of the course. Instructors and TAs monitor the tab between sessions.
Office hours can be set up on request.

**You'll need a free GitHub account to post.** If you don't have one, sign up at
**[github.com/signup](https://github.com/signup)** — it takes two minutes, a
personal email address is fine, and you do not need to pay for anything. This is
also a genuinely useful thing to own as a scientist: it's where nearly all
bioinformatics software lives, and where you'll eventually put your own analysis
code.

Two issue templates are set up for you:

- **👋 Introduce yourself** — please do this before Session 1 (step 4 above).
- **❓ Ask a question** — for anything that isn't working.

No question is too basic here. Several of you told us on the interest form that
you have no programming experience at all, which is exactly who this course is
built for.

> One exception: **account, password and MFA problems** belong to Scientific
> Computing at **hpchelp@hpc.mssm.edu**, not to us. Anything about the course
> project itself — folder access, the `acc_BiNGS_bulk` account, the `BINGS_1`
> reservation — is ours.

## Schedule

| #  | Date         | Time         | Location   | Topic                                          | 
|----|--------------|--------------|------------|------------------------------------------------|
| 01 | Wed · Sep 9  | 4:00-5:30 PM | Hess 8-101 | Introduction to bulk RNA-seq & BiNGS           | 
| Optional | Wed · Sep 16 | 3:30-4:00 PM | Hess 5-101 | If you need help with logging into msmc-green wifi, HPC setup, logging into minerva, opening OnDemand, getting ready to run code during session 2, feel free to come in half an hour early and we will be there to help you!       | 
| 02 | Wed · Sep 16 | 4:00-5:30 PM | Hess 5-101 | Introduction to HPC & R in biostatistics       | 
| 03 | Wed · Sep 23 | 4:00-5:30 PM | Hess 5-101 | QC metrics & tools                             | 
| 04 | Wed · Sep 30 | 4:00-5:30 PM | Hess 5-101 | Preprocessing                                  |
| 05 | Wed · Oct 7  | 4:00-5:30 PM | Hess 5-101 | Differential expression analysis               |
| 06 | Wed · Oct 14 | 4:00-5:30 PM | James 12   | Functional enrichment analysis                 | 
| 07 | Wed · Oct 21 | 4:00-5:30 PM | James 12   | Interpreting results · public datasets · single-cell intro  |
| 08 | Wed · Oct 28 | 4:00-5:30 PM | Hess 5-101 | Integration: ATAC-seq / ChIP-seq & epigenetics |
| 09 | Wed · Nov 4  | 4:00-5:30 PM | Hess 5-101 | Bring your own data · full pipeline            | 
| 10 | Wed · Nov 11 | 4:00-5:30 PM | Hess 8-101 | Bring your own data · full pipeline            |

## Syllabus

### Session 1 — Introduction to bulk RNA-seq & BiNGS
- **Style:** Lecture
- **Content:** Overview of bulk RNA-seq: project design, library preparation, sequencing, alignment, and analysis. Participants learn what questions bulk RNA-seq is best suited to answer and when to reach for it, plus a brief tour of adjacent transcriptomic modalities (splicing analysis, long-read RNA-seq). Covers the most important outputs of the analysis and its key takeaways. The session also introduces the course project and the BiNGS Core.

### Session 2 — Introduction to HPC & R in biostatistics
- **Style:** Lecture + Hands-on
- **Content:** Combined tour of Mount Sinai's HPC (Minerva) and R for RNA-seq. On the HPC side: navigating the Minerva file system, basic Linux commands, submitting batch jobs (LSF), and what Singularity containers are for and why the pipeline uses them. On the R side: loading packages, running commands interactively vs from scripts, and preparing the packages the DEG and functional-analysis sessions will need. Participants finish the session with a working environment for the rest of the course.
- **Tools:** Linux, LSF, Singularity containers, OnDemand Code Server, OnDemand RStudio Server, R

### Session 3 — QC metrics & tools
- **Style:** Lecture
- **Content:** What "quality" means for a bulk RNA-seq experiment and how we measure it. Walkthrough of the QC metrics that matter at each stage — raw-read quality, adapter/base bias, alignment rate, duplication, gene-body coverage, strandedness, rRNA contamination, sample correlation, PCA, variance patterns — and the tools that emit them. Participants learn how to read a MultiQC report and what should trigger dropping or reprocessing a sample before differential-expression analysis.
- **Tools:** FastQC, MultiQC, Qualimap, samtools, deepTools, PCA / correlation heatmaps

### Session 4 — Preprocessing
- **Style:** Lecture + Hands-on
- **Content:** Raw FASTQs → gene- and transcript-level counts. Participants run the full preprocessing chain on chr5-subset reads: read QC, adapter/quality trimming, alignment to the reference genome, alignment-level QC, alignment-based counting, and lightweight-mapping quantification. The session ends with a MultiQC report aggregating every tool run and a comparison of gene-level (featureCounts) vs transcript-level (Salmon) count matrices, both valid DESeq2 inputs.
- **Tools:** FastQC, TrimGalore, STAR, subread/featureCounts, samtools, deepTools bamCoverage, Qualimap, Salmon, MultiQC

### Session 5 — Differential expression analysis
- **Style:** Lecture + Hands-on
- **Content:** From count matrix to a ranked list of differentially expressed genes. Participants define the sample-condition design, import Salmon quantifications, normalize with DESeq2, fit the model, and extract DEGs with shrunken log2 fold changes. Outputs include the DEG table, a top-DEG expression heatmap, a volcano plot, an MA plot, per-gene expression boxplots, and a genome-wide karyogram of DEG log2FC — the standard figure panel for a DEG story.
- **Tools:** tximport, DESeq2, apeglm, ggplot2, ggrepel, ComplexHeatmap, circlize, patchwork

### Session 6 — Functional enrichment analysis
- **Style:** Lecture + Hands-on
- **Content:** From a DEG list to biology. Participants run both over-representation analysis (ORA) on up- and down-regulated DEGs and Gene Set Enrichment Analysis (GSEA) on the ranked full gene list against MSigDB Hallmark gene sets, then interpret the results as dotplots, cnetplots, running-enrichment plots, and leading-edge heatmaps for the top pathways (MYC targets, interferon-gamma response). The session closes by connecting the enrichment story back to the biological question introduced in Session 1.
- **Tools:** clusterProfiler, msigdbr, enrichplot, ggplot2, Enrichr

### Session 7 — Presenting & interpreting · single-cell intro
- **Style:** Lecture
- **Content:** How to tell the story once the analysis is done: which QC metrics belong in a figure, how to justify sample removal from PCA, framing the biological question and background, presenting DEGs and highlighting the most relevant functional-analysis output, and writing the methods section for the whole pipeline. Also covers where to find public bulk RNA-seq datasets (GEO, SRA, TCGA/cBioPortal, UCSC Genome Browser) to complement or benchmark your own, and closes with a brief look at how the same experimental question would be tackled with single-cell RNA-seq.
- **Tools:** GEO, SRA, TCGA / cBioPortal

### Session 8 — Integration: ATAC-seq / ChIP-seq & epigenetics
- **Style:** Lecture + Hands-on
- **Content:** How chromatin-accessibility (ATAC-seq) and TF/histone-mark (ChIP-seq, CUT&RUN) data can be layered onto RNA-seq results to explain *why* genes are differentially expressed. Participants work with matched ATAC and ChIP peak sets on the course dataset, cross-reference peaks with DEGs, visualize coverage at DEG loci, and annotate peak-gene associations. The session covers when to reach for each assay and how to design an experiment that supports integrated analysis from the start.
- **Tools:** deepTools, ChIPseeker, bedtools, IGV, UCSC Genome Browser

### Sessions 9 & 10 — Bring your own data · full pipeline
- **Style:** Hands-on
- **Content:** Two open sessions where participants run the full production BiNGS bulk RNA-seq pipeline on data they bring themselves. TAs and instructors work through it alongside each participant end-to-end: sample-metadata setup, preprocessing, DEG analysis, functional enrichment, and (where relevant) epigenetic integration. Goal: leave the course with a reproducible analysis of your own data, not just the course dataset.
- **Tools:** the full BiNGS bulk RNA-seq pipeline (all tools from Sessions 4–8)

## Working on the course

Launch [OnDemand Code Server](minerva/README.md#ondemand-code-server) (bash
sessions) or [OnDemand RStudio Server](minerva/README.md#ondemand-rstudio-server)
(R sessions), and set **Reservation ID** to **`BINGS_1`**.

**`BINGS_1` is live continuously from Wed Sep 2, 3:55 PM to Wed Nov 12,
5:30 PM** — the whole course, not just Wednesday class hours. Use it any
time you're working on course material, including evenings and weekends.
It puts your job on the two nodes we've reserved for the class, so it
usually starts immediately instead of queuing behind the rest of the
institution. Full details: [Using the reservation](minerva/README.md#using-the-reservation).

