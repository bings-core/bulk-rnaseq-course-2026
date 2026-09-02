#!/bin/bash
# =============================================================================
#  BiNGS Bulk RNA-seq Course (Fall 2026) - one-time Minerva access check
#
#  Run this once, before Session 1, from a Minerva login node:
#
#      bash /sc/arion/projects/BiNGS_bulk/shared/access_check.sh
#
#  It checks five things, creates your personal course folder, and submits a
#  short test job to the classroom reservation. Takes about 30 seconds.
#  It changes nothing except creating your own folder, and it is safe to
#  run more than once.
# =============================================================================

PROJECT_DIR="/sc/arion/projects/BiNGS_bulk"
MY_DIR="${PROJECT_DIR}/${USER}"
ACCOUNT="acc_BiNGS_bulk"
RESERVATION="BINGS_1"
FAILED=0

pass() { printf '  \033[32mOK\033[0m       %s\n' "$1"; }
fail() { printf '  \033[31mFAILED\033[0m   %s\n' "$1"; FAILED=1; }
info() { printf '           %s\n' "$1"; }

echo
echo "============================================================"
echo "  BiNGS Bulk RNA-seq Course - Minerva access check"
echo "============================================================"
echo

# --- 1. Who am I -------------------------------------------------------------
echo "[1/5] Your Minerva account"
if [ -n "$USER" ]; then
    pass "You are logged in as: $USER"
else
    fail "Could not determine your username."
fi
echo

# --- 2. Group membership -----------------------------------------------------
echo "[2/5] Course group membership"
if id -nG 2>/dev/null | tr ' ' '\n' | grep -qx "BiNGS_bulk"; then
    pass "You are in the BiNGS_bulk group."
else
    fail "You are NOT in the BiNGS_bulk group."
    info "-> Open an issue at github.com/bings-core/bulk-rnaseq-course-2026/issues"
    info "   This one is ours to fix, not yours."
fi
echo

# --- 3. Personal folder ------------------------------------------------------
echo "[3/5] Your personal course folder"
if mkdir -p "$MY_DIR" 2>/dev/null && [ -w "$MY_DIR" ]; then
    pass "Ready: $MY_DIR"
    echo "Welcome $USER - created $(date)" > "${MY_DIR}/welcome.txt"
else
    fail "Could not create or write to $MY_DIR"
    info "-> Usually means check 2 failed. Open an issue."
fi
echo

# --- 4 & 5. Test job on the reservation --------------------------------------
echo "[4/5] Submitting a test job (this takes ~30 seconds, please wait)"
if [ "$FAILED" -eq 1 ]; then
    info "Skipped - fix the failures above first."
    echo
else
    cd "$MY_DIR" || exit 1
    JOB_OUT="${MY_DIR}/access_check_job.out"
    JOB_LSF="${MY_DIR}/.access_check_job.lsf"
    rm -f "$JOB_OUT"

    cat > "$JOB_LSF" <<'LSFJOB'
#!/bin/bash
#BSUB -J course_access_check
#BSUB -P acc_BiNGS_bulk
#BSUB -q premium
#BSUB -U BINGS_1
#BSUB -n 1
#BSUB -W 00:05
#BSUB -R rusage[mem=2000]

echo "CHECK_NODE=$(hostname)"
module purge
module load R/4.2.0
Rscript -e 'cat("CHECK_R=", R.version.string, "\n", sep="")'
LSFJOB

    bsub -K -o "$JOB_OUT" -e /dev/null < "$JOB_LSF" > /dev/null 2>&1

    # Only look at what the job printed, not LSF's echo of the script itself.
    RESULT=$(sed -n '/The output (if any) follows:/,$p' "$JOB_OUT" 2>/dev/null)

    if [ -f "$JOB_OUT" ] && grep -q "Successfully completed" "$JOB_OUT"; then
        pass "Job ran under project account $ACCOUNT."
        echo
        echo "[5/5] Classroom reservation"
        NODE=$(printf '%s\n' "$RESULT" | grep -m1 '^CHECK_NODE=' | cut -d= -f2 | tr -d '[:space:]')
        RVER=$(printf '%s\n' "$RESULT" | grep -m1 '^CHECK_R='    | cut -d= -f2-)
        case "$NODE" in
            lc07e04|lc07e05) pass "Job landed on reserved node $NODE - reservation $RESERVATION works." ;;
            "")              fail "Job ran but its output could not be read. See $JOB_OUT" ;;
            *)               pass "Job ran on $NODE (not a reserved node, but it ran)." ;;
        esac
        [ -n "$RVER" ] && info "R on the compute node: $RVER"
        rm -f "$JOB_LSF"
    else
        fail "The test job did not complete."
        info "-> Look at $JOB_OUT for the reason, and paste it into an issue at"
        info "   github.com/bings-core/bulk-rnaseq-course-2026/issues"
        if [ -f "$JOB_OUT" ]; then
            ERRLINE=$(grep -m1 -i 'not a member\|not valid\|error' "$JOB_OUT")
            [ -n "$ERRLINE" ] && info "   (first error line: $ERRLINE)"
        fi
    fi
    echo
fi

# --- Summary -----------------------------------------------------------------
echo "------------------------------------------------------------"
if [ "$FAILED" -eq 0 ]; then
    echo "  ALL CHECKS PASSED"
    echo
    echo "  Your course folder:  $MY_DIR"
    echo "  You are ready for Session 1. Nothing else to do."
else
    echo "  SOME CHECKS FAILED"
    echo
    echo "  Please open an issue with the output above:"
    echo "  github.com/bings-core/bulk-rnaseq-course-2026/issues"
    echo "  We would much rather fix this now than on September 9."
fi
echo "------------------------------------------------------------"
echo
echo "  You can now type 'exit' to log out of Minerva."
echo
