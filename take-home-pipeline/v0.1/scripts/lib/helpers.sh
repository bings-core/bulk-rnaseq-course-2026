# Shared helpers for every take-home script. Source it, do not run it.

# Load config.sh from the project root (the folder that holds config.sh) and work from there.
load_config() {
  # Run from the project root (where config.sh is); otherwise use the folder above the calling script.
  local root="$PWD"
  if [ ! -f "$root/config.sh" ]; then
    root="$(cd "$(dirname "${BASH_SOURCE[1]:-$0}")/.." && pwd)"
  fi
  if [ ! -f "$root/config.sh" ]; then
    echo "ERROR: config.sh not found in $root. Copy config.example.sh to config.sh and edit it." >&2
    return 1
  fi
  cd "$root" || return 1
  set -a; . "$root/config.sh"; set +a      # export every setting so the R scripts see them
  mkdir -p logs
}

# Print one line per sample: sample_id fastq_1 fastq_2(or -) strandedness.
# Skips rows not in SAMPLE_FILTER when that is set.
read_samples() {
  awk -F, -v filt="$SAMPLE_FILTER" 'NR>1 && $1!="" {
    if (filt != "" && index(" " filt " ", " " $1 " ") == 0) next
    print $1, $4, ($5==""?"-":$5), ($6==""?"auto":$6)
  }' "$SAMPLES_CSV"
}

# Species-specific reference paths. Fails clearly if a reference is missing.
species_refs() {
  REF="$RESOURCE_DIR/$SPECIES"
  STAR_INDEX="$REF/star/2.7.5b/index"
  SALMON_INDEX="$REF/salmon/1.2.1/index"
  GTF="$REF/gtf/annotation.gtf"
  for p in "$STAR_INDEX" "$SALMON_INDEX" "$GTF"; do
    if [ ! -e "$p" ]; then
      echo "ERROR: reference $p is missing. Check RESOURCE_DIR and SPECIES in config.sh." >&2
      return 1
    fi
  done
}

# Run a command, log it, and stop the script if it fails. The log line goes to stderr so $(submit ...) gets only the job ID.
run() {
  echo "[$(date +%H:%M:%S)] $*" >&2
  "$@" || { echo "FAILED: $*" >&2; return 1; }
}

# submit NAME CORES MEM_GB [DEPENDENCY] -- COMMAND...
#   With LSF_ACCOUNT set: submits one LSF job (named NAME) and prints its job ID.
#   Without it: runs the command here, so the same scripts work interactively.
#   DEPENDENCY is an LSF -w expression built from job IDs, e.g. "done(278498211) && done(278498213)".
#   Prints the LSF job ID (empty when run here).
submit() {
  local name="$1" cores="$2" mem="$3"; shift 3
  local dep=""
  if [ "$1" != "--" ]; then dep="$1"; shift; fi
  shift   # the --
  if [ -z "$LSF_ACCOUNT" ]; then
    run "$@"; return
  fi
  local args=(-P "$LSF_ACCOUNT" -q "$LSF_QUEUE" -n "$cores"
              -R "span[hosts=1]" -R "rusage[mem=$((mem*1024))]" -M "$((mem*1024))"
              -J "$name" -o "logs/$name.%J.out" -e "logs/$name.%J.err")
  [ -n "$LSF_WALLTIME" ] && args+=(-W "$LSF_WALLTIME")
  [ -n "$dep" ] && args+=(-w "$dep")
  local out jid
  out=$(bsub "${args[@]}" "$@") || { echo "bsub failed for $name" >&2; return 1; }
  echo "$out" >&2
  jid=$(echo "$out" | grep -o '<[0-9]*>' | tr -d '<>' | head -1)
  echo "$jid"
}
