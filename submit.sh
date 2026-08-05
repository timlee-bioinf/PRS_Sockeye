#!/usr/bin/env bash
# Orchestrator: prepare score files (login node), then submit the SLURM DAG
#   step1 extract (array) -> step2 merge -> step3 score
# Run via submit.local.sh (which sets your account/paths), or set the env vars
# yourself and run `bash submit.sh`. No personal info lives in this file.
set -euo pipefail

SUBMIT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG="${CONFIG:-$SUBMIT_DIR/config.sh}"
source "$CONFIG"

RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
INPUTS="$BASE_OUT/inputs"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$INPUTS" "$LOG_DIR"

# Duplicate guard: refuse if a pipeline run is already queued/running.
JOB_RE='^(grs_extract|grs_merge|grs_score)$'
if squeue -u "$USER" -h -o "%j %T" 2>/dev/null \
   | awk -v re="$JOB_RE" '$1 ~ re && ($2=="R"||$2=="PD"){f=1} END{exit(f?0:1)}'; then
  echo "[ERROR] a grs_* pipeline job is already running/pending. Check: squeue -u $USER"
  exit 1
fi

# --- prep on the login node: build snp_list + score file(s) before the array ---
echo "[PREP] preparing inputs from: $SNP_INPUT (mode=$SCORE_MODE)"
module load gcc/9.4.0 r/4.4.0 2>/dev/null || true
Rscript "$SUBMIT_DIR/scripts/prepare_copa_score_files.R" \
  --input "$SNP_INPUT" --mode "$SCORE_MODE" --outdir "$INPUTS"
[[ -f "$INPUTS/snp_list.txt" ]] || { echo "[ERROR] prep failed: no snp_list.txt"; exit 1; }
echo "[PREP] $(wc -l < "$INPUTS/snp_list.txt") SNPs in snp_list.txt"

# Account is taken from SBATCH_ACCOUNT in the environment (set in submit.local.sh).
# Optional email notifications: set MAIL_USER (and MAIL_TYPE) in submit.local.sh.
# Kept here, not in the committed slurm files, so no personal email is committed.
MAIL_ARGS=()
if [[ -n "${MAIL_USER:-}" ]]; then
  MAIL_ARGS=(--mail-user="$MAIL_USER" --mail-type="${MAIL_TYPE:-ALL}")
fi
sub() {
  sbatch --parsable --chdir="$BASE_OUT" \
    "${MAIL_ARGS[@]+"${MAIL_ARGS[@]}"}" \
    --export=ALL,CONFIG="$CONFIG",RUN_ID="$RUN_ID" "$@"
}

STEP1=$(sub --output="$LOG_DIR/extract_%A_%a.out" --error="$LOG_DIR/extract_%A_%a.err" \
            "$SUBMIT_DIR/slurm/step1_extract.slurm")
STEP2=$(sub --dependency=afterok:"$STEP1" \
            --output="$LOG_DIR/merge_%j.out" --error="$LOG_DIR/merge_%j.err" \
            "$SUBMIT_DIR/slurm/step2_merge.slurm")
STEP3=$(sub --dependency=afterok:"$STEP2" \
            --output="$LOG_DIR/score_%j.out" --error="$LOG_DIR/score_%j.err" \
            "$SUBMIT_DIR/slurm/step3_score.slurm")
STEP4=$(sub --dependency=afterok:"$STEP2" \
            --output="$LOG_DIR/report_%j.out" --error="$LOG_DIR/report_%j.err" \
            "$SUBMIT_DIR/slurm/step4_report.slurm")

echo "[DONE] RUN_ID=$RUN_ID"
echo "  step1(extract)=$STEP1  step2(merge)=$STEP2  step3(score)=$STEP3  step4(report)=$STEP4"
echo "  outputs: $BASE_OUT"
echo "  logs:    $LOG_DIR"
