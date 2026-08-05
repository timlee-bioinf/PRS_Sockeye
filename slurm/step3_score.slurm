#!/bin/bash -l
#SBATCH --job-name=grs_score
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/score_%j.out
#SBATCH --error=logs/score_%j.err

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
: "${RUN_ID:?RUN_ID not set}"

BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
INPUTS="$BASE_OUT/inputs"
TASK2_OUT="$BASE_OUT/task2_merge"
TASK3_OUT="$BASE_OUT/task3_score"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$TASK3_OUT" "$LOG_DIR"
STEP_LOG="$LOG_DIR/score.log"; : > "$STEP_LOG"

PFILE="$TASK2_OUT/all_snps"
PVAR="$PFILE.pvar"
[[ -f "$PVAR" ]] || { echo "[ERROR] merged pvar missing: $PVAR"; exit 1; }

# Optional R2 filter (applied at scoring time, not extraction). R2_THRESH=0 -> keep all.
EXTRACT_OPT=()
if awk "BEGIN{exit !($R2_THRESH > 0)}"; then
  KEEP="$TASK3_OUT/keep_r2_${R2_THRESH}.txt"
  awk -v th="$R2_THRESH" '
    /^##/      { next }
    $1=="#CHROM"{ for(i=1;i<=NF;i++) if($i=="INFO") ic=i; next }
    /^#/       { next }
    {
      info=$(ic); r="";
      n=split(info, a, ";");
      for (i=1; i<=n; i++) if (a[i] ~ /^R2=/) r=substr(a[i], 4);
      if (r != "" && r+0 >= th) print $3
    }' "$PVAR" > "$KEEP"
  echo "[INFO] R2>=$R2_THRESH keeps $(wc -l < "$KEEP") of $(grep -vc '^#' "$PVAR") variants" | tee -a "$STEP_LOG"
  EXTRACT_OPT=(--extract "$KEEP")
else
  echo "[INFO] R2_THRESH=$R2_THRESH -> no R2 filter (keep all)" | tee -a "$STEP_LOG"
fi

# no-mean-imputation: small cohort (<50) blocks PLINK's auto allele-frequency
# imputation; imputed dosages have ~no missingness, so score observed alleles only.
run_score() {
  local sfile="$1" out="$2" label="$3"
  [[ -f "$sfile" ]] || { echo "[ERROR] missing score file: $sfile"; exit 1; }
  echo "[INFO] scoring: $label" | tee -a "$STEP_LOG"
  "$PLINK2" \
    --pfile "$PFILE" \
    "${EXTRACT_OPT[@]+"${EXTRACT_OPT[@]}"}" \
    --score "$sfile" 1 2 3 header cols=+scoresums no-mean-imputation \
    --out "$out" >> "$STEP_LOG" 2>&1
}

case "$SCORE_MODE" in
  both)       DO_W=1; DO_U=1 ;;
  weighted)   DO_W=1; DO_U=0 ;;
  unweighted) DO_W=0; DO_U=1 ;;
  *) echo "[ERROR] invalid SCORE_MODE: $SCORE_MODE (both|weighted|unweighted)"; exit 1 ;;
esac

[[ "$DO_U" -eq 1 ]] && run_score "$INPUTS/copa_unweighted_score.tsv" "$TASK3_OUT/grs_unweighted" "unweighted"
[[ "$DO_W" -eq 1 ]] && run_score "$INPUTS/copa_weighted_score.tsv"   "$TASK3_OUT/grs_weighted"   "weighted"

echo "[DONE] scores in $TASK3_OUT" | tee -a "$STEP_LOG"
ls -1 "$TASK3_OUT"/*.sscore 2>/dev/null | tee -a "$STEP_LOG" || true
