#!/bin/bash -l
#SBATCH --job-name=grs_merge
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
: "${RUN_ID:?RUN_ID not set}"

BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
TASK1_OUT="$BASE_OUT/task1_extract"
TASK2_OUT="$BASE_OUT/task2_merge"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$TASK2_OUT" "$LOG_DIR"
STEP_LOG="$LOG_DIR/merge.log"; : > "$STEP_LOG"

# Collect the per-chromosome pgens that actually have variants.
PGEN_LIST="$TASK2_OUT/pgen_prefixes.txt"; : > "$PGEN_LIST"
for chr in {1..22}; do
  PREFIX="$TASK1_OUT/chr${chr}"
  if [[ -f "${PREFIX}.pgen" ]]; then
    echo "$PREFIX" >> "$PGEN_LIST"
  else
    echo "[WARN] missing ${PREFIX}.pgen - skip chr${chr}" | tee -a "$STEP_LOG"
  fi
done
[[ -s "$PGEN_LIST" ]] || { echo "[ERROR] no per-chr pgen found in $TASK1_OUT" | tee -a "$STEP_LOG"; exit 1; }

N=$(wc -l < "$PGEN_LIST")
echo "[INFO] merging $N chromosome pgen(s)" | tee -a "$STEP_LOG"

"$PLINK2" \
  --pmerge-list "$PGEN_LIST" \
  --make-pgen \
  --out "$TASK2_OUT/all_snps" >> "$STEP_LOG" 2>&1

MERGED=$(grep -vc '^#' "$TASK2_OUT/all_snps.pvar")
echo "[DONE] merged -> $TASK2_OUT/all_snps ($MERGED variants)" | tee -a "$STEP_LOG"
