#!/bin/bash -l
#SBATCH --job-name=grs_ancestry
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/ancestry_%j.out
#SBATCH --error=logs/ancestry_%j.err

# Only submitted by submit.sh when RUN_ANCESTRY=1. Projects this run's merged
# genotypes onto the cached reference PCA (see scripts/setup_ancestry_reference.sh),
# classifies each sample's most-similar reference population (Mahalanobis
# distance), and writes one --keep list per ancestry for step3_score.slurm.

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
: "${RUN_ID:?RUN_ID not set}"
SUBMIT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

[[ "$RUN_ANCESTRY" == "1" ]] || { echo "[ERROR] step2b_ancestry.slurm run with RUN_ANCESTRY!=1 - should not have been submitted"; exit 1; }

module load gcc/9.4.0 r/4.4.0 2>/dev/null || true

BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
TASK2_OUT="$BASE_OUT/task2_merge"
TASK2B_OUT="$BASE_OUT/task2b_ancestry"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$TASK2B_OUT" "$LOG_DIR"
STEP_LOG="$LOG_DIR/ancestry.log"; : > "$STEP_LOG"

PFILE="$TASK2_OUT/all_snps"
[[ -f "${PFILE}.pvar" ]] || { echo "[ERROR] merged pvar missing: ${PFILE}.pvar"; exit 1; }

REF_PCA="$ANCESTRY_REF_CACHE/ref_pca"
ALLELE_FILE="${REF_PCA}.eigenvec.allele"
FREQ_FILE="${REF_PCA}.afreq"
PRUNE_LIST="$ANCESTRY_REF_CACHE/prune.prune.in"
POP_MODELS="$ANCESTRY_REF_CACHE/pop_models.rds"
for f in "$ALLELE_FILE" "$FREQ_FILE" "$PRUNE_LIST" "$POP_MODELS"; do
  [[ -f "$f" ]] || { echo "[ERROR] missing $f - run scripts/setup_ancestry_reference.sh first"; exit 1; }
done
ANC_MAP="$(ancestry_map)" || exit 1
ANCESTRY_LIST=$(cut -f1 <<< "$ANC_MAP" | paste -sd, -)
[[ -n "$ANCESTRY_LIST" ]] || { echo "[ERROR] no ancestries found in $WEIGHTS_DIR"; exit 1; }
MIN_PCS=$(( ANCESTRY_N_POPCOMP > ANCESTRY_N_NORM ? ANCESTRY_N_POPCOMP : ANCESTRY_N_NORM ))

# .eigenvec.allele column layout was not verified against a live plink2
# binary - detect ID/A1/PC1 column positions from its own header rather than
# hardcoding numbers, and fail loudly (rather than silently misscoring) if
# they can't be found.
HEADER=$(head -1 "$ALLELE_FILE")
ID_COL=$(awk -v h="$HEADER" 'BEGIN{n=split(h,a,/[ \t]+/); for(i=1;i<=n;i++) if(a[i]=="ID") print i}')
A1_COL=$(awk -v h="$HEADER" 'BEGIN{n=split(h,a,/[ \t]+/); for(i=1;i<=n;i++) if(a[i]=="A1") print i}')
PC1_COL=$(awk -v h="$HEADER" 'BEGIN{n=split(h,a,/[ \t]+/); for(i=1;i<=n;i++) if(a[i]=="PC1") print i}')
LAST_COL=$(awk -v h="$HEADER" 'BEGIN{print split(h,a,/[ \t]+/)}')
if [[ -z "$ID_COL" || -z "$A1_COL" || -z "$PC1_COL" ]]; then
  echo "[ERROR] could not find ID/A1/PC1 columns in $ALLELE_FILE header:" | tee -a "$STEP_LOG"
  echo "  $HEADER" | tee -a "$STEP_LOG"
  echo "  ID_COL=$ID_COL A1_COL=$A1_COL PC1_COL=$PC1_COL" | tee -a "$STEP_LOG"
  exit 1
fi
echo "[INFO] .eigenvec.allele columns: ID=$ID_COL A1=$A1_COL PC1=$PC1_COL last=$LAST_COL" | tee -a "$STEP_LOG"

echo "=== [ancestry] projection start $(date -Is) ===" | tee -a "$STEP_LOG"
# --read-freq: variance-standardize target genotypes using the REFERENCE
# panel's own allele frequencies (cached at setup time), not frequencies
# re-estimated from whatever cohort is being scored this run - required for
# the projection to land in the same PC space the classifier/normalization
# models were fit in.
"$PLINK2" \
  --pfile "$PFILE" \
  --extract "$PRUNE_LIST" \
  --read-freq "$FREQ_FILE" \
  --score "$ALLELE_FILE" "$ID_COL" "$A1_COL" header-read no-mean-imputation variance-standardize \
  --score-col-nums "${PC1_COL}-${LAST_COL}" \
  --out "$TASK2B_OUT/target_proj" >> "$STEP_LOG" 2>&1

echo "[INFO] classifying ancestry" | tee -a "$STEP_LOG"
Rscript "$SUBMIT_DIR/scripts/classify_ancestry.R" \
  --proj-sscore "$TASK2B_OUT/target_proj.sscore" \
  --pop-models "$POP_MODELS" \
  --pval-thresh "$ANCESTRY_PVAL_THRESH" \
  --ancestries "$ANCESTRY_LIST" \
  --min-pcs "$MIN_PCS" \
  --out-calls "$TASK2B_OUT/ancestry_calls.tsv" \
  --out-keep-dir "$TASK2B_OUT" >> "$STEP_LOG" 2>&1

echo "[DONE] ancestry calls in $TASK2B_OUT/ancestry_calls.tsv" | tee -a "$STEP_LOG"
