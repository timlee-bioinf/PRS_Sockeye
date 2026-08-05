#!/bin/bash -l
#SBATCH --job-name=grs_extract
#SBATCH --array=1-22
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/extract_%A_%a.out
#SBATCH --error=logs/extract_%A_%a.err

set -euo pipefail

# CONFIG and RUN_ID are exported by submit.sh.
CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
: "${RUN_ID:?RUN_ID not set}"

CHR="${SLURM_ARRAY_TASK_ID}"
BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
INPUTS="$BASE_OUT/inputs"
TASK1_OUT="$BASE_OUT/task1_extract"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$TASK1_OUT" "$LOG_DIR"

SNP_LIST="$INPUTS/snp_list.txt"
LOGFILE="$LOG_DIR/extract_chr${CHR}.log"
MASTER_LOG="$LOG_DIR/extract_all.log"
LOCKFILE="$LOG_DIR/.extract.lock"

log_master() { { flock -x 200; echo "$*"; } 200>"$LOCKFILE" >> "$MASTER_LOG"; }

# Preflight on the first array task only.
if [[ "$CHR" -eq 1 ]]; then
  command -v "$PLINK2" >/dev/null 2>&1 || [[ -x "$PLINK2" ]] || { echo "[ERROR] PLINK2 not found: $PLINK2"; exit 1; }
  [[ -f "$SNP_LIST" ]]   || { echo "[ERROR] SNP list not found: $SNP_LIST (run submit.sh prep first)"; exit 1; }
  [[ -d "$IMPUTE_DIR" ]] || { echo "[ERROR] IMPUTE_DIR not found: $IMPUTE_DIR"; exit 1; }
fi

IN_VCF="$IMPUTE_DIR/chr_${CHR}/chr${CHR}.dose.vcf.gz"
OUT_PREFIX="$TASK1_OUT/chr${CHR}"

if [[ ! -f "$IN_VCF" ]]; then
  echo "[chr$CHR] missing VCF: $IN_VCF - skip" | tee -a "$LOGFILE"
  log_master "[chr$CHR] SKIP (missing vcf)"
  exit 0
fi

echo "=== [chr$CHR] extract start $(date -Is) ===" | tee -a "$LOGFILE"

# Pass the full rsID list to every chromosome: each per-chr VCF only contains
# its own variants, so only the matching rsIDs are extracted (no chr/pos needed).
# Split multiallelic records can share an rsID; keep the first copy.
set +e
"$PLINK2" \
  --vcf "$IN_VCF" dosage=DS \
  --extract "$SNP_LIST" \
  --rm-dup force-first \
  --make-pgen \
  --out "$OUT_PREFIX" >>"$LOGFILE" 2>&1
RC=$?
set -e

if [[ $RC -ne 0 ]]; then
  # 0 variants on this chromosome is expected, not a failure.
  if grep -qiE "No variants remaining|0 variants remaining|Error.*--extract" "$LOGFILE"; then
    echo "[chr$CHR] 0 variants after --extract - skip" | tee -a "$LOGFILE"
    log_master "[chr$CHR] SKIP (0 variants)"
    rm -f "${OUT_PREFIX}.pgen" "${OUT_PREFIX}.pvar" "${OUT_PREFIX}.psam"
    exit 0
  fi
  echo "[chr$CHR] PLINK2 failed (rc=$RC). See $LOGFILE" | tee -a "$LOGFILE"
  log_master "[chr$CHR] FAIL rc=$RC"
  exit $RC
fi

COUNT=0
[[ -f "${OUT_PREFIX}.pvar" ]] && COUNT=$(grep -vc '^#' "${OUT_PREFIX}.pvar" || echo 0)
echo "[chr$CHR] extracted $COUNT variants" | tee -a "$LOGFILE"
if [[ "$COUNT" -eq 0 ]]; then
  rm -f "${OUT_PREFIX}.pgen" "${OUT_PREFIX}.pvar" "${OUT_PREFIX}.psam"
fi
log_master "[chr$CHR] DONE variants=$COUNT"
