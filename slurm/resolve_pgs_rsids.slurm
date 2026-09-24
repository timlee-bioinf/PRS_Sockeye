#!/bin/bash -l
#SBATCH --job-name=grs_resolve_rsids
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# One-time conversion (per trait, NOT part of the per-run DAG): turn
# position-only weights files into rsID weights files the pipeline can use.
#
# Input : $WEIGHTS_DIR/pgs_raw/<ANC>.txt.gz   (e.g. PGS Catalog *_hmPOS_GRCh38 files)
# Output: $WEIGHTS_DIR/<ANC>.tsv               rsID, effect_allele, other_allele,
#                                              effect_weight, chr, pos
#         $WEIGHTS_DIR/pgs_raw/<ANC>.unmatched.tsv  what was dropped, and why
#         $WEIGHTS_DIR/pgs_raw/resolve.log
#
# Each weights variant is looked up by GRCh38 chr + pos + alleles in the imputed
# VCFs (IMPUTE_DIR) and gets that record's rsID. Variants absent from the
# imputed data, or present but with ID ".", are dropped and listed.
#
# Submit with: ./submit.local.sh rsids

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
SUBMIT_DIR="${SUBMIT_DIR:?SUBMIT_DIR not exported}"

module load gcc/9.4.0 r/4.4.0 2>/dev/null || true

RAW_DIR="$WEIGHTS_DIR/pgs_raw"
LOG="$RAW_DIR/resolve.log"; : > "$LOG"
TMP="$RAW_DIR/.resolve_tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
log() { echo "$*" | tee -a "$LOG"; }

command -v "$PLINK2" >/dev/null 2>&1 || [[ -x "$PLINK2" ]] || { log "[ERROR] PLINK2 not found: $PLINK2"; exit 1; }
[[ -d "$IMPUTE_DIR" ]] || { log "[ERROR] IMPUTE_DIR not found: $IMPUTE_DIR"; exit 1; }

shopt -s nullglob
FILES=("$RAW_DIR"/*.txt.gz "$RAW_DIR"/*.txt "$RAW_DIR"/*.tsv.gz "$RAW_DIR"/*.tsv)
shopt -u nullglob
[[ ${#FILES[@]} -gt 0 ]] || { log "[ERROR] no weights files in $RAW_DIR (expected e.g. EUR.txt.gz)"; exit 1; }

log "=== [1/3] collecting positions from ${#FILES[@]} file(s) $(date -Is) ==="
Rscript "$SUBMIT_DIR/scripts/resolve_pgs_rsids.R" ranges --out-bed "$TMP/ranges.bed" "${FILES[@]}" 2>&1 | tee -a "$LOG"

log "=== [2/3] looking up positions in the imputed VCFs $(date -Is) ==="
: > "$TMP/imputed.pvar"
for CHR in $(cut -f1 "$TMP/ranges.bed" | sort -u); do
  IN_VCF="$IMPUTE_DIR/chr_${CHR}/chr${CHR}.dose.vcf.gz"
  if [[ ! -f "$IN_VCF" ]]; then
    log "[chr$CHR] missing VCF: $IN_VCF - its variants will be unmatched"
    continue
  fi
  set +e
  "$PLINK2" --vcf "$IN_VCF" \
    --extract bed1 "$TMP/ranges.bed" \
    --make-pgen \
    --threads "${SLURM_CPUS_PER_TASK:-1}" \
    --out "$TMP/chr${CHR}" > "$TMP/chr${CHR}.plink.log" 2>&1
  RC=$?
  set -e
  cat "$TMP/chr${CHR}.plink.log" >> "$LOG"
  if [[ $RC -ne 0 ]]; then
    if grep -qiE "No variants remaining|0 variants remaining" "$TMP/chr${CHR}.plink.log"; then
      log "[chr$CHR] 0 positions found"
      continue
    fi
    log "[ERROR] plink2 failed on chr$CHR (rc=$RC) - see $LOG"; exit $RC
  fi
  grep -v '^#' "$TMP/chr${CHR}.pvar" | cut -f1-5 >> "$TMP/imputed.pvar"
  log "[chr$CHR] $(grep -vc '^#' "$TMP/chr${CHR}.pvar") imputed record(s) at weights positions"
  rm -f "$TMP/chr${CHR}".{pgen,psam}
done

log "=== [3/3] writing rsID weights files $(date -Is) ==="
for f in "${FILES[@]}"; do
  base="$(basename "$f")"; label="${base%.gz}"; label="${label%.*}"
  Rscript "$SUBMIT_DIR/scripts/resolve_pgs_rsids.R" join \
    --pgs "$f" --pvar "$TMP/imputed.pvar" \
    --out "$WEIGHTS_DIR/${label}.tsv" \
    --unmatched "$RAW_DIR/${label}.unmatched.tsv" 2>&1 | tee -a "$LOG"
done

log "=== DONE $(date -Is) - rsID weights files in $WEIGHTS_DIR ==="
