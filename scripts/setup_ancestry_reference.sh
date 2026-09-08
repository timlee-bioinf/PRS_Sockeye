#!/usr/bin/env bash
# One-time setup for multi-ancestry mode: builds the ancestry reference cache
# (LD-pruned marker set, reference PCA basis + projection loadings, population
# classifier models, per-ancestry empirical + PC-regression normalization
# models) from the reference panel (ANCESTRY_REF_PFILE) and the ancestry
# SNP/weight lists (ANCESTRY_MAP).
#
# This is NOT part of the per-run submit.sh DAG - run it once (on the login
# node; wrap in salloc/sbatch yourself if the reference panel is large enough
# that LD-pruning + PCA need more than a login-node allocation) whenever the
# reference panel, ANCESTRY_MAP, or the ANCESTRY_N_* settings in config.sh
# change. step2b_ancestry.slurm reads this cache on every pipeline run.
#
# Usage: bash scripts/setup_ancestry_reference.sh
set -euo pipefail

SUBMIT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CONFIG="${CONFIG:-$SUBMIT_DIR/config.sh}"
source "$CONFIG"

[[ "$RUN_ANCESTRY" == "1" ]] || { echo "[ERROR] RUN_ANCESTRY=1 required (set it before sourcing config, e.g. in submit.local.sh)"; exit 1; }
[[ -f "$ANCESTRY_MAP" ]] || { echo "[ERROR] ANCESTRY_MAP not found: $ANCESTRY_MAP"; exit 1; }
[[ -f "${ANCESTRY_REF_PFILE}.pvar" || -f "${ANCESTRY_REF_PFILE}.pvar.zst" ]] || {
  echo "[ERROR] reference panel not found: ${ANCESTRY_REF_PFILE}.pvar[.zst]"; exit 1; }
command -v "$PLINK2" >/dev/null 2>&1 || [[ -x "$PLINK2" ]] || { echo "[ERROR] PLINK2 not found: $PLINK2"; exit 1; }

module load gcc/9.4.0 r/4.4.0 2>/dev/null || true

mkdir -p "$ANCESTRY_REF_CACHE"
LOG="$ANCESTRY_REF_CACHE/setup.log"; : > "$LOG"
log() { echo "$*" | tee -a "$LOG"; }

PLINK_EXCLUDE_OPT=()
R_EXCLUDE_OPT=()
if [[ -n "$ANCESTRY_REF_EXCLUDE" ]]; then
  [[ -f "$ANCESTRY_REF_EXCLUDE" ]] || { echo "[ERROR] ANCESTRY_REF_EXCLUDE not found: $ANCESTRY_REF_EXCLUDE"; exit 1; }
  PLINK_EXCLUDE_OPT=(--remove "$ANCESTRY_REF_EXCLUDE")
  R_EXCLUDE_OPT=(--exclude "$ANCESTRY_REF_EXCLUDE")
fi

log "=== [1/4] LD-pruning reference panel $(date -Is) ==="
"$PLINK2" --pfile "$ANCESTRY_REF_PFILE" "${PLINK_EXCLUDE_OPT[@]+"${PLINK_EXCLUDE_OPT[@]}"}" \
  --autosome --maf 0.05 --geno 0.05 \
  --indep-pairwise 200 50 0.1 \
  --out "$ANCESTRY_REF_CACHE/prune" >> "$LOG" 2>&1
log "  $(wc -l < "$ANCESTRY_REF_CACHE/prune.prune.in") ancestry-informative markers -> $ANCESTRY_REF_CACHE/prune.prune.in"

log "=== [2/4] Computing reference PCA ($ANCESTRY_N_PCS PCs) $(date -Is) ==="
# --freq alongside --pca caches the reference panel's own allele frequencies
# (ref_pca.afreq) so step2b_ancestry.slurm's projection --score can pass
# --read-freq and variance-standardize target genotypes against the SAME
# frequencies the PCA basis was built on, instead of re-estimating them from
# whatever (possibly small) cohort is being scored that run.
"$PLINK2" --pfile "$ANCESTRY_REF_PFILE" "${PLINK_EXCLUDE_OPT[@]+"${PLINK_EXCLUDE_OPT[@]}"}" \
  --extract "$ANCESTRY_REF_CACHE/prune.prune.in" \
  --pca allele-wts "$ANCESTRY_N_PCS" \
  --freq \
  --out "$ANCESTRY_REF_CACHE/ref_pca" >> "$LOG" 2>&1
log "  -> $ANCESTRY_REF_CACHE/ref_pca.eigenvec (+ .eigenvec.allele for projection, .afreq for --read-freq)"

log "=== [3/4] Fitting population classifier models $(date -Is) ==="
Rscript "$SUBMIT_DIR/scripts/build_pop_models.R" \
  --eigenvec "$ANCESTRY_REF_CACHE/ref_pca.eigenvec" \
  --psam "${ANCESTRY_REF_PFILE}.psam" \
  --label-col "$ANCESTRY_REF_LABEL_COL" \
  --npcs "$ANCESTRY_N_POPCOMP" \
  --out "$ANCESTRY_REF_CACHE/pop_models.rds" \
  --unrelated-out "$ANCESTRY_REF_CACHE/ref_unrelated.txt" \
  "${R_EXCLUDE_OPT[@]+"${R_EXCLUDE_OPT[@]}"}" >> "$LOG" 2>&1
log "  -> $ANCESTRY_REF_CACHE/pop_models.rds"

log "=== [4/4] Per-ancestry score files + normalization models $(date -Is) ==="
case "$SCORE_MODE" in
  both)       MODES=(weighted unweighted) ;;
  weighted)   MODES=(weighted) ;;
  unweighted) MODES=(unweighted) ;;
  *) echo "[ERROR] invalid SCORE_MODE: $SCORE_MODE"; exit 1 ;;
esac

while IFS=$'\t' read -r ANC SNP_IN; do
  [[ -z "$ANC" || "$ANC" == \#* ]] && continue
  log "--- ancestry: $ANC ($SNP_IN) ---"
  OUTDIR="$ANCESTRY_REF_CACHE/ancestry/$ANC"
  mkdir -p "$OUTDIR"
  Rscript "$SUBMIT_DIR/scripts/prepare_copa_score_files.R" \
    --input "$SNP_IN" --mode "$SCORE_MODE" --outdir "$OUTDIR" >> "$LOG" 2>&1

  for MODE in "${MODES[@]}"; do
    SFILE="$OUTDIR/copa_${MODE}_score.tsv"
    [[ -f "$SFILE" ]] || { echo "[ERROR] missing $SFILE"; exit 1; }
    "$PLINK2" --pfile "$ANCESTRY_REF_PFILE" "${PLINK_EXCLUDE_OPT[@]+"${PLINK_EXCLUDE_OPT[@]}"}" \
      --extract "$OUTDIR/snp_list.txt" \
      --score "$SFILE" 1 2 3 header cols=+scoresums no-mean-imputation \
      --out "$OUTDIR/ref_score_${MODE}" >> "$LOG" 2>&1
    Rscript "$SUBMIT_DIR/scripts/build_ancestry_norm_models.R" \
      --sscore "$OUTDIR/ref_score_${MODE}.sscore" \
      --eigenvec "$ANCESTRY_REF_CACHE/ref_pca.eigenvec" \
      --psam "${ANCESTRY_REF_PFILE}.psam" \
      --label-col "$ANCESTRY_REF_LABEL_COL" \
      --ancestry "$ANC" \
      --npcs "$ANCESTRY_N_NORM" \
      --out "$OUTDIR/norm_models_${MODE}.rds" \
      "${R_EXCLUDE_OPT[@]+"${R_EXCLUDE_OPT[@]}"}" >> "$LOG" 2>&1
  done
done < <(tr -d '\r' < "$ANCESTRY_MAP"; printf '\n')

log "=== DONE $(date -Is) - cache ready at $ANCESTRY_REF_CACHE ==="
