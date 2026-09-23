#!/usr/bin/env bash
# Orchestrator: prepare score files (login node), then submit the SLURM DAG
#   step1 extract (array) -> step2 merge -> step3 score
# Run via submit.local.sh (which sets your account/paths - see
# submit.local.sh.example), or set the env vars yourself and run `bash submit.sh`.
# No personal info lives in this file.
set -euo pipefail

SUBMIT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG="${CONFIG:-$SUBMIT_DIR/config.sh}"
source "$CONFIG"

RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
INPUTS="$BASE_OUT/inputs"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$INPUTS" "$LOG_DIR"

# Prune old run directories (RETENTION_DAYS=0 disables). Only touches this
# run's siblings under RUN_BASE, never the one just created above.
if [[ "${RETENTION_DAYS:-0}" -gt 0 ]]; then
  find "$RUN_BASE" -maxdepth 1 -mindepth 1 -type d -name 'SNP_extract_*' \
    -mtime "+${RETENTION_DAYS}" -exec rm -rf -- {} +
fi

# Duplicate guard: refuse if a pipeline run is already queued/running.
JOB_RE='^(grs_extract|grs_merge|grs_ancestry|grs_score)$'
if squeue -u "$USER" -h -o "%j %T" 2>/dev/null \
   | awk -v re="$JOB_RE" '$1 ~ re && ($2=="R"||$2=="PD"){f=1} END{exit(f?0:1)}'; then
  echo "[ERROR] a grs_* pipeline job is already running/pending. Check: squeue -u $USER"
  exit 1
fi

# --- prep on the login node: build snp_list + score file(s) before the array ---
echo "[PREP] preparing inputs (trait=$TRAIT, mode=$SCORE_MODE, RUN_ANCESTRY=$RUN_ANCESTRY)"
module load gcc/9.4.0 r/4.4.0 2>/dev/null || true

if [[ "$RUN_ANCESTRY" == "1" ]]; then
  # Multi-ancestry: build each ancestry's own score files (reusing
  # prepare_copa_score_files.R unchanged, once per weights file in
  # WEIGHTS_DIR), then union their SNPs (+ the cached ancestry-PCA marker set)
  # into a single snp_list.txt / copa_score_variants_qc.tsv so
  # step1_extract.slurm and step4_report.slurm need no changes at all.
  ANC_MAP="$(ancestry_map)" || exit 1
  [[ -f "$ANCESTRY_REF_CACHE/prune.prune.in" && -f "$ANCESTRY_REF_CACHE/pop_models.rds" \
     && -f "$ANCESTRY_REF_CACHE/ref_pca.eigenvec.allele" && -f "$ANCESTRY_REF_CACHE/ref_pca.afreq" ]] || {
    echo "[ERROR] ancestry reference cache missing/incomplete at $ANCESTRY_REF_CACHE - run scripts/setup_ancestry_reference.sh first"
    exit 1
  }
  # Cache is built per-trait and per-mode (see setup_ancestry_reference.sh) and
  # can go stale if WEIGHTS_DIR gains a file or SCORE_MODE changes without re-running setup -
  # check every norm-model file this run will actually need BEFORE submitting
  # the (expensive) DAG, not after step3 has already scored everything.
  case "$SCORE_MODE" in
    both)       CHECK_MODES=(weighted unweighted) ;;
    weighted)   CHECK_MODES=(weighted) ;;
    unweighted) CHECK_MODES=(unweighted) ;;
    *) echo "[ERROR] invalid SCORE_MODE: $SCORE_MODE"; exit 1 ;;
  esac
  MISSING_CACHE=0
  while IFS=$'\t' read -r ANC SNP_IN; do
    [[ -z "$ANC" || "$ANC" == \#* ]] && continue
    for MODE in "${CHECK_MODES[@]}"; do
      [[ -f "$ANCESTRY_NORM_DIR/$ANC/norm_models_${MODE}.rds" ]] || {
        echo "[ERROR] missing $ANCESTRY_NORM_DIR/$ANC/norm_models_${MODE}.rds" \
             "- re-run scripts/setup_ancestry_reference.sh (new trait, new weights file, or SCORE_MODE changed, since it was last built?)"
        MISSING_CACHE=1
      }
    done
  done <<< "$ANC_MAP"
  [[ "$MISSING_CACHE" -eq 0 ]] || exit 1

  : > "$INPUTS/snp_list.txt"
  : > "$INPUTS/copa_score_variants_qc.tsv"
  QC_HEADER_WRITTEN=0
  while IFS=$'\t' read -r ANC SNP_IN; do
    [[ -z "$ANC" || "$ANC" == \#* ]] && continue
    echo "[PREP] ancestry=$ANC input=$SNP_IN"
    ANC_OUT="$INPUTS/ancestry/$ANC"; mkdir -p "$ANC_OUT"
    Rscript "$SUBMIT_DIR/scripts/prepare_copa_score_files.R" \
      --input "$SNP_IN" --mode "$SCORE_MODE" --outdir "$ANC_OUT"
    [[ -f "$ANC_OUT/snp_list.txt" ]] || { echo "[ERROR] prep failed for ancestry $ANC"; exit 1; }
    cat "$ANC_OUT/snp_list.txt" >> "$INPUTS/snp_list.txt"
    if [[ "$QC_HEADER_WRITTEN" -eq 0 ]]; then
      cat "$ANC_OUT/copa_score_variants_qc.tsv" >> "$INPUTS/copa_score_variants_qc.tsv"
      QC_HEADER_WRITTEN=1
    else
      tail -n +2 "$ANC_OUT/copa_score_variants_qc.tsv" >> "$INPUTS/copa_score_variants_qc.tsv"
    fi
  done <<< "$ANC_MAP"

  # dedup: a SNP shared by two ancestry lists would otherwise double-count in
  # step4's coverage report (keep header + first occurrence per SNP).
  awk -F'\t' 'NR==1 || !seen[$1]++' "$INPUTS/copa_score_variants_qc.tsv" > "$INPUTS/copa_score_variants_qc.tsv.tmp"
  mv "$INPUTS/copa_score_variants_qc.tsv.tmp" "$INPUTS/copa_score_variants_qc.tsv"

  cat "$ANCESTRY_REF_CACHE/prune.prune.in" >> "$INPUTS/snp_list.txt"
  sort -u -o "$INPUTS/snp_list.txt" "$INPUTS/snp_list.txt"
  echo "[PREP] $(wc -l < "$INPUTS/snp_list.txt") SNPs in union snp_list.txt (PRS lists + ancestry markers)"
else
  Rscript "$SUBMIT_DIR/scripts/prepare_copa_score_files.R" \
    --input "$SNP_INPUT" --mode "$SCORE_MODE" --outdir "$INPUTS"
  [[ -f "$INPUTS/snp_list.txt" ]] || { echo "[ERROR] prep failed: no snp_list.txt"; exit 1; }
  echo "[PREP] $(wc -l < "$INPUTS/snp_list.txt") SNPs in snp_list.txt"
fi

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

STEP3_DEP="$STEP2"
STEP2B=""
if [[ "$RUN_ANCESTRY" == "1" ]]; then
  STEP2B=$(sub --dependency=afterok:"$STEP2" \
              --output="$LOG_DIR/ancestry_%j.out" --error="$LOG_DIR/ancestry_%j.err" \
              "$SUBMIT_DIR/slurm/step2b_ancestry.slurm")
  STEP3_DEP="$STEP2B"
fi

STEP3=$(sub --dependency=afterok:"$STEP3_DEP" \
            --output="$LOG_DIR/score_%j.out" --error="$LOG_DIR/score_%j.err" \
            "$SUBMIT_DIR/slurm/step3_score.slurm")
STEP4=$(sub --dependency=afterok:"$STEP2" \
            --output="$LOG_DIR/report_%j.out" --error="$LOG_DIR/report_%j.err" \
            "$SUBMIT_DIR/slurm/step4_report.slurm")

echo "[DONE] RUN_ID=$RUN_ID"
if [[ -n "$STEP2B" ]]; then
  echo "  step1(extract)=$STEP1  step2(merge)=$STEP2  step2b(ancestry)=$STEP2B  step3(score)=$STEP3  step4(report)=$STEP4"
else
  echo "  step1(extract)=$STEP1  step2(merge)=$STEP2  step3(score)=$STEP3  step4(report)=$STEP4"
fi
echo "  outputs: $BASE_OUT"
echo "  logs:    $LOG_DIR"
