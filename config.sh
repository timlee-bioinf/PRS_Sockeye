#!/usr/bin/env bash
# Shared configuration for the GRS pipeline.
# This is the ONE file to edit: set your data, weights, run directory, and
# SLURM account below, then run `bash submit.sh`. Every script (submit.sh,
# the slurm steps, setup_ancestry_reference.sh) sources this file.
set -euo pipefail

CONFIG_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$CONFIG_DIR"

# --- SLURM ---
# Allocation to charge jobs to. Leave empty to use your cluster default.
SBATCH_ACCOUNT=""
# Optional email notifications. Leave MAIL_USER empty to disable.
MAIL_USER=""
MAIL_TYPE="ALL"
if [[ -n "$SBATCH_ACCOUNT" ]]; then export SBATCH_ACCOUNT; fi

# --- inputs ---
# Expected layout under DATA_ROOT:
#   imputation/chr_<N>/chr<N>.dose.vcf.gz   genotypes
#   weights/<TRAIT>/<ANC>.{tsv,csv,xlsx}     one weights file per ancestry
#   ancestry_ref/ref.{pgen,pvar,psam}       reference panel (multi-ancestry only)
#   ancestry_ref_cache/                     built by setup_ancestry_reference.sh
DATA_ROOT="/path/to/your/grs_data"
# Per-chromosome imputed dose VCFs: $IMPUTE_DIR/chr_<N>/chr<N>.dose.vcf.gz
IMPUTE_DIR="$DATA_ROOT/imputation"
# Trait being scored = name of its folder under weights/.
TRAIT="copd"
WEIGHTS_DIR="$DATA_ROOT/weights/$TRAIT"
# Weights files: .xlsx / .csv / .tsv with at least an rsID column and a
# risk/effect allele column; a weight column is only needed for weighted scoring.
#   RUN_ANCESTRY=0: SNP_INPUT is the ONE file used for every sample (assumes
#                   the whole cohort is that ancestry).
#   RUN_ANCESTRY=1: SNP_INPUT is ignored; every weights file in WEIGHTS_DIR is
#                   used, and its file name (minus extension) is its ancestry
#                   label, e.g. EUR.tsv -> EUR. Labels MUST exactly match the
#                   reference panel's ANCESTRY_REF_LABEL_COL values.
SNP_INPUT="$WEIGHTS_DIR/EUR.tsv"

# --- tools ---
PLINK2="plink2"   # binary name on PATH, or an absolute path

# --- options ---
SCORE_MODE="both"   # both | weighted | unweighted
R2_THRESH="0"       # 0 = keep all; e.g. 0.8 to filter at scoring time

# --- output ---
RUN_BASE="$PROJECT_ROOT/run_output/$TRAIT"
# Auto-prune SNP_extract_* run directories older than this many days on each
# new submit.sh invocation. 0 disables pruning.
RETENTION_DAYS="30"

# --- multi-ancestry (optional; default off, classic single-list behavior unchanged) ---
# RUN_ANCESTRY=1 turns on: (a) per-ancestry weights files from WEIGHTS_DIR instead
# of one SNP_INPUT, (b) ancestry inference (PCA projection + population
# classification) in step2b, (c) ancestry-routed scoring + discrete/continuous
# ancestry-normalized scores in step3.
# See scripts/setup_ancestry_reference.sh for the required cache build (once per trait).
RUN_ANCESTRY="0"
# Reference ancestry panel (e.g. your "g2k" data): PLINK2 pgen/pvar/psam, path
# WITHOUT extension. Must carry rsIDs in its .pvar (matched by rsID, same as
# the rest of this pipeline - no build/liftover requirement as long as both
# sides use rsIDs).
ANCESTRY_REF_PFILE="$DATA_ROOT/ancestry_ref/ref"
# Column in the reference .psam holding each sample's population label.
ANCESTRY_REF_LABEL_COL="SuperPop"
# Optional: sample-ID file (one per line) of related/duplicate reference
# samples to exclude before computing PCA + reference score distributions
# (e.g. a *.king.cutoff.out.id file). Leave empty to use all reference samples.
ANCESTRY_REF_EXCLUDE=""
# Setup output cache (see scripts/setup_ancestry_reference.sh). Not tied to
# RUN_ID - reused by every pipeline run. The PCA/classifier part is shared by
# all traits; the per-ancestry normalization models live under traits/<TRAIT>/.
ANCESTRY_REF_CACHE="$DATA_ROOT/ancestry_ref_cache"
ANCESTRY_NORM_DIR="$ANCESTRY_REF_CACHE/traits/$TRAIT"
ANCESTRY_N_PCS="10"            # PCs computed for the reference PCA basis
ANCESTRY_N_POPCOMP="5"         # top PCs used for population classification (Mahalanobis)
ANCESTRY_N_NORM="4"            # top PCs used for continuous PC-regression normalization (Z_norm1)
# Mahalanobis chi-sq p-value floor below which a call is flagged LowConfidence
# (samples are still assigned + scored - this never rejects, only flags).
ANCESTRY_PVAL_THRESH="1e-10"

# Prints "<ancestry_label><TAB><weights file>" for every weights file in
# WEIGHTS_DIR (label = file name minus extension). Fails if the folder is
# missing/empty or two files share a label (e.g. EUR.tsv + EUR.xlsx).
# Use as: ANC_MAP="$(ancestry_map)" || exit 1
ancestry_map() {
  [[ -d "$WEIGHTS_DIR" ]] || { echo "[ERROR] WEIGHTS_DIR not found: $WEIGHTS_DIR" >&2; return 1; }
  local f base label out="" seen=" "
  for f in "$WEIGHTS_DIR"/*.xlsx "$WEIGHTS_DIR"/*.csv "$WEIGHTS_DIR"/*.tsv; do
    [[ -f "$f" ]] || continue
    base="$(basename "$f")"
    [[ "$base" == [~.]* ]] && continue   # skip Excel lock files / hidden files
    label="${base%.*}"
    [[ "$seen" == *" $label "* ]] && { echo "[ERROR] two weights files for ancestry '$label' in $WEIGHTS_DIR" >&2; return 1; }
    seen+="$label "
    out+="${label}"$'\t'"${f}"$'\n'
  done
  [[ -n "$out" ]] || { echo "[ERROR] no .xlsx/.csv/.tsv weights files in $WEIGHTS_DIR" >&2; return 1; }
  printf '%s' "$out"
}

if [[ "${DEBUG_CONFIG:-0}" == "1" ]]; then
  for v in PROJECT_ROOT SBATCH_ACCOUNT MAIL_USER DATA_ROOT IMPUTE_DIR TRAIT WEIGHTS_DIR SNP_INPUT PLINK2 SCORE_MODE R2_THRESH \
           RUN_BASE RETENTION_DAYS RUN_ANCESTRY ANCESTRY_REF_PFILE ANCESTRY_REF_LABEL_COL ANCESTRY_REF_EXCLUDE \
           ANCESTRY_REF_CACHE ANCESTRY_NORM_DIR ANCESTRY_N_PCS ANCESTRY_N_POPCOMP ANCESTRY_N_NORM ANCESTRY_PVAL_THRESH; do
    echo "[CONFIG] $v=${!v}"
  done
  if [[ "$RUN_ANCESTRY" == "1" ]]; then
    ancestry_map | sed 's/^/[CONFIG] ancestry: /' || true
  fi
fi
