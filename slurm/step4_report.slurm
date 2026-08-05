#!/bin/bash -l
#SBATCH --job-name=grs_report
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/report_%j.out
#SBATCH --error=logs/report_%j.err

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"
: "${RUN_ID:?RUN_ID not set}"

BASE_OUT="$RUN_BASE/SNP_extract_${RUN_ID}"
INPUTS="$BASE_OUT/inputs"
TASK2_OUT="$BASE_OUT/task2_merge"
TASK4_OUT="$BASE_OUT/task4_report"
LOG_DIR="$BASE_OUT/logs"
mkdir -p "$TASK4_OUT" "$LOG_DIR"
STEP_LOG="$LOG_DIR/report.log"; : > "$STEP_LOG"

QC="$INPUTS/copa_score_variants_qc.tsv"          # expected SNPs: SNP CHR POS A1 [WEIGHT]
PVAR="$TASK2_OUT/all_snps.pvar"                   # variants actually present
PFILE="$TASK2_OUT/all_snps"
[[ -f "$QC"   ]] || { echo "[ERROR] missing $QC";   exit 1; }
[[ -f "$PVAR" ]] || { echo "[ERROR] missing $PVAR"; exit 1; }

MISS="$TASK4_OUT/missing_snps.txt"
BYCHR="$TASK4_OUT/coverage_by_chr.tsv"
: > "$MISS"

# ---- expected (from QC, by chr) vs found (in merged pvar), per chromosome ----
awk -v MISS="$MISS" '
  FNR==NR { if ($1 !~ /^#/) found[$3]=1; next }   # pvar: collect found IDs (col 3); skip ## and #CHROM
  FNR==1  { next }                                 # qc header
  {
    snp=$1; chr=$2;
    ex[chr]++;
    if (snp in found) fd[chr]++;
    else print snp "\t" chr >> MISS
  }
  END { for (c in ex) print c "\t" ex[c] "\t" (fd[c]+0) "\t" (ex[c]-fd[c]) }
' "$PVAR" "$QC" | sort -k1,1n > "$TASK4_OUT/_chr_rows.tsv"

{
  printf "CHR\texpected\tfound\tmissing\n"
  cat "$TASK4_OUT/_chr_rows.tsv"
  awk -F'\t' '{e+=$2; f+=$3; m+=$4} END{printf "TOTAL\t%d\t%d\t%d\n", e, f, m}' "$TASK4_OUT/_chr_rows.tsv"
} > "$BYCHR"
rm -f "$TASK4_OUT/_chr_rows.tsv"

echo "[INFO] coverage by chromosome:" | tee -a "$STEP_LOG"
cat "$BYCHR" | tee -a "$STEP_LOG"
echo "[INFO] $(wc -l < "$MISS") SNP(s) not found in the imputed data -> $MISS" | tee -a "$STEP_LOG"

# ---- per-sample completeness: do all samples have these SNPs? ----
"$PLINK2" --pfile "$PFILE" --missing --out "$TASK4_OUT/all_snps" >> "$STEP_LOG" 2>&1 || true
SMISS="$TASK4_OUT/all_snps.smiss"
if [[ -f "$SMISS" ]]; then
  awk '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="F_MISS") fc=i; next }
    { n++; if ($fc+0 > 0) miss++; if ($fc+0 > max) max=$fc+0 }
    END { printf "[INFO] samples=%d  with_any_missing_genotype=%d  max_F_MISS=%g\n", n, miss+0, max+0 }
  ' "$SMISS" | tee -a "$STEP_LOG"
else
  echo "[WARN] no .smiss produced; skipped sample-completeness summary" | tee -a "$STEP_LOG"
fi

echo "[DONE] coverage report in $TASK4_OUT" | tee -a "$STEP_LOG"
