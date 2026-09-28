#!/bin/bash -l
#SBATCH --job-name=grs_ref_rsids
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# One-time conversion (NOT part of the per-run DAG): give the reference panel
# rsID variant IDs, so it matches the rsID-based weights and imputed data.
#
# pgsc_calc's reference panels (e.g. GRCh38_HGDP+1kGP_ALL) use positional IDs
# (1:10410:CCCTAA:C) and keep the rsID in INFO as gnomad_ID=rs.... This writes
# a new .pvar.zst whose ID is that rsID where present, else the original
# positional ID. The .pgen/.psam don't store variant IDs, so they are reused
# as-is (the copy step links them under the new prefix - no 13 GB copy).
#
# One rsID per variant: if several records share an rsID (split multi-allelic
# sites), it goes to the first SNV record (else the first record); the others
# keep their positional IDs, so every ID stays unique.
#
# Input : $ANCESTRY_REF_PFILE.{pvar.zst|pvar}   (read-only)
# Output: $REF_RSID_WORK/ref_rsid.pvar.zst + ref_rsid.log   (on /scratch)
# `./submit.local.sh refids` waits for this job, then copies the .pvar.zst next
# to the original as <prefix>_rsid.pvar.zst and links the .pgen/.psam.

set -euo pipefail

CONFIG="${CONFIG:?CONFIG not exported}"
source "$CONFIG"

REF="$ANCESTRY_REF_PFILE"
WORK="$REF_RSID_WORK"
mkdir -p "$WORK"
rm -f "$WORK/ref_rsid.pvar.zst"
LOG="$WORK/ref_rsid.log"; : > "$LOG"
TMP="$WORK/tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
log() { echo "$*" | tee -a "$LOG"; }

command -v "$PLINK2" >/dev/null 2>&1 || [[ -x "$PLINK2" ]] || { log "[ERROR] PLINK2 not found: $PLINK2"; exit 1; }
if   [[ -f "${REF}.pvar.zst" ]]; then PVAR="${REF}.pvar.zst"
elif [[ -f "${REF}.pvar"     ]]; then PVAR="${REF}.pvar"
else log "[ERROR] reference .pvar(.zst) not found for prefix: $REF"; exit 1; fi

log "=== [1/3] reading rsIDs from INFO/gnomad_ID in $PVAR $(date -Is) ==="
# Stream the .pvar (plink2 decompresses .zst itself - no zstd module needed).
# For each record with gnomad_ID=rs..., emit: rsID, not-SNV flag, record #, old ID.
# The INFO column position is read from the #CHROM header line.
read_pvar() { if [[ "$PVAR" == *.zst ]]; then "$PLINK2" --zst-decompress "$PVAR"; else cat "$PVAR"; fi; }
read_pvar | awk -F'\t' -v OFS='\t' -v nfile="$TMP/n_total" '
  /^##/ { next }
  /^#CHROM/ { for (i = 1; i <= NF; i++) { if ($i == "ID") id = i; if ($i == "REF") r = i; if ($i == "ALT") a = i; if ($i == "INFO") inf = i }
              if (!id || !r || !a || !inf) { print "[ERROR] missing ID/REF/ALT/INFO column in header" > "/dev/stderr"; exit 1 }
              next }
  { n++
    if (match($inf, /gnomad_ID=rs[0-9]+/)) {
      rs = substr($inf, RSTART + 10, RLENGTH - 10)
      snv = (length($r) == 1 && length($a) == 1) ? 0 : 1
      print rs, snv, n, $id
    } }
  END { print n + 0 > nfile }' > "$TMP/cand.tsv"
N_TOTAL=$(cat "$TMP/n_total")
N_CAND=$(wc -l < "$TMP/cand.tsv")
log "  $N_TOTAL variants, $N_CAND with a gnomad_ID rsID"

log "=== [2/3] picking one record per rsID $(date -Is) ==="
# Sort by rsID, then SNVs first, then file order; keep the first per rsID.
# Output: old ID <TAB> new ID (for plink2 --update-name).
LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -S 4G -T "$TMP" "$TMP/cand.tsv" \
  | awk -F'\t' -v OFS='\t' '$1 != prev { print $4, $1; prev = $1 }' > "$TMP/rename.tsv"
N_RENAME=$(wc -l < "$TMP/rename.tsv")
log "  $N_RENAME variants get an rsID ($((N_CAND - N_RENAME)) duplicate-rsID records keep positional IDs)"

log "=== [3/3] writing the rsID .pvar.zst $(date -Is) ==="
"$PLINK2" --pfile "$REF" \
  --update-name "$TMP/rename.tsv" 2 1 \
  --make-just-pvar zs \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --memory 28000 \
  --out "$WORK/ref_rsid" >> "$LOG" 2>&1

[[ -f "$WORK/ref_rsid.pvar.zst" ]] || { log "[ERROR] plink2 did not write $WORK/ref_rsid.pvar.zst - see $LOG"; exit 1; }
log "=== DONE $(date -Is) - $WORK/ref_rsid.pvar.zst ($N_RENAME / $N_TOTAL variants with rsIDs) ==="
