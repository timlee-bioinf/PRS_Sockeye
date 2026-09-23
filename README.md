# GRS / PRS Calculation Pipeline

Compute **weighted and/or unweighted genetic risk scores (GRS/PRS)** for a
cohort from per-chromosome imputed dosage VCFs, using PLINK2.

The pipeline is generic: give it any candidate SNP list (rsID + risk allele,
optionally a weight) and it extracts those variants from the imputed genotypes,
merges them, and scores every sample. Weighting is optional.

## How it works

`submit.sh` runs a SLURM **DAG** (directed acyclic graph): a set of jobs chained
by `--dependency=afterok`, so each step starts only after the one(s) it depends
on finish successfully. One `bash submit.sh` submits the whole thing.

```text
        prep (login node: build snp_list + score files)
                          |
        step1  extract  (array: 22 chromosomes in parallel)
                          |  afterok (all 22)
        step2  merge    (--pmerge-list -> all_snps)
                          |  afterok
              +-----------+-----------+
              v                       v
        step3  score            step4  report     (run in parallel)
        (--score -> PRS)        (coverage + missingness)
```

1. **prep** (login node) - `scripts/prepare_copa_score_files.R` turns the SNP
   input into a plain rsID list (for extraction) and the score file(s).
2. **extract** (array, 1 task per chromosome) - pull the candidate rsIDs from
   each `chr<N>.dose.vcf.gz` with PLINK2 `--extract --make-pgen`.
3. **merge** - combine the per-chromosome filesets with `--pmerge-list`.
4. **score** - PLINK2 `--score` -> per-sample GRS (weighted and/or unweighted),
   with an optional R2 quality filter applied here.
5. **report** - per-chromosome coverage (expected vs found) and per-sample
   completeness.

`afterok` means a step is cancelled if anything it depends on fails; an array
task that finds 0 variants on a chromosome exits successfully (it is a skip,
not a failure), so it does not block the merge.

## Repository layout

```text
config.sh                          # EDIT THIS: account, data/weights paths, run dir, options
submit.sh                          # orchestrator: prep + submit the DAG
scripts/
  prepare_copa_score_files.R       # SNP input -> snp_list.txt + score file(s)
  analyze_copa_grs.R               # OPTIONAL, standalone association models (not in the DAG)
slurm/
  step1_extract.slurm              # array 1-22, extract by rsID
  step2_merge.slurm                # --pmerge-list
  step3_score.slurm                # optional R2 filter + --score
  step4_report.slurm               # coverage (expected vs found) + sample completeness
```

## Requirements

- **PLINK 2.0** - a `plink2` binary on `PATH` (or set `PLINK2` to its full path).
- **R** with `readxl`, `readr`, `dplyr` (install once into your personal library):
  ```bash
  Rscript -e 'install.packages(c("readxl","readr","dplyr"), repos="https://cloud.r-project.org")'
  ```
  `submit.sh` loads R via `module load gcc/9.4.0 r/4.4.0` - adjust for your cluster.

## Inputs

Expected layout under `DATA_ROOT` (set in `config.sh`):

```text
DATA_ROOT/
├── imputation/                      ← IMPUTE_DIR
│   └── chr_<N>/chr<N>.dose.vcf.gz   (N = 1..22)
├── weights/
│   └── <TRAIT>/                     ← WEIGHTS_DIR (one folder per trait, e.g. copd/)
│       ├── EUR.tsv                  one weights file per ancestry; the file name
│       ├── EAS.tsv                  (minus extension) is the ancestry label
│       ├── AFR.tsv
│       └── AMR.tsv
├── ancestry_ref/ref.{pgen,pvar,psam}  ← ANCESTRY_REF_PFILE (multi-ancestry only)
└── ancestry_ref_cache/                built by setup_ancestry_reference.sh
```

- **Single ancestry** (`RUN_ANCESTRY=0`): `SNP_INPUT` points at **one** weights
  file (e.g. `$WEIGHTS_DIR/EUR.tsv`), used for every sample.
- **Multi-ancestry** (`RUN_ANCESTRY=1`): every weights file in `WEIGHTS_DIR`
  is used; `SNP_INPUT` is ignored.

- **Genotypes**: per-chromosome imputed dosage VCFs at
  `$IMPUTE_DIR/chr_<N>/chr<N>.dose.vcf.gz` (e.g. Michigan/TOPMed output). The
  variant ID column must be **rsID** (matching is build-independent - no liftOver
  needed even if the data is GRCh38 and the weights are b37).
- **Weights files** (`SNP_INPUT`, or each file in `WEIGHTS_DIR`): an `.xlsx`
  (e.g. the GWAS weight workbook) or a `.csv`/`.tsv` candidate list. Columns
  are auto-detected:
  - rsID - required (`Marker Name`, `rsID`, `SNP`, `ID`, ...)
  - risk/effect allele - required (`Risk allele`, `effect_allele`, `A1`, ...)
  - weight - optional (`Weight`, `Beta`); only needed for weighted scoring.

  For the supplementary workbook, the script reads sheet
  `Supp. Table 28 RiskScoreWeights` (skipping its 2 header rows).

## Running

Edit `config.sh` with your values, then run:

```bash
# config.sh sets:
#   SBATCH_ACCOUNT, MAIL_USER, DATA_ROOT, IMPUTE_DIR, TRAIT, SNP_INPUT, PLINK2,
#   SCORE_MODE (both|weighted|unweighted), R2_THRESH, RUN_BASE, RUN_ANCESTRY
bash submit.sh
```

`config.sh` exports `SBATCH_ACCOUNT` (when set) so Slurm charges that
allocation. Leave it empty to use your cluster default.

### Options

- `SCORE_MODE` - `both` (default), `weighted`, or `unweighted`. A GRS does not
  require weights; `unweighted` just needs rsID + risk allele.
- `R2_THRESH` - imputation-quality filter applied **at scoring** (`0` = keep all;
  e.g. `0.8`). Extraction keeps every matched variant, so you can re-score at a
  different threshold, or compare thresholds, without re-extracting.

## Outputs

Under `$RUN_BASE/SNP_extract_<RUN_ID>/` (`RUN_BASE` defaults to
`run_output/<TRAIT>/`, so each trait's runs are kept apart):

```text
inputs/snp_list.txt                # rsIDs used for extraction
inputs/copa_*_score.tsv            # score file(s)
task1_extract/chr<N>.{pgen,pvar,psam}
task2_merge/all_snps.{pgen,pvar,psam}
task3_score/grs_unweighted.sscore  # per-sample PRS (if mode includes it)
task3_score/grs_weighted.sscore
task4_report/coverage_by_chr.tsv   # expected vs found, per chromosome
task4_report/missing_snps.txt      # candidate SNPs not present in the data
task4_report/all_snps.{smiss,vmiss}# per-sample / per-variant missingness
logs/                              # per-step logs
```

**The PRS** is in each `.sscore`: `#IID` is the sample, `SCORE1_SUM` is its
score (weighted sum, or risk-allele count in the unweighted file):

```text
#IID                NMISS_ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  SCORE1_AVG  SCORE1_SUM
sample_001          556              278.3                    0.00123     0.6841
sample_002          556              279.0                    0.00115     0.6390
```

**The coverage table** (`coverage_by_chr.tsv`) shows how many candidate SNPs
were actually found per chromosome:

```text
CHR  expected  found  missing
1    33        32     1
2    28        28     0
...
TOTAL 279      27x    x
```

`missing_snps.txt` lists the specific rsIDs not found, and the step4 log prints
a sample-completeness line (e.g. `samples=31 with_any_missing_genotype=0
max_F_MISS=0`) confirming whether every sample has these SNPs.

## Multi-ancestry (optional)

Adapts the ancestry-inference idea from
[pgscatalog/pgsc_calc](https://github.com/pgscatalog/pgsc_calc) into this
PLINK2/bash/SLURM pipeline, without adopting Nextflow, containers, or its
Python package. Key difference from pgsc_calc: pgsc_calc scores every sample
with **one** PGS weight file and only adjusts the result post-hoc; this mode
instead **routes each sample to its own ancestry-matched SNP/weight list**
(what you already have), and additionally reports pgsc_calc-style
normalization columns.

Off by default (`RUN_ANCESTRY=0`) - every file/step above is unchanged unless
you opt in. Turning it on replaces the single `SNP_INPUT` with every
per-ancestry weights file in `WEIGHTS_DIR` and adds one pipeline step:

```text
        prep (per ancestry: build snp_list + score files; union -> step1)
                          |
        step1  extract  (unchanged - snp_list.txt is now bigger)
                          |
        step2  merge    (unchanged)
                          |
              +-----------+-----------+
              v                       v
        step2b ancestry          step4  report   (unchanged)
        (PCA projection +
         classification)
              |
        step3  score
        (per-ancestry --keep + --score,
         then empirical + PCA-regression
         normalization)
```

### Setup (not part of the per-run DAG)

1. A reference ancestry panel in **PLINK2 pgen/pvar/psam** format, with rsIDs
   in its `.pvar` (matched by rsID, same as everywhere else in this pipeline -
   no genome-build/liftover requirement) and a population-label column in its
   `.psam` (e.g. `SuperPop` with values like `EUR`/`EAS`/`SAS`/`AFR`/`AMR`).
   pgsc_calc's own prebuilt 1000G / HGDP+1kGP reference panels are in exactly
   this format.
2. One weights file per ancestry in `weights/<TRAIT>/`, named by ancestry
   label (`EUR.tsv`, `EAS.tsv`, `AFR.tsv`, `AMR.tsv`, ...). **The file name
   (minus extension) must exactly match** a value in the reference panel's
   label column - e.g. with 1000 Genomes super-populations, Latino/admixed
   American is `AMR`.
3. Set `RUN_ANCESTRY=1`, `TRAIT`, and the `ANCESTRY_*` variables in
   `config.sh` (reference panel path, label column, PC counts - see
   comments in `config.sh`).
4. Run `bash scripts/setup_ancestry_reference.sh` once per trait. The first
   time, it LD-prunes the reference panel, computes its PCA and fits a
   Mahalanobis population classifier (shared by all traits, cached in
   `ANCESTRY_REF_CACHE`). Then, per ancestry, it scores the reference panel
   with that ancestry's weights file to build empirical + PC-regression
   normalization models (cached per trait in
   `ANCESTRY_REF_CACHE/traits/<TRAIT>/`). For a new trait it only redoes the
   per-trait part; re-run it when a trait's weights files or `SCORE_MODE`
   change, and use `REBUILD_REF=1` if the reference panel,
   `ANCESTRY_REF_EXCLUDE`, or the PC-count settings change.

### New per-run step: `step2b_ancestry.slurm`

Projects this run's samples onto the cached reference PCA (`plink2 --score`
against the PCA allele-weights, a bare-PLINK2 analog of pgsc_calc's
FRAPOSA/OADP projection), then classifies each sample's most-similar
reference population via Mahalanobis distance (mirrors pgsc_calc's
`--ancestry_method Mahalanobis`) and writes one `--keep` list per ancestry.
Never rejects a sample - low-confidence calls are flagged (`LowConfidence`
column), not dropped; decide filtering downstream if you want it.

### `step3_score.slurm` output: `task3_score/grs_combined.tsv`

Per sample, per `SCORE_MODE`:
- `ancestry_assigned`, `LowConfidence`, `PC1..PCk` - the classification result
  ("continuous ancestry").
- `Score_<mode>` - raw score from the sample's own ancestry-matched weight file.
- `Z_MostSimilarPop_<mode>` / `percentile_MostSimilarPop_<mode>` - empirical
  Z-score/percentile vs. the reference panel's same-population score
  distribution ("discrete ancestry" normalization, pgsc_calc-style).
- `Z_norm1_<mode>` - continuous PCA-regression residual Z-score (score
  regressed on top PCs within the matched population; Khera et al.
  2019-style, like pgsc_calc's `Z_norm1`). Fit **per ancestry bucket** here
  (not across the whole reference panel like pgsc_calc) since different
  ancestries are scored on different weight-file scales.
- pgsc_calc's variance-adjusted `Z_norm2` is **not** implemented (scope-cut
  for a bare-minimum port).

### Caveats - please verify on first run

- **Untested end-to-end**: written without cluster/PLINK2/live-data access,
  then checked with an independent adversarial code review (bash correctness,
  R/statistics correctness, cross-file interface consistency) that found and
  fixed 10 concrete bugs before this was ever run for real - see git history
  for what changed. Still worth a careful look at `logs/ancestry.log` and
  `grs_combined.tsv` on your first real run:
  - The exact column layout `plink2 --pca allele-wts` writes to
    `.eigenvec.allele` was not verified against a live binary -
    `step2b_ancestry.slurm` and `classify_ancestry.R` detect columns from the
    file's own header rather than hardcoding positions, and fail loudly
    (dumping the header) rather than silently miscomputing if detection
    fails.
  - Ancestries are discovered from the `.xlsx`/`.csv`/`.tsv` files in
    `WEIGHTS_DIR` (other files are ignored). Run
    `DEBUG_CONFIG=1 bash -c 'source config.sh'` to see which labels were
    found, and check each shows up in `$INPUTS/ancestry/<label>/` after prep
    and in `task3_score/<label>/` after scoring.
- `ANCESTRY_REF_EXCLUDE` (related/duplicate reference samples to drop, e.g. a
  `*.king.cutoff.out.id` file) is optional but recommended if your reference
  panel includes related individuals - they bias the PCA/classifier otherwise.
  It's applied via `plink2 --remove` (the primary mechanism) to every
  reference-panel plink2 call; the R scripts also filter by it defensively,
  but since they read plink2's already-filtered output, that R-side filter
  is normally a no-op.
- A sample with no reference population resembling it at all is still
  assigned to its nearest match (flagged `LowConfidence`), matching pgsc_calc's
  own default behavior - there's no reject/"unassigned" category.
- `setup_ancestry_reference.sh` builds normalization models only for the
  `SCORE_MODE` in effect when you run it. `submit.sh`'s preflight now checks
  that every ancestry/mode this run needs actually exists in the cache for the
  current `TRAIT` before submitting the DAG - re-run
  `setup_ancestry_reference.sh` if you add a trait or weights file, or change
  `SCORE_MODE`.
- A monomorphic/degenerate weight file for one ancestry (zero variance in the
  reference panel's scores for that ancestry) makes `setup_ancestry_reference.sh`
  fail loudly for that ancestry rather than silently emitting `Inf`/`NaN`
  normalized scores.

## Notes

- **Coverage**: expect fewer variants than the input list - some SNPs may be
  absent or (if `R2_THRESH>0`) filtered out. The merge and score logs report the
  counts. If coverage is very low, the rsIDs may have drifted across dbSNP builds
  (fallback: match by `chr:pos`, which would require build-aligned positions).
- **Small cohorts**: `--score` uses `no-mean-imputation` (imputed dosages have
  ~no missingness, and PLINK refuses to auto-impute allele frequencies for
  fewer than 50 samples).
- **Strand-ambiguous SNPs** (A/T, C/G) carry orientation risk under `--score`;
  worth checking if your list contains many.
- `scripts/analyze_copa_grs.R` is an optional, standalone association step
  (GRS-only models for a small cohort) - not part of the DAG.
