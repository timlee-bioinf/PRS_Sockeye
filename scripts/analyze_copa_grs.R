suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

pheno_path <- get_arg("--pheno", "COPA_Pheno.xlsx")
weighted_path <- get_arg("--weighted", "results/copa_weighted_grs.sscore")
unweighted_path <- get_arg("--unweighted", "results/copa_unweighted_grs.sscore")
outdir <- get_arg("--outdir", "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clean_plink_names <- function(x) {
  sub("^#", "", x)
}

read_sscore <- function(path, score_name) {
  x <- read_table(path, show_col_types = FALSE)
  names(x) <- clean_plink_names(names(x))
  if (!"IID" %in% names(x)) stop("No IID column found in ", path)
  score_col <- if ("SCORE1_SUM" %in% names(x)) "SCORE1_SUM" else "SCORE1_AVG"
  x %>%
    transmute(
      short_id = str_extract(IID, "COPA-[0-9]+"),
      "{score_name}" := as.numeric(.data[[score_col]]),
      "{paste0(score_name, '_z')}" := as.numeric(scale(as.numeric(.data[[score_col]])))
    )
}

pheno <- read_excel(pheno_path) %>%
  mutate(
    short_id = str_extract(as.character(ID), "COPA-?[0-9]+"),
    FEV1_FVC = as.numeric(`Post FEV1/FVC`),
    age = as.numeric(Age),
    sex = factor(Sex),
    COPD = as.integer(COPD)
  )

weighted <- read_sscore(weighted_path, "GRS_weighted")
unweighted <- read_sscore(unweighted_path, "GRS_unweighted")

dat <- pheno %>%
  left_join(weighted, by = "short_id") %>%
  left_join(unweighted, by = "short_id")

if (sum(!is.na(dat$GRS_weighted)) == 0) {
  stop("No samples matched between phenotype ID and PLINK IIDs - check ID formats")
}

write_csv(dat, file.path(outdir, "copa_pheno_grs_merged.csv"))

# Small cohort (n=31): GRS-only models, no covariates, no PCs (avoids overfitting).
fit_lm <- function(score) {
  lm(as.formula(paste("FEV1_FVC ~", score)), data = dat)
}

fit_glm <- function(score) {
  glm(as.formula(paste("COPD ~", score)), data = dat, family = binomial())
}

models <- list(
  fev1_fvc_weighted = fit_lm("GRS_weighted_z"),
  fev1_fvc_unweighted = fit_lm("GRS_unweighted_z"),
  copd_weighted = fit_glm("GRS_weighted_z"),
  copd_unweighted = fit_glm("GRS_unweighted_z")
)

model_tables <- bind_rows(lapply(names(models), function(nm) {
  tab <- as.data.frame(coef(summary(models[[nm]])))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab$model <- nm
  tab %>% select(model, term, everything())
}))

write_csv(model_tables, file.path(outdir, "copa_grs_model_terms.csv"))

sink(file.path(outdir, "copa_grs_model_summary.txt"))
cat("Merged samples with weighted GRS:", sum(!is.na(dat$GRS_weighted)), "\n")
cat("Merged samples with unweighted GRS:", sum(!is.na(dat$GRS_unweighted)), "\n")
cat("Models: GRS-only (no covariates, no PCs) - small-cohort setting\n\n")
for (nm in names(models)) {
  cat("\n====================\n", nm, "\n====================\n", sep = "")
  print(summary(models[[nm]]))
}
sink()
