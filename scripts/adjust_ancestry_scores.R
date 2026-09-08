#!/usr/bin/env Rscript
# Per-run (end of step3_score.slurm): combine each ancestry's ancestry-routed
# raw score with pgsc_calc-style ancestry normalization:
#   SCORE_raw               - raw score, computed with the sample's own
#                              ancestry-matched weight file
#   Z_MostSimilarPop /
#   percentile_MostSimilarPop - empirical Z-score / percentile vs. the
#                              reference panel's same-population, same-weight-
#                              file score distribution
#   Z_norm1                  - continuous PCA-regression residual Z-score
#                              (score regressed on top PCs within the same
#                              reference population; Khera et al. 2019-style),
#                              alongside the sample's own PCs ("continuous
#                              ancestry"). This mirrors pgsc_calc's Z_norm1 but
#                              is fit per ancestry bucket rather than across
#                              the whole reference panel, since (unlike
#                              pgsc_calc) different ancestries here are scored
#                              with different weight files on different scales.
#   pgsc_calc's variance-adjusted Z_norm2 is NOT implemented (out of scope for
#   a bare-minimum port); Z_norm1 covers mean-normalization across ancestry.
#
# Usage:
#   Rscript adjust_ancestry_scores.R --calls ancestry_calls.tsv
#       --ancestries EUR,EAS,SAS,AFR --score-mode both
#       --score-dir task3_score/ --norm-models-dir ancestry_ref_cache/ancestry
#       --out task3_score/grs_combined.tsv

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

calls_path      <- get_arg("--calls");           if (is.null(calls_path))      stop("Provide --calls")
ancestries      <- strsplit(get_arg("--ancestries"), ",")[[1]]
score_mode      <- get_arg("--score-mode", "both")
score_dir       <- get_arg("--score-dir");       if (is.null(score_dir))       stop("Provide --score-dir")
norm_models_dir <- get_arg("--norm-models-dir"); if (is.null(norm_models_dir)) stop("Provide --norm-models-dir")
out_path        <- get_arg("--out");             if (is.null(out_path))        stop("Provide --out")

modes <- switch(score_mode,
  both       = c("weighted", "unweighted"),
  weighted   = "weighted",
  unweighted = "unweighted",
  stop("--score-mode must be one of: both, weighted, unweighted"))

calls <- read.delim(calls_path, check.names = FALSE, stringsAsFactors = FALSE)
pc_cols_all <- grep("^PC[0-9]+$", names(calls), value = TRUE)
base <- calls %>% select(IID, ancestry_assigned = MostSimilarPop, LowConfidence, all_of(pc_cols_all))

read_sscore <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(x) <- sub("^#", "", names(x))
  score_col <- if ("SCORE1_SUM" %in% names(x)) "SCORE1_SUM" else if ("SCORE1_AVG" %in% names(x)) "SCORE1_AVG" else {
    stop("No SCORE1_SUM/SCORE1_AVG column in ", path, ". Columns present: ", paste(names(x), collapse = ", "))
  }
  x %>% transmute(IID, Score = as.numeric(.data[[score_col]]))
}

for (mode in modes) {
  rows <- list()
  for (a in ancestries) {
    sscore_path <- file.path(score_dir, a, paste0("grs_", mode, ".sscore"))
    if (!file.exists(sscore_path)) {
      message(sprintf("[%s/%s] no sscore file (0 samples assigned or scoring skipped) - skip", a, mode))
      next
    }
    sc <- read_sscore(sscore_path)
    if (nrow(sc) == 0) next

    model_path <- file.path(norm_models_dir, a, paste0("norm_models_", mode, ".rds"))
    if (!file.exists(model_path)) {
      stop("Missing normalization model: ", model_path,
           " - did setup_ancestry_reference.sh build norm models for ancestry '", a, "' mode '", mode, "'?")
    }
    m <- readRDS(model_path)

    sc <- sc %>% mutate(
      Z_MostSimilarPop = (Score - m$mean) / m$sd,
      percentile_MostSimilarPop = 100 * findInterval(Score, m$sorted_scores) / length(m$sorted_scores)
    )

    pcs_for_fit <- base %>% filter(IID %in% sc$IID) %>% select(IID, all_of(m$pc_cols))
    sc <- sc %>% inner_join(pcs_for_fit, by = "IID")
    predicted <- predict(m$fit, newdata = sc)
    sc$Z_norm1 <- (sc$Score - predicted) / m$resid_sd
    sc <- sc %>% select(IID, Score, Z_MostSimilarPop, percentile_MostSimilarPop, Z_norm1)
    names(sc)[-1] <- paste0(names(sc)[-1], "_", mode)
    rows[[a]] <- sc
  }
  if (length(rows) == 0) {
    message(sprintf("[%s] no samples scored in any ancestry - skipping this mode's columns", mode))
    next
  }
  combined_mode <- bind_rows(rows)
  base <- base %>% left_join(combined_mode, by = "IID")
}

write.table(base, out_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
message(sprintf("Wrote combined ancestry-adjusted scores for %d samples -> %s", nrow(base), out_path))
