#!/usr/bin/env Rscript
# One-time setup helper (called by setup_ancestry_reference.sh), run once per
# ancestry: fits the empirical (mean/sd/percentile) and continuous
# (PC-regression) score-normalization models against the reference panel's own
# samples for that population, using that ancestry's own weight file.
#
# Usage:
#   Rscript build_ancestry_norm_models.R --sscore ref_score.sscore
#       --eigenvec ref_pca.eigenvec --psam ref.psam --label-col SuperPop
#       --ancestry EUR --npcs 4 --out norm_models.rds [--exclude related_ids.txt]

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

# Sample-ID file reader tolerant of both a bare one-ID-per-line format and
# plink2's native "#FID<TAB>IID" 2-column format (e.g. a *.king.cutoff.out.id
# file) - takes the LAST whitespace-delimited field of each non-comment line,
# which is IID either way.
read_id_list <- function(path) {
  lines <- trimws(readLines(path))
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  vapply(strsplit(lines, "[ \t]+"), function(x) x[length(x)], character(1))
}

sscore_path  <- get_arg("--sscore");   if (is.null(sscore_path))  stop("Provide --sscore")
eigenvec_path<- get_arg("--eigenvec"); if (is.null(eigenvec_path))stop("Provide --eigenvec")
psam_path    <- get_arg("--psam");     if (is.null(psam_path))    stop("Provide --psam")
label_col    <- get_arg("--label-col", "SuperPop")
ancestry     <- get_arg("--ancestry"); if (is.null(ancestry))     stop("Provide --ancestry")
npcs         <- as.integer(get_arg("--npcs", "4"))
out_path     <- get_arg("--out");      if (is.null(out_path))     stop("Provide --out")
exclude_path <- get_arg("--exclude", NA)

read_plink_table <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(x) <- sub("^#", "", names(x))
  x
}

sc <- read_plink_table(sscore_path)
if (!"IID" %in% names(sc)) stop("No IID column in ", sscore_path)
score_col <- if ("SCORE1_SUM" %in% names(sc)) "SCORE1_SUM" else if ("SCORE1_AVG" %in% names(sc)) "SCORE1_AVG" else {
  stop("No SCORE1_SUM/SCORE1_AVG column in ", sscore_path, ". Columns present: ", paste(names(sc), collapse = ", "))
}
sc <- sc %>% transmute(IID, Score = as.numeric(.data[[score_col]]))

eig <- read_plink_table(eigenvec_path)
if (!"IID" %in% names(eig)) stop("No IID column in ", eigenvec_path)
pc_cols <- grep("^PC[0-9]+$", names(eig), value = TRUE)
if (length(pc_cols) < npcs) {
  stop("Requested --npcs ", npcs, " but ", eigenvec_path, " only has ", length(pc_cols), " PC columns")
}
pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))][seq_len(npcs)]
eig <- eig %>% select(IID, all_of(pc_cols))

psam <- read_plink_table(psam_path)
if (!"IID" %in% names(psam)) stop("No IID column in ", psam_path)
if (!label_col %in% names(psam)) stop("Label column '", label_col, "' not found in ", psam_path)
psam <- psam %>% transmute(IID, Pop = .data[[label_col]])

dat <- sc %>% inner_join(eig, by = "IID") %>% inner_join(psam, by = "IID")

if (!is.na(exclude_path)) {
  excl <- read_id_list(exclude_path)
  dat <- dat %>% filter(!IID %in% excl)
}

dat <- dat %>% filter(Pop == ancestry)
n_before_na <- nrow(dat)
dat <- dat %>% filter(!is.na(Score))
if (n_before_na - nrow(dat) > 0) {
  message(sprintf("[%s] dropped %d reference sample(s) with NA score (e.g. 0 non-missing genotypes in this ancestry's weight file)",
                   ancestry, n_before_na - nrow(dat)))
}
if (nrow(dat) <= npcs + 1) {
  stop("Only ", nrow(dat), " reference samples labeled '", ancestry, "' with a non-NA score after joining score+PCs+labels - ",
       "need more than --npcs (", npcs, ") + 1 to fit a regression. Check --label-col / --ancestry spelling.")
}

score_sd <- sd(dat$Score)
if (!is.finite(score_sd) || score_sd <= 0) {
  stop("Reference scores for ancestry '", ancestry, "' have zero/non-finite variance (sd=", score_sd, ") - ",
       "likely a monomorphic or degenerate weight file for this ancestry+mode. Cannot build an empirical/",
       "regression normalization model from a constant score.")
}

fit <- lm(as.formula(paste("Score ~", paste(pc_cols, collapse = " + "))), data = dat)
resid_sd <- sigma(fit)
if (!is.finite(resid_sd) || resid_sd <= 0) {
  stop("PC-regression residual variance for ancestry '", ancestry, "' is zero/non-finite (resid_sd=", resid_sd,
       ") - the regression fit degenerately (e.g. too few distinct score values). Cannot compute Z_norm1.")
}

model <- list(
  ancestry = ancestry,
  n = nrow(dat),
  mean = mean(dat$Score),
  sd = score_sd,
  sorted_scores = sort(dat$Score),
  pc_cols = pc_cols,
  fit = fit,
  resid_sd = resid_sd
)
saveRDS(model, out_path)
message(sprintf("[%s] n=%d empirical mean=%.4g sd=%.4g | regression on %s, resid_sd=%.4g -> %s",
                 ancestry, model$n, model$mean, model$sd, paste(pc_cols, collapse = ","), model$resid_sd, out_path))
