#!/usr/bin/env Rscript
# Per-run (step2b_ancestry.slurm): classify each target sample's most-similar
# reference population from its projected PCs (Mahalanobis distance, mirroring
# pgsc_calc's alternative --ancestry_method Mahalanobis), and write one --keep
# list per ancestry for step3_score.slurm to score against.
#
# NOTE ON COLUMN DETECTION: plink2's exact output column naming for a
# multi-column `--score ... --score-col-nums a-b header-read` projection was
# not verified against a live plink2 binary when this script was written. It
# tries two plausible naming schemes (see detect_pc_cols()) and errors with the
# actual column names if neither matches - if that happens, inspect the
# projection .sscore header and adjust detect_pc_cols() accordingly.
#
# Usage:
#   Rscript classify_ancestry.R --proj-sscore target_proj.sscore
#       --pop-models pop_models.rds --pval-thresh 1e-10
#       --ancestries EUR,EAS,SAS,AFR --min-pcs 5
#       --out-calls ancestry_calls.tsv --out-keep-dir task2b_ancestry/

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

proj_path    <- get_arg("--proj-sscore"); if (is.null(proj_path))    stop("Provide --proj-sscore")
pop_models_p <- get_arg("--pop-models");  if (is.null(pop_models_p)) stop("Provide --pop-models")
pval_thresh  <- as.numeric(get_arg("--pval-thresh", "1e-10"))
ancestries   <- strsplit(get_arg("--ancestries"), ",")[[1]]
if (length(ancestries) == 0) stop("Provide --ancestries as a comma-separated list")
min_pcs      <- as.integer(get_arg("--min-pcs")); if (is.na(min_pcs)) stop("Provide --min-pcs")
out_calls    <- get_arg("--out-calls");   if (is.null(out_calls))    stop("Provide --out-calls")
out_keep_dir <- get_arg("--out-keep-dir"); if (is.null(out_keep_dir)) stop("Provide --out-keep-dir")
dir.create(out_keep_dir, showWarnings = FALSE, recursive = TRUE)

detect_pc_cols <- function(nms, min_pcs) {
  m <- grep("^PC([0-9]+)(_AVG|_SUM)?$", nms, value = TRUE)
  if (length(m) >= min_pcs) {
    ord <- as.integer(sub("^PC([0-9]+).*$", "\\1", m))
    return(m[order(ord)])
  }
  m2 <- grep("^SCORE([0-9]+)_(AVG|SUM)$", nms, value = TRUE)
  if (length(m2) >= min_pcs) {
    ord <- as.integer(sub("^SCORE([0-9]+)_.*$", "\\1", m2))
    return(m2[order(ord)])
  }
  stop("Could not identify >= ", min_pcs, " PC/score columns in ", proj_path, "\n",
       "Columns present: ", paste(nms, collapse = ", "), "\n",
       "plink2's --score output column naming may differ from what this script expects - ",
       "inspect the header above and adjust detect_pc_cols() in classify_ancestry.R.")
}

proj <- read.delim(proj_path, check.names = FALSE, stringsAsFactors = FALSE)
names(proj) <- sub("^#", "", names(proj))
if (!"IID" %in% names(proj)) stop("No IID column in ", proj_path)
pc_cols <- detect_pc_cols(names(proj), min_pcs)
message("Using projection columns (in PC order): ", paste(pc_cols, collapse = ", "))

pcs <- as.matrix(proj[, pc_cols, drop = FALSE])
rownames(pcs) <- proj$IID
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))  # normalize names regardless of source naming

pm <- readRDS(pop_models_p)
npcs_classify <- pm$npcs
if (ncol(pcs) < npcs_classify) {
  stop("pop_models.rds was built with npcs=", npcs_classify,
       " but only ", ncol(pcs), " PC columns were detected in ", proj_path)
}
classify_pcs <- pcs[, seq_len(npcs_classify), drop = FALSE]

pop_names <- names(pm$models)
pvals <- sapply(pop_names, function(p) {
  mdl <- pm$models[[p]]
  d2 <- mahalanobis(classify_pcs, center = mdl$mean, cov = mdl$cov)
  pchisq(d2, df = npcs_classify, lower.tail = FALSE)
})
if (is.null(dim(pvals))) pvals <- matrix(pvals, nrow = 1, dimnames = list(NULL, pop_names))  # 1-sample edge case
colnames(pvals) <- pop_names

most_similar <- pop_names[max.col(pvals, ties.method = "first")]
max_p <- apply(pvals, 1, max)
low_conf <- max_p < pval_thresh

calls <- data.frame(IID = rownames(pcs), MostSimilarPop = most_similar, LowConfidence = low_conf,
                     check.names = FALSE)
colnames(pvals) <- paste0("p_", pop_names)
calls <- cbind(calls, pvals, pcs)
write.table(calls, out_calls, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Classified %d samples -> %s", nrow(calls), out_calls))
print(table(calls$MostSimilarPop, useNA = "ifany"))
if (any(low_conf)) message(sprintf("%d sample(s) flagged LowConfidence (max p < %g)", sum(low_conf), pval_thresh))

unmapped <- setdiff(unique(most_similar), ancestries)
if (length(unmapped) > 0) {
  message("[WARN] samples were assigned to reference population(s) with no matching entry in ",
          "ANCESTRY_MAP: ", paste(unmapped, collapse = ", "),
          " - those samples have no ancestry-specific score file and will be skipped in step3.")
}

for (a in ancestries) {
  ids <- calls$IID[calls$MostSimilarPop == a]
  writeLines(ids, file.path(out_keep_dir, paste0("keep_", a, ".txt")))
  message(sprintf("  keep_%s.txt: %d sample(s)", a, length(ids)))
}
