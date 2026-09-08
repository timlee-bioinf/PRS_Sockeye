#!/usr/bin/env Rscript
# One-time setup helper (called by setup_ancestry_reference.sh): fit a
# per-population Mahalanobis model (mean vector + covariance) on the reference
# panel's own PCA, for classifying target samples later in classify_ancestry.R.
#
# Usage:
#   Rscript build_pop_models.R --eigenvec ref_pca.eigenvec --psam ref.psam
#       --label-col SuperPop --npcs 5 --out pop_models.rds
#       [--exclude related_ids.txt] [--unrelated-out ref_unrelated.txt]

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
# which is IID either way. NOTE: --eigenvec is expected to already come from a
# plink2 run that applied this same exclude file via --remove (see
# setup_ancestry_reference.sh), so this filter is normally a no-op safety net,
# not the primary exclusion mechanism - kept correct anyway in case these
# scripts are ever invoked with an --eigenvec that wasn't pre-filtered.
read_id_list <- function(path) {
  lines <- trimws(readLines(path))
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  vapply(strsplit(lines, "[ \t]+"), function(x) x[length(x)], character(1))
}

eigenvec_path <- get_arg("--eigenvec"); if (is.null(eigenvec_path)) stop("Provide --eigenvec")
psam_path     <- get_arg("--psam");     if (is.null(psam_path))     stop("Provide --psam")
label_col     <- get_arg("--label-col", "SuperPop")
npcs          <- as.integer(get_arg("--npcs", "5"))
out_path      <- get_arg("--out");      if (is.null(out_path))      stop("Provide --out")
exclude_path  <- get_arg("--exclude", NA)
unrelated_out <- get_arg("--unrelated-out", NA)

read_plink_table <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(x) <- sub("^#", "", names(x))
  x
}

eig <- read_plink_table(eigenvec_path)
iid_col_eig <- if ("IID" %in% names(eig)) "IID" else stop("No IID column in ", eigenvec_path)
pc_cols <- grep("^PC[0-9]+$", names(eig), value = TRUE)
if (length(pc_cols) < npcs) {
  stop("Requested --npcs ", npcs, " but ", eigenvec_path, " only has ", length(pc_cols), " PC columns")
}
pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))][seq_len(npcs)]

psam <- read_plink_table(psam_path)
iid_col_psam <- if ("IID" %in% names(psam)) "IID" else stop("No IID column in ", psam_path)
if (!label_col %in% names(psam)) {
  stop("Label column '", label_col, "' not found in ", psam_path,
       ". Columns present: ", paste(names(psam), collapse = ", "))
}

dat <- eig %>%
  select(IID = all_of(iid_col_eig), all_of(pc_cols)) %>%
  inner_join(psam %>% select(IID = all_of(iid_col_psam), Pop = all_of(label_col)), by = "IID")

if (!is.na(exclude_path)) {
  excl <- read_id_list(exclude_path)
  before <- nrow(dat)
  dat <- dat %>% filter(!IID %in% excl)
  message(sprintf("Excluded %d related/duplicate reference samples (%d -> %d)",
                   before - nrow(dat), before, nrow(dat)))
}

dat <- dat %>% filter(!is.na(Pop), Pop != "")
if (nrow(dat) == 0) stop("No reference samples left after joining PCs to labels - check --label-col and ID matching")

if (!is.na(unrelated_out)) writeLines(dat$IID, unrelated_out)

pops <- unique(dat$Pop)
models <- list()
for (p in pops) {
  sub <- dat %>% filter(Pop == p) %>% select(all_of(pc_cols)) %>% as.matrix()
  if (nrow(sub) <= npcs) {
    stop("Population '", p, "' has only ", nrow(sub), " reference samples, need > --npcs (",
         npcs, ") to estimate a covariance matrix. Reduce --npcs or check this population is real.")
  }
  models[[p]] <- list(mean = colMeans(sub), cov = cov(sub), n = nrow(sub))
}

saveRDS(list(pc_cols = pc_cols, npcs = npcs, models = models), out_path)
message(sprintf("Fit population models for %d population(s) on %d PCs -> %s",
                 length(models), npcs, out_path))
for (p in names(models)) message(sprintf("  %-10s n=%d", p, models[[p]]$n))
