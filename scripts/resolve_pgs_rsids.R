#!/usr/bin/env Rscript
# Helper for slurm/resolve_pgs_rsids.slurm: turn position-only weights files
# (e.g. PGS Catalog scoring files with an empty hm_rsID column) into rsID
# weights files, using the rsIDs in the imputed data itself.
#
#   ranges : union of every input file's GRCh38 positions -> plink2 bed1 file
#   join   : match one input file to the extracted imputed variants by
#            chr + pos + alleles, write rsID / effect_allele / other_allele /
#            effect_weight / chr / pos, and list what could not be matched.
#
# Positions: hm_chr/hm_pos (PGS Catalog harmonized, requires #HmPOS_build=GRCh38),
# else chr_name/chr_position when the file's #genome_build is GRCh38.
#
# Usage:
#   Rscript resolve_pgs_rsids.R ranges --out-bed ranges.bed <file> [<file> ...]
#   Rscript resolve_pgs_rsids.R join --pgs <file> --pvar imputed.pvar
#       --out <ANC>.tsv --unmatched <ANC>.unmatched.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: resolve_pgs_rsids.R ranges|join ...")
cmd <- args[[1]]
args <- args[-1]
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

# Read a PGS-style file: leading "#" lines are metadata (key=value), then a
# tab-separated table. Returns a tibble CHR, POS, EA, OA, WEIGHT (+ row number).
read_pgs <- function(path) {
  head_lines <- readLines(path, n = 1000, warn = FALSE)
  n_meta <- match(FALSE, startsWith(head_lines, "#"), nomatch = length(head_lines) + 1) - 1
  meta_lines <- sub("^#+", "", head_lines[seq_len(n_meta)])
  kv <- strsplit(meta_lines[grepl("=", meta_lines)], "=", fixed = TRUE)
  meta <- setNames(vapply(kv, function(x) paste(x[-1], collapse = "="), ""),
                   vapply(kv, `[`, "", 1))

  raw <- read_tsv(path, skip = n_meta, col_types = cols(.default = col_character()),
                  show_col_types = FALSE)
  need <- function(col) {
    if (!col %in% names(raw)) stop(path, ": no '", col, "' column. Columns present: ",
                                   paste(names(raw), collapse = ", "))
    raw[[col]]
  }

  has_hm <- all(c("hm_chr", "hm_pos") %in% names(raw)) && any(!is.na(raw$hm_pos) & raw$hm_pos != "")
  if (has_hm) {
    hm_build <- meta[["HmPOS_build"]]
    if (is.null(hm_build) || is.na(hm_build) || hm_build != "GRCh38") {
      stop(path, ": harmonized positions are ", hm_build, ", but the imputed data is GRCh38. ",
           "Download the *_hmPOS_GRCh38.txt.gz version of this score.")
    }
    chr <- raw$hm_chr; pos <- raw$hm_pos; src <- "hm_chr/hm_pos (GRCh38)"
  } else {
    build <- meta[["genome_build"]]
    if (is.null(build) || is.na(build) || !build %in% c("GRCh38", "hg38")) {
      stop(path, ": no harmonized GRCh38 positions and #genome_build=", build,
           " - download the *_hmPOS_GRCh38.txt.gz version of this score.")
    }
    chr <- need("chr_name"); pos <- need("chr_position"); src <- "chr_name/chr_position (GRCh38)"
  }

  oa <- if ("other_allele" %in% names(raw)) raw$other_allele
        else if ("hm_inferOtherAllele" %in% names(raw)) raw$hm_inferOtherAllele
        else NA_character_

  message(sprintf("%s: %d variants, positions from %s, weight_type=%s",
                  basename(path), nrow(raw), src,
                  if (is.null(meta[["weight_type"]])) "unspecified" else meta[["weight_type"]]))
  tibble(
    row    = seq_len(nrow(raw)),
    CHR    = sub("^chr", "", chr, ignore.case = TRUE),
    POS    = suppressWarnings(as.integer(pos)),
    EA     = toupper(need("effect_allele")),
    OA     = toupper(oa),
    WEIGHT = suppressWarnings(as.numeric(need("effect_weight")))
  )
}

if (cmd == "ranges") {
  out_bed <- get_arg("--out-bed"); if (is.null(out_bed)) stop("Provide --out-bed")
  files <- args[!args %in% c("--out-bed", out_bed)]
  if (length(files) == 0) stop("No input files given")
  pos <- bind_rows(lapply(files, read_pgs)) %>%
    filter(!is.na(POS), !is.na(CHR), CHR != "") %>%
    distinct(CHR, POS) %>%
    arrange(CHR, POS)
  # bed1 = 1-based, inclusive start/end
  write_tsv(pos %>% transmute(CHR, START = POS, END = POS), out_bed, col_names = FALSE)
  message(sprintf("Wrote %d unique positions -> %s", nrow(pos), out_bed))

} else if (cmd == "join") {
  pgs_path  <- get_arg("--pgs");       if (is.null(pgs_path))  stop("Provide --pgs")
  pvar_path <- get_arg("--pvar");      if (is.null(pvar_path)) stop("Provide --pvar")
  out_path  <- get_arg("--out");       if (is.null(out_path))  stop("Provide --out")
  unm_path  <- get_arg("--unmatched"); if (is.null(unm_path))  stop("Provide --unmatched")

  pgs <- read_pgs(pgs_path)
  pvar <- read_tsv(pvar_path, col_names = c("CHR", "POS", "ID", "REF", "ALT"),
                   col_types = "ciccc", show_col_types = FALSE) %>%
    mutate(CHR = sub("^chr", "", CHR, ignore.case = TRUE), REF = toupper(REF), ALT = toupper(ALT))

  # Allele-aware match: the effect allele must be REF or ALT at that position,
  # and (when given) the other allele must be the remaining one.
  cand <- pgs %>%
    filter(!is.na(POS)) %>%
    inner_join(pvar, by = c("CHR", "POS"), relationship = "many-to-many") %>%
    filter((EA == REF & (is.na(OA) | OA == "" | OA == ALT)) |
           (EA == ALT & (is.na(OA) | OA == "" | OA == REF)))

  # A variant with >1 candidate record (only possible without an other allele)
  # is ambiguous - drop it rather than guess.
  n_cand <- cand %>% count(row, name = "n")
  matched <- cand %>%
    semi_join(n_cand %>% filter(n == 1), by = "row") %>%
    filter(ID != ".", !is.na(ID), !is.na(WEIGHT)) %>%
    distinct(ID, .keep_all = TRUE)

  reason <- pgs %>%
    left_join(n_cand, by = "row") %>%
    left_join(cand %>% group_by(row) %>% summarise(ID = first(ID), .groups = "drop"), by = "row") %>%
    mutate(reason = case_when(
      is.na(POS)            ~ "no position (liftover failed?)",
      is.na(WEIGHT)         ~ "no weight",
      is.na(n)              ~ "position/alleles not in imputed data",
      n > 1                 ~ "ambiguous: several records match",
      ID == "." | is.na(ID) ~ "in imputed data but has no rsID",
      TRUE                  ~ "duplicate rsID"
    )) %>%
    filter(!row %in% matched$row) %>%
    select(CHR, POS, EA, OA, WEIGHT, reason)

  write_tsv(matched %>% transmute(rsID = ID, effect_allele = EA, other_allele = OA,
                                  effect_weight = WEIGHT, chr = CHR, pos = POS),
            out_path)
  write_tsv(reason, unm_path)
  message(sprintf("%s: matched %d / %d variants -> %s", basename(pgs_path),
                  nrow(matched), nrow(pgs), out_path))
  if (nrow(reason) > 0) {
    message(sprintf("  %d unmatched (see %s):", nrow(reason), unm_path))
    print(table(reason$reason))
  }

} else {
  stop("Unknown command '", cmd, "' - use ranges or join")
}
