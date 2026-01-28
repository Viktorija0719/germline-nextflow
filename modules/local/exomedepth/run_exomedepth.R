#!/usr/bin/env Rscript
# run_exomedepth.R
#
# ExomeDepth CNV calling from BAMs using a BED target list (Workflow B style).
# - No enforced QC thresholds (writes summaries only).
# - Optional chrX calling with sex-matched reference pools inferred from samtools idxstats.
#
# Example:
#   Rscript run_exomedepth.R \
#     --fasta /path/ref.fa \
#     --bed /path/targets.bed \
#     --bam_dir /path/bams \
#     --samplesheet /path/samplesheet.csv \
#     --outdir /path/out_exomedepth \
#     --idxstats_dir /path/results/qc/idxstats \
#     --include_x TRUE
#
# Notes:
# - --bed_coord bed0 assumes standard BED (0-based start, half-open). Use --bed_coord one if 1-based.
# - BED chromosomes are stored WITHOUT 'chr' prefix to avoid 'chrchr1' when include.chr=TRUE.
# - include.chr is auto-detected from BAM header unless explicitly set.
# - chrX calling requires BED to include X intervals.

suppressPackageStartupMessages({
  library(ExomeDepth)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
})

# ---------------------------
# Utils / CLI parsing
# ---------------------------
parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!grepl("^--", a)) stop("Unexpected argument: ", a)
    a <- sub("^--", "", a)

    if (grepl("=", a, fixed = TRUE)) {
      kv <- strsplit(a, "=", fixed = TRUE)[[1]]
      key <- kv[1]; val <- kv[2]
    } else {
      key <- a
      if (i == length(argv)) stop("Missing value for --", key)
      val <- argv[i + 1]
      i <- i + 1
    }
    out[[key]] <- val
    i <- i + 1
  }
  out
}

as_logical <- function(x) {
  if (is.null(x)) return(NA)
  x <- tolower(trimws(x))
  if (x %in% c("true","t","1","yes","y")) return(TRUE)
  if (x %in% c("false","f","0","no","n")) return(FALSE)
  NA
}

stop_if_missing <- function(args, keys) {
  miss <- keys[!keys %in% names(args)]
  if (length(miss) > 0) stop("Missing required arguments: ", paste0("--", miss, collapse = ", "))
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

# Create a robust reference-set row (even if choice is empty)
make_ref_rows <- function(sample, set_name, reference_choice) {
  if (length(reference_choice) == 0) {
    return(data.frame(sample = sample, set = set_name, reference_sample = NA_character_, stringsAsFactors = FALSE))
  }
  data.frame(sample = sample, set = set_name, reference_sample = reference_choice, stringsAsFactors = FALSE)
}

# Add sample/set columns safely to a (possibly 0-row) data.frame
stamp_calls <- function(df, sample, set_name) {
  if (is.null(df)) return(NULL)
  if (!is.data.frame(df)) return(NULL)
  df$sample <- rep(sample, nrow(df))
  df$set    <- rep(set_name, nrow(df))
  df
}

# ---------------------------
# Input readers
# ---------------------------
read_samplesheet_samples <- function(samplesheet_csv) {
  ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% names(ss)) {
    stop("samplesheet.csv must contain a column named 'sample'. Found: ", paste(names(ss), collapse = ", "))
  }
  samples <- unique(ss$sample)
  samples <- samples[!is.na(samples) & nzchar(samples)]
  if (length(samples) < 2) stop("Need at least 2 samples for ExomeDepth. Found: ", length(samples))
  samples
}

read_bed_targets <- function(bed_path, bed_coord = c("bed0", "one")) {
  bed_coord <- match.arg(bed_coord)

  lines <- readLines(bed_path, warn = FALSE)
  lines <- lines[!grepl("^(track|browser)\\b", lines)]
  if (length(lines) == 0) stop("BED appears empty after removing header lines: ", bed_path)

  con <- textConnection(lines)
  on.exit(close(con), add = TRUE)

  bed <- tryCatch(
    read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = ""),
    error = function(e) stop("Failed to read BED: ", bed_path, "\n", conditionMessage(e))
  )
  if (ncol(bed) < 3) stop("BED must have at least 3 columns (chrom, start, end). Found: ", ncol(bed))

  chrom_raw <- as.character(bed[[1]])
  start_in <- suppressWarnings(as.integer(bed[[2]]))
  end_in   <- suppressWarnings(as.integer(bed[[3]]))
  if (anyNA(start_in) || anyNA(end_in)) stop("BED start/end columns contain non-integer values.")

  start1 <- if (bed_coord == "bed0") start_in + 1L else start_in
  end1 <- end_in

  nm <- if (ncol(bed) >= 4) as.character(bed[[4]]) else paste0("region_", seq_len(nrow(bed)))
  nm[is.na(nm) | !nzchar(nm)] <- paste0("region_", which(is.na(nm) | !nzchar(nm)))

  # CRITICAL: store WITHOUT 'chr' prefix; include.chr will add it if needed.
  chrom <- strip_chr(chrom_raw)

  df <- data.frame(
    chromosome = chrom,
    start = start1,
    end = end1,
    name = make.unique(nm),
    stringsAsFactors = FALSE
  )

  bad <- which(df$start > df$end)
  if (length(bad) > 0) stop("BED has intervals with start > end. Example row: ", bad[1])

  df
}

# ---------------------------
# BAM discovery + include.chr
# ---------------------------
find_bams_for_samples <- function(samples, bam_dir) {
  all_bams <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
  if (length(all_bams) == 0) stop("No .bam files found in bam_dir: ", bam_dir)

  bam_map <- setNames(rep(NA_character_, length(samples)), samples)

  for (s in samples) {
    hits <- all_bams[grepl(paste0("^", s, ".*\\.bam$"), basename(all_bams))]
    if (length(hits) == 0) hits <- all_bams[grepl(s, basename(all_bams), fixed = TRUE)]
    if (length(hits) == 0) stop("No BAM found for sample '", s, "' in ", bam_dir)

    exact <- hits[basename(hits) == paste0(s, ".bam")]
    if (length(exact) == 1) {
      bam_map[[s]] <- exact
    } else {
      hits_sorted <- hits[order(nchar(basename(hits)))]
      if (length(hits_sorted) > 1) {
        message("Multiple BAMs matched sample '", s, "'. Using: ", basename(hits_sorted[1]),
                "  (candidates: ", paste(basename(hits_sorted), collapse = ", "), ")")
      }
      bam_map[[s]] <- hits_sorted[1]
    }
  }

  # Validate BAM + index presence
  for (s in samples) {
    bam <- bam_map[[s]]
    if (!file.exists(bam)) stop("BAM missing on disk: ", bam)
    bai1 <- paste0(bam, ".bai")
    bai2 <- sub("\\.bam$", ".bai", bam)
    bai3 <- paste0(bam, ".bam.bai")
    if (!(file.exists(bai1) || file.exists(bai2) || file.exists(bai3))) {
      stop("Missing BAM index (.bai) for sample '", s, "': ", bam,
           "\nExpected one of: ", bai1, " OR ", bai2, " OR ", bai3)
    }
  }

  bam_map
}

detect_include_chr_from_bam <- function(bam_path) {
  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    hdr <- Rsamtools::scanBamHeader(bam_path)
    targets <- hdr[[1]]$targets
    if (!is.null(names(targets)) && length(names(targets)) > 0) {
      return(any(grepl("^chr", names(targets))))
    }
  }
  FALSE
}

# ---------------------------
# Option B: Infer sex from idxstats
# ---------------------------
read_idxstats <- function(path) {
  df <- read.table(path, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  colnames(df) <- c("contig","length","mapped","unmapped")
  df$contig <- strip_chr(df$contig)
  df
}

find_idxstats_for_sample <- function(sample, idxstats_dir) {
  exact <- file.path(idxstats_dir, paste0(sample, ".idxstats"))
  if (file.exists(exact)) return(exact)

  hits <- list.files(idxstats_dir, pattern = "\\.idxstats$", full.names = TRUE)
  hits <- hits[grepl(sample, basename(hits), fixed = TRUE)]
  if (length(hits) == 0) return(NA_character_)

  hits <- hits[order(nchar(basename(hits)))]
  if (length(hits) > 1) {
    message("Multiple idxstats matched sample '", sample, "'. Using: ", basename(hits[1]),
            "  (candidates: ", paste(basename(hits), collapse = ", "), ")")
  }
  hits[1]
}

infer_sex_from_idxstats_dir <- function(samples, idxstats_dir) {
  if (is.null(idxstats_dir) || is.na(idxstats_dir) || !nzchar(idxstats_dir)) {
    stop("idxstats_dir is not set.")
  }
  if (!dir.exists(idxstats_dir)) stop("idxstats_dir does not exist: ", idxstats_dir)

  res <- lapply(samples, function(sample) {
    f <- find_idxstats_for_sample(sample, idxstats_dir)
    if (is.na(f) || !file.exists(f)) {
      return(data.frame(sample=sample, idxstats_file=NA_character_, x_auto_ratio=NA_real_,
                        y_mapped=NA_real_, sex="U", stringsAsFactors=FALSE))
    }

    df <- read_idxstats(f)

    get_mapped <- function(ctg) {
      v <- df$mapped[df$contig == ctg]
      if (length(v) == 0) 0 else sum(v, na.rm = TRUE)
    }
    get_len <- function(ctg) {
      v <- df$length[df$contig == ctg]
      if (length(v) == 0) NA_real_ else sum(v, na.rm = TRUE)
    }

    x_m <- get_mapped("X")
    y_m <- get_mapped("Y")

    auto_ctg <- as.character(1:22)
    auto_m <- sum(df$mapped[df$contig %in% auto_ctg], na.rm = TRUE)
    auto_l <- sum(df$length[df$contig %in% auto_ctg], na.rm = TRUE)

    x_l <- get_len("X")

    x_norm <- if (!is.na(x_l) && x_l > 0) x_m / x_l else NA_real_
    auto_norm <- if (!is.na(auto_l) && auto_l > 0) auto_m / auto_l else NA_real_
    ratio <- x_norm / (auto_norm + 1e-12)

    sex <- "U"
    if (!is.na(y_m) && y_m > 0) {
      sex <- "M"
    } else if (!is.na(ratio)) {
      if (ratio < 0.75) sex <- "M"
      if (ratio > 0.85) sex <- "F"
    }

    data.frame(sample=sample, idxstats_file=f, x_auto_ratio=ratio, y_mapped=y_m, sex=sex,
               stringsAsFactors=FALSE)
  })

  do.call(rbind, res)
}

# ---------------------------
# getBamCounts conversion (robust)
# ---------------------------
get_counts_objects <- function(targets_df, bam_files, include_chr, fasta) {
  counts_obj <- ExomeDepth::getBamCounts(
    bed.frame = targets_df,
    bam.files = bam_files,
    include.chr = include_chr,
    referenceFasta = fasta
  )

  if (inherits(counts_obj, "GRanges")) {
    counts_gr <- counts_obj
    counts_df <- as.data.frame(counts_gr)

    chrom_vec <- as.character(GenomeInfoDb::seqnames(counts_gr))
    if (length(chrom_vec) != nrow(counts_df)) {
      stop("Internal error: seqnames length (", length(chrom_vec),
           ") != nrow(counts_df) (", nrow(counts_df), ").")
    }
    counts_df$chromosome <- chrom_vec

    if (!"names" %in% names(counts_df)) {
      nm <- NULL
      if ("name" %in% names(mcols(counts_gr))) nm <- mcols(counts_gr)$name
      if (is.null(nm) && "names" %in% names(mcols(counts_gr))) nm <- mcols(counts_gr)$names
      if (is.null(nm)) nm <- paste0("region_", seq_len(nrow(counts_df)))
      counts_df$names <- as.character(nm)
    }

    return(list(counts_df = counts_df))
  }

  if (is.data.frame(counts_obj)) {
    counts_df <- counts_obj
    if ("seqnames" %in% names(counts_df)) {
      counts_df$chromosome <- as.character(counts_df$seqnames)
    } else if (!"chromosome" %in% names(counts_df)) {
      stop("getBamCounts returned a data.frame without seqnames/chromosome columns.")
    }
    if (!"names" %in% names(counts_df)) counts_df$names <- paste0("region_", seq_len(nrow(counts_df)))
    return(list(counts_df = counts_df))
  }

  stop("Unexpected getBamCounts output class: ", paste(class(counts_obj), collapse = ", "))
}

# ---------------------------
# ExomeDepth runner for subset rows (safe)
# ---------------------------
run_exomedepth_subset_safe <- function(ExomeCount_mat, counts_df, row_idx,
                                      test_sample, ref_candidates,
                                      n_bins_reduced, transition_prob) {
  if (!any(row_idx)) {
    return(list(calls = NULL, reference_choice = character(0), cor_test_vs_reference = NA_real_))
  }

  tryCatch({
    test_counts <- ExomeCount_mat[row_idx, test_sample]
    ref_mat <- ExomeCount_mat[row_idx, ref_candidates, drop = FALSE]

    bin_length <- (counts_df$end[row_idx] - counts_df$start[row_idx]) / 1000

    choice <- ExomeDepth::select.reference.set(
      test.counts = test_counts,
      reference.counts = as.matrix(ref_mat),
      bin.length = bin_length,
      n.bins.reduced = n_bins_reduced
    )

    ref_selected <- apply(
      X = ExomeCount_mat[row_idx, choice$reference.choice, drop = FALSE],
      MAR = 1,
      FUN = sum
    )

    cor_val <- suppressWarnings(cor(test_counts, ref_selected, use = "pairwise.complete.obs"))

    ed <- new(
      "ExomeDepth",
      test = test_counts,
      reference = ref_selected,
      formula = "cbind(test, reference) ~ 1"
    )

    ed <- ExomeDepth::CallCNVs(
      x = ed,
      transition.probability = transition_prob,
      chromosome = counts_df$chromosome[row_idx],
      start = counts_df$start[row_idx],
      end = counts_df$end[row_idx],
      name = counts_df$names[row_idx]
    )

    list(
      calls = ed@CNV.calls,
      reference_choice = choice$reference.choice,
      cor_test_vs_reference = cor_val
    )
  }, error = function(e) {
    message("ExomeDepth failed for sample ", test_sample, " on subset: ", conditionMessage(e))
    list(calls = NULL, reference_choice = character(0), cor_test_vs_reference = NA_real_)
  })
}

# ---------------------------
# Main
# ---------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))
stop_if_missing(args, c("fasta", "bed", "bam_dir", "samplesheet", "outdir"))

fasta <- normalizePath(args$fasta, mustWork = TRUE)
bed_path <- normalizePath(args$bed, mustWork = TRUE)
bam_dir <- normalizePath(args$bam_dir, mustWork = TRUE)
samplesheet <- normalizePath(args$samplesheet, mustWork = TRUE)
outdir <- args$outdir

bed_coord <- if (!is.null(args$bed_coord)) args$bed_coord else "bed0"
include_chr_arg <- if (!is.null(args$include_chr)) tolower(args$include_chr) else "auto"
transition_prob <- if (!is.null(args$transition_prob)) as.numeric(args$transition_prob) else 1e-4
n_bins_reduced  <- if (!is.null(args$n_bins_reduced)) as.integer(args$n_bins_reduced) else 10000L

idxstats_dir <- if (!is.null(args$idxstats_dir)) args$idxstats_dir else NA_character_
include_x <- if (!is.null(args$include_x)) as_logical(args$include_x) else FALSE
if (is.na(include_x)) stop("--include_x must be TRUE or FALSE")

safe_mkdir(outdir)
safe_mkdir(file.path(outdir, "calls"))

samples <- read_samplesheet_samples(samplesheet)
bam_map <- find_bams_for_samples(samples, bam_dir)

bam_files <- unname(bam_map[samples])
names(bam_files) <- samples

include_chr <- if (include_chr_arg == "auto") {
  detect_include_chr_from_bam(bam_files[[1]])
} else {
  v <- as_logical(include_chr_arg)
  if (is.na(v)) stop("--include_chr must be one of: auto, TRUE, FALSE")
  v
}

targets_df <- read_bed_targets(bed_path, bed_coord = bed_coord)

# Persist manifests
write.table(
  data.frame(sample = samples, bam = bam_files, stringsAsFactors = FALSE),
  file = file.path(outdir, "bam_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(
    fasta = fasta,
    bed = bed_path,
    bam_dir = bam_dir,
    samplesheet = samplesheet,
    idxstats_dir = ifelse(is.na(idxstats_dir), "", idxstats_dir),
    bed_coord = bed_coord,
    include_chr = include_chr,
    include_x = include_x,
    transition_prob = transition_prob,
    n_bins_reduced = n_bins_reduced,
    stringsAsFactors = FALSE
  ),
  file = file.path(outdir, "run_config.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Count reads
message("Counting reads with ExomeDepth::getBamCounts ...")
counts <- get_counts_objects(
  targets_df = targets_df,
  bam_files = bam_files,
  include_chr = include_chr,
  fasta = fasta
)
counts_df <- counts$counts_df

# Identify numeric count columns and rename to sample names
known_cols <- c("seqnames","start","end","width","strand","GC","names","chromosome")
candidate <- setdiff(names(counts_df), known_cols)
candidate <- candidate[sapply(counts_df[candidate], is.numeric)]
if (length(candidate) != length(samples)) {
  stop(
    "Could not reliably identify count columns from getBamCounts output.\n",
    "Expected ", length(samples), " numeric count columns, found ", length(candidate), ".\n",
    "Columns seen: ", paste(names(counts_df), collapse = ", ")
  )
}
names(counts_df)[match(candidate, names(counts_df))] <- samples
ExomeCount_mat <- as.matrix(counts_df[, samples, drop = FALSE])

# Subset indices for autosomes / X using stripped chromosome labels
chrom_naked <- strip_chr(counts_df$chromosome)
auto_rows <- chrom_naked %in% as.character(1:22)
x_rows    <- chrom_naked %in% "X"

if (!any(auto_rows)) stop("No autosome targets detected (1-22) in BED/counts.")
if (isTRUE(include_x) && !any(x_rows)) {
  stop("include_x=TRUE but no chrX targets detected in BED/counts. Ensure BED includes X intervals.")
}

# Option B: infer sex from idxstats (for chrX reference pools)
sex_df <- data.frame(sample = samples, sex = "U", stringsAsFactors = FALSE)
if (isTRUE(include_x)) {
  if (!is.na(idxstats_dir) && nzchar(idxstats_dir)) {
    idxstats_dir <- normalizePath(idxstats_dir, mustWork = TRUE)
    sex_full <- infer_sex_from_idxstats_dir(samples, idxstats_dir)
    sex_df <- sex_full[, c("sample", "sex")]
    write.table(sex_full, file = file.path(outdir, "inferred_sex_from_idxstats.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    message("include_x=TRUE but --idxstats_dir not provided; sex='U' and chrX uses all-sample references.")
    write.table(sex_df, file = file.path(outdir, "inferred_sex_from_idxstats.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
} else {
  write.table(sex_df, file = file.path(outdir, "inferred_sex_from_idxstats.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
sex_map <- setNames(sex_df$sex, sex_df$sample)

# Summary collectors
summary_rows <- vector("list", length(samples))
ref_rows <- list()

for (i in seq_along(samples)) {
  s <- samples[i]
  message("\n==============================")
  message("Sample: ", s, " (", i, "/", length(samples), ")")
  message("==============================")

  ref_all <- setdiff(samples, s)

  # AUTOSOMES
  auto_res <- run_exomedepth_subset_safe(
    ExomeCount_mat = ExomeCount_mat,
    counts_df = counts_df,
    row_idx = auto_rows,
    test_sample = s,
    ref_candidates = ref_all,
    n_bins_reduced = n_bins_reduced,
    transition_prob = transition_prob
  )

  calls_auto <- auto_res$calls
  if (!is.null(calls_auto) && nrow(calls_auto) > 0 && "BF" %in% names(calls_auto)) {
    calls_auto <- calls_auto[order(calls_auto$BF, decreasing = TRUE), , drop = FALSE]
  }
  calls_auto <- stamp_calls(calls_auto, s, "autosomes")

  # chrX
  calls_x <- NULL
  x_res <- list(reference_choice = character(0), cor_test_vs_reference = NA_real_, calls = NULL)

  if (isTRUE(include_x)) {
    sx <- sex_map[[s]]
    if (!is.null(sx) && sx %in% c("M", "F")) {
      ref_same_sex <- setdiff(samples[sex_map[samples] == sx], s)
      if (length(ref_same_sex) < 1) ref_same_sex <- ref_all
    } else {
      ref_same_sex <- ref_all
    }

    x_res <- run_exomedepth_subset_safe(
      ExomeCount_mat = ExomeCount_mat,
      counts_df = counts_df,
      row_idx = x_rows,
      test_sample = s,
      ref_candidates = ref_same_sex,
      n_bins_reduced = n_bins_reduced,
      transition_prob = transition_prob
    )

    calls_x <- x_res$calls
    if (!is.null(calls_x) && nrow(calls_x) > 0 && "BF" %in% names(calls_x)) {
      calls_x <- calls_x[order(calls_x$BF, decreasing = TRUE), , drop = FALSE]
    }
    calls_x <- stamp_calls(calls_x, s, "chrX")
  }

  # Combine calls without rbind problems
  calls_all <- NULL
  if (!is.null(calls_auto) && nrow(calls_auto) > 0) calls_all <- calls_auto
  if (!is.null(calls_x) && nrow(calls_x) > 0) {
    if (is.null(calls_all)) calls_all <- calls_x else calls_all <- rbind(calls_all, calls_x)
  }
  if (is.null(calls_all)) calls_all <- data.frame()

  out_csv <- file.path(outdir, "calls", paste0(s, ".exomedepth.cnv.csv"))
  write.csv(calls_all, file = out_csv, row.names = FALSE)

  # Summary
  n_calls_auto <- if (!is.null(calls_auto)) nrow(calls_auto) else 0L
  n_calls_x    <- if (!is.null(calls_x)) nrow(calls_x) else 0L

  summary_rows[[i]] <- data.frame(
    sample = s,
    bam = bam_files[[s]],
    inferred_sex = sex_map[[s]],
    include_chr = include_chr,
    n_ref_autosomes = length(auto_res$reference_choice),
    cor_autosomes = as.numeric(auto_res$cor_test_vs_reference),
    n_calls_autosomes = n_calls_auto,
    n_ref_chrX = if (isTRUE(include_x)) length(x_res$reference_choice) else NA_integer_,
    cor_chrX = if (isTRUE(include_x)) as.numeric(x_res$cor_test_vs_reference) else NA_real_,
    n_calls_chrX = if (isTRUE(include_x)) n_calls_x else NA_integer_,
    stringsAsFactors = FALSE
  )

  # Reference set tables (always safe)
  ref_rows[[length(ref_rows) + 1]] <- make_ref_rows(s, "autosomes", auto_res$reference_choice)
  if (isTRUE(include_x)) {
    ref_rows[[length(ref_rows) + 1]] <- make_ref_rows(s, "chrX", x_res$reference_choice)
  }
}

summary_df <- do.call(rbind, summary_rows)
ref_df <- do.call(rbind, ref_rows)

write.table(summary_df, file = file.path(outdir, "exomedepth_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ref_df, file = file.path(outdir, "reference_sets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\nDone.")
message("Calls written to: ", file.path(outdir, "calls"))
message("Summary: ", file.path(outdir, "exomedepth_summary.tsv"))
message("Reference sets: ", file.path(outdir, "reference_sets.tsv"))
message("Sex inference (idxstats): ", file.path(outdir, "inferred_sex_from_idxstats.tsv"))
