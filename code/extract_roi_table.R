# =============================================================================
# Gate 3 data extraction: build the per-ROI table the size-spectrum + trait engine
# needs, from raw IFCB data. RUN ON THE MACHINE THAT HAS THE RAW DATA.
#
# Because raw-data access is transient, this captures EVERYTHING in one pass:
#   - per-ROI biovolume (um3) + carbon (pg) + cell_count/cell_count_resolved (chains)
#       via iRfcb::ifcb_extract_biovolumes() on the feature-chain-counts branch
#   - per-ROI morphometric/shape features (the classification-free traits)
#       via iRfcb::ifcb_read_features(), feature_version = 4 (SLIM set, ~30 cols).
#       v4 is already slim (no texture/Wedge/Ring/HOG) so we keep ALL columns. v4
#       adds SurfaceArea + RepresentativeWidth -> SA:V is a DIRECT measurement
#       (SurfaceArea/Biovolume), and max/minFeretDiameter give elongation.
#   - per-ROI class + classification confidence score from the production *_class.h5
#
# Sample volume is NOT stored here — it lives in data.txt (SMVOL, the analyzed mL,
# ~5 mL) and is joined later by SMPNO. A per-sample summary table is written for QC.
#
# Output: one .txt per year (one row per ROI) + a per-sample summary + a metadata
# sidecar recording micron_factor, iRfcb branch/commit, classifier, counts.
#
# Deps: iRfcb (feature-chain-counts branch), data.table, dplyr, hdf5r. arrow optional.
# =============================================================================

suppressPackageStartupMessages({
  library(iRfcb); library(data.table); library(dplyr); library(hdf5r)
})

# ---- config (EDIT PATHS for the raw-data machine) ---------------------------
cfg <- list(
  years = c(2016, 2022, 2023, 2024, 2025, 2026),

  # raw .adc/.hdr/.roi per year  (e.g. .../ifcb/data/data/2022)
  raw_root      = "Z:/data/data/%d",
  # extracted v4 SLIM feature files (_fea_v4.csv); created by ifcb_extract_features if absent
  features_root = "Z:/data/features/v4/%d",
  # blob files folder (needed by ifcb_extract_features)
  blobs_root    = "Z:/data/blobs/v4/%d",
  # production classifier output: one D..._class.h5 per sample. Folder names vary
  # (e.g. class2016_v3, class2024_v1, ...) so resolve by GLOB per year rather than
  # assume a suffix. If the glob matches >1 dir the script errors -> use class_dir.
  classified_root = "Z:/data/classified/SMHI-NIVA-SYKE-SAMS-SZN-ResNet50-V6_1",
  class_glob      = "class%d_v3*",                 # %d -> year
  class_dir       = list(),                      # optional per-year override, e.g.
                                                 #   list(`2016` = ".../class2016_v3")

  out_dir       = "Z:/data/out/roi_level_blacklist",          # gitignored; per-year per-ROI .txt
  out_format    = "txt",                      # "txt" (tab, gz) or "parquet"

  micron_factor   = 1/2.77,                   # CURRENT SMHI factor (was 1/3.4)
  feature_version = 4,                         # v4 SLIM features (NOT v2)
  use_cell_counts = TRUE,
  drop_zero_volume = TRUE,
  diatom_equation = "auto",
  carbon_conversion = "cell",
  diatom_include = c("Actinocyclus_spp", "Actinocyclus", "Navicula-like", "Navicula", "cf_Proboscia_rhizosolenia"),
  run_feature_extraction = FALSE,             # TRUE only if feature files not yet built (SLOW)
  parallel = TRUE, n_cores = NULL,

  # v4 is already slim, so keep ALL feature columns (NULL = keep everything except
  # the roi_number join key). v4 cols (confirmed against real data, 2026 sample):
  #   Area Biovolume BoundingBox_xwidth BoundingBox_ywidth ConvexArea ConvexPerimeter
  #   Eccentricity EquivDiameter Extent MajorAxisLength MinorAxisLength Orientation
  #   Perimeter RepresentativeWidth Solidity SurfaceArea maxFeretDiameter
  #   minFeretDiameter numBlobs summed{Area,Biovolume,ConvexArea,ConvexPerimeter,
  #   MajorAxisLength,MinorAxisLength,Perimeter,SurfaceArea} Area_over_Perimeter
  #   Area_over_PerimeterSquared summedConvexPerimeter_over_Perimeter
  # NB the feature 'Biovolume'/'SurfaceArea' are in PIXEL units; biovolume_um3 (from
  # ifcb_extract_biovolumes, micron-scaled) is authoritative for size. For SA:V keep
  # SA and V in the SAME unit system (both pixel, or scale SA by micron_factor^2).
  morpho_cols = NULL,

  # Sample blacklist (tab-separated: sample, reason). Entries are PREFIXES of the
  # sample ID: a full ID (D20240607T102611_IFCB134) drops one sample, a day
  # (D20231205) drops that day, a shorter stem (D2026051) drops every sample whose
  # ID starts with it. Blacklisted samples are removed BEFORE extraction.
  blacklist_file = "C:/R/ifcb-data-pipeline/data/sample_blacklist.tsv"
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----------------------------------------------------------------
strip_sample <- function(x) sub("_(fea.*|class|multiblob).*$", "", tools::file_path_sans_ext(basename(x)))

# read the blacklist; tolerate CRLF, stray whitespace, blank rows and duplicates.
# Entries that are not sample-ID prefixes (e.g. raw-folder names like
# "IFCBdata_11_march_2023_...") cannot be matched here and are reported.
read_blacklist <- function(path) {
  if (!file.exists(path)) stop("blacklist not found: ", path)
  bl <- fread(path, sep = "\t", colClasses = "character", encoding = "UTF-8", fill = TRUE)
  setnames(bl, 1:2, c("entry", "reason"))
  bl[, `:=`(entry = trimws(entry), reason = trimws(reason))]
  bl <- unique(bl[nzchar(entry)], by = "entry")
  bl[, kind := fcase(grepl("^D\\d{8}T\\d{6}_IFCB\\d+$", entry), "sample",
                     grepl("^D\\d{8}$", entry),                     "day",
                     grepl("^D\\d{1,7}$", entry),                   "stem",
                     default = "unmatchable")]
  message(sprintf("blacklist: %d entries (%s)", nrow(bl),
                  paste(bl[, .N, by = kind][, sprintf("%s %d", kind, N)], collapse = ", ")))
  if (any(bl$kind == "stem"))
    message("  short stems match every sample that starts with them: ",
            paste(bl[kind == "stem", entry], collapse = ", "))
  if (any(bl$kind == "unmatchable"))
    warning("blacklist entries that are not sample-ID prefixes (exclude these at folder level): ",
            paste(bl[kind == "unmatchable", entry], collapse = "; "), call. = FALSE)
  bl[kind != "unmatchable"]
}

# first blacklist entry each sample starts with (NA = keep); returns entry + reason
match_blacklist <- function(samples, bl) {
  idx <- rep(NA_integer_, length(samples))
  for (i in seq_len(nrow(bl))) idx[is.na(idx) & startsWith(samples, bl$entry[i])] <- i
  data.table(sample = samples, entry = bl$entry[idx], reason = bl$reason[idx])
}

# resolve the class .h5 folder for a year: explicit override wins; else glob under
# classified_root. Errors on ambiguity so an old dev folder can't be picked silently.
resolve_class_dir <- function(year) {
  ov <- cfg$class_dir[[as.character(year)]]
  if (!is.null(ov)) return(ov)
  hits <- Sys.glob(file.path(cfg$classified_root, sprintf(cfg$class_glob, year)))
  hits <- hits[dir.exists(hits)]
  if (length(hits) == 0) stop("no class folder for ", year, " under ", cfg$classified_root)
  if (length(hits) > 1)
    stop("ambiguous class folders for ", year, " (", paste(basename(hits), collapse = ", "),
         ") -> set cfg$class_dir[['", year, "']] to the production one")
  hits
}

# per-ROI class + confidence from a production _class.h5. Column names follow the
# iRfcb::ifcb_classify_sample() convention: class_name (threshold-adjusted, incl
# "unclassified"), class_name_auto (no threshold), score. model_name attached for
# provenance (verify it is the intended PRODUCTION classifier, not an old dev build).
read_h5_scores <- function(h5_path) {
  h <- H5File$new(h5_path, mode = "r")
  on.exit(h$close_all())
  # read first existing dataset among candidate names (tolerates h5 name generations)
  rd <- function(names) {
    for (nm in names) { v <- tryCatch(h[[nm]]$read(), error = function(e) NULL)
                        if (!is.null(v)) return(v) }
    NULL
  }
  roi <- rd("roi_numbers")
  if (is.null(roi) || length(roi) == 0) return(NULL)          # empty sample -> skip
  n    <- length(roi)
  thr  <- rd(c("class_name", "class_labels_above_threshold")) # CANONICAL: threshold-adjusted (built data.txt)
  auto <- rd(c("class_name_auto", "class_labels_auto"))       # argmax, no threshold
  scr  <- rd("output_scores")
  model <- rd(c("model_name", "classifier_name", "classifierName"))

  # winning probability per ROI, robust to a 1-ROI sample (scr read as a vector)
  # and to either [n_class x n_roi] or [n_roi x n_class] orientation.
  if (is.null(dim(scr))) {
    score <- max(scr)                                         # single ROI -> vector of class scores
  } else {
    roi_dim <- if (ncol(scr) == n) 2L else 1L                 # axis whose length == n ROIs
    score <- apply(scr, roi_dim, max)
  }
  data.table(
    sample          = strip_sample(h5_path),
    roi_number      = as.integer(roi),
    class_name      = thr,
    class_name_auto = auto,
    score           = score,
    model_name      = if (is.null(model)) NA_character_ else model[1]
  )
}

# ---- per-year extraction ----------------------------------------------------
process_year <- function(year) {
  raw    <- sprintf(cfg$raw_root, year)
  feats  <- sprintf(cfg$features_root, year)
  blobs  <- sprintf(cfg$blobs_root, year)
  classd <- resolve_class_dir(year)
  message("\n==== year ", year, " ====")
  message("  class folder: ", classd, "  (VERIFY this is the production classifier)")

  # 0) feature extraction (SLOW; only if needed) — writes _fea_v2.csv into `feats`
  if (cfg$run_feature_extraction) {
    message("  extracting features (slow) ...")
    iRfcb::ifcb_extract_features(data_folder = raw, features_folder = feats,
                                 blobs_folder = blobs, parallel = cfg$parallel,
                                 n_cores = cfg$n_cores, overwrite = FALSE, verbose = TRUE)
  }
  stopifnot(dir.exists(feats))

  # 0b) drop blacklisted samples before anything is read. Same listing rule iRfcb
  # applies to a folder; the class .h5 files are matched to these samples inside
  # ifcb_extract_biovolumes(), and the h5 score list is filtered below as well.
  ff <- list.files(feats, pattern = "D.*\\.csv", full.names = TRUE, recursive = TRUE)
  ff_bl <- match_blacklist(strip_sample(ff), blacklist)
  excluded <- unique(ff_bl[!is.na(entry)], by = "sample")
  ff <- ff[is.na(ff_bl$entry)]
  message(sprintf("  blacklist: %d samples excluded, %d feature files kept",
                  nrow(excluded), length(ff)))
  if (length(ff) == 0) stop("no feature files left for ", year, " after the blacklist")

  # 1) biovolume + carbon + chain counts (one row per ROI)
  message("  ifcb_extract_biovolumes() ...")
  bv <- iRfcb::ifcb_extract_biovolumes(
    feature_files = ff, class_files = classd,
    micron_factor = cfg$micron_factor, use_cell_counts = cfg$use_cell_counts,
    drop_zero_volume = cfg$drop_zero_volume,
    diatom_equation = cfg$diatom_equation,
    carbon_conversion = cfg$carbon_conversion,
    diatom_include = cfg$diatom_include) %>%
    as.data.table()
  bv[, roi_number := as.integer(roi_number)]

  # 2) curated morphometric features
  message("  ifcb_read_features() + select morpho subset ...")
  fl <- iRfcb::ifcb_read_features(ff, feature_version = cfg$feature_version)
  feat <- rbindlist(lapply(names(fl), function(nm) {
    d <- as.data.table(fl[[nm]]); d[, sample := strip_sample(nm)]
    keep <- if (is.null(cfg$morpho_cols)) setdiff(names(d), c("sample", "roi_number"))
            else intersect(cfg$morpho_cols, names(d))
    d[, c("sample", "roi_number", keep), with = FALSE]
  }), use.names = TRUE, fill = TRUE)
  feat[, roi_number := as.integer(roi_number)]

  # 3) class confidence score from production .h5
  message("  reading *_class.h5 scores ...")
  h5s <- list.files(classd, pattern = "_class\\.h5$", full.names = TRUE)
  h5s <- h5s[is.na(match_blacklist(strip_sample(h5s), blacklist)$entry)]
  scores <- rbindlist(lapply(h5s, read_h5_scores), use.names = TRUE, fill = TRUE)

  # 4) assemble (left-join on sample + roi_number; biovolume table is the spine).
  # NB bv already carries a `class` from ifcb_extract_biovolumes(class_files=); we ALSO
  # attach class_thresholded from the h5 — prefer the latter as canonical (it matches
  # data.txt). ALL ROIs are kept (incl. "unclassified") — the size spectrum is
  # classification-free; class only gates the per-group/trait analyses.
  roi <- bv %>%
    merge(feat,   by = c("sample", "roi_number"), all.x = TRUE) %>%
    merge(scores[, .(sample, roi_number, class_name, class_name_auto, score, model_name)],
          by = c("sample", "roi_number"), all.x = TRUE)
  setorder(roi, sample, roi_number)
  leaked <- match_blacklist(unique(roi$sample), blacklist)[!is.na(entry)]
  if (nrow(leaked) > 0) stop("blacklisted samples in the ", year, " table: ",
                             paste(head(leaked$sample, 10), collapse = ", "))

  # 5) write per-ROI table
  out <- file.path(cfg$out_dir, sprintf("roi_%d.%s", year,
                   ifelse(cfg$out_format == "parquet", "parquet", "txt.gz")))
  if (cfg$out_format == "parquet" && requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(roi, out)
  } else {
    fwrite(roi, out, sep = "\t")   # .gz inferred from extension
  }
  message(sprintf("  wrote %s  (%s ROIs, %s samples)", out,
                  format(nrow(roi), big.mark = ","), uniqueN(roi$sample)))

  # 6) per-sample summary for QC (n_roi, n_cells, total biovolume; volume joined later)
  smry <- roi[, .(n_roi = .N,
                  n_cells = sum(cell_count_resolved %||% cell_count, na.rm = TRUE),
                  total_biovol_um3 = sum(biovolume_um3, na.rm = TRUE),
                  median_score = median(score, na.rm = TRUE)), by = sample]
  list(excluded = excluded[, .(year = year, sample, entry, reason)],
       roi_meta = data.table(year = year, n_roi = nrow(roi), n_samples = uniqueN(roi$sample),
                             n_blacklisted = nrow(excluded),
                             class_dir = basename(classd),
                             model_name = scores$model_name[1] %||% NA), summary = smry)
}

`%||%` <- function(a, b) if (is.null(a) || all(is.na(a))) b else a

# ---- run --------------------------------------------------------------------
blacklist <- read_blacklist(cfg$blacklist_file)
results <- lapply(cfg$years, function(y) tryCatch(process_year(y),
                  error = function(e) { message("  YEAR ", y, " FAILED: ", conditionMessage(e)); NULL }))
results <- Filter(Negate(is.null), results)

# per-sample summary across years
smry_all <- rbindlist(lapply(results, `[[`, "summary"), use.names = TRUE, fill = TRUE)
fwrite(smry_all, file.path(cfg$out_dir, "sample_summary.txt"), sep = "\t")

# samples removed by the blacklist, with the entry that matched them
excl_all <- rbindlist(lapply(results, `[[`, "excluded"), use.names = TRUE, fill = TRUE)
fwrite(excl_all, file.path(cfg$out_dir, "excluded_samples.txt"), sep = "\t")

# ---- reproducibility sidecar ------------------------------------------------
desc <- tryCatch(packageDescription("iRfcb"), error = function(e) NULL)
meta <- c(
  sprintf("extracted_on: %s", Sys.Date()),
  sprintf("iRfcb_version: %s", as.character(utils::packageVersion("iRfcb"))),
  sprintf("iRfcb_remote_ref: %s", desc$RemoteRef %||% "NA"),
  sprintf("iRfcb_remote_sha: %s", desc$RemoteSha %||% "NA"),
  sprintf("micron_factor: %s (= %.5f)", "1/2.77", cfg$micron_factor),
  sprintf("feature_version: %s", cfg$feature_version),
  sprintf("use_cell_counts: %s", cfg$use_cell_counts),
  sprintf("biovolume_units: um3 (per-ROI/per-chain; divide by cell_count for per-cell)"),
  sprintf("morpho_cols: %s", paste(cfg$morpho_cols, collapse = ", ")),
  sprintf("years_done: %s", paste(sapply(results, function(r) r$roi_meta$year), collapse = ", ")),
  sprintf("classifiers (VERIFY = production, not dev): %s",
          paste(sapply(results, function(r) sprintf("%s -> %s [%s]",
                r$roi_meta$year, r$roi_meta$model_name, r$roi_meta$class_dir)), collapse = "; ")),
  sprintf("blacklist: %s (md5 %s, modified %s, %d usable entries) -> %d samples excluded, listed in excluded_samples.txt",
          normalizePath(cfg$blacklist_file), unname(tools::md5sum(cfg$blacklist_file)),
          format(file.mtime(cfg$blacklist_file), "%Y-%m-%d %H:%M"), nrow(blacklist), nrow(excl_all)),
  sprintf("note: sample volume (ml analyzed) NOT here -> join data.txt SMVOL by SMPNO")
)
writeLines(meta, file.path(cfg$out_dir, "extraction_metadata.txt"))

cat("\n================ EXTRACTION DONE ================\n")
print(rbindlist(lapply(results, `[[`, "roi_meta"), fill = TRUE))
cat("\nper-ROI tables + sample_summary.txt + extraction_metadata.txt in", cfg$out_dir, "\n")
cat("join key to delivered data: sample == SMPNO ; volume from data.txt SMVOL\n")
cat("================================================\n")
