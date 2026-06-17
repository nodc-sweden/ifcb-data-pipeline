#!/usr/bin/env Rscript
#
# Reproduce the merged training dataset from the file index.
# Reads V6/file_index.csv and copies images from their source
# locations into the destination directory structure.
#
# Usage:
#   Rscript 02_copy_from_index.R [dest_dir]
#
# If dest_dir is not provided, defaults to V6/SMHI-NIVA-SYKE-SAMS-SZN_reproduced/
#
# SMHI files (IFCB110/IFCB134) are expected to already exist in
# the destination directory or a separate SMHI source directory,
# since they are not part of original_data.

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- file.path(getwd(), "V6")
original_dir <- file.path(base_dir, "original_data")

dest_dir <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(base_dir, "SMHI-NIVA-SYKE-SAMS-SZN2")
}

# Optional: path to SMHI source images (set to NA to skip SMHI files)
smhi_source_dir <- if (length(args) >= 2) args[2] else NA

index_file <- file.path(base_dir, "file_index.csv")
if (!file.exists(index_file)) {
  stop("file_index.csv not found. Run 01_create_file_index.R first.")
}

cat(sprintf("Reading index from %s\n", index_file))
index <- fread(index_file)

# For files with multiple source class matches, use the first match only
index <- index[, .SD[1], by = .(filename, dest_class)]

cat(sprintf("  %d files to copy\n", nrow(index)))
cat(sprintf("  Destination: %s\n", dest_dir))

# --- Build source paths ---
index[source_owner %in% c("NIVA", "SAMS", "SYKE", "SZN"),
      source_path := file.path(original_dir, source_owner, source_relpath)]

if (!is.na(smhi_source_dir)) {
  index[source_owner == "SMHI",
        source_path := file.path(smhi_source_dir, filename)]
}

# --- Copy files ---
external <- index[!is.na(source_path)]
skipped_smhi <- index[source_owner == "SMHI" & is.na(source_path)]

if (nrow(skipped_smhi) > 0) {
  cat(sprintf("\nSkipping %d SMHI files (no source directory provided).\n",
              nrow(skipped_smhi)))
  cat("  To include SMHI files, pass the SMHI source directory as\n")
  cat("  the second argument.\n")
}

cat(sprintf("\nCopying %d files...\n", nrow(external)))

# Create all destination directories
dest_classes <- unique(index$dest_class)
for (dc in dest_classes) {
  dir.create(file.path(dest_dir, dc), recursive = TRUE, showWarnings = FALSE)
}

n_copied <- 0
n_missing <- 0
errors <- character(0)

for (i in seq_len(nrow(external))) {
  src <- external$source_path[i]
  dst <- file.path(dest_dir, external$dest_class[i], external$filename[i])

  if (!file.exists(src)) {
    n_missing <- n_missing + 1
    if (n_missing <= 10) {
      errors <- c(errors, sprintf("  Missing: %s", src))
    }
    next
  }

  file.copy(src, dst, overwrite = FALSE)
  n_copied <- n_copied + 1

  if (n_copied %% 10000 == 0) {
    cat(sprintf("  %d / %d copied\n", n_copied, nrow(external)))
  }
}

cat(sprintf("\nDone. Copied %d files.\n", n_copied))
if (n_missing > 0) {
  cat(sprintf("Missing source files: %d\n", n_missing))
  cat(paste(errors, collapse = "\n"), "\n")
  if (n_missing > 10) {
    cat(sprintf("  ... and %d more\n", n_missing - 10))
  }
}
if (nrow(skipped_smhi) > 0) {
  cat(sprintf("Skipped SMHI files: %d\n", nrow(skipped_smhi)))
}
