# delete_duplicate_images.R
#
# For each duplicated image listed in duplicate_files.csv, this script:
#   1. Checks which subfolder the image lives in under the ClassiPyR_2026-03-11 reference folder
#   2. Deletes the copy in SMHI-NIVA-SYKE-SAMS-SZN that does NOT match that subfolder
#
# Run in dry-run mode first (dry_run = TRUE) to preview deletions before committing.

library(dplyr)
library(readr)
library(purrr)

# ── Configuration ──────────────────────────────────────────────────────────────

csv_path <- "Z:/data/training_libraries/V6/duplicate_files.csv"   # path to your duplicate_files.csv

reference_root <- "Z:/data/png_images/ClassiPyR_2026-03-12"
target_root    <- "Z:/data/training_libraries/V6/SMHI-NIVA-SYKE-SAMS-SZN2"

dry_run <- TRUE   # Set to FALSE to actually delete files

# ── Load and reshape duplicates ────────────────────────────────────────────────

dupes <- read_csv(csv_path, show_col_types = FALSE)

# Each filename appears twice with two candidate dest_class values
# Group by filename so we have both candidates together
dupes_grouped <- dupes |>
  group_by(filename) |>
  summarise(candidates = list(dest_class), .groups = "drop") |>
  filter(lengths(candidates) == 2)   # only handle pairs (skip if >2 somehow)

cat(sprintf("Found %d duplicated filenames to process.\n\n", nrow(dupes_grouped)))

# ── Resolve correct class via ClassiPyR reference folder ──────────────────────

results <- dupes_grouped |>
  rowwise() |>
  mutate(
    # Check which candidate subfolder contains this file in the reference tree
    found_in = {
      hits <- candidates[
        file.exists(file.path(reference_root, candidates, filename))
      ]
      if (length(hits) == 1) hits[[1]]
      else if (length(hits) == 0) NA_character_
      else paste(hits, collapse = "|")   # both found — flag it
    },
    # The incorrect class is the other one
    delete_class = {
      if (is.na(found_in) || grepl("|", found_in, fixed = TRUE)) NA_character_
      else setdiff(candidates, found_in)[[1]]
    },
    delete_path = if_else(
      !is.na(delete_class),
      file.path(target_root, delete_class, filename),
      NA_character_
    )
  ) |>
  ungroup()

# ── Report unresolved cases ────────────────────────────────────────────────────

not_in_reference <- results |> filter(is.na(found_in))
ambiguous        <- results |> filter(grepl("|", found_in %||% "", fixed = TRUE))

if (nrow(not_in_reference) > 0) {
  cat(sprintf(
    "WARNING: %d file(s) not found in either candidate subfolder of the reference tree.\n",
    nrow(not_in_reference)
  ))
  cat("  These will be skipped:\n")
  walk(not_in_reference$filename, ~ cat("  -", .x, "\n"))
  cat("\n")
}

if (nrow(ambiguous) > 0) {
  cat(sprintf(
    "WARNING: %d file(s) found in BOTH candidate subfolders of the reference tree.\n",
    nrow(ambiguous)
  ))
  cat("  These will be skipped (manual review needed):\n")
  walk(ambiguous$filename, ~ cat("  -", .x, "\n"))
  cat("\n")
}

# ── Delete (or preview) ────────────────────────────────────────────────────────

to_delete <- results |>
  filter(!is.na(delete_path)) |>
  mutate(target_exists = file.exists(delete_path))

missing_targets <- to_delete |> filter(!target_exists)
if (nrow(missing_targets) > 0) {
  cat(sprintf(
    "NOTE: %d file(s) already absent from target (nothing to delete):\n",
    nrow(missing_targets)
  ))
  walk(missing_targets$delete_path, ~ cat("  -", .x, "\n"))
  cat("\n")
}

actionable <- to_delete |> filter(target_exists)

if (dry_run) {
  cat(sprintf("DRY RUN — %d file(s) would be deleted:\n", nrow(actionable)))
  walk(actionable$delete_path, ~ cat("  DELETE:", .x, "\n"))
  cat("\nSet dry_run <- FALSE to perform actual deletion.\n")
} else {
  cat(sprintf("Deleting %d file(s)...\n", nrow(actionable)))
  deleted  <- 0
  failed   <- 0
  
  for (path in actionable$delete_path) {
    ok <- tryCatch({
      file.remove(path)
    }, error = function(e) FALSE)
    
    if (isTRUE(ok)) {
      cat("  DELETED:", path, "\n")
      deleted <- deleted + 1
    } else {
      cat("  FAILED: ", path, "\n")
      failed <- failed + 1
    }
  }
  
  cat(sprintf("\nDone. Deleted: %d  |  Failed: %d\n", deleted, failed))
}
