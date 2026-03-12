#!/usr/bin/env Rscript
#
# Create an index of files in the merged training dataset
# (V6/SMHI-NIVA-SYKE-SAMS-SZN/) mapping each image back to its
# original source dataset (SMHI / NIVA / SAMS / SYKE / SZN).
#
# Output: V6/file_index.csv

library(data.table)

base_dir <- file.path(getwd(), "V6")
merged_dir <- file.path(base_dir, "SMHI-NIVA-SYKE-SAMS-SZN")
original_dir <- file.path(base_dir, "original_data")
sources <- c("NIVA", "SAMS", "SYKE", "SZN")

# --- 1. Index all files in original_data sources ---
cat("Indexing original_data sources...\n")
source_index <- rbindlist(lapply(sources, function(src) {
  src_dir <- file.path(original_dir, src)
  files <- list.files(src_dir, pattern = "\\.png$", recursive = TRUE,
                      full.names = FALSE)
  if (length(files) == 0) return(data.table())
  data.table(
    filename = basename(files),
    source_owner = src,
    source_class = dirname(files),
    source_relpath = files
  )
}))

# Some sources have nested subdirectories (e.g. SYKE had
# Syke-plankton_IFCB_2025_02/ClassName). Use only the final directory
# component as the class name.
source_index[, source_class := basename(source_class)]

cat(sprintf("  Found %d files across %d sources\n",
            nrow(source_index), length(sources)))

# --- 2. List all png files in the merged directory ---
cat("Listing merged directory...\n")
merged_files <- list.files(merged_dir, pattern = "\\.png$",
                           recursive = TRUE, full.names = FALSE)
merged_dt <- data.table(
  filename = basename(merged_files),
  dest_class = dirname(merged_files)
)
cat(sprintf("  Found %d files in merged directory\n", nrow(merged_dt)))

# --- 3. Match merged files to original sources ---
cat("Matching files to sources...\n")

# For files appearing in multiple source classes within the same source,
# keep all matches so the provenance is fully documented.
index <- merge(merged_dt, source_index, by = "filename", all.x = TRUE)

# Files not found in any original source are SMHI (IFCB110 / IFCB134)
index[is.na(source_owner), source_owner := "SMHI"]

# --- 4. Summary ---
cat("\nSummary of matched files per source:\n")
summary_dt <- index[, .(
  n_files = .N,
  n_unique_files = uniqueN(filename)
), by = source_owner]
print(summary_dt)

# Check for files with multiple source matches
multi_match <- index[, .N, by = filename][N > 1]
if (nrow(multi_match) > 0) {
  cat(sprintf("\nNote: %d filenames matched multiple source classes.\n",
              nrow(multi_match)))
  cat("These appear as multiple rows in the index.\n")
}

# Warn about unmatched non-SMHI files (files that should have been
# found in original_data but weren't)
smhi_no_ifcb <- index[source_owner == "SMHI" &
                       !grepl("IFCB110|IFCB134", filename)]
if (nrow(smhi_no_ifcb) > 0) {
  cat(sprintf("\nWarning: %d files classified as SMHI do not contain ",
              nrow(smhi_no_ifcb)))
  cat("IFCB110 or IFCB134 in their filename.\n")
  cat("First few:\n")
  print(head(smhi_no_ifcb[, .(filename, dest_class)], 10))
}

# --- 5. Check for duplicates in dest (same file in multiple classes) ---
dest_dupes <- merged_dt[, .N, by = filename][N > 1]
if (nrow(dest_dupes) > 0) {
  cat(sprintf("\nWARNING: %d filenames appear in multiple destination classes!\n",
              nrow(dest_dupes)))
  dupes_detail <- merged_dt[filename %in% dest_dupes$filename
                            ][order(filename)]
  dupes_file <- file.path(base_dir, "duplicate_files.csv")
  fwrite(dupes_detail, dupes_file)
  cat(sprintf("  Written to %s\n", dupes_file))
  cat("  First examples:\n")
  print(head(dupes_detail, 20))
} else {
  cat("\nNo duplicate filenames across destination classes.\n")
}

# --- 6. Class name translation table ---
# For each source owner, show which source_class mapped to which dest_class
# and how many files were involved.
cat("\nCreating class name translation table...\n")
class_map <- index[, .(n_files = .N),
                   by = .(source_owner, source_class, dest_class)]
setorder(class_map, source_owner, source_class, dest_class)
class_map_file <- file.path(base_dir, "class_name_translation.csv")
fwrite(class_map, class_map_file)
cat(sprintf("  Written to %s (%d mappings)\n", class_map_file, nrow(class_map)))

# Print a per-owner summary
for (owner in sort(unique(class_map$source_owner))) {
  owner_map <- class_map[source_owner == owner]
  cat(sprintf("\n  %s: %d source classes -> %d dest classes (%d files)\n",
              owner,
              uniqueN(owner_map$source_class),
              uniqueN(owner_map$dest_class),
              sum(owner_map$n_files)))
}

# --- 7. Write output ---
output_file <- file.path(base_dir, "file_index.csv")
fwrite(index, output_file)
cat(sprintf("\nIndex written to %s (%d rows)\n", output_file, nrow(index)))
