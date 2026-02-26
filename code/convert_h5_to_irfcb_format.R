#!/usr/bin/env Rscript
#
# Convert GPU-format classification H5 files to Gradio format.
#
# Differences handled:
#   classifierName (scalar)          -> classifier_name (1D)
#   class_labels_auto                -> class_name_auto
#   class_labels_above_threshold     -> class_name
#   class labels with _NNN suffix    -> suffix stripped (underscores kept)
#
# Usage:
#   Rscript convert_gpu_to_gradio_h5.R /path/to/folder
#
# Overwrites each file in place. Back up first if needed.

library(hdf5r)

# Set the folder containing GPU-format H5 files
folder <- "/path/to/classified"

if (!dir.exists(folder)) {
  stop("Folder not found: ", folder)
}

h5_files <- list.files(folder, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) == 0) {
  stop("No .h5 files found in: ", folder)
}

strip_suffix <- function(x) {
  gsub("_\\d{3}$", "", x)
}

for (h5_path in h5_files) {
  message("Converting: ", basename(h5_path))
  
  h5 <- H5File$new(h5_path, mode = "r")
  datasets <- h5$names
  
  # Check if this is GPU format (has classifierName)
  if (!"classifierName" %in% datasets) {
    message("  Skipping (not GPU format)")
    h5$close_all()
    next
  }
  
  # Read all data
  classifier_name <- h5[["classifierName"]]$read()
  class_labels <- strip_suffix(h5[["class_labels"]]$read())
  class_labels_auto <- strip_suffix(h5[["class_labels_auto"]]$read())
  class_labels_above_threshold <- strip_suffix(h5[["class_labels_above_threshold"]]$read())
  output_scores <- h5[["output_scores"]]$read()
  roi_numbers <- h5[["roi_numbers"]]$read()
  thresholds <- h5[["thresholds"]]$read()
  h5$close_all()
  
  # Write in Gradio format
  h5 <- H5File$new(h5_path, mode = "w")
  h5[["classifier_name"]] <- classifier_name
  h5[["class_labels"]] <- class_labels
  h5[["class_name_auto"]] <- class_labels_auto
  h5[["class_name"]] <- class_labels_above_threshold
  h5[["output_scores"]] <- output_scores
  h5[["roi_numbers"]] <- roi_numbers
  h5[["thresholds"]] <- thresholds
  h5$close_all()
  
  message("  Done (", length(roi_numbers), " ROIs, ", length(class_labels), " classes)")
}

message("\nConverted ", length(h5_files), " file(s).")
