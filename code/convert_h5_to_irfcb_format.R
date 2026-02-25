#' Convert legacy _class.h5 files to iRfcb ifcb_save_classification() format
#'
#' Renames HDF5 datasets from the old naming convention to the new one:
#'   classifierName            -> classifier_name
#'   class_labels_auto         -> class_name_auto
#'   class_labels_above_threshold -> class_name
#'
#' Datasets that already use the new names are left untouched.
#' A backup of each original file is created as *_class_backup.h5.

library(hdf5r)

convert_h5_to_irfcb <- function(input_dir,
                                output_dir = input_dir,
                                backup = TRUE,
                                verbose = TRUE) {

  h5_files <- list.files(input_dir,
                         pattern = "_class\\.h5$",
                         full.names = TRUE,
                         recursive = TRUE)

  if (length(h5_files) == 0) {
    message("No *_class.h5 files found in: ", input_dir)
    return(invisible(character(0)))
  }

  if (verbose) message("Found ", length(h5_files), " file(s) to convert.")

  # Field mapping: old name -> new name
  rename_map <- list(
    classifierName            = "classifier_name",
    class_labels_auto         = "class_name_auto",
    class_labels_above_threshold = "class_name"
  )

  converted <- character(0)

  for (h5_path in h5_files) {
    if (verbose) message("Processing: ", basename(h5_path))

    # Determine output path
    rel_path <- sub(paste0("^", normalizePath(input_dir, winslash = "/"), "/?"),
                    "", normalizePath(h5_path, winslash = "/"))
    out_path <- file.path(output_dir, rel_path)
    out_dir <- dirname(out_path)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Read all datasets from the source file
    src <- H5File$new(h5_path, mode = "r")
    on.exit(src$close_all(), add = TRUE)

    dataset_names <- names(src)

    # Check if already in new format
    already_new <- all(c("classifier_name", "class_name_auto", "class_name") %in% dataset_names)
    has_old <- any(names(rename_map) %in% dataset_names)

    if (already_new && !has_old) {
      if (verbose) message("  Already in iRfcb format, skipping.")
      src$close_all()
      next
    }

    # Read all datasets into memory
    data <- list()
    for (ds_name in dataset_names) {
      data[[ds_name]] <- src[[ds_name]]$read()
    }
    src$close_all()

    # Create backup if writing in-place
    if (backup && normalizePath(h5_path) == normalizePath(out_path)) {
      backup_path <- sub("_class\\.h5$", "_class_backup.h5", h5_path)
      file.copy(h5_path, backup_path, overwrite = TRUE)
      if (verbose) message("  Backup: ", basename(backup_path))
    }

    # Apply renames
    for (old_name in names(rename_map)) {
      new_name <- rename_map[[old_name]]
      if (old_name %in% names(data) && !(new_name %in% names(data))) {
        data[[new_name]] <- data[[old_name]]
        data[[old_name]] <- NULL
        if (verbose) message("  Renamed: ", old_name, " -> ", new_name)
      }
    }

    # Strip trailing _NNN suffixes from class label/name strings
    truncate <- iRfcb:::truncate_folder_name
    label_fields <- c("class_labels", "class_name_auto", "class_name")
    for (field in label_fields) {
      if (field %in% names(data) && is.character(data[[field]])) {
        data[[field]] <- vapply(data[[field]], truncate, character(1),
                                USE.NAMES = FALSE)
        if (verbose) message("  Truncated trailing _NNN in: ", field)
      }
    }

    # Write new file
    dst <- H5File$new(out_path, mode = "w")

    for (ds_name in names(data)) {
      dst[[ds_name]] <- data[[ds_name]]
    }

    dst$close_all()
    converted <- c(converted, out_path)
    if (verbose) message("  Saved: ", out_path)
  }

  if (verbose) message("Converted ", length(converted), " of ", length(h5_files), " file(s).")
  invisible(converted)
}

# --- Run conversion ---
# Edit the path below to point to your classified .h5 files
input_dir <- "classified/niva_smhi_baltic/V5/"

convert_h5_to_irfcb(input_dir)
