# This script copies new files to the work directory, extract v4 blobs and features, and runs PSD quality control
# The script can be scheduled to run on a regular basis

# Do not use the default venv for this script
Sys.setenv(USE_IRFCB_PYTHON = FALSE)

# -------------------------------
# Load libraries
# -------------------------------

library(readr)
library(fs)
library(reticulate)
library(iRfcb)
library(tools)
library(dplyr)
library(stringr)

# -------------------------------
# Configuration
# -------------------------------

year <- as.numeric(format(Sys.Date(), "%Y"))
ifcb_path <- Sys.getenv("ifcb_path")
ifcb_base_path <- Sys.getenv("ifcb_base_path")
repo_dir <- "C:/R/ifcb-data-pipeline/"

feature_folder <- file.path(ifcb_path, "features", "v4", year)
blobs_folder <- file.path(ifcb_path, "blobs", "v4", year)
delivery_dir <- file.path(ifcb_base_path, "delivery", year)
raw_folder <- file.path(ifcb_path, "data", year)
psd_plot_folder <- file.path(ifcb_path, "psd", "figures", year)

blacklist_file <- file.path(repo_dir, "data", "sample_blacklist.tsv")
log_file <- file.path(repo_dir, "large_bins_log.txt")
psd_data_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_data.csv"))
psd_fits_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_fits.csv"))
psd_flags_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_flags.csv"))

# File size threshold: 1 GB
size_threshold <- 1073741824  # bytes

# -------------------------------
# Copy new delivery files
# -------------------------------

# Read sample blacklist
blacklist <- read_tsv(blacklist_file, show_col_types = FALSE, progress = FALSE)$sample

# Convert blacklist patterns to a single regex for efficient matching
blacklist_regex <- paste(blacklist, collapse = "|")

# Get list of folders in source
source_folders <- dir_ls(delivery_dir, type = "directory")

# Process each folder
for (src_folder in source_folders) {

  folder_name <- path_file(src_folder)
  dest_folder <- path(raw_folder, folder_name)

  # Skip folder if already exists in destination
  if (dir_exists(dest_folder)) {
    message("Skipping existing folder: ", folder_name)
    next
  }

  # List all files recursively inside the source folder
  files <- dir_ls(src_folder, recurse = TRUE, type = "file")

  # --------------------------------------------------------
  # Filter out blacklisted files
  # --------------------------------------------------------
  if (!is.null(blacklist_regex)) {
    files <- files[!grepl(blacklist_regex, files)]
  }

  # --------------------------------------------------------
  # Detect large files and log them
  # --------------------------------------------------------
  sizes <- file_info(files)$size
  large_files_idx <- which(sizes > size_threshold)

  if (length(large_files_idx) > 0) {

    large_files <- files[large_files_idx]
    # Extract bin names (file without extension)
    large_bins <- file_path_sans_ext(basename(large_files))

    log_entry <- paste(
      Sys.time(),
      "- Skipped bins:",
      paste(large_bins, collapse = ", ")
    )

    cat(log_entry, file = log_file, sep = "\n", append = TRUE)

    warning(sprintf(
      "Large files detected (>1 GB). These were skipped. See log file for details: %s",
      log_file
    ))

    # Remove large files from list of files to copy
    files <- files[-large_files_idx]
  }

  if (length(files) > 0) {

    message("Copying new folder: ", folder_name)
    dir_create(dest_folder)

    # --------------------------------------------------------
    # Copy remaining files preserving folder structure
    # --------------------------------------------------------
    for (file in files) {
      rel_path <- path_rel(file, start = src_folder)
      dest_path <- path(dest_folder, rel_path)

      dir_create(path_dir(dest_path))
      file_copy(file, dest_path, overwrite = TRUE)
    }
  }
}

# -------------------------------
# v4 Feature extraction
# -------------------------------

# Use venv
use_virtualenv(file.path(repo_dir, "venv"))

# Add root repo to sys.path
py_run_string("import sys; sys.path.append('C:/R/ifcb-data-pipeline/code/python/ifcb-features')")

# Import function to extract all features
extract_slim_features <- import("extract_slim_features")

# List already extracted feature files
feature_files <- list.files(feature_folder, ".csv")
feature_bins <- sub("^([^_]*_[^_]*)_.*$", "\\1", feature_files)

# List roi files
roi_files <- list.files(raw_folder, pattern = "\\.roi$", recursive = TRUE, full.names = TRUE)
roi_bins <- tools::file_path_sans_ext(basename(roi_files))

# Skip already processed bins
bins_to_process <- roi_bins[!roi_bins %in% feature_bins]

# Only extract features if there a any new files to process
if (length(bins_to_process) > 0) {
  # Extract v4 features
  extract_slim_features$extract_and_save_all_features(as.character(raw_folder),
                                                      as.character(feature_folder),
                                                      as.character(blobs_folder),
                                                      bins = as.list(bins_to_process))
}

# -------------------------------
# PSD QC
# -------------------------------

# Read stored PSD results

if (file.exists(psd_fits_file)) {
  psd_fits <- read_csv(psd_fits_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique_quiet")

  psd_fits <- psd_fits %>%
    rename(sample = dplyr::any_of("...1"))
  
  processed_bins <- as.character(psd_fits$sample)
} else {
  psd_fits <- data.frame()
  processed_bins <- as.character()
}

if (file.exists(psd_data_file)) {
  psd_data <- read_csv(psd_data_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique_quiet", locale = locale(encoding = "UTF-8"))
  
  psd_data <- psd_data %>%
    rename(sample = dplyr::any_of("...1"))
  
  colnames(psd_data) <- gsub("\u03BC", "u", colnames(psd_data))
} else {
  psd_data <- data.frame()
}

if (file.exists(psd_flags_file)) {
  psd_flags <- read_csv(psd_flags_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique_quiet")
} else {
  psd_flags <- data.frame()
}

# Keep only bins whose timestamp is NOT in processed_bins
bins_to_psd <- roi_bins[!roi_bins %in% processed_bins]

# Calculate the particle size distribution (PSD) using IFCB data from the specified folders
psd <- ifcb_psd(feature_folder = feature_folder,
                hdr_folder = raw_folder,
                bins = bins_to_psd,
                save_data = FALSE,
                output_file = NULL,
                plot_folder = psd_plot_folder,
                use_marker = FALSE,
                start_fit = 15,
                r_sqr = 0.5,
                beads = 10 ** 20,
                bubbles = 110,
                incomplete = c(1500, 3),
                missing_cells = 0.5,
                biomass = 3000,
                bloom = 10,
                humidity = 75,
                micron_factor = 1/2.77,
                fea_v = 4,
                use_plot_subfolders = FALSE)

# Bind data
psd_data <- bind_rows(psd_data, psd$data)
psd_fits <- bind_rows(psd_fits, psd$fits)
psd_flags <- bind_rows(psd_flags, psd$flags)

# Store PSD results
write_csv(psd_flags, psd_flags_file, progress = FALSE)
write_csv(psd_fits, psd_fits_file, progress = FALSE)
write_csv(psd_data, psd_data_file, progress = FALSE)

# -------------------------------
# Log successful completion
# -------------------------------

completion_log <- file.path(repo_dir, "pipeline_run_log.txt")
log_entry <- paste(Sys.time(), "- Pipeline completed successfully")
cat(log_entry, file = completion_log, sep = "\n", append = TRUE)
