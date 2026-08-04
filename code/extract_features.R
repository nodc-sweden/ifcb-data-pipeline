# Do not use the default venv for this script
Sys.setenv(USE_IRFCB_PYTHON = FALSE)

# -------------------------------
# Load libraries
# -------------------------------

suppressPackageStartupMessages({
  library(iRfcb, quietly = TRUE)
  library(reticulate, quietly = TRUE)
})

# -------------------------------
# Configuration
# -------------------------------

year <- 2026
ifcb_path <- Sys.getenv("ifcb_path")
repo_dir <- "C:/R/ifcb-data-pipeline"

feature_folder <- file.path(ifcb_path, "features_1_1_1", "v4", year)
blobs_folder <- file.path(ifcb_path, "blobs", "v4", year)
raw_folder <- file.path(ifcb_path, "data", year)

# Use venv
use_virtualenv(file.path(repo_dir, ".virtualenvs", "features"))

# -------------------------------
# Update metadata
# -------------------------------

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
  
  cat("Extracting features...\n")
  
  # Extract v4 features
  ifcb_extract_features(data_folder = as.character(raw_folder),
                        features_folder = as.character(feature_folder),
                        blobs_folder = as.character(blobs_folder),
                        bins = bins_to_process,
                        parallel = TRUE,
                        n_cores = 3,
                        feature_tag = "fea")
}
