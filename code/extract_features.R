### Superseded by iRfcb::ifcb_extract_features()



# Do not use the default venv for this script
Sys.setenv(USE_IRFCB_PYTHON = FALSE)

library(reticulate)
library(iRfcb)

year <- 2025
ifcb_path <- Sys.getenv("ifcb_path")

# Define paths to features and .roi files
feature_folder <- file.path(ifcb_path, "features", "v4", year)
raw_folder <- file.path(ifcb_path, "data", year)
blobs_folder <- file.path(ifcb_path, "blobs_v4", year)

# List required packages for pyifcb and ifcb_features
additional_packages <- c("pyyaml", "smbprotocol", "pysmb", "rectpack", "Pillow", "requests", 
                         "h5py", "scikit-image", "pyfftw", "scikit-learn")

# Install and use python virtual environment
ifcb_py_install("C:/R/ifcb-data-pipeline/venv", packages = additional_packages)

# Install custom phasepack dependency
py_install("git+https://github.com/WHOIGit/phasepack@v1.6.1", 
           envname = "./venv")

# Install pyifcb dependency
py_install("git+https://github.com/joefutrelle/pyifcb", 
           envname = "./venv", 
           pip_options = "--no-deps")

# Install ifcb-features
py_install("git+https://github.com/WHOIGit/ifcb-features", 
           envname = "./venv", 
           pip_options = "--no-deps")

# Add root repo to sys.path
py_run_string("import sys; sys.path.append('C:/R/ifcb-data-pipeline/src/python/ifcb-features')")

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

# Get file info including size (bytes)
roi_info <- file.info(roi_files)

# Find files larger than 1 GB (1 GB = 1,073,741,824 bytes)
large_files <- rownames(roi_info[roi_info$size > 1073741824, ])
large_bins <- tools::file_path_sans_ext(basename(large_files))

# Store large files in a logfile for later inspection
log_file <- file.path("data", "large_bins_log.txt")

if (length(large_bins) > 0) {
  # Log entry text
  log_entry <- paste(Sys.time(), "- Skipped bins:", paste(large_bins, collapse = ", "))
  
  # Append to log file
  cat(log_entry, file = log_file, sep = "\n", append = TRUE)
  
  # Warning message pointing to logfile
  warning(sprintf(
    "Large files detected; the script will not extract features for these bins. See log file for details: %s",
    log_file
  ))
}

# Skip large files (bubbles, most likely)
bins_to_process <- bins_to_process[!bins_to_process %in% large_bins]

# Only extract features if there a any new files to process
if (length(bins_to_process) > 0) {
  # Extract v4 features
  extract_slim_features$extract_and_save_all_features(as.character(raw_folder),
                                                      as.character(feature_folder),
                                                      as.character(blobs_folder),
                                                      bins = as.list(bins_to_process))
}
