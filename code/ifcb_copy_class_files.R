# This script merges .mat class files from the Baltic and Skagerrak-Kattegat classifier into a merged folder that is used by the IFCB dashboard
# It also counts the annotations and stores the results at havgem

# -------------------------------
# Load libraries
# -------------------------------

suppressPackageStartupMessages({
  library(readr, quietly = TRUE)
  library(fs, quietly = TRUE)
  library(iRfcb, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
})

# -------------------------------
# Configuration
# -------------------------------

year <- as.numeric(format(Sys.Date(), "%Y"))
ifcb_path <- Sys.getenv("ifcb_path")
havgem_path <- Sys.getenv("havgem_path")
ferrybox_path <- Sys.getenv("ferrybox_path")
raw_folder <- file.path(ifcb_path, "data", year)
class_score_version <- 1

# Define path for the generic class folders
class_folder_sk <- file.path(ifcb_path, "classified", "Skagerrak-Kattegat")
class_folder_baltic <- file.path(ifcb_path, "classified", "Baltic")

# List all class versions
class_subfolders_sk <- list.dirs(class_folder_sk, full.names = TRUE, recursive = FALSE)

# Find all subfolders from the selected year
class_subfolders_sk <- class_subfolders_sk[grepl(as.character(year), class_subfolders_sk)]

# Extract the version numbers from the subfolder names
versions_sk <- sub(paste0("class", year , "_", "Skagerrak-Kattegat", "_v"), "", basename(class_subfolders_sk))
versions_sk <- as.numeric(versions_sk)

# Get the subfolder with the highest version number
max_version_sk <- max(versions_sk, na.rm = TRUE)
class_folder_sk <- class_subfolders_sk[which(versions_sk == max_version_sk)]

# List all class versions
class_subfolders_baltic <- list.dirs(class_folder_baltic, full.names = TRUE, recursive = FALSE)

# Find all subfolders from the selected year
class_subfolders_baltic <- class_subfolders_baltic[grepl(as.character(year), class_subfolders_baltic)]

# Extract the version numbers from the subfolder names
versions_baltic <- sub(paste0("class", year , "_", "Baltic", "_v"), "", basename(class_subfolders_baltic))
versions_baltic <- as.numeric(versions_baltic)

# Get the subfolder with the highest version number
max_version_baltic <- max(versions_baltic, na.rm = TRUE)
class_folder_baltic <- class_subfolders_baltic[which(versions_baltic == max_version_baltic)]

# Define the path for the combined classification folder for the specific year and class score version
class_folder <- file.path(ifcb_path, 
                          "classified", 
                          "Baltic-Skagerrak-Kattegat-dashboard", 
                          paste0("class", year, "_v", class_score_version))

# -------------------------------
# Identify files
# -------------------------------

# List already merged class files
combined_class_files <- list.files(class_folder, ".mat")

# Download sample metadata
metadata <- ifcb_download_dashboard_metadata("https://ifcb-dashboard-utv.smhi.se/", "RV_Svea")

# Skip samples with missing position
hdr_data <- metadata %>%
  bind_cols(ifcb_convert_filenames(metadata$pid)) %>%
  filter(!latitude == -999) %>%
  filter(!longitude == -999) %>%
  filter(year == .env$year)

# -------------------------------
# Copy files
# -------------------------------

if (nrow(hdr_data) > 0) {
  
  # Determine if positions are in the Baltic basin
  hdr_data <- hdr_data %>%
    mutate(in_baltic = ifcb_is_in_basin(
      latitude,
      longitude))
  
  # To be removed
  available_bins <- gsub(paste0("_class_v", class_score_version, ".mat"), "", combined_class_files)
  bins_to_remove <- available_bins[!available_bins %in% hdr_data$pid]
  files_to_remove <- paste0(bins_to_remove, "_class_v", class_score_version, ".mat")
  
  # Delete files no longer relevant
  unlink(file.path(class_folder, files_to_remove))
  
  # Filter samples in the Baltic Sea
  baltic_samples <- hdr_data %>%
    filter(in_baltic)
  
  # Filter samples outside the Baltic Sea (e.g., Skagerrak-Kattegat)
  sk_samples <- hdr_data %>%
    filter(!in_baltic)
  
  # Retrieve class files for Skagerrak-Kattegat for the current iteration
  sk_class_files <- data.frame(files = list.files(class_folder_sk, 
                                                  pattern = ".mat", 
                                                  full.names = TRUE)) %>%
    mutate(pid = gsub(paste0("_class_v", class_score_version, ".mat"), "", basename(files)))
  
  # Retrieve class files for Baltic for the current iteration
  baltic_class_files <- data.frame(files = list.files(class_folder_baltic, 
                                                      pattern = ".mat",
                                                      full.names = TRUE)) %>%
    mutate(pid = gsub(paste0("_class_v", class_score_version, ".mat"), "", basename(files)))
  
  # Filter Skagerrak-Kattegat class files to include only those matching pids in sk_samples
  sk_files_to_copy <- sk_class_files %>%
    filter(pid %in% sk_samples$sample)
  
  # Filter Baltic class files to include only those matching pids in baltic_samples
  baltic_files_to_copy <- baltic_class_files %>%
    filter(pid %in% baltic_samples$sample)
  
  # Create the directory if it doesn’t exist, allowing for nested folder creation
  if (!dir.exists(class_folder)) {
    dir.create(class_folder, recursive = TRUE)
  }
  
  # Copy selected Skagerrak-Kattegat files to the combined classification folder
  copy_sk <- file.copy(sk_files_to_copy$files,
                       class_folder,
                       overwrite = TRUE)
  
  cat("Copied", sum(copy_sk), "Skagerrak-Kattegat files to folder", class_folder, "\n")
  
  
  # Copy selected Baltic files to the combined classification folder
  copy_baltic <- file.copy(baltic_files_to_copy$files,
                           class_folder,
                           overwrite = TRUE)
  
  cat("Copied", sum(copy_baltic), "Baltic files to folder", class_folder, "\n")
}

# -------------------------------
# Count manually annotated images
# -------------------------------

# Define the different classifiers
classifiers <- c("Baltic", "Skagerrak-Kattegat")

# Count the manual annotations and save as .txt
for (classifier in classifiers) {
  class_counts <- ifcb_count_mat_annotations(file.path(ifcb_path, "manual", classifier),
                                             file.path(ifcb_path, "config", paste0("class2use_", classifier, ".mat")),
                                             skip_class = 1)
  
  class_counts <- arrange(class_counts, desc(n))
  
  out_path <- file.path(havgem_path, "annotation",  year(Sys.Date()), Sys.Date())
  
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }
  
  write_tsv(class_counts,
            file.path(out_path, paste0("class_counts_", classifier, "_", Sys.Date(), ".txt")),
            na = "")
  
}
