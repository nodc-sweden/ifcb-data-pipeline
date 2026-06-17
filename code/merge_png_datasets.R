library(tidyverse)
library(fs)
library(stringr)

ifcb_path <- Sys.getenv("ifcb_path")

# Input folders
baltic_dir <- file.path(ifcb_path, "png_images/Baltic/2025-11-04")
westcoasts_dir <- file.path(ifcb_path, "png_images/Skagerrak-Kattegat/2025-11-04")
tangesund_dir <- file.path(ifcb_path, "png_images/Tangesund/2025-11-04")
syke_uto <- file.path(ifcb_path, "png_images/SYKE/phytoplankton_Uto_2021_labeled")
syke_labeled <- file.path(ifcb_path, "png_images/SYKE/labeled_20201020")
niva <- file.path(ifcb_path, "ifcb139/Classifier_rev/Classifier_rev")
med <- file.path(ifcb_path, "png_images/MedPlanktonSet/IFCB_images/IFCB_images_selected")
sams <- file.path(ifcb_path, "sams/png_images/extracted_images_selected")

# Output folder (adjust if you want a different destination)
merged_dir <- file.path(ifcb_path, "png_images/CNN_merged/2025-11-04")
dir_create(merged_dir)

# Function to copy and remove trailing numbers in folder names
process_and_copy <- function(source_dir, target_dir) {
  folders <- dir_ls(source_dir, type = "directory")
  
  for (folder in folders) {
    folder_name <- path_file(folder)
    
    # Remove trailing underscore + digits, e.g. "_096", "_12"
    clean_name <- str_replace(folder_name, "_\\d+$", "")
    
    # Create target folder
    target_folder <- path(target_dir, clean_name)
    dir_create(target_folder)
    
    # Copy files into merged folder
    file_copy(dir_ls(folder, glob = "*.png"), target_folder, overwrite = TRUE)
  }
}

# Apply to both folders
process_and_copy(baltic_dir, merged_dir)
process_and_copy(westcoasts_dir, merged_dir)
process_and_copy(tangesund_dir, merged_dir)
process_and_copy(syke_uto, merged_dir)
process_and_copy(syke_labeled, merged_dir)
process_and_copy(niva, merged_dir)
process_and_copy(med, merged_dir)
process_and_copy(sams, merged_dir)


# List all subdirectories in merged_dir
subdirs <- list.dirs(merged_dir, full.names = TRUE, recursive = FALSE)

# Define valid image file extensions
image_extensions <- c("jpg", "jpeg", "png", "tif", "tiff", "bmp", "gif")

# Count images in each subdirectory
image_counts <- lapply(subdirs, function(subdir) {
  # List all files in the subdirectory
  files <- list.files(subdir, full.names = TRUE)
  
  # Filter for image files
  image_files <- files[tolower(tools::file_ext(files)) %in% image_extensions]
  
  # Return count
  length(image_files)
}) %>%
  setNames(basename(subdirs)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subdirectory") %>%
  dplyr::rename(image_count = V1)

# Print result
print(image_counts)
