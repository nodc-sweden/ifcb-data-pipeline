library(iRfcb)
library(tidyverse)
library(zip)

# Edit the .Renviron by usethis::edit_r_environ("project")
ifcb_path <- Sys.getenv("ifcb_path")

# Define results folders
output_folder <- file.path(ifcb_path, "shark")

# Read sample blacklist
blacklist <- read_tsv("data/sample_blacklist.tsv", col_types = cols())

# Proceed with PSD calculation if cache doesn't exist or files have changed
psd_result_files <- file.path(ifcb_path, "psd")

# Read psd results
psd <- read_csv(file.path(psd_result_files, "psd_2024_flags.csv"), col_types = cols())

# Find samples with bubbles
bubbles <- psd %>%
  filter(flag == "Bubbles")

# Combine blacklisted and bubble samples
files_to_skip <- c(blacklist$sample, bubbles$file)

# List all roi files
roi_files <- list.files(file.path(ifcb_path, "data", "2024"), pattern = "D.*\\.roi", full.names = TRUE, recursive = TRUE)

# Extract only the relevant parts of the file paths (filenames)
roi_filenames <- basename(roi_files)

# Identify which files should be skipped
to_keep <- !roi_filenames %in% files_to_skip

# Filter out the files to skip
filtered_roi_files <- roi_files[to_keep]

# Define output directory
output_dir <- file.path(ifcb_path, "extracted_images", "2024")

# Extract each sample
for (file in filtered_roi_files) {
  tryCatch(
    {
      ifcb_extract_pngs(file, output_dir)
    },
    error = function(e) {
      message("Skipping file due to error: ", file, " - ", e$message)
    }
  )
}

png_files <- list.files(file.path(ifcb_path, "extracted_images", "2024"), pattern = "D.*\\.png", full.names = TRUE, recursive = TRUE)

# Extract object_id (remove .png extension)
object_id <- tools::file_path_sans_ext(basename(png_files))

# Extract relative path
relative_path <- gsub(paste0("^", file.path(ifcb_path, "extracted_images", "2024"), "/"), "", png_files)

# Create a data frame
index_df <- data.frame(object_id = object_id, path = relative_path, stringsAsFactors = FALSE)

# Write to CSV
write.csv(index_df, file.path(ifcb_path, "extracted_images", "2024", "index.csv"), row.names = FALSE, quote = FALSE)

# Define paths
extracted_images_path <- file.path(ifcb_path, "extracted_images", "2024")
zip_file <- file.path(ifcb_path, "extracted_images", "ifcb_extracted_2024.zip")
index_file <- file.path(extracted_images_path, "index.csv")

# Use existing list of images + index file
files_to_zip <- c(png_files, index_file)

# Create ZIP archive (keeping only relative paths)
zip::zip(zipfile = zip_file, files = files_to_zip, root = extracted_images_path, mode = "chdir") # Root did not work as expected

message("ZIP archive created successfully!")