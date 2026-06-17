library(tidyverse)

zip_files <- list.files(path = "Z:/data/whoi_plankton/png_images/zip-archives", pattern = "\\.zip$", full.names = TRUE)
png_path <- file.path(ifcb_path, "whoi_plankton", "png_images", "extracted_images")

ifcb_path <- Sys.getenv("ifcb_path")

# Destination folder
dest_dir <- file.path(ifcb_path, "whoi_plankton", "png_images", "extracted_images")

# Loop through each zip file and extract
for (zip_file in zip_files) {

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Extract the zip file
  unzip(zip_file, exdir = dest_dir)
  
  message("Extracted: ", zip_file, " -> ", dest_dir)
}

png_files <- list.files(path = png_path, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)

png_df <- data.frame(year = basename(dirname(dirname(png_files))),
                     class = basename(dirname(png_files)),
                     image = basename(png_files))

write_tsv(png_df, file.path(png_path, "images.txt"), na = "")
