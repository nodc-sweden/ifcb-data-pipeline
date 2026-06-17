# Load required library
library(dplyr)

# Define the main directory
merged_dir <- "Z:/data/png_images/CNN_merged/2025-11-11"

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

image_threshold <- 50 # Number of images to include in classifier

above_threshold <- image_counts %>%
  filter(image_count > image_threshold) %>%
  nrow()

# Print result
print(above_threshold)
