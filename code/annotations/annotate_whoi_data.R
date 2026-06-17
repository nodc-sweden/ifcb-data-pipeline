library(iRfcb)
library(tidyverse)
# ifcb_py_install()

# Define year
year <- "2014"

# Define paths
ifcb_path <- Sys.getenv("ifcb_path")
class2use_file <- file.path(ifcb_path, "whoi_plankton", "config", "class2use_whoi.mat")
png_path <- file.path(ifcb_path, "whoi_plankton", "png_images", "extracted_images", year)
manual_folder <- file.path(ifcb_path, "whoi_plankton", "manual")
data_folder <- file.path(ifcb_path, "whoi_plankton", "data")

# List all png files
all_png_images <- list.files(png_path, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)

# Get sample info
sample_info <- ifcb_convert_filenames(basename(all_png_images))

# Extract info to build new filenames
sample_info <- sample_info %>%
  mutate(
    formatted_ifcb_number = sprintf("IFCB%03d", as.integer(gsub("IFCB", "", ifcb_number))),
    formatted_roi = sprintf("%05d", roi),  # Ensure roi is a five-digit number
    new_name = paste0("I", 
                      format(date, "%Y%m%d"), 
                      "T", 
                      format(as.POSIXct(time, format="%H:%M:%S"), "%H%M%S"), 
                      "_",
                      formatted_ifcb_number,
                      "_",
                      formatted_roi,
                      ".png")
  ) %>%
  select(-formatted_ifcb_number, -formatted_roi)

# Create a df with new names
rename_df <- data.frame(all_png_images, new_name = paste(dirname(all_png_images), sample_info$new_name, sep = "/"))

# Rename png files to match ifcb files
file.rename(rename_df$all_png_images, rename_df$new_name)

# List all classes from the current year
classes <- c("unclassified", basename(list.dirs(png_path, recursive = FALSE)))

# Create a new class2use file          
ifcb_create_class2use(classes, class2use_file)

# List the adc files
adcfiles <- list.files(data_folder, pattern = "adc$", full.names = TRUE, recursive = TRUE)

for (i in 2:length(classes)) {
  # Select class to annotate
  selected_class <- classes[i]
  
  # Find path to pngs
  class_png_path <- file.path(png_path, selected_class)
  
  # List png filenames
  png_images <- list.files(class_png_path, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)
  
  # Annotate batch
  ifcb_annotate_batch(
    png_images,
    selected_class,
    manual_folder,
    adcfiles,
    class2use_file,
    manual_output = NULL,
    manual_recursive = FALSE,
    unclassified_id = 1,
    do_compression = TRUE
  )
}

# Extract images to verify correct annotation
ifcb_extract_annotated_images(
  manual_folder,
  class2use_file,
  data_folder,
  file.path(ifcb_path, "whoi_plankton", "extracted_images"),
  skip_class = 1,
  old_adc = FALSE)
