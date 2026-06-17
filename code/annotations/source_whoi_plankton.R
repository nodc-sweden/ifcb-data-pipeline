library(iRfcb)

ifcb_py_install("venv")

# Define years
years <- 2006:2014

Sys.time()

# Define paths
ifcb_path <- Sys.getenv("ifcb_path")
raw_folder <- file.path(ifcb_path, "whoi_plankton", "data")
png_folder <- file.path(ifcb_path, "whoi_plankton", "png_images", "extracted_images")
png_extraction_folder <- file.path(ifcb_path, "whoi_plankton", "png_images", "reextracted_images")
class2use_file <- file.path(ifcb_path, "whoi_plankton", "config", "class2use_whoi.mat")
manual_folder <- file.path(ifcb_path, "whoi_plankton", "manual")
blobs_folder <- file.path(ifcb_path, "whoi_plankton", "blobs")

# Download WHOI-Plankton for all years
ifcb_prepare_whoi_plankton(years, 
                           png_folder, 
                           raw_folder, 
                           manual_folder, 
                           class2use_file, 
                           download_blobs = TRUE, 
                           extract_images = FALSE, 
                           blobs_folder = blobs_folder)

Sys.time()

# Extract images for verification
ifcb_extract_annotated_images(manual_folder, 
                              class2use_file, 
                              raw_folder, 
                              png_extraction_folder, 
                              skip_class = c("unclassified", "mix"))

Sys.time()

