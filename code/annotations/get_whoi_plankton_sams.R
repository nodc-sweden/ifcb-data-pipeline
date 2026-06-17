#year <- "2014"
ifcb_path <- Sys.getenv("ifcb_path")
classes = NULL
#dashboard_url = "https://ifcb-data.whoi.edu/mvco/"
# dashboard_url = "https://ifcb-data.sams.ac.uk/Scalloway/"
dashboard_url = "https://ifcb-data.sams.ac.uk/data/"
#dashboard_url = "https://ifcb-data.sams.ac.uk/MI_RV_Tom_Crean/"

#png_path <- file.path(ifcb_path, "whoi_plankton", "png_images", year)
png_path <- file.path(ifcb_path, "sams", "png_images", "extracted_images")

#data_path <- file.path(ifcb_path, "whoi_plankton", "raw", year)
data_path <- file.path(ifcb_path, "sams", "raw")

convert_filenames <- FALSE

library(tidyverse)
library(curl)
library(iRfcb)
library(digest)

png_files <- list.files(path = png_path, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)

image_df <- data.frame(folder = folder_name <- basename(dirname(png_files)),
                       image = png_files)

if (!is.null(classes)) {
  selected_images <- dplyr::filter(image_df, folder %in% classes)
} else {
  selected_images <- image_df
}

# Input string
ifcb_string <- tools::file_path_sans_ext(basename(selected_images$image))

# # Extract components using regex
# ifcb_parts <- str_match(ifcb_string, "^(IFCB\\d+)_(\\d{4})_(\\d{3})_(\\d{6})_(\\d+)$")
# 
# # Convert day of year to date
# date_value <- as.Date(paste0(ifcb_parts[,3], "-01-01")) + as.integer(ifcb_parts[,4]) - 1
# 
# # Construct dataframe
# sample_df <- tibble(
#   sample = paste(ifcb_parts[,2], ifcb_parts[,3], ifcb_parts[,4], ifcb_parts[,5], sep = "_"),
#   timestamp = str_glue("{date_value} {substr(ifcb_parts[,5], 1, 2)}:{substr(ifcb_parts[,5], 3, 4)}:{substr(ifcb_parts[,5], 5, 6)}"),
#   date = date_value,
#   year = as.integer(ifcb_parts[,3]),
#   month = month(date_value),
#   day = day(date_value),
#   time = str_glue("{substr(ifcb_parts[,5], 1, 2)}:{substr(ifcb_parts[,5], 3, 4)}:{substr(ifcb_parts[,5], 5, 6)}"),
#   ifcb_number = ifcb_parts[,2],
#   roi = as.integer(ifcb_parts[,6])
# )

sample_df <- ifcb_convert_filenames(ifcb_string)

# URL of the file
urls <- unique(paste0(dashboard_url, sample_df$sample))

# File extensions to download
extensions <- c("zip")

download_ifcb_files <- function(base_url, extensions, dest_dir, convert_filenames = convert_filenames) {
  for (ext in extensions) {
    file_url <- paste0(base_url, ".", ext)  # Construct full file URL
    
    # Extract the IFCB string
    ifcb_string <- tools::file_path_sans_ext(basename(base_url))
    
    # Initialize date_object as NULL
    date_object <- NULL
    
    # Check if the string matches the first format (IFCB1_2014_188_222013)
    if (grepl("^IFCB\\d+_\\d{4}_\\d{3}_\\d{6}$", ifcb_string)) {

      # Extract components using regex
      ifcb_parts <- str_match(ifcb_string, "^(IFCB\\d+)_(\\d{4})_(\\d{3})_(\\d{6})$")
      
      # Convert day of year to date
      date_object <- format(as.Date(paste0(ifcb_parts[,3], "-01-01")) + as.integer(ifcb_parts[,4]) - 1, "D%Y%m%d")
    } else if (grepl("^D\\d{8}T\\d{6}_IFCB\\d+$", ifcb_string)) {

      # Extract components using regex
      ifcb_parts <- str_match(ifcb_string, "^D(\\d{8})T(\\d{6})_IFCB(\\d+)$")
      
      # Extract date (YYYYMMDD) from the match
      date_object <- paste0("D", ifcb_parts[,2])
      
    } else {
      message("Unknown format")
      date_object <- NULL
    }
    
    # Filename creation logic
    if (convert_filenames && !is.null(date_object)) {
      # Extract the file name part (before the extension)
      ifcb_string <- tools::file_path_sans_ext(basename(file_url))
      
      # Extract components using regex
      ifcb_parts <- str_match(ifcb_string, "^(IFCB\\d+)_(\\d{4})_(\\d{3})_(\\d{6})$")
      
      # Extract year, day of year, and time
      year <- ifcb_parts[,3]
      day_of_year <- as.integer(ifcb_parts[,4])
      time <- ifcb_parts[,5]
      ifcb_number <- sub("IFCB(\\d+)", sprintf("IFCB%03d", as.integer(sub("IFCB", "", ifcb_parts[,2]))), ifcb_parts[,2])
      
      # Convert day of year to date
      date_object <- paste0("D", format(as.Date(paste0(year, "-01-01")) + day_of_year - 1, "%Y%m%d"))
      
      # Construct the filename in the desired format: DYYYYMMDD_THHMMSS_IFCB<ifcb_number>_YYYY_DDD_HHMMSS.ext
      filename <- paste0(date_object, "T", substr(time, 1, 2), substr(time, 3, 4), substr(time, 5, 6), "_", ifcb_number, ".", ext)
    } else {
      filename <- basename(file_url)
    }
    
    # Set the destination file path
    if (!is.null(date_object)) {
      destfile <- file.path(dest_dir, date_object, filename)  # Use date_object as directory
    } else {
      destfile <- file.path(dest_dir, filename)  # If no date_object, just the filename
    }
    
    # Create the destination folder if it doesn't exist
    dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)
    
    # Skip download if the file already exists
    if (file.exists(destfile)) {
      message("File already exists, skipping: ", destfile)
      next
    }
    
    # Set up error handling with curl
    tryCatch({
      handle <- curl::new_handle(low_speed_time = 30, low_speed_limit = 1000)
      curl::handle_setopt(handle, failonerror = TRUE)
      
      # Check if file exists on the server by fetching the file (do not download yet)
      con <- curl::curl_fetch_memory(file_url, handle = handle)
      if (con$status_code == 200) {
        # File exists, so download
        curl::curl_download(file_url, destfile)
        message("Download successful: ", destfile)
      } else {
        stop("File not found on server.")
      }
    }, error = function(e) {
      warning("Download failed for ", file_url, " - ", conditionMessage(e))
    })
  }
}

# Loop over all base URLs and download the files
for (url in urls) {
  download_ifcb_files(url, extensions, data_path, convert_filenames)
}

raw_path <- file.path(ifcb_path, "sams", "raw")

zip_files <- list.files(path = raw_path, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)

# Loop over zip files and extract
for (zip_file in zip_files) {
  # Get folder name by removing .zip extension
  extract_folder <- tools::file_path_sans_ext(zip_file)
  
  # Create folder if it doesn't exist
  if (!dir.exists(extract_folder)) {
    dir.create(extract_folder, recursive = TRUE)
  }
  
  # Extract zip file into the created folder
  unzip(zip_file, exdir = extract_folder)
}

# List all extracted files from previous step
unzipped_files <- list.files(raw_path, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)

# Get file sizes for unzipped files
unzipped_df <- data.frame(
  filepath = unzipped_files,
  sample = basename(dirname(unzipped_files)),
  size = fs::file_size(unzipped_files)
)

# Get file sizes for png_files
png_df <- data.frame(
  filepath = png_files,
  sample = str_extract(basename(png_files), "^D\\d{8}T\\d{6}_IFCB\\d+"),
  size = fs::file_size(png_files)
)

# Join by size to find matches
matched_files <- inner_join(png_df, unzipped_df, by = c("sample", "size"), suffix = c("_png", "_unzipped"))

# Extract filename components to compare patterns
matched_files <- matched_files %>%
  mutate(
    base_png = str_remove(basename(filepath_png), "_\\d+\\.png$"),
    base_unzipped = str_remove(basename(filepath_unzipped), "_\\d+\\.png$"),
    portal_id = basename(filepath_png),
    dashboard_id = basename(filepath_unzipped)
  ) %>%
  filter(base_png == base_unzipped)  # Ensure filenames match except for last digits

# multiple_matches <- matched_files %>%
#   filter(duplicated(filepath_png) | duplicated(filepath_png, fromLast = TRUE))
# 
# unique_matches <- matched_files %>%
#   filter(!duplicated(filepath_png))

# Function to compute hash of a file
compute_hash <- function(filepath) {
  digest::digest(file = filepath, algo = "md5", serialize = FALSE)
}

# Compute hashes for PNG files
hash_df <- matched_files %>%
  mutate(hash_png = sapply(filepath_png, compute_hash),
         hash_unzipeed = sapply(filepath_unzipped, compute_hash))

# Find all images with identical hashes
single_matches <- hash_df %>%
  filter(hash_png == hash_unzipeed) %>%
  mutate(class = basename(dirname(filepath_png)))

# Some images are not found in Dashboard, they are listed here
missing_pngs <- png_df$filepath[!png_df$filepath %in% single_matches$filepath_png]

# Save translation list
write_tsv(single_matches,
          file.path(ifcb_path, "sams", "ifcb_portal_to_dashboard_ids.txt"),
          na = "")
