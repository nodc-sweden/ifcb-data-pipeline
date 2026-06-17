library(tidyverse)
library(curl)
library(iRfcb)

years <- 2014:2014
ifcb_path <- Sys.getenv("ifcb_path")
classes = NULL
dashboard_url = "https://ifcb-data.whoi.edu/data/"
data_folder <- file.path(ifcb_path, "whoi_plankton", "data")
png_folder <- file.path(ifcb_path, "whoi_plankton", "png_images", "extracted_images")
convert_filenames <- TRUE
convert_adc <- TRUE

for (year in years) {

  png_path <- file.path(png_folder, year)
  data_path <- file.path(data_folder, year)
  
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
  
  sample_df <- ifcb_convert_filenames(ifcb_string)
  
  # URL of the file
  urls <- unique(paste0(dashboard_url, sample_df$sample))
  
  # File extensions to download
  extensions <- c("roi", "hdr", "adc")
  
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
        date_object <- ifcb_parts[,2]
        
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
        ifcb_number <- ifcb_parts[,2]
        
        # Convert day of year to date
        if (convert_adc) {
          first_letter <- "D"
        } else {
          first_letter <- "I"
        }
        date_object <- paste0(first_letter, format(as.Date(paste0(year, "-01-01")) + day_of_year - 1, "%Y%m%d"))
        
        # Construct the filename in the desired format: DYYYYMMDD_THHMMSS_IFCB<ifcb_number>_YYYY_DDD_HHMMSS.ext
        filename <- paste0(date_object, "T", substr(time, 1, 2), substr(time, 3, 4), substr(time, 5, 6), "_", ifcb_number, ".", ext)
      } else {
        filename <- basename(file_url)
      }
      
      # Set the destination file path
      if (!is.null(date_object)) {
        if (convert_filenames & !convert_adc) {
          date_object <- sub("I", "D", date_object)
        }
        
        destfile <- file.path(dest_dir, date_object, filename)  # Use date_object with a D as directory
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
        handle <- curl::new_handle()
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
      
      if (ext == "adc" & convert_adc) {
        # Read adc file
        adc_file <- read.csv(destfile, header = FALSE)
        
        # Insert four empty columns after column 7
        adc_file <- cbind(adc_file[, 1:7], matrix(0, nrow = nrow(adc_file), ncol = 4), adc_file[, 8:ncol(adc_file)])
        
        # Rename columns to maintain numeric order (optional)
        colnames(adc_file) <- paste0("V", seq_len(ncol(adc_file)))
        
        # Store new file
        write.table(adc_file, file(destfile, "wb"), sep = ",", row.names = FALSE, col.names = FALSE, na = "", eol = "\n")
        
        # Print message
        message("Adjusted ADC file to new format: ", destfile)
      }
    }
  }
  
  # Loop over all base URLs and download the files
  for (url in urls) {
    download_ifcb_files(url, extensions, data_path, convert_filenames)
  }
}

