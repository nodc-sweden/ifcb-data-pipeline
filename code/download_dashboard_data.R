### Superseded by `algaware` package


library(tidyverse)
library(sf)
library(iRfcb)
library(SHARK4R)

# Load helper functions
source(file.path("src", "R", "helpers.R"))

# Define the URL to the IFCB Dashboard
dashboard_url <- "https://ifcb-dashboard-utv.smhi.se"

# Get metadata from IFCB Dashboard
svea_metadata <- ifcb_get_dashboard_metadata(dashboard_url, "RV_Svea")

# Identify available cruise numbers
cruise_numbers <- unique(svea_metadata$cruise)

# Select the cruise of interest
selected_cruise <- last(cruise_numbers)
selected_cruise <- "SVEA_2025_022"

# Filter metada for the selected cruise
latest_cruise_metadata <- svea_metadata %>%
  filter(cruise == selected_cruise)

# Read station list
station_list <- SHARK4R:::load_station_bundle(verbose = FALSE)

# Define AlgAware stations
algaware_stations <- read_tsv(file.path("resources", 
                                        "algaware_stations.txt"),
                              col_types = cols(),
                              progress = FALSE)

# Filter AlgAware station data
algaware_station_data <- station_list %>%
  filter(STATION_NAME %in% algaware_stations$STATION_NAME)

# Convert both datasets to sf objects
# assuming lat/long are in decimal degrees (WGS84)
stations_sf <- algaware_station_data %>%
  st_as_sf(coords = c("LONGITUDE_WGS84_SWEREF99_DD", "LATITUDE_WGS84_SWEREF99_DD"), crs = 4326)

metadata_sf <- latest_cruise_metadata %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Transform to a projected coordinate system (for meter-based distance)
# SWEREF99 TM (EPSG:3006) is good for Sweden
stations_sf <- st_transform(stations_sf, 3006)
metadata_sf <- st_transform(metadata_sf, 3006)

# Create station buffers using OUT_OF_BOUNDS_RADIUS ----
station_buffers <- st_buffer(stations_sf, dist = algaware_station_data$OUT_OF_BOUNDS_RADIUS)

# Spatial join: find all metadata points inside any station buffer ----
metadata_joined <- st_join(
  metadata_sf,
  station_buffers %>% select(STATION_NAME),
  join = st_within
)

# Convert back to tibble
algaware_metadata <- metadata_joined %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  filter(!is.na(STATION_NAME))

# Download all raw files
ifcb_download_dashboard_data(dashboard_url = dashboard_url,
                             samples = algaware_metadata$pid,
                             file_types = c("roi", "adc", "hdr"),
                             dest_dir = file.path("data", "raw"))

# Download all available feature files
ifcb_download_dashboard_data(dashboard_url = dashboard_url,
                             samples = algaware_metadata$pid,
                             file_types = "features",
                             dest_dir = file.path("data", "features"))

# List all roi files
roi_files <- list.files(file.path("data", "raw"), 
                        ".roi", 
                        recursive = TRUE, 
                        full.names = TRUE)

# Create progress bar
cat("Extracting .png images...\n")
pb <- utils::txtProgressBar(min = 0, max = length(roi_files), style = 3)

for (i in seq_along(roi_files)) {
  roi_path <- roi_files[i]
  
  # Extract into that directory
  ifcb_extract_pngs(roi_path,
                    "data/png",
                    verbose = FALSE)
  
  # Update progress bar
  utils::setTxtProgressBar(pb, i)
}

# Close progress bar
close(pb)

# Store metadata file for later use
write_tsv(algaware_metadata,
          file.path("data", "metadata", paste0(selected_cruise, "_metadata.txt")),
          progress = FALSE)
