# 

# Use the default venv for this script
Sys.setenv(USE_IRFCB_PYTHON = TRUE)

# -------------------------------
# Load libraries
# -------------------------------

suppressPackageStartupMessages({
  library(SHARK4R, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(fs, quietly = TRUE)
  library(reticulate, quietly = TRUE)
  library(iRfcb, quietly = TRUE)
  library(tools, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(lubridate, quietly = TRUE)
  library(curl, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(tidyr, quietly = TRUE)
})

# -------------------------------
# Configuration
# -------------------------------

remove_flagged_data <- c("bubbles", "incomplete", "near land")

year <- as.numeric(format(Sys.Date(), "%Y"))
class_score_version <- 1

ifcb_path <- Sys.getenv("ifcb_path")
ifcb_base_path <- Sys.getenv("ifcb_base_path")
havgem_path <- Sys.getenv("havgem_path")
ferrybox_path <- Sys.getenv("ferrybox_path")
repo_dir <- Sys.getenv("repo_path")
smtp_server <- Sys.getenv("SMTP_SERVER")
at_email <- Sys.getenv("at_email")
bk_email <- Sys.getenv("bk_email")
ifcb_email <- Sys.getenv("ifcb_email")

output_folder <- file.path(havgem_path, "data_to_shark", "auto_processed")

# Define data deliviery path
data_delivery_path <- file.path(output_folder, 
                                paste0("SHARK_PlanktonImaging_", 
                                       as.character(year), 
                                       "_SMHI"))

# Define paths for processed, received, and correspondence data folders for current year
processed_data <- file.path(data_delivery_path, "processed_data")
received_data <- file.path(data_delivery_path, "received_data")
correspondence <- file.path(data_delivery_path, "correspondence")

feature_folder <- file.path(ifcb_path, "features", year)
raw_folder <- file.path(ifcb_path, "data", year)
psd_plot_folder <- file.path(ifcb_path, "psd", "figures", year)
metadata_folder <- file.path(ifcb_path, "ifcbdb_metadata")

blacklist_file <- file.path(repo_dir, "data", "sample_blacklist.tsv")
psd_data_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_data.csv"))
psd_fits_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_fits.csv"))
psd_flags_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_flags.csv"))

manual_folder_baltic <- file.path(ifcb_path, "manual", "Baltic")
manual_folder_westcoast <- file.path(ifcb_path, "manual", "Skagerrak-Kattegat")

# Define path to class2use files
class2use_file_baltic <- file.path(ifcb_path, "config", "class2use_Baltic.mat")
class2use_file_westcoast <- file.path(ifcb_path, "config", "class2use_Skagerrak-Kattegat.mat")

# Define the path for the combined classification folder for the specific year and class score version
class_folder <- file.path(ifcb_path, 
                          "classified", 
                          "Baltic-Skagerrak-Kattegat-dashboard", 
                          paste0("class", year, "_v", class_score_version))

# Define the URL of a EEA shapefile and the target directory
url <- "https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-2/gis-data/eea-coastline-polygon/at_download/file"
zip_file <- file.path(repo_dir, "data/shapefiles/EEA_Coastline_Polygon_Shape.zip")
extracted_dir <- file.path(repo_dir, "data/shapefiles/EEA_Coastline_Polygon_Shape")
shapefile_path <- file.path(extracted_dir, "EEA_Coastline_20170228.shp")

# Check if the shapefile already exists
if (!file.exists(shapefile_path)) {
  
  # Create the target directory if it doesn't exist
  if (!dir.exists(extracted_dir)) {
    dir.create(extracted_dir, recursive = TRUE)
  }
  
  # Download the ZIP file
  download.file(url, zip_file, mode = "wb")
  
  # Extract the ZIP file
  unzip(zip_file, exdir = extracted_dir)
  
  # Optionally, remove the downloaded ZIP file after extraction
  file.remove(zip_file)
  
  message("File downloaded and extracted.")
} else {
  message("Shapefile already exists. No download needed.")
}

# Taxa that should be forced to be regarded as diatoms, as they are not according to WoRMS
diatom_include = c("Actinocyclus_spp", "Actinocyclus", "Navicula-like", "Navicula")

# -------------------------------
# Get data
# -------------------------------

# List already merged class files
combined_class_files <- list.files(class_folder, ".mat", full.names = TRUE)
class_bins <- gsub(paste0("_class_v", class_score_version, ".mat"), "", basename(combined_class_files))

# Download sample metadata
metadata <- ifcb_download_dashboard_metadata("https://ifcb-dashboard-utv.smhi.se/", "RV_Svea")

# Filter classified bins
metadata_selected <- metadata %>%
  filter(pid %in% class_bins)

# Read sample blacklist
blacklist <- read_tsv(blacklist_file, progress = FALSE, col_types = cols())

# Create a regex pattern by joining all blacklist samples with `|`
blacklist_pattern <- paste(blacklist$sample, collapse = "|")

if (file.exists(file.path(processed_data, "data.txt"))) {
  # Read already delivered data
  delivered_data <- read_tsv(file.path(processed_data, "data.txt"), progress = FALSE, col_types = cols())
  
  # Identify the last 50 samples that have been classified
  valid_class_bins <- class_bins[!grepl(blacklist_pattern, class_bins)]
  last_50_samples <- sort(valid_class_bins, decreasing = TRUE)[1:50]
  
  # Check if any of the last 50 bins are already delivered
  all_delivered <- any(last_50_samples %in% delivered_data$SMPNO)
  
  # Stop if all data are already delivered
  if (all_delivered) {
    message("No rows in data frame, stopping script.")
    quit(save = "no", status = 0)
  }
}

cat("Calculating distance from land...\n")

# Add geographical information to metadata
metadata_selected <- metadata_selected %>%
  mutate(in_baltic = ifelse(latitude == -999, "missing_position", iRfcb::ifcb_is_in_basin(latitude, longitude)),
         near_land = iRfcb::ifcb_is_near_land(
           latitude,
           longitude,
           distance = 500, # Släggö is about 600 m from Land
           shape = shapefile_path))

# Read PSD flags
psd_flags <- read_csv(psd_flags_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique_quiet")

# Select ferrybox parameters to retrieve
ferrybox_parameters <- c("70", "80070", "8063", "88063", "8165", "88165", "8173", 
                         "88173", "8166","88166", "8172", "88172", "8174", "88174", 
                         "8177", "88177", "8179", "88179", "8181", "88181", "8190", 
                         "88190", "8191", "88191")

# Extract timestamps to get ferrybox data from
timestamps <- as.POSIXct(metadata_selected$sample_time)

cat("Getting ferrybox data...\n")

# Retrieve ferrybox data
ferrybox_data <- ifcb_get_ferrybox_data(timestamps, ferrybox_path, parameters = ferrybox_parameters,
                                        max_time_diff_min = 5)

# Set as NA if QC is not 1 (=Good)
ferrybox_data <- ferrybox_data %>%
  mutate(`80070` = ifelse(`70` == 1, `80070`, NA),
         `8063` = ifelse(`88063` == 1, `8063`, NA),
         `8165` = ifelse(`88165` == 1, `8165`, NA),
         `8173` = ifelse(`88173` == 1, `8173`, NA),
         `8166` = ifelse(`88166` == 1, `8166`, NA),
         `8172` = ifelse(`88172` == 1, `8172`, NA),
         `8174` = ifelse(`88174` == 1, `8174`, NA),
         `8177` = ifelse(`88177` == 1, `8177`, NA),
         `8179` = ifelse(`88179` == 1, `8179`, NA),
         `8181` = ifelse(`88181` == 1, `8181`, NA),
         `8190` = ifelse(`88190` == 1, `8190`, NA),
         `8191` = ifelse(`88191` == 1, `8191`, NA),
  )

# Add ferrybox data to metadata
metadata_selected <- bind_cols(metadata_selected, ferrybox_data)

cat("Calculating biovolume and biomass data...\n")

# Summarize biovolume and biomass data
biovolumes <- ifcb_summarize_biovolumes(feature_folder = feature_folder,
                                        mat_files = combined_class_files,
                                        hdr_folder = raw_folder,
                                        micron_factor = 1 / 3.4,
                                        diatom_class = "Bacillariophyceae",
                                        diatom_include = diatom_include,
                                        use_python = FALSE,
                                        feature_version = 2,
                                        verbose = TRUE,
                                        drop_zero_volume = TRUE)

# List manual Baltic files
manual_files_baltic <- list.files(manual_folder_baltic, ".mat", recursive = FALSE, full.names = TRUE)
manual_files_baltic <- manual_files_baltic[grepl(paste0("/D", year), manual_files_baltic)]

# Summarize data
if (length(manual_files_baltic) > 0) {
  cat("Calculating biovolume and biomass data for manual data...\n")
  
  manual_biovolumes_baltic <- ifcb_summarize_biovolumes(
    feature_folder = feature_folder,
    mat_files = manual_files_baltic,
    class2use_file = class2use_file_baltic,
    hdr_folder = raw_folder,
    mat_recursive = FALSE,
    use_python = FALSE,
    verbose = TRUE,
    drop_zero_volume = FALSE
  )
} else {
  manual_biovolumes_baltic <- data.frame(sample = character(), class = character())
}

# List manual Skagerrak-Kattegat files
manual_files_westcoast <- list.files(manual_folder_westcoast, ".mat", recursive = FALSE, full.names = TRUE)
manual_files_westcoast <- manual_files_westcoast[grepl(paste0("/D", year), manual_files_westcoast)]

# Summarize data
if (length(manual_files_baltic) > 0) {
  manual_biovolumes_westcoast <- ifcb_summarize_biovolumes(
    feature_folder = feature_folder,
    mat_files = manual_files_westcoast,
    class2use_file = class2use_file_westcoast,
    hdr_folder = raw_folder,
    mat_recursive = FALSE,
    use_python = FALSE,
    verbose = TRUE,
    drop_zero_volume = FALSE
  )
} else {
  manual_biovolumes_westcoast <- data.frame(sample = character(), class = character())
}

# Bind manual data together
manual_biovolumes <- bind_rows(manual_biovolumes_baltic, manual_biovolumes_westcoast)

# Find samples in the Baltic sea
baltic_samples <- metadata_selected %>%
  filter(in_baltic) %>%
  pull(pid)

# Find samples NOT in the baltic sea
westcoast_samples <- metadata_selected %>%
  filter(!in_baltic) %>%
  pull(pid)

# Extract the classifier name from the first file in class_files
classifier_baltic <- ifcb_get_mat_variable(combined_class_files[grepl(baltic_samples[1], combined_class_files)], "classifierName")[1]
classifier_westcoast <- ifcb_get_mat_variable(combined_class_files[grepl(westcoast_samples[1], combined_class_files)], "classifierName")[1]

# Construct the class_score_file name based on threshold and classifier
if (grepl("mat", classifier_baltic)) {
  class_score_file_baltic <- gsub(".mat", "_opt.csv", classifier_baltic)
} else {
  class_score_file_baltic <- paste0(classifier_baltic, "_opt.csv")
}

# Remap path to cross-platform path
if (!ifcb_path == "Z:/data") {
  
  # Replace the backslashes with forward slashes for consistency
  path_string <- gsub("\\\\", "/", class_score_file_baltic)
  
  # Replace volume path
  path_string <- gsub("Z:/data", ifcb_path, path_string)
  
  # Split the path into components
  path_components <- strsplit(path_string, "/")[[1]]
  
  # Reassemble the path using file.path
  class_score_file_baltic <- do.call(file.path, as.list(path_components))
}

# Check if the class_score_file_baltic exists
if (file.exists(class_score_file_baltic)) {
  # Read the class scores from the CSV file
  class_scores_baltic <- read_csv(class_score_file_baltic,
                                  show_col_types = FALSE)
} else {
  # Create a data frame with unique class names and NA values for PR, PD, and PM
  class_scores_baltic <- data.frame(class = unique(biovolume_data$class),
                                    precision = NA,
                                    detection_probability = NA,
                                    miss_probability = NA)
}

# Calculate F1-score
class_scores_baltic <- class_scores_baltic %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(f1 = 2 * (precision * detection_probability) / (precision + detection_probability),
         in_baltic = TRUE) 

# Construct the class_score_file name based on threshold and classifier
if (grepl("mat", classifier_westcoast)) {
  class_score_file_westcoast <- gsub(".mat", "_opt.csv", classifier_westcoast)
} else {
  class_score_file_westcoast <- paste0(classifier_westcoast, "_opt.csv")
}

# Remap path to cross-platform path
if (!ifcb_path == "Z:/data") {
  
  # Replace the backslashes with forward slashes for consistency
  path_string <- gsub("\\\\", "/", class_score_file_westcoast)
  
  # Replace volume path
  path_string <- gsub("Z:/data", ifcb_path, path_string)
  
  # Split the path into components
  path_components <- strsplit(path_string, "/")[[1]]
  
  # Reassemble the path using file.path
  class_score_file_westcoast <- do.call(file.path, as.list(path_components))
}

# Check if the class_score_file_westcoast exists
if (file.exists(class_score_file_westcoast)) {
  # Read the class scores from the CSV file
  class_scores_westcoast <- read_csv(class_score_file_westcoast,
                                  show_col_types = FALSE)
} else {
  # Create a data frame with unique class names and NA values for PR, PD, and PM
  class_scores_westcoast <- data.frame(class = unique(biovolume_data$class),
                                    precision = NA,
                                    detection_probability = NA,
                                    miss_probability = NA)
}

# Calculate F1-score
class_scores_westcoast <- class_scores_westcoast %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(f1 = 2 * (precision * detection_probability) / (precision + detection_probability),
         in_baltic = FALSE) 

# Bind data together in a single df
class_scores <- bind_rows(class_scores_westcoast, class_scores_baltic)

# List all .mat files in the class directory
mat_files_baltic <- list.files(path = manual_folder_baltic, pattern = "\\.mat$", full.names = TRUE)
mat_files_westcoast <- list.files(path = manual_folder_westcoast, pattern = "\\.mat$", full.names = TRUE)

# Get analysis date from manual files
mat_files <- c(mat_files_baltic, mat_files_westcoast)

manual_analysis_dates <- data.frame(sample = NULL, analysis_date = NULL)
for (mat in mat_files) {
  file_info <- file.info(mat)$ctime
  
  temp <- data.frame(sample = sub(".mat$", "", basename(mat)), analysis_date = as.Date(file_info))
  
  manual_analysis_dates <- rbind(manual_analysis_dates, temp)
}

# Get analysis date from autoclass files
autoclass_analysis_dates <- data.frame(sample = NULL, analysis_date = NULL)
for (mat in combined_class_files) {
  file_info <- file.info(mat)$ctime
  
  temp <- data.frame(sample = sub(paste0("_class_v", class_score_version, ".mat$"), "", basename(mat)), analysis_date = as.Date(file_info))
  
  autoclass_analysis_dates <- rbind(autoclass_analysis_dates, temp)
}

# -------------------------------
# Correct taxa names
# -------------------------------

# Extract unique taxa names from biovolumes and manual_biovolumes
taxa_names <- unique(biovolumes$class)
manual_taxa_names <- unique(manual_biovolumes$class)

# Combine manual and class names
taxa_names <- c(taxa_names, manual_taxa_names)

# Clean taxa_names by substituting specific patterns with spaces or empty strings
taxa_names_clean <- gsub("_", " ", taxa_names)
taxa_names_clean <- gsub(" single cell", "", taxa_names_clean)
taxa_names_clean <- gsub(" single", "", taxa_names_clean)
taxa_names_clean <- gsub(" chain", "", taxa_names_clean)
taxa_names_clean <- gsub(" coil", "", taxa_names_clean)
taxa_names_clean <- gsub(" filament", "", taxa_names_clean)
taxa_names_clean <- gsub(" pair", "", taxa_names_clean)
taxa_names_clean <- gsub("-like", "", taxa_names_clean)
taxa_names_clean <- gsub(" like", "", taxa_names_clean)
taxa_names_clean <- gsub(" bundle", "", taxa_names_clean)
taxa_names_clean <- gsub(" larger than 30", "", taxa_names_clean)
taxa_names_clean <- gsub(" larger than 30unidentified", "", taxa_names_clean)
taxa_names_clean <- gsub(" smaller than 30unidentified", "", taxa_names_clean)
taxa_names_clean <- gsub(" smaller than 30", "", taxa_names_clean)

# Remove species flags from class names
taxa_names_clean <- gsub("\\<cf\\>", "", taxa_names_clean)
taxa_names_clean <- gsub("\\<spp\\>", "", taxa_names_clean)
taxa_names_clean <- gsub("\\<sp\\>", "", taxa_names_clean)
taxa_names_clean <- gsub(" group", "", taxa_names_clean)
taxa_names_clean <- gsub("  ", " ", taxa_names_clean)

# Turn f to f. for forma
taxa_names_clean <- gsub("\\bf\\b", "f.", taxa_names_clean)

# Add "/" for multiple names with capital letters
# e.g. Snowella_Woronichinia to Snowella/Woronichinia
taxa_names_clean <- gsub(" ([A-Z])", "/\\1", taxa_names_clean)
taxa_names_clean <- gsub(" ([A-Z])", "/\\1", taxa_names_clean)

# Use the first name of combined classes, 
# e.g. Nodularia_spumigena_coil,Nodularia_spumigena_filament to Nodularia spumigena
taxa_names_clean <- sapply(strsplit(taxa_names_clean, ","), `[`, 1)

# Remove any whitespace
taxa_names_clean <- trimws(taxa_names_clean)

# Retrieve worms records with retry mechanism
worms_records <- ifcb_match_taxa_names(taxa_names_clean,
                                       max_retries = 5,
                                       sleep_time = 60,
                                       marine_only = FALSE,
                                       verbose = FALSE)

# Extract Aphia IDs and class names from WoRMS
aphia_id <- worms_records$AphiaID

# Select relevant information from WoRMS
worms_df <- bind_rows(worms_records) %>%
  select(AphiaID, scientificname, rank, kingdom, phylum, class, order, family, genus, parentNameUsageID) %>%
  rename(worms_kingdom = kingdom, worms_phylum = phylum, worms_class = class, worms_order = order,
         worms_family = family, worms_genus = genus) %>%
  distinct()

# Create class_names data frame with taxa information
class_names <- data.frame(class = taxa_names,
                          class_clean = taxa_names_clean,
                          aphia_id,
                          sflag = ifelse(grepl("-like", taxa_names) | 
                                           grepl("_cf_", taxa_names)| 
                                           grepl("_like", taxa_names),
                                         "CF", NA),
                          is_diatom = ifcb_is_diatom(taxa_names_clean,
                                                     diatom_include = diatom_include)) %>%
  mutate(sflag = ifelse(grepl("\\<spp\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "SPP"), 
                        sflag)) %>%
  mutate(sflag = ifelse(grepl("\\<group\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "GRP"), 
                        sflag)) %>%
  mutate(sflag = ifelse(grepl("\\<sp\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "SP"), 
                        sflag)) %>%
  mutate(sflag = str_trim(sflag)) %>%
  mutate(trophic_type = ifcb_get_trophic_type(class_clean)) %>%
  # left_join(class_scores, by = "class") %>%
  left_join(worms_df, by = c("aphia_id" = "AphiaID")) %>%
  mutate(
    species_flag_taxon_name = case_when(
      is.na(sflag) ~ NA,
      sflag == "SPP" ~ paste(class_clean, "spp."),
      sflag == "GRP" ~ paste(class_clean, "group"),
      sflag == "CF" & rank == "Species" ~ sub(" ", " cf. ", class_clean), 
      sflag == "CF" & rank == "Genus" ~ paste0(worms_order, ": ", worms_family, " cf. ", worms_genus),
      TRUE ~ class_clean # Keeps class_clean as default
    )
  ) %>%
  distinct() %>%
  select(-rank, -worms_kingdom, -worms_phylum, 
         -worms_order, -worms_family, -worms_genus)

# Find all taxa classified with cf.
cf_taxa <- class_names %>%
  filter(sflag == "CF")

parent_ids <- unique(cf_taxa$parentNameUsageID)[!is.na(unique(cf_taxa$parentNameUsageID))]

# Extract parent records from cf. classes
parent_records <- get_worms_records(parent_ids) %>%
  select(AphiaID, scientificname) %>%
  rename(parentNameUsageID = AphiaID,
         parentName = scientificname)

# Replace scientific names for cf. taxa with the parent names
class_names <- class_names %>%
  left_join(parent_records, by = "parentNameUsageID") %>%
  mutate(scientificname_merge = parentName) %>%
  mutate(scientificname_merge = coalesce(scientificname_merge, scientificname)) %>%
  mutate(scientificname_merge = coalesce(scientificname_merge, class_clean)) %>%
  mutate(aphia_id_merge = ifelse(is.na(parentName), NA, parentNameUsageID)) %>%
  mutate(aphia_id_merge = coalesce(aphia_id_merge, aphia_id)) %>%
  mutate(sflag = ifelse(sflag == "CF", NA, sflag)) %>%
  select(-aphia_id, -scientificname, -parentName, -parentNameUsageID) %>%
  rename(aphia_id = aphia_id_merge,
         scientificname = scientificname_merge) %>%
  arrange(class_clean)

# -------------------------------
# Join data
# -------------------------------

autoclass_joined <- biovolumes %>%
  left_join(ifcb_convert_filenames(metadata_selected$pid), by ="sample") %>%
  left_join(select(metadata_selected, -ml_analyzed, -timestamp), by = c("sample" = "pid")) %>%
  left_join(psd_flags, by ="sample") %>%
  left_join(class_names, by = "class") %>%
  left_join(class_scores, by = c("class", "in_baltic")) %>%
  left_join(autoclass_analysis_dates, by = "sample") %>%
  mutate(platform = "IFCB",
         verification = "PredictedByMachine",
         verified_by = NA,
         classifier_created_by = "Anders Torstensson",
         # classifier_used = paste0("SMHI-", params$classifier, " v.", max_versions[[i]]),
         metoa = "IMA-SW",
         associated_media = NA)

manual_joined <- manual_biovolumes %>%
  left_join(ifcb_convert_filenames(metadata_selected$pid), by ="sample") %>%
  left_join(select(metadata_selected, -ml_analyzed, -timestamp), by = c("sample" = "pid")) %>%
  left_join(psd_flags, by ="sample") %>%
  left_join(class_names, by = "class") %>%
  left_join(manual_analysis_dates, by = "sample") %>%
  mutate(platform = "IFCB",
         verification = "ValidatedByHuman",
         verified_by = "Ann-Turi Skjevik",
         classifier_created_by = NA,
         classifier_used = NA,
         metoa = "IMA",
         associated_media = ifelse(in_baltic, "https://ecotaxa.obs-vlfr.fr/prj/822", "https://ecotaxa.obs-vlfr.fr/prj/14392"))

# Combine human verified data (manual annotations) with automatic classification
if (nrow(manual_joined) > 0) {
  data_aggregated <- bind_rows(autoclass_joined, manual_joined)
} else {
  data_aggregated <- autoclass_joined
}


# -------------------------------
# Remove bad data
# -------------------------------

# Remove samples before and after a bubble run, as these samples often contain poor images
if ("bubbles" %in% remove_flagged_data) {
  
  # Identify the indices where the flag contains "bubbles"
  bubbles_indices <- which(grepl("bubbles", tolower(data_aggregated$flag), fixed = TRUE))
  
  # Find the indices before and after the "bubbles" samples
  before_indices <- bubbles_indices - 1
  after_indices <- bubbles_indices + 1
  
  # Remove indices that are out of bounds
  before_indices <- before_indices[before_indices > 0]
  after_indices <- after_indices[after_indices <= nrow(data_aggregated)]
  
  # Extract the samples before and after the "bubbles" samples
  samples_before <- data_aggregated$sample[before_indices]
  samples_after <- data_aggregated$sample[after_indices]
  
  # Combine the results into a data frame for easier viewing
  adjacent_bubble_samples <- data.frame(
    BubblesSample = data_aggregated$sample[bubbles_indices],
    SampleBefore = samples_before,
    SampleAfter = samples_after
  )
  
  # Remove data from samples before and after bubble runs
  data_aggregated <- data_aggregated %>%
    filter(!sample %in% adjacent_bubble_samples$SampleBefore) %>%
    filter(!sample %in% adjacent_bubble_samples$SampleAfter)
}

# Filter flagged data
data_aggregated <- data_aggregated[!sapply(data_aggregated$flag, function(x) {
  any(sapply(remove_flagged_data, function(y) grepl(y, tolower(x), fixed = TRUE)))
}), ]

# Filter blacklisted data
data_aggregated <- data_aggregated %>%
  filter(!str_detect(sample, blacklist_pattern))

# Function to check if a column is entirely NA
is_all_na <- function(x) {
  all(is.na(x))
}

# Extract unclassified data and ferrybox data
unclassified <- data_aggregated %>%
  filter(class_clean == "unclassified") %>%
  select(-worms_class, -class_clean, -scientificname, -species_flag_taxon_name, -aphia_id, -sflag, -trophic_type, -is_diatom, -biovolume_mm3, -carbon_ug, 
         -class, -precision, -detection_probability, -miss_probability, -f1, -carbon_ug_per_liter) %>%
  rename(unclassified_counts = counts,
         unclassified_abundance = counts_per_liter,
         unclassified_volume = biovolume_mm3_per_liter,
         ph = `70`,
         chl = `8063`,
         cdom = `8165`,
         pc = `8173`,
         pe = `8166`,
         waterflow = `8172`,
         turb = `8174`,
         pco2 = `8177`,
         watertemp = `8179`,
         salinity = `8181`,
         oxygen_saturation = `8190`,
         oxygen_ml_l = `8191`) %>%
  mutate(across(all_of(intersect(c("ph", "chl", "cdom", "pc", "pe", 
                                   "waterflow", "turb", "pco2", "oxygen_saturation"), 
                                 names(.)[!sapply(., is_all_na)])),
                ~ na_if(.x, -999))) %>% # Replace -999 values with NA
  mutate(across(all_of(intersect(c("ph", "pco2"), 
                                 names(.)[!sapply(., is_all_na)])),
                ~ na_if(.x, 0))) # Replace 0 values with NA

# Remove unclassified from main data frame and bind with unclassified and ferrybox data
data_aggregated <- data_aggregated %>%
  filter(!class_clean == "unclassified") %>%
  rename(ph = `70`,
         chl = `8063`,
         cdom = `8165`,
         pc = `8173`,
         pe = `8166`,
         waterflow = `8172`,
         turb = `8174`,
         pco2 = `8177`,
         watertemp = `8179`,
         salinity = `8181`,
         oxygen_saturation = `8190`,
         oxygen_ml_l = `8191`) %>%
  bind_rows(unclassified) %>%
  arrange(desc(verification), sample, class_clean)

# -------------------------------
# Plot data
# -------------------------------

# Aggregate biovoume and image count data
sample_sum <- data_aggregated %>%
  group_by(sample, timestamp) %>%
  summarise(
    biovolume_mm3 = sum(biovolume_mm3, na.rm = TRUE),
    n_images = unique(n_images, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(vol_to_image_ratio = biovolume_mm3/n_images*1000) %>%
  arrange(desc(vol_to_image_ratio))

# Transform to long
sample_long <- sample_sum %>%
  pivot_longer(
    cols = c(biovolume_mm3, n_images),
    names_to = "variable",
    values_to = "value"
  )

# Plot data
qc_plot <- ggplot(sample_long, aes(x = timestamp, y = value)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 1.4, alpha = 0.8) +
  facet_wrap(
    ~ variable,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(c(
      biovolume_mm3 = expression(paste("Total biovolume (mm"^3, ")")),
      n_images = "Number of images"
    ))
  ) +
  labs(x = "Time", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

if (!dir.exists(received_data)) {
  dir.create(received_data, recursive = TRUE)
}

# Store results
ggsave(file.path(received_data, "qc_plot.png"),
       bg = "white")

# Store data
write_tsv(sample_sum, file.path(received_data, "sample_sums.txt"), 
          progress = FALSE, na = "")

# Store taxonomic info
write_tsv(class_names, file.path(received_data, "class_names.txt"), 
          progress = FALSE, na = "")

# -------------------------------
# Map column names
# -------------------------------

# Retrieve column names for Shark database integration
shark_col <- ifcb_get_shark_colnames()

# Create a data frame with empty rows matching the length of data
shark_col[nrow(shark_col) + seq_len(nrow(data_aggregated)), ] <- NA

station_name <- paste0("RV_FB_", data_aggregated$sample)
depth <- 4

# Create shark_df by mapping relevant data from 'data' to Shark database columns
shark_df <- shark_col %>%
  mutate(MYEAR = data_aggregated$year,
         STATN = station_name,
         SAMPLING_PLATFORM = data_aggregated$platform,
         PROJ = ifelse(year < 2025, "IFCB, DTO, JERICO", "IFCB"),
         ORDERER = "SMHI",
         SHIPC = "77SE",
         CRUISE_NO = data_aggregated$cruise,
         SDATE = data_aggregated$date,
         DATE_TIME = as.character(format(data_aggregated$timestamp, "%Y%m%d%H%M%S")),
         TIMEZONE = "UTC",
         STIME = data_aggregated$time,
         LATIT = data_aggregated$latitude,
         LONGI = data_aggregated$longitude,
         POSYS = "GPS",
         PDMET = "MIX",
         METFP = "NON",
         IFCBNO	= data_aggregated$ifcb_number,
         MPROG = "PROJ",
         MNDEP = depth,
         MXDEP = depth,
         SLABO = "SMHI",
         ACKR_SMP = "N",
         SMTYP = "FerryBox",
         SMVOL = round(data_aggregated$ml_analyzed, 3), # VOLUME
         SMPNO = data_aggregated$sample, # SAMPLE NAME
         LATNM = data_aggregated$scientificname, # SPECIES
         SFLAG = data_aggregated$sflag, # SP or SPP
         LATNM_SFLAG = data_aggregated$species_flag_taxon_name,
         TRPHY = data_aggregated$trophic_type,
         APHIA_ID = data_aggregated$aphia_id,
         IMAGE_VERIFICATION = data_aggregated$verification,
         VERIFIED_BY = data_aggregated$verified_by,
         CLASS_NAME = data_aggregated$class,
         CLASS_F1 = round(data_aggregated$f1*100, 3), # In percent
         COUNT = data_aggregated$counts, # COUNTS per SAMPLE
         COEFF = round(1000/data_aggregated$ml_analyzed, 1),
         ABUND = round(data_aggregated$counts_per_liter, 1), #COUNTS PER LITER
         QFLAG = data_aggregated$flag,
         C_CONC = signif(data_aggregated$carbon_ug_per_liter, 6), # CARBON PER LITER
         BIOVOL = signif(data_aggregated$biovolume_mm3_per_liter, 6), # BIOVOLUME PER LITER
         METOA = data_aggregated$metoa,
         ASSOCIATED_MEDIA = data_aggregated$associated_media,
         CLASSPROG = paste("https://github.com/hsosik/ifcb-analysis/releases/tag/v2.0"),
         ALABO = "SMHI",
         ACKR_ANA = "N",
         ANADATE = data_aggregated$analysis_date,
         METDC = paste("https://github.com/hsosik/ifcb-analysis",
                       "https://github.com/kudelalab/PSD",
                       "https://github.com/EuropeanIFCBGroup/iRfcb",
                       sep = ", "), # METHOD
         TRAINING_SET = "https://doi.org/10.17044/scilifelab.25883455.v4",
         CLASSIFIER_USED = basename(gsub("\\\\", "/", iconv(data_aggregated$classifier, from = "Windows-1252", to = "UTF-8"))),
         MANUAL_QC_DATE = as.Date(Sys.Date()),
         PRE_FILTER_SIZE ="150", # unit um
         UNCLASSIFIED_COUNTS = round(data_aggregated$unclassified_counts, 1),
         UNCLASSIFIED_ABUNDANCE = signif(data_aggregated$unclassified_abundance, 6),
         UNCLASSIFIED_VOLUME = signif(data_aggregated$unclassified_volume, 6),
         PH_FB = NA, # round(data_aggregated$ph, 3), # Wait until method is validated
         CHL_FB = round(data_aggregated$chl, 4),
         CDOM_FB = round(data_aggregated$cdom, 4),
         PHYC_FB = data_aggregated$pc,
         PHER_FB = round(data_aggregated$pe, 4),
         WATERFLOW_FB = round(data_aggregated$waterflow, 4),
         TURB_FB = round(data_aggregated$turb, 4),
         PCO2_FB = data_aggregated$pco2,
         TEMP_FB = round(data_aggregated$watertemp, 4),
         PSAL_FB = round(data_aggregated$salinity, 4),
         OSAT_FB = round(data_aggregated$oxygen_saturation, 4),
         DOXY_FB = round(data_aggregated$oxygen_ml_l, 4)
  ) 

# -------------------------------
# Data delivery
# -------------------------------

# Create directories if they do not exist
if (!dir.exists(processed_data)) {
  dir.create(processed_data, recursive = TRUE)
}

if (!dir.exists(received_data)) {
  dir.create(received_data, recursive = TRUE)
}

if (!dir.exists(correspondence)) {
  dir.create(correspondence, recursive = TRUE)
}

# Save shark_df data to a tab-delimited file in processed_data folder for current year
write_tsv(shark_df, file = file.path(processed_data, "data.txt"), na = "", progress = FALSE) # Save as tab-delimited file

# Save shark_df data to a tab-delimited file in received_data folder with dynamic filename
filename <- paste0("shark_data_", 
                   as.character(year), "_", 
                   Sys.Date(), ".txt")
write_tsv(shark_df, file = file.path(received_data, filename), na = "", progress = FALSE) # Save as tab-delimited file

# Define delivery note content as a character vector for current year
delivery_note_content <- c(
  paste("provtagningsår:", year),
  "datatyp: Plankton Imaging",
  "rapporterande institut: SMHI",
  paste("rapporteringsdatum:", Sys.Date()),
  "kontaktperson: Anders Torstensson",
  "format: PlanktonImaging:IFCB_SMHI",
  "data kontrollerad av: Leverantör",
  "övervakningsprogram: PROJ",
  "beställare: SMHI",
  "projekt: IFCB, DTO, JERICO",
  "kommentarer:",
  "status: test"
)

# Write delivery note content to a .txt file in processed_data folder for current year
writeLines(delivery_note_content, file.path(processed_data, "delivery_note.txt"))

# -------------------------------
# Send E-mail notification
# -------------------------------

# Make path clickable on Windows
data_delivery_path <- gsub("/", "\\\\", data_delivery_path)

# Dynamic emails
names <- c("Anders Torstensson", "Bengt Karlson")
emails <- c(at_email, bk_email)
to_header <- paste0('"', names, '" <', emails, '>')  # no spaces added
to_header <- paste(to_header, collapse = ",")        # no spaces after comma

# Build message like working format
message <- c(
  paste0('From: "IFCB" <ifcb.u@smhi.se>
To:', to_header, '
Subject: Dataleverans IFCB

Hej,

Ny data:'),
  data_delivery_path,
  '

Mvh
ifcb.u'
)

# Deliver mail to recipients
send_mail(
  mail_from = ifcb_email,
  mail_rcpt = emails,
  message = message,
  smtp_server = smtp_server,
  use_ssl = "force",
  verbose = FALSE
)
