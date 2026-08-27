# Tångesund (IFCB110, 2016) cell-count data delivery to SHARK.
#
# This is the Tångesund counterpart of ifcb_delivery_data_to_shark_cellcounts.R.
# Tångesund was sampled at a fixed location (Hakefjorden) but at different depths,
# so the metadata is handled differently than the R/V Svea ferrybox dataset:
#   - Sample metadata is downloaded from the "Tangesund" dashboard dataset, which
#     already carries the fixed position and per-sample depth.
#   - There is no ferrybox or SVEPA cruise data for 2016, so those steps are dropped
#     and the *_FB columns are reported as NA.
#   - Station/metadata fields differ: STATN = "Tangesund", SHIPC = NA, PDMET = "SEP",
#     SMTYP = "IFCB", MNDEP/MXDEP from the per-sample depth, SMPNO = full pid.
# Classification follows the pytorch V6 cell-count pipeline, identical to the Svea script.

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
  library(ClassiPyR, quietly = TRUE)
  library(hdf5r, quietly = TRUE)
  library(jsonlite, quietly = TRUE)
  library(DBI, quietly = TRUE)
})

source("code/utils/clean_taxa_fn.R")

# -------------------------------
# Configuration
# -------------------------------

internal_use <- TRUE
cnn_model <- "SMHI-NIVA-SYKE-SAMS-SZN-ResNet50-V6_1" # niva_smhi_baltic or NA for MATLAB
threshold_file <- "V6_1/V6_1-resnet50_dataset_squarepad_augmented_b64_lr0.0001_e30_thresholds_and_metrics.json"
classlist_file <- "V6_1-resnet50_dataset_squarepad_augmented_b64_lr0.0001_e30_classes.txt"
instrument_number <- "IFCB110"
dashboard_dataset <- "Tangesund"

# Tångesund was only sampled in 2016 (IFCB110)
year <- 2016

remove_flagged_data <- c("bubbles", "incomplete")

ifcb_path <- Sys.getenv("ifcb_path")
ifcb_base_path <- Sys.getenv("ifcb_base_path")
havgem_path <- Sys.getenv("havgem_path")
repo_dir <- Sys.getenv("repo_path")
smtp_server <- Sys.getenv("SMTP_SERVER")
at_email <- Sys.getenv("at_email")
bk_email <- Sys.getenv("bk_email")
ifcb_email <- Sys.getenv("ifcb_email")

if (internal_use) {
  output_folder <- file.path(havgem_path, "data_internal", "cell_counts_new")
} else {
  output_folder <- file.path(havgem_path, "data_to_shark")
}

output_folder <- file.path(output_folder, cnn_model)

# Define data deliviery path
data_delivery_path <- file.path(output_folder,
                                paste0("SHARK_PlanktonImaging_",
                                       as.character(year),
                                       "_Tangesund_SMHI"))

# Define paths for processed, received, and correspondence data folders for current year
processed_data <- file.path(data_delivery_path, "processed_data")
received_data <- file.path(data_delivery_path, "received_data")
correspondence <- file.path(data_delivery_path, "correspondence")

feature_folder <- file.path(ifcb_path, "features", "v4", year)
raw_folder <- file.path(ifcb_path, "data", year)
psd_plot_folder <- file.path(ifcb_path, "psd", "figures", year)
metadata_folder <- file.path(ifcb_path, "ifcbdb_metadata")

blacklist_file <- file.path(repo_dir, "data", "sample_blacklist.tsv")
psd_data_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_data.csv"))
psd_fits_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_fits.csv"))
psd_flags_file <- file.path(ifcb_path, "psd", paste0("psd_", year, "_flags.csv"))
threshold_file <- file.path(ifcb_path, "pytorch", "models", threshold_file)
classlist_file <- file.path(ifcb_path, "pytorch", "models", classlist_file)

sqlite_path <- file.path(tools::R_user_dir("ClassiPyR", "data"), "annotations.sqlite")

class_folder <- file.path(ifcb_path, "classified", cnn_model)

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
diatom_include = c("Actinocyclus_spp", "Actinocyclus", "Navicula-like", "Navicula", "cf_Proboscia_rhizosolenia")

# -------------------------------
# Get data
# -------------------------------

# List already merged class files
combined_class_files <- list.files(class_folder, ".h5", full.names = TRUE, recursive = TRUE)
combined_class_files <- combined_class_files[grepl(paste0("D", year), basename(combined_class_files))]
class_bins <- gsub("_class.h5", "", basename(combined_class_files))

# Download sample metadata (Tångesund dataset carries fixed position and per-sample depth)
metadata <- ifcb_download_dashboard_metadata("https://ifcb-dashboard-utv.smhi.se/", dashboard_dataset)

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

valid_idx <- metadata_selected$latitude != -999 & metadata_selected$longitude != -999

# Preallocate
metadata_selected$near_land <- NA

# Only run function on valid positions
metadata_selected$near_land[valid_idx] <- as.character(
  iRfcb::ifcb_is_near_land(
    latitudes = metadata_selected$latitude[valid_idx],
    longitudes = metadata_selected$longitude[valid_idx],
    distance = 0,
    shape = shapefile_path
  )
)

# Mark the invalid ones
metadata_selected$near_land[!valid_idx] <- NA

# Read PSD flags
psd_flags <- read_csv(psd_flags_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique_quiet") %>%
  dplyr::select(-dplyr::any_of("...1"))

# No ferrybox or SVEPA cruise data exist for Tångesund 2016; those steps are skipped
# and the *_FB parameters are reported as NA in the SHARK output.

cat("Calculating biovolume and biomass data...\n")

biovolumes <- ifcb_summarize_biovolumes(feature_folder = feature_folder,
                                        class_files = combined_class_files,
                                        hdr_folder = raw_folder,
                                        micron_factor = 1 / 2.77,
                                        diatom_class = "Bacillariophyceae",
                                        diatom_include = diatom_include,
                                        use_python = FALSE,
                                        feature_version = 4,
                                        verbose = TRUE,
                                        drop_zero_volume = TRUE,
                                        use_cell_counts = TRUE,          # report cell counts
                                        diatom_equation = "auto",
                                        carbon_conversion = "cell",
                                        single_cell_values = c(-1, 0))   # default: uncounted/zero ROI = 1 cell

# Read manual annotations from SQlite database
con <- dbConnect(RSQLite::SQLite(), sqlite_path)
annotations <- dbReadTable(con, "annotations")
dbDisconnect(con)

annotations <- annotations %>%
  filter(grepl(paste0("D", year), sample_name)) %>%
  filter(grepl(paste0(instrument_number), sample_name)) %>%
  mutate(image_name = paste0(sample_name, "_", sprintf("%05d", roi_number)))

if (nrow(annotations) > 0) {
  manual_biovolumes <- ifcb_summarize_biovolumes(
    feature_folder = feature_folder,
    custom_images = annotations$image_name,
    custom_classes = annotations$class_name,
    hdr_folder = raw_folder,
    use_python = FALSE,
    verbose = TRUE,
    drop_zero_volume = FALSE,
    use_cell_counts = FALSE,
    diatom_equation = "auto",
    carbon_conversion = "roi"
  )

  # Manual annotations are single-image verifications (no chain counting);
  # treat each ROI as one cell so columns align with the autoclass data.
  manual_biovolumes$cell_counts <- manual_biovolumes$counts
  manual_biovolumes$cell_counts_per_liter <- manual_biovolumes$counts_per_liter
}

# Find samples
all_samples <- metadata_selected %>%
  pull(pid)

# raw <-readLines(classlist_file)

# cnn_class_names <- iRfcb:::truncate_folder_name(gsub("^'|'$", "", raw))

class_scores <- fromJSON(threshold_file)$class_metrics |>
  bind_rows(.id = "class")

manual_analysis_dates <- annotations %>%
  select(sample_name, annotator, timestamp) %>%
  distinct() %>%
  arrange(sample_name, timestamp) %>%
  group_by(sample_name) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(analysis_date = as.Date(timestamp)) %>%
  select(-timestamp) %>%
  rename(sample = sample_name)

# Get analysis date from autoclass files
autoclass_analysis_dates <- data.frame(sample = NULL, analysis_date = NULL)
for (mat in combined_class_files) {
  file_info <- file.info(mat)$ctime

  temp <- data.frame(sample = sub("_class.h5", "", basename(mat)), analysis_date = as.Date(file_info))

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

class_names <- clean_taxa(taxa_names) %>%
  mutate(trophic_type = ifcb_get_trophic_type(cleaned_name)) %>%
  rename(class_clean = cleaned_name,
         scientificname = reported_name,
         class = class_name,
         sflag = flag) %>%
  distinct()

# -------------------------------
# Join data
# -------------------------------

autoclass_joined <- biovolumes %>%
  left_join(ifcb_convert_filenames(metadata_selected$pid), by ="sample") %>%
  left_join(select(metadata_selected, -ml_analyzed), by = c("sample" = "pid")) %>%
  left_join(psd_flags, by ="sample") %>%
  left_join(class_names, by = "class") %>%
  left_join(class_scores, by = "class") %>%
  left_join(autoclass_analysis_dates, by = "sample") %>%
  mutate(platform = "IFCB",
         verification = "PredictedByMachine",
         verified_by = NA,
         classifier_created_by = "Anders Torstensson",
         # classifier_used = paste0("SMHI-", params$classifier, " v.", max_versions[[i]]),
         metoa = "IMA-SW",
         associated_media = NA,
         class_prog = "https://github.com/nodc-sweden/ifcb-pytorch-classify/releases/tag/v0.2.0")

autoclass_aggregated <- autoclass_joined %>%
  # Group by species where known, else by classifier class, so classes lacking a
  # WoRMS scientificname (Debris, Beads, Air_bubbles, Scrippsiella_group, mixed,
  # unclassified) each keep their own row instead of collapsing into one NA bin.
  mutate(agg_key = coalesce(scientificname, class)) %>%
  group_by(sample, agg_key) %>%
  summarise(
    # Sum the numeric measurements
    counts = sum(counts, na.rm = TRUE),
    cell_counts = sum(cell_counts, na.rm = TRUE),
    biovolume_mm3 = sum(biovolume_mm3, na.rm = TRUE),
    carbon_ug = sum(carbon_ug, na.rm = TRUE),
    ml_analyzed = first(ml_analyzed),  # Should be identical per sample
    counts_per_liter = sum(counts_per_liter, na.rm = TRUE),
    cell_counts_per_liter = sum(cell_counts_per_liter, na.rm = TRUE),
    biovolume_mm3_per_liter = sum(biovolume_mm3_per_liter, na.rm = TRUE),
    carbon_ug_per_liter = sum(carbon_ug_per_liter, na.rm = TRUE),

    # Concatenate these fields
    class = paste(class, collapse = ", "),
    # class_id = paste(class_id, collapse = ", "),
    threshold = paste(threshold, collapse = ", "),
    precision = paste(precision, collapse = ", "),
    recall = paste(recall, collapse = ",  "),
    support = paste(support, collapse = ", "),
    f1 = paste(round(f1, 3)*100, collapse = ", "),

    # Keep first value of all other columns (assumed identical)
    across(
      -c(counts, cell_counts, biovolume_mm3, carbon_ug, ml_analyzed,
         counts_per_liter, cell_counts_per_liter, biovolume_mm3_per_liter, carbon_ug_per_liter,
         class, precision, threshold, f1, recall, support),
      first
    ),
    .groups = "drop"
  ) %>%
  select(-agg_key)  # drop temporary grouping key

manual_joined <- manual_biovolumes %>%
  left_join(ifcb_convert_filenames(metadata_selected$pid), by ="sample") %>%
  left_join(select(metadata_selected, -ml_analyzed), by = c("sample" = "pid")) %>%
  left_join(psd_flags, by ="sample") %>%
  left_join(class_names, by = "class") %>%
  left_join(manual_analysis_dates, by = "sample") %>%
  mutate(platform = "IFCB",
         verification = "ValidatedByHuman",
         verified_by = annotator,
         classifier_created_by = NA,
         classifier_used = NA,
         metoa = "IMA",
         associated_media = "https://ecotaxa.obs-vlfr.fr/prj/899",
         class_prog = "https://github.com/EuropeanIFCBGroup/ClassiPyR/releases/tag/v0.2.1")

# Aggregate based on scientific name
manual_aggregated <- manual_joined %>%
  # Group by species where known, else by classifier class (see autoclass note).
  mutate(agg_key = coalesce(scientificname, class)) %>%
  group_by(sample, agg_key) %>%
  summarise(
    # Sum the numeric measurements
    counts = sum(counts, na.rm = TRUE),
    cell_counts = sum(cell_counts, na.rm = TRUE),
    biovolume_mm3 = sum(biovolume_mm3, na.rm = TRUE),
    carbon_ug = sum(carbon_ug, na.rm = TRUE),
    ml_analyzed = first(ml_analyzed),  # Should be identical per sample
    counts_per_liter = sum(counts_per_liter, na.rm = TRUE),
    cell_counts_per_liter = sum(cell_counts_per_liter, na.rm = TRUE),
    biovolume_mm3_per_liter = sum(biovolume_mm3_per_liter, na.rm = TRUE),
    carbon_ug_per_liter = sum(carbon_ug_per_liter, na.rm = TRUE),

    # Concatenate these fields
    class = paste(class, collapse = ","),

    # Keep first value of all other columns (assumed identical)
    across(
      -c(counts, cell_counts, biovolume_mm3, carbon_ug, ml_analyzed,
         counts_per_liter, cell_counts_per_liter, biovolume_mm3_per_liter, carbon_ug_per_liter,
         class),
      first
    ),
    .groups = "drop"
  ) %>%
  select(-agg_key)  # drop temporary grouping key

# Combine human verified data (manual annotations) with automatic classification
if (nrow(manual_aggregated) > 0) {
  data_aggregated <- bind_rows(autoclass_aggregated, manual_aggregated)
} else {
  data_aggregated <- autoclass_aggregated
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

# Extract unclassified data (no ferrybox columns exist for Tångesund)
unclassified <- data_aggregated %>%
  filter(class_clean == "unclassified") %>%
  select(
    -worms_class, -class_clean, -scientificname, -species_flag_taxon_name,
    -aphia_id, -sflag, -trophic_type, -biovolume_mm3,
    -carbon_ug, -carbon_ug_per_liter, -counts, -counts_per_liter,
    -any_of(c("precision", "detection_probability", "miss_probability",
              "f1", "class_id", "threshold", "recall", "support"))
  ) %>%
  rename(unclassified_counts = cell_counts,
         unclassified_abundance = cell_counts_per_liter,
         unclassified_volume = biovolume_mm3_per_liter)

# Remove unclassified from main data frame and bind with unclassified data
data_aggregated <- data_aggregated %>%
  filter(!class_clean == "unclassified") %>%
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

# Tångesund: fixed station sampled at different depths (per-sample depth from metadata)
station_name <- "Tangesund"
depth <- data_aggregated$depth

# Create shark_df by mapping relevant data from 'data' to Shark database columns
shark_df <- shark_col %>%
  mutate(MYEAR = data_aggregated$year,
         STATN = station_name,
         SAMPLING_PLATFORM = data_aggregated$platform,
         PROJ = ifelse(year < 2025, "IFCB, DTO, JERICO", "IFCB"),
         ORDERER = "SMHI",
         SHIPC = NA,
         CRUISE_NO = data_aggregated$cruise,
         SDATE = data_aggregated$date,
         DATE_TIME = as.character(format(data_aggregated$timestamp, "%Y%m%d%H%M%S")),
         TIMEZONE = "UTC",
         STIME = data_aggregated$time,
         LATIT = round(data_aggregated$latitude, 4),
         LONGI = round(data_aggregated$longitude, 4),
         POSYS = "GPS",
         PDMET = "SEP",
         METFP = "NON",
         IFCBNO	= data_aggregated$ifcb_number,
         MPROG = "PROJ",
         MNDEP = depth,
         MXDEP = depth,
         SLABO = "SMHI",
         ACKR_SMP = "N",
         SMTYP = "IFCB",
         SMVOL = round(data_aggregated$ml_analyzed, 3), # VOLUME
         SMPNO = data_aggregated$sample, # SAMPLE NAME (full pid)
         LATNM = data_aggregated$scientificname, # SPECIES
         SFLAG = data_aggregated$sflag, # SP or SPP
         LATNM_SFLAG = data_aggregated$species_flag_taxon_name,
         TRPHY = data_aggregated$trophic_type,
         APHIA_ID = data_aggregated$aphia_id,
         IMAGE_VERIFICATION = data_aggregated$verification,
         VERIFIED_BY = data_aggregated$verified_by,
         CLASS_NAME = data_aggregated$class,
         CLASS_F1 = data_aggregated$f1, # In percent
         COUNT = round(data_aggregated$cell_counts), # CELL COUNTS per SAMPLE
         IMAGE_COUNT = data_aggregated$counts, # NUMBER OF ROIs (images); mean chain length = COUNT / IMAGE_COUNT
         COEFF = round(1000/data_aggregated$ml_analyzed, 1),
         ABUND = round(data_aggregated$cell_counts_per_liter, 1), # CELLS PER LITER
         QFLAG = data_aggregated$flag,
         C_CONC = signif(data_aggregated$carbon_ug_per_liter, 6), # CARBON PER LITER
         BIOVOL = signif(data_aggregated$biovolume_mm3_per_liter, 6), # BIOVOLUME PER LITER
         METOA = data_aggregated$metoa,
         ASSOCIATED_MEDIA = data_aggregated$associated_media,
         CLASSPROG = data_aggregated$class_prog,
         ALABO = "SMHI",
         ACKR_ANA = "N",
         ANADATE = data_aggregated$analysis_date,
         METDC = paste("https://github.com/EuropeanIFCBGroup/ClassiPyR",
                       "https://github.com/nodc-sweden/ifcb-pytorch-classify",
                       "https://github.com/kudelalab/PSD",
                       "https://github.com/WHOIGit/ifcb-features/releases/tag/v1.0.0",
                       "https://github.com/EuropeanIFCBGroup/iRfcb",
                       sep = ", "), # METHOD
         TRAINING_SET = "https://doi.org/10.17044/scilifelab.25883455.v6",
         CLASSIFIER_USED = basename(gsub("\\\\", "/", iconv(data_aggregated$classifier, from = "Windows-1252", to = "UTF-8"))),
         MANUAL_QC_DATE = as.Date(Sys.Date()),
         PRE_FILTER_SIZE ="150", # unit um
         UNCLASSIFIED_COUNTS = round(data_aggregated$unclassified_counts, 1),
         UNCLASSIFIED_ABUNDANCE = signif(data_aggregated$unclassified_abundance, 6),
         UNCLASSIFIED_VOLUME = signif(data_aggregated$unclassified_volume, 6),
         # No ferrybox data for Tångesund 2016 - report *_FB parameters as NA
         PH_FB = NA,
         CHL_FB = NA,
         CDOM_FB = NA,
         PHYC_FB = NA,
         PHER_FB = NA,
         WATERFLOW_FB = NA,
         TURB_FB = NA,
         PCO2_FB = NA,
         TEMP_FB = NA,
         PSAL_FB = NA,
         OSAT_FB = NA,
         DOXY_FB = NA
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
                   as.character(year), "_Tangesund_",
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
  "kommentarer: Tångesund (IFCB110), fixed station sampled at different depths",
  "status: test"
)

# Write delivery note content to a .txt file in processed_data folder for current year
writeLines(delivery_note_content, file.path(processed_data, "delivery_note.txt"))

# -------------------------------
# Send E-mail notification
# -------------------------------
# 
# # Make path clickable on Windows
# data_delivery_path <- gsub("/", "\\\\", data_delivery_path)
# 
# # Dynamic emails
# names <- c("Anders Torstensson", "Bengt Karlson")
# emails <- c(at_email, bk_email)
# to_header <- paste0('"', names, '" <', emails, '>')  # no spaces added
# to_header <- paste(to_header, collapse = ",")        # no spaces after comma
# 
# # Build message like working format
# message <- c(
#   paste0('From: "IFCB" <ifcb.u@smhi.se>
# To:', to_header, '
# Subject: Dataleverans IFCB Tångesund
# 
# Hej,
# 
# Ny data:'),
#   data_delivery_path,
#   '
# 
# Mvh
# ifcb.u'
# )
# 
# # Deliver mail to recipients
# send_mail(
#   mail_from = ifcb_email,
#   mail_rcpt = emails,
#   message = message,
#   smtp_server = smtp_server,
#   use_ssl = "force",
#   verbose = FALSE
# )
