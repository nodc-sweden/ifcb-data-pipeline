library(iRfcb)
library(dplyr)

ifcb_path <- Sys.getenv("ifcb_path")

# Download sample metadata
metadata <- ifcb_download_dashboard_metadata("https://ifcb-dashboard-utv.smhi.se/", "RV_Svea")

west_coast_samples <- metadata %>%
  filter(tag1 == "skagerrak_kattegat")

classifier <- "niva_smhi"
year <- 2023

raw_folder <- file.path(ifcb_path, "data", year)
manual_folder <- file.path(ifcb_path, "manual", classifier)
class_folder <- file.path(ifcb_path, "classified", classifier, paste0("class", year, "_", classifier, "_v1"))

class_files <- list.files(class_folder, recursive = TRUE)
manual_files <- list.files(manual_folder, pattern = ".mat", recursive = FALSE)

# Remove already annotated samples
class_files <- class_files[!gsub("_class_v1", "", class_files) %in% manual_files]

# Optionally, remove baltic or west coast samples
class_files <- class_files[gsub("_class_v1.mat", "", class_files) %in% west_coast_samples$pid]

class_files_samples <- sample(class_files, 20)
class_files_samples <- gsub("_class_v1.mat", "", class_files_samples)

for (sample in class_files_samples) {
  ifcb_extract_classified_images(sample,
                                 class_folder,
                                 raw_folder,
                                 file.path(ifcb_path, "classified_images", classifier, Sys.Date(), sample))
}
