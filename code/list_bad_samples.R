library(tidyverse)

ifcb_path <- Sys.getenv("ifcb_path")

files <- list.files(file.path(ifcb_path, "data/2025"), recursive = FALSE)

# Define the start and end timestamps
start_timestamp <- "D20250316"
end_timestamp <- "D20250411"

# Extract the timestamps from the file paths
timestamps <- sub(".*(D\\d{8}T\\d{6}_IFCB\\d+).*", "\\1", files)

# Filter the files that fall within the desired range
filtered_files <- files[timestamps >= start_timestamp & timestamps <= end_timestamp]

# Print the filtered files
files <- unique(tools::file_path_sans_ext(basename(filtered_files)))

# Save as tsv
write_tsv(as_tibble(files), "output/bad_files.tsv")
