library(iRfcb)

year <- "2025"

ifcb_path <- Sys.getenv("ifcb_path")

data_folder <- file.path(ifcb_path, "data", year)

ext_dir <- file.path(ifcb_path, "extracted_images", year)

files <- list.files(data_folder, pattern = ".roi", recursive = TRUE, full.names = TRUE)
existing_dirs <- list.dirs(ext_dir, full.names = FALSE, recursive = FALSE)

missing_files <- files[!vapply(files, function(f)
  any(grepl(basename(tools::file_path_sans_ext(f)), existing_dirs, fixed = TRUE)), logical(1))]

for (file in missing_files) {
  ifcb_extract_pngs(file, ext_dir)
}
