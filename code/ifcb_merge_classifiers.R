library(iRfcb)

# Setup python environment
ifcb_py_install(envname = ".virtualenvs/iRfcb")

ifcb_path <- Sys.getenv("ifcb_path")

# Name the new, merged dataset (folder and class2use names)
merged_dataset_name <- "test"

# Define paths to class2use files
class2use_file_base <- file.path(ifcb_path, "config", "class2use_Skagerrak-Kattegat.mat")
class2use_file_additions <- file.path(ifcb_path, "config", "class2use_Tångesund.mat")
class2use_file_combined <- file.path(ifcb_path, "config", paste0("class2use_", merged_dataset_name, ".mat"))

# Define paths to manual folders
manual_folder_base <- file.path(ifcb_path, "manual", "Skagerrak-Kattegat")
manual_folder_additions <- file.path(ifcb_path, "manual", "Tångesund")
manual_folder_combined <- file.path(ifcb_path, "manual", merged_dataset_name)

# Merge the datasets
ifcb_merge_manual(class2use_file_base,
                  class2use_file_additions,
                  class2use_file_combined,
                  manual_folder_base,
                  manual_folder_additions,
                  manual_folder_combined)

# Extract the images to verify the result
ifcb_extract_annotated_images(manual_folder_combined,
                              class2use_file_combined,
                              file.path(ifcb_path, "data"),
                              file.path(ifcb_path, "png_images", merged_dataset_name),
                              skip_class = "unclassified",
                              verbose = FALSE,
                              overwrite = TRUE)
