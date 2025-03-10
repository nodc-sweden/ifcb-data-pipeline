library(iRfcb)

# Setup python environment
ifcb_py_install(envname = ".virtualenvs/iRfcb")

# Get path to data from env
ifcb_path <- Sys.getenv("ifcb_path")

# Name the new, merged dataset (folder and class2use names)
smhi_merged_dataset_name <- "Skagerrak-Kattegat-Tangesund"
smhi_niva_dataset_name <- "niva_smhi"

# Define paths to class2use files
class2use_file_base <- file.path(ifcb_path, "config", "class2use_Skagerrak-Kattegat.mat")
class2use_file_tangesund <- file.path(ifcb_path, "config", "class2use_Tangesund.mat")
class2use_file_niva <- file.path(ifcb_path, "ifcb139", "classifier", "class2use_20250310.mat")
class2use_file_combined_smhi <- file.path(ifcb_path, "config", paste0("class2use_", smhi_merged_dataset_name, ".mat"))
class2use_file_combined_smhi_niva <- file.path(ifcb_path, "config", paste0("class2use_", smhi_niva_dataset_name, ".mat"))

# Define paths to manual folders
manual_folder_base <- file.path(ifcb_path, "manual", "Skagerrak-Kattegat")
manual_folder_tangesund <- file.path(ifcb_path, "manual", "Tangesund")
manual_folder_niva <- file.path(ifcb_path, "ifcb139", "manual")
manual_folder_smhi <- file.path(ifcb_path, "manual", smhi_merged_dataset_name)
manual_folder_smhi_niva <- file.path(ifcb_path, "manual", smhi_niva_dataset_name)

# Merge Tångesund and Skagerrak-Kattegat datasets
ifcb_merge_manual(class2use_file_base,
                  class2use_file_tangesund,
                  class2use_file_combined_smhi,
                  manual_folder_base,
                  manual_folder_tangesund,
                  manual_folder_smhi)

# Merge the SMHI (Skagerrak-Kattegat-Tångesund) and NIVA datasets
ifcb_merge_manual(class2use_file_combined_smhi,
                  class2use_file_niva,
                  class2use_file_combined_smhi_niva,
                  manual_folder_smhi,
                  manual_folder_niva,
                  manual_folder_smhi_niva)

# Extract the images to verify the result
ifcb_extract_annotated_images(manual_folder_smhi_niva,
                              class2use_file_combined_smhi_niva,
                              c(file.path(ifcb_path, "data"), file.path(ifcb_path, "ifcb139", "raw")),
                              file.path(ifcb_path, "png_images", smhi_niva_dataset_name, Sys.Date()),
                              skip_class = "unclassified",
                              verbose = TRUE,
                              overwrite = TRUE)
