library(iRfcb)

ifcb_path <- Sys.getenv("ifcb_path")

png_images <- list.files("data/Karins_suggestions/Tripos_muelleri")
class <- "Tripos_muelleri"
manual_folder <- "data/manual"
adc_folder <- file.path(ifcb_path, "data")
class2use_file <- file.path(ifcb_path, "config", "class2use_Tångesund.mat")

ifcb_py_install()

ifcb_write_manual_file(png_images,
                       class,
                       manual_folder,
                       adc_folder,
                       class2use_file)
