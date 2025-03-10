library(iRfcb)
library(tidyverse)

ifcb_path <- Sys.getenv("ifcb_path")

classlist_smhi <- as.data.frame(ifcb_get_mat_variable(file.path(ifcb_path, "config", "class2use_Skagerrak-Kattegat.mat")))
classlist_tangesund <- as.data.frame(ifcb_get_mat_variable(file.path(ifcb_path, "config", "class2use_Tangesund.mat")))
classlist_niva <- as.data.frame(ifcb_get_mat_variable(file.path(ifcb_path, "NIVA IFCB", "classifier", "class2use_20250310.mat")))

names(classlist_smhi) <- "class"
names(classlist_niva) <- "class"
names(classlist_tangesund) <- "class"

classlist_smhi$smhi <- "x"
classlist_tangesund$tangesund <- "x"
classlist_niva$niva <- "x"

classlist_join <- full_join(classlist_smhi, classlist_niva, by = "class")
classlist_join <- arrange(classlist_join, class)

classlist_west_coast <- full_join(classlist_smhi, classlist_tangesund, by = "class")
classlist_west_coast <- arrange(classlist_west_coast, class)
