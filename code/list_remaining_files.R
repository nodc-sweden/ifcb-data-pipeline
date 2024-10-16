library(tidyverse)
library(iRfcb)

list <- read.table("output/export_899_20241013_1759_selected_images_cryptomonadales.txt",
                   header = TRUE)

list2 <- read.table("output/export_899_20241013_1759_selected_images_eutreptiella.txt",
                   header = TRUE)

list <- rbind(list, list2)

all_files <- list.files("images/export_899_20241013_1759")

unclassified <- all_files %>%
  as.data.frame() %>%
  rename(image_filename = ".") %>%
  filter(!image_filename %in% list$image_filename)

write_tsv(unclassified, "output/export_899_20241013_1759_selected_images_unclassified.txt")
