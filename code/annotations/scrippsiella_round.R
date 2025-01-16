library(tidyverse)

image_metadata <- read_tsv("/Users/anderstorstensson/Dropbox/R/ifcb-data-pipeline/images/Tångesund/Scrippsiella_group/ecotaxa_Scrippsiella_group.tsv")

round <- image_metadata %>%
  filter(object_eccentricity < 0.35) %>%
  filter(!img_file_name == "[t]") %>%
  rename(image_filename = img_file_name) %>%
  select(image_filename)

write_tsv(round,
          "output/round_scrippsiella.txt")
