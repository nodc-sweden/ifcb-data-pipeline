base_path <- "Z:/data/cnn/training_data/V3_niva_smhi_only/_filtered_dataset"

Sys.sleep(3600)

# Get all subdirectories (including the base folder itself)
subdirs <- list.dirs(path = base_path, recursive = TRUE)

# Count PNGs in each subdirectory
png_counts <- sapply(subdirs, function(d) {
  length(list.files(d, pattern = "\\.png$", ignore.case = TRUE))
})

# Combine results into a data frame
result <- data.frame(
  class_name = basename(subdirs),
  png_count = png_counts,
  row.names = NULL
)

# Extract unique taxa names from biovolumes and manual_biovolumes
taxa_names <- unique(result$class_name)

taxa_names_clean <- iRfcb:::truncate_folder_name(taxa_names)

# Special cases
taxa_names_clean <- gsub("Gymnodiniales_S", "Gymnodiniales", taxa_names_clean)
taxa_names_clean <- gsub("Gymnodiniales_CS", "Gymnodiniales", taxa_names_clean)
taxa_names_clean <- gsub("Gymnodiniales_CC", "Gymnodiniales", taxa_names_clean)

taxa_names_clean <- gsub("Ciliophora_S", "Ciliophora", taxa_names_clean)
taxa_names_clean <- gsub("cf_Proboscia_rhizosolenia", "Proboscia_Rhizosolenia", taxa_names_clean)

# Clean taxa_names by substituting specific patterns with spaces or empty strings
taxa_names_clean <- gsub("_", " ", taxa_names_clean)
taxa_names_clean <- gsub(" single cell", "", taxa_names_clean)
taxa_names_clean <- gsub(" single", "", taxa_names_clean)
taxa_names_clean <- gsub(" chain", "", taxa_names_clean)
taxa_names_clean <- gsub(" coil", "", taxa_names_clean)
taxa_names_clean <- gsub("-coiled", "", taxa_names_clean)
taxa_names_clean <- gsub(" filament", "", taxa_names_clean)
taxa_names_clean <- gsub(" pair", "", taxa_names_clean)
taxa_names_clean <- gsub("-like", "", taxa_names_clean)
taxa_names_clean <- gsub(" like", "", taxa_names_clean)
taxa_names_clean <- gsub(" bundle", "", taxa_names_clean)
taxa_names_clean <- gsub(" larger than 30", "", taxa_names_clean)
taxa_names_clean <- gsub(" larger than 30unidentified", "", taxa_names_clean)
taxa_names_clean <- gsub(" than 30", "", taxa_names_clean)
taxa_names_clean <- gsub(" smaller than 30unidentified", "", taxa_names_clean)
taxa_names_clean <- gsub(" smaller than 30", "", taxa_names_clean)
taxa_names_clean <- gsub(" smaller", "", taxa_names_clean)
taxa_names_clean <- gsub(" elliptical", "", taxa_names_clean)
taxa_names_clean <- gsub(" thin", "", taxa_names_clean)
taxa_names_clean <- gsub(" small", "", taxa_names_clean)
taxa_names_clean <- gsub(" sideview", "", taxa_names_clean)
taxa_names_clean <- gsub(" bottomview", "", taxa_names_clean)
taxa_names_clean <- gsub(" heterotrof", "", taxa_names_clean)
taxa_names_clean <- gsub(" heterotropic", "", taxa_names_clean)
taxa_names_clean <- gsub(" large", "", taxa_names_clean)
taxa_names_clean <- gsub(" thick", "", taxa_names_clean)
taxa_names_clean <- gsub(" sp01", "", taxa_names_clean)
taxa_names_clean <- gsub(" sp02", "", taxa_names_clean)
taxa_names_clean <- gsub(" sp1", "", taxa_names_clean)
taxa_names_clean <- gsub(" sp2", "", taxa_names_clean)

# Remove species flags from class names
taxa_names_clean <- gsub("\\<cf\\>", "", taxa_names_clean)
taxa_names_clean <- gsub("\\<spp\\>", "", taxa_names_clean)
taxa_names_clean <- gsub("\\<sp\\>", "", taxa_names_clean)
taxa_names_clean <- gsub(" group", "", taxa_names_clean)
taxa_names_clean <- gsub("  ", " ", taxa_names_clean)

# Turn f to f. for forma
taxa_names_clean <- gsub("\\bf\\b", "f.", taxa_names_clean)

# Add "/" for multiple names with capital letters
# e.g. Snowella_Woronichinia to Snowella/Woronichinia
taxa_names_clean <- gsub(" ([A-Z])", "/\\1", taxa_names_clean)
taxa_names_clean <- gsub(" ([A-Z])", "/\\1", taxa_names_clean)

# Use the first name of combined classes, 
# e.g. Nodularia_spumigena_coil,Nodularia_spumigena_filament to Nodularia spumigena
taxa_names_clean <- sapply(strsplit(taxa_names_clean, ","), `[`, 1)

# Remove any whitespace
taxa_names_clean <- trimws(taxa_names_clean)

# Retrieve worms records with retry mechanism
worms_records <- match_worms_taxa(taxa_names_clean,
                                  max_retries = 5,
                                  sleep_time = 60,
                                  marine_only = FALSE,
                                  verbose = FALSE)

# Extract Aphia IDs and class names from WoRMS
aphia_id <- worms_records$AphiaID

# Select relevant information from WoRMS
worms_df <- bind_rows(worms_records) %>%
  select(AphiaID, scientificname, rank, kingdom, phylum, class, order, family, genus, parentNameUsageID) %>%
  rename(worms_kingdom = kingdom, worms_phylum = phylum, worms_class = class, worms_order = order,
         worms_family = family, worms_genus = genus) %>%
  distinct()

# Create class_names data frame with taxa information
class_names <- data.frame(class = taxa_names,
                          class_clean = taxa_names_clean,
                          aphia_id,
                          sflag = ifelse(grepl("-like", taxa_names) | 
                                           grepl("_cf_", taxa_names)| 
                                           grepl("_like", taxa_names),
                                         "CF", NA)) %>%
  mutate(sflag = ifelse(grepl("\\<spp\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "SPP"), 
                        sflag)) %>%
  mutate(sflag = ifelse(grepl("\\<group\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "GRP"), 
                        sflag)) %>%
  mutate(sflag = ifelse(grepl("\\<sp\\>", gsub("_", " ", taxa_names)), 
                        paste(ifelse(is.na(sflag), "", sflag), "SP"), 
                        sflag)) %>%
  mutate(sflag = str_trim(sflag)) %>%
  mutate(trophic_type = ifcb_get_trophic_type(class_clean)) %>%
  # left_join(class_scores, by = "class") %>%
  left_join(worms_df, by = c("aphia_id" = "AphiaID")) %>%
  mutate(
    species_flag_taxon_name = case_when(
      is.na(sflag) ~ NA,
      sflag == "SPP" ~ paste(class_clean, "spp."),
      sflag == "GRP" ~ paste(class_clean, "group"),
      sflag == "CF" & rank == "Species" ~ sub(" ", " cf. ", class_clean), 
      sflag == "CF" & rank == "Genus" ~ paste0(worms_order, ": ", worms_family, " cf. ", worms_genus),
      TRUE ~ class_clean # Keeps class_clean as default
    )
  ) %>%
  distinct() %>%
  select(-rank, -worms_kingdom, -worms_phylum, 
         -worms_order, -worms_family, -worms_genus)

# Find all taxa classified with cf.
cf_taxa <- class_names %>%
  filter(sflag == "CF")

parent_ids <- unique(cf_taxa$parentNameUsageID)[!is.na(unique(cf_taxa$parentNameUsageID))]

# Extract parent records from cf. classes
parent_records <- get_worms_records(parent_ids) %>%
  select(AphiaID, scientificname) %>%
  rename(parentNameUsageID = AphiaID,
         parentName = scientificname)

# Replace scientific names for cf. taxa with the parent names
class_names <- class_names %>%
  left_join(parent_records, by = "parentNameUsageID") %>%
  mutate(scientificname_merge = parentName) %>%
  mutate(scientificname_merge = coalesce(scientificname_merge, scientificname)) %>%
  mutate(scientificname_merge = coalesce(scientificname_merge, class_clean)) %>%
  mutate(aphia_id_merge = ifelse(is.na(parentName), NA, parentNameUsageID)) %>%
  mutate(aphia_id_merge = coalesce(aphia_id_merge, aphia_id)) %>%
  mutate(sflag = ifelse(sflag == "CF", NA, sflag)) %>%
  select(-aphia_id, -scientificname, -parentName, -parentNameUsageID) %>%
  rename(aphia_id = aphia_id_merge,
         scientificname = scientificname_merge) %>%
  arrange(class_clean)

output <- result %>%
  left_join(class_names, by = c("class_name" = "class")) %>%
  filter(!class_name == "_filtered_dataset")

# Save as .txt (tab-delimited)
write.table(output, file = "png_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
