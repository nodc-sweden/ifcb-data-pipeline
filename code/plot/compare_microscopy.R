library(SHARK4R)
library(tidyverse)
library(geosphere)

# Define paths
ifcb_path <- Sys.getenv("ifcb_path")
baltic_path <- file.path("multiyear", "SHARK_IFCB_2022_2024_Baltic_SMHI")
west_coast_path <- file.path("multiyear", "SHARK_IFCB_2022_2024_Skagerrak-Kattegat_SMHI")

# Read data
data_baltic <- read_tsv(file.path(ifcb_path, "shark", baltic_path, "processed_data", "data.txt"))
data_west_coast <- read_tsv(file.path(ifcb_path, "shark", west_coast_path, "processed_data", "data.txt"))

# Bind data together and remove manual annotations
ifcb_data <- rbind(data_baltic, data_west_coast) %>%
  filter(IMAGE_VERIFICATION == "PredictedByMachine")

# Download Phytoplankton data from SHARK
shark_data <- get_shark_data(tableView = "sharkdata_phytoplankton",
                             toYear = 2024,
                             fromYear = 2022,
                             dataTypes = "Phytoplankton")


# Identify disitinct IFCB samples
ifcb_samples <- ifcb_data %>%
  select(LATIT, LONGI, SDATE, SMPNO) %>%
  distinct()

# Identify disitinct Microscopy samples
microscopy_samples <- shark_data %>%
  select(sample_date, sample_longitude_dd, sample_latitude_dd) %>%
  distinct()

# Perform an inner join to match IFCB and microscopy samples by date
combined_samples <- ifcb_samples %>%
  inner_join(microscopy_samples, by = c("SDATE" = "sample_date"))

# Calculate the distance between each pair of joined samples
combined_samples <- combined_samples %>%
  rowwise() %>%
  mutate(
    distance_km = distHaversine(
      c(LONGI, LATIT),  # IFCB sample coordinates
      c(sample_longitude_dd, sample_latitude_dd)  # Microscopy sample coordinates
    ) / 1000  # Convert meters to kilometers
  ) %>%
  ungroup() %>%
  group_by(sample_longitude_dd, sample_latitude_dd) %>%
  mutate(shark_id = uuid::UUIDgenerate()) %>%
  ungroup()

# Filter to keep only pairs within 5 km
nearby_samples <- combined_samples %>%
  filter(distance_km <= 5) %>%
  select(sample_longitude_dd, sample_latitude_dd, SDATE, shark_id) %>%
  mutate(near_ifcb = TRUE) %>%
  rename(sample_date = SDATE) %>%
  distinct()

# Filter to keep only pairs within 5 km
nearby_ifcb_samples <- combined_samples %>%
  filter(distance_km <= 5) %>%
  select(SMPNO, distance_km, shark_id) %>%
  mutate(near_microscopy = TRUE) %>%
  distinct()

# Join microscopy data and filter samples that are sampled near IFCB samples
filtered_shark_data <- shark_data %>%
  left_join(nearby_samples, by = c("sample_longitude_dd", "sample_latitude_dd", "sample_date")) %>%
  filter(near_ifcb) %>%
  filter(scientific_name %in% ifcb_data$LATNM) %>%
  filter(!station_name == "FARSTAVIKEN") %>%  # Remove the Calluna station
  filter(!visit_year == 2022)

# Aggregate size classes and transform measurements into wide format
shark_data_wide <- filtered_shark_data %>%
  select(shark_id, visit_year, sample_date, sample_latitude_dd, sample_longitude_dd, scientific_name, size_class, parameter, value) %>%
  group_by(shark_id, visit_year, sample_date, sample_latitude_dd, sample_longitude_dd, scientific_name, parameter) %>%
  summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  select(-`# counted`)

# Join IFCB data and samples that are sampled near microscopy samples
filtered_ifcb_data <- ifcb_data %>%
  left_join(nearby_ifcb_samples, by = "SMPNO") %>%
  filter(near_microscopy) %>%
  filter(!MYEAR == 2022) %>%
  filter(LATNM %in% filtered_shark_data$scientific_name)

# Calculate average values for multiple IFCB samples near a microscopy sample
average_ifcb_data <- filtered_ifcb_data %>%
  group_by(shark_id, LATNM) %>%
  summarise(LATIT = mean(LATIT, na.rm = TRUE),
            LONGI = mean(LONGI, na.rm = TRUE),
            ABUND = mean(ABUND, na.rm = TRUE),
            BIOVOL = mean(BIOVOL, na.rm = TRUE),
            C_CONC = mean(C_CONC, na.rm = TRUE)) %>%
  ungroup()
  
# Join IFCB and microscopy data in a single dataframe
joined_data <- average_ifcb_data %>%
  rename(scientific_name = LATNM) %>%
  left_join(shark_data_wide) %>%
  filter(scientific_name %in% c("Aphanizomenon", "Dolichospermum", "Nodularia spumigena")) %>%
  filter(!is.na(sample_date))

# Scatter plot with linear regression line, faceted by scientific_name
regression_plot <- ggplot(joined_data, aes(x = `Carbon concentration`, y = C_CONC)) +
  # Scatter points with reduced transparency and appropriate size
  geom_point(alpha = 0.7, size = 2, color = "darkblue") +
  # Linear regression line with confidence interval
  # geom_smooth(method = "lm", color = "red", fill = "lightpink", se = TRUE, linewidth = 0.8) +
  # Facet by scientific name, wrapping nicely for many panels
  facet_wrap(~ scientific_name, scales = "free", ncol = 4) +
  # Labels with improved units
  labs(
    x = "Microscopy-Derived Carbon Concentration (µg C/L)",
    y = "IFCB-Derived Carbon Concentration (µg C/L)",
    # title = "Comparison of Microscopy vs. IFCB Carbon Concentrations",
    # subtitle = "Faceted by Taxonomic Groups"
  ) +
  # Aesthetic improvements for minimal yet elegant theme
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),  # Remove minor gridlines for clarity
    panel.grid.major = element_line(color = "gray90"),  # Subtle gridlines
    strip.text = element_text(size = 12, face = "bold", margin = margin(b = 5)),  # Emphasize facet labels
    axis.text = element_text(size = 12, color = "black"),  # Enhance axis text readability
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles for emphasis
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered plot title
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),  # Centered subtitle
    legend.position = "none"  # Remove legend if unnecessary
  )

# Save the plot as a PNG
ggsave(filename = "plots/microscopy_comparison_cyanbacteria.png", 
       plot = regression_plot, 
       width = 12, height = 10, dpi = 300,
       bg = "white")
