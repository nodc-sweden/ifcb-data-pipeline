library(SHARK4R)
library(tidyverse)
library(geosphere)
library(ggOceanMaps)
library(iRfcb)
library(ggh4x)

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
  mutate(in_baltic = ifcb_is_in_basin(sample_latitude_dd, sample_longitude_dd)) %>%
  filter(near_ifcb) %>%
  filter(scientific_name %in% ifcb_data$LATNM) %>%
  filter(!station_name == "FARSTAVIKEN") %>%  # Remove the Calluna station
  filter(!visit_year == 2022) %>%
  filter(in_baltic)
# %>%
#   filter(!sample_project_name_sv == "PROJ Utsjöprovtagningar Algprognos och Algtoxiner")

# Aggregate size classes and transform measurements into wide format
shark_data_wide <- filtered_shark_data %>%
  select(shark_id, visit_year, sample_date, sampler_type_code, sample_latitude_dd, sample_longitude_dd, scientific_name, size_class, parameter, value) %>%
  group_by(shark_id, visit_year, sample_date, sampler_type_code, sample_latitude_dd, sample_longitude_dd, scientific_name, parameter) %>%
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
  full_join(shark_data_wide) %>%
  group_by(shark_id) %>%
  fill(visit_year, sample_date, sample_latitude_dd, sample_longitude_dd, .direction = "downup") %>%
  ungroup() %>%
  filter(scientific_name %in% c("Aphanizomenon", "Dolichospermum", "Nodularia spumigena")) %>%
  filter(!is.na(sample_date)) %>%
  mutate(sampler_type_code = if_else(is.na(sampler_type_code), "NSK", sampler_type_code)) %>%
  mutate(LATIT = if_else(is.na(LATIT), sample_latitude_dd, LATIT),
         LONGI = if_else(is.na(LONGI), sample_longitude_dd, LONGI))

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

# Restructure data for sample_type
joined_data_long <- joined_data %>%
  select(LATIT, LONGI, scientific_name, sampler_type_code, sample_date, C_CONC, `Carbon concentration`) %>%
  pivot_longer(
    cols = c(C_CONC, `Carbon concentration`),
    names_to = "sample_type",
    values_to = "carbon_concentration"
  ) %>%
  mutate(
    sample_type = case_when(
      sample_type == "C_CONC" ~ "IFCB",
      sample_type == "Carbon concentration" ~ "Microscopy"
    ),
    month = format(sample_date, "%B"),  # Extract month name
    month = factor(
      month,
      levels = c("January", "February", "March", "April", "May", "June",
                 "July", "August", "September", "October", "November", "December")
    )  # Set custom order for months
  ) %>%
  filter(!month %in% c("March", "October", "September")) %>%
  mutate(sample_depth = if_else(sampler_type_code == "NSK", "Surface", "Integrated 0-10 or 0-20 m")) %>%
  mutate(sample_depth = if_else(sample_type == "IFCB", "Surface", sample_depth)) %>%
  distinct() %>%
  filter(!is.na(carbon_concentration))

# Convert the data to an sf object
joined_data_sf <- joined_data_long %>%
  st_as_sf(coords = c("LONGI", "LATIT"), crs = 4326)  # Use WGS84 CRS (EPSG:4326)

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(joined_data_long$LONGI) - 1, max(joined_data_long$LONGI) + 1, min(joined_data_long$LATIT) - 1, max(joined_data_long$LATIT) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

# Create the faceted map with modified order
faceted_map <- baltic_sea_map +
  geom_sf(
    data = joined_data_sf,
    aes(
      size = carbon_concentration,
      color = sample_depth,         # Color for the border
      fill = sample_type           # Fill for the interior of the points
    ),
    alpha = 0.7,                   # Adjust the transparency for both fill and border
    shape = 21                     # Set shape to make the points have borders
  ) +
  facet_nested(scientific_name ~ month + sample_type,       # Facet columns by month and sample_type
               nest_line = element_line(colour = "red"),
  ) +
  scale_fill_manual(
    values = c("IFCB" = "#1f77b4",    # Soft blue (for IFCB)
               "Microscopy" = "#ff7f0e"),  # Orange (for Microscopy)
    guide = "none"
  ) +
  scale_color_manual(
    values = c("Integrated 0-10 or 0-20 m" = "#2ca02c",    # Green (for Hose)
               "Surface" = "#d62728")  # Red (for Surface)
  ) +
  labs(
    size = "Carbon Concentration (µg C/L)",
    fill = "Sample Type",
    color = "Sample Depth"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),  # Customize facet labels
    axis.title = element_text(size = 12),                # Axis titles
    axis.text = element_text(size = 10),                 # Axis text
    legend.title = element_text(size = 12),              # Legend title
    legend.text = element_text(size = 10),               # Legend text
    # panel.spacing = unit(1, "lines"),                    # Adjust spacing between facets
    legend.position = "bottom",                             # Move legend to the top
    legend.title.position = "top",                    # Adjust legend title position
    legend.direction = "horizontal"                      # Horizontal legend layout
  )

# Save the plot as a PNG
ggsave(filename = "plots/microscopy_comparison_maps.png", 
       plot = faceted_map, 
       width = 10, height = 8, dpi = 300,
       bg = "white")
