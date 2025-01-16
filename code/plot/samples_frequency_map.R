library(tidyverse)
library(sf)
library(ggOceanMaps)
library(cowplot) 

# Define paths
ifcb_path <- Sys.getenv("ifcb_path")
baltic_path <- file.path("multiyear", "SHARK_IFCB_2022_2024_Baltic_SMHI")
west_coast_path <- file.path("multiyear", "SHARK_IFCB_2022_2024_Skagerrak-Kattegat_SMHI")

# Read data
data_baltic <- read_tsv(file.path(ifcb_path, "shark", baltic_path, "processed_data", "data.txt"))
data_west_coast <- read_tsv(file.path(ifcb_path, "shark", west_coast_path, "processed_data", "data.txt"))

# Select position for the baltic dataset
samples_baltic <- data_baltic %>%
  select(LATIT, LONGI, SMPNO, SDATE) %>%
  distinct()

# Select position for the Skagerrak-Kattegat dataset
samples_west_coast <- data_west_coast %>%
  select(LATIT, LONGI, SMPNO, SDATE) %>%
  distinct()

# Bind samples 
samples <- rbind(samples_baltic, samples_west_coast)

# Convert to an sf object
samples_sf <- st_as_sf(samples, coords = c("LONGI", "LATIT"), crs = 4326) %>% 
  mutate(sample = TRUE)

# Define grid
grid <- st_make_grid(samples_sf, cellsize = c(0.3, 0.25), what = "polygons")
grid_sf <- st_sf(grid_id = seq_along(grid), geometry = grid)

# Perform spatial join and count frequencies
grid_counts <- st_join(grid_sf, samples_sf, join = st_intersects) %>%
  filter(sample) %>%  # Keep grid_ids with a decimal point
  group_by(grid_id) %>%
  summarise(frequency = n(), .groups = "drop")  # Count the number of matches

# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(samples$LONGI) - 1, max(samples$LONGI) + 1, min(samples$LATIT) - 1, max(samples$LATIT) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

# Visualize the gridded data
map <- baltic_sea_map +
  geom_sf(data = grid_counts, aes(fill = frequency), colour = NA) +  # No border for grid cells
  scale_fill_viridis_c(option = "plasma", na.value = "white") +     # Set NA to white
  # borders("world", xlim = c(5, 25), ylim = c(54, 66), fill = "gray80", colour = "black") +
  # coord_sf(xlim = c(5, 25), ylim = c(54, 66), expand = FALSE) +
  labs(
    # title = "Gridded Sampling Frequency in the Baltic Sea",
    fill = "Sampling frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.title.position = "bottom",
        legend.direction = "horizontal")

# Reorder layers
map <- reorder_layers(map)


### Sampling frequency histogram

# Set locale to English
Sys.setlocale("LC_TIME", "en_US.UTF-8")

# Extract month from SDATE and count sampling frequencies
monthly_counts <- samples %>%
  mutate(month = month(SDATE, label = TRUE, abbr = TRUE)) %>% # Extract month as a labeled factor
  count(month, name = "frequency")                           # Count the occurrences per month

month_palette <- c(
  "Dec" = "#4682B4",  # Steel Blue
  "Jan" = "#B0E0E6",  # Powder Blue
  "Feb" = "#87CEFA",  # Light Sky Blue
  "Mar" = "#98FB98",  # Pale Green
  "Apr" = "#00FA9A",  # Medium Spring Green
  "May" = "#7CFC00",  # Lawn Green
  "Jun" = "#FFD700",  # Gold
  "Jul" = "#FFA500",  # Orange
  "Aug" = "#FF6347",  # Tomato
  "Sep" = "#FF8C00",  # Dark Orange
  "Oct" = "#CD853F",  # Peru
  "Nov" = "#8B4513"   # Saddle Brown
)

# Convert to factor to ensure ordering
monthly_counts <- monthly_counts %>%
  mutate(month = factor(month, levels = month.abb))

# Add colors to ggplot
histogram <- ggplot(monthly_counts, aes(x = month, y = frequency, fill = month)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = month_palette) +
  labs(x = "", y = "Sampling Frequency") +
  theme_minimal()

# Combine the two plots into a panel
combined_panel <- plot_grid(map, histogram, ncol = 2, 
                            # labels = c("A", "B"), 
                            rel_widths = c(1,1), rel_heights = c(1,1), align = "hv", scale = c(1,.73))

# Save the plot as a PNG
ggsave(filename = "plots/sampling_frequency_map.png", 
       plot = map, 
       width = 12, height = 10, dpi = 300,
       bg = "white")

# Save the plot as a PNG
ggsave(filename = "plots/sampling_frequency_panel.png", 
       plot = combined_panel, 
       width = 12, height = 10, dpi = 300,
       bg = "white")
