library(tidyverse)
library(worrms)
library(ggOceanMaps)
library(maps)


tångesund <- read_tsv("data/SHARK_IFCB_2016_Tångesund_SMHI/processed_data/data.txt")
baltic <- read_tsv("data/SHARK_IFCB_2022_2024_Baltic_SMHI/processed_data/data.txt")
sk <- read_tsv("data/SHARK_IFCB_2022_2024_Skagerrak-Kattegat_SMHI/processed_data/data.txt")

all_data <- rbind(tångesund, baltic, sk)
  
aphia_id <- unique(all_data$APHIA_ID)

aphia_id <- aphia_id[!is.na(aphia_id)]

worms_records <- wm_record(aphia_id)

all_data <- all_data %>%
  left_join(worms_records, by = c("APHIA_ID" = "AphiaID"))

# Count the number of rows for each class
class_counts <- all_data %>%
  mutate(class = ifelse(CLASS_NAME == "Rhizosolenia_Pseudosolenia", "Bacillariophyceae", class)) %>%
  mutate(class = ifelse(CLASS_NAME == "Cylindrotheca_closterium_Nitzschia_longissima", "Bacillariophyceae", class)) %>%
  mutate(class = ifelse(CLASS_NAME == "Leptocylindrus_danicus_Leptocylindrus_minimus", "Bacillariophyceae", class)) %>%
  mutate(class = ifelse(CLASS_NAME == "Gyrosigma_Pleurosigma", "Bacillariophyceae", class)) %>%
  mutate(class = ifelse(CLASS_NAME == "Heterocapsa_Azadinium", "Dinophyceae", class)) %>%
  mutate(class = ifelse(CLASS_NAME == "Snowella_Woronichinia", "Cyanophyceae", class)) %>%
  filter(!is.na(class)) %>%
  group_by(class) %>%
  summarise(count = n())

# Create the bar graph
occurences <- ggplot(class_counts, aes(x = reorder(class, count), y = count, fill = class)) +
  geom_bar(stat = "identity", color = "black") +  # Black outlines for clarity
  coord_flip() +  # Flip the coordinates for horizontal bars
  # scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +  # Use a 'plasma' color palette
  labs(x = "", y = "Number of Occurences") +
  theme_minimal(base_size = 16) +  # Increase base font size for poster readability
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis labels
    axis.title = element_text(size = 18, face = "bold"),  # Increase axis title size and bold
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Center and enlarge plot title
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Add margins around the plot
    legend.position = "none",  # Remove the legend to declutter the graph
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # Add space at the top of the bars

ggsave("plots/occurences.png", occurences, device = "png")

ferry <- all_data %>%
  filter(!MYEAR == 2016)

length(unique(ferry$SMPNO))

# Extract unique LATIT and LONGI points
unique_points <- all_data %>%
  select(LATIT, LONGI) %>%
  distinct()

# Get the range of latitude and longitude for zooming
lat_range <- range(unique_points$LATIT, na.rm = TRUE)
long_range <- range(unique_points$LONGI, na.rm = TRUE)

# Create a world map using the maps package
world_map <- map_data("world")

# Plot the map and zoom into the area of interest
ggplot() +
  # Draw the world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "black") +
  # Overlay the unique data points
  geom_point(data = unique_points, aes(x = LONGI, y = LATIT), color = "red", size = 1.5, alpha = 0.7) +
  labs(x = "Longitude", y = "Latitude", title = "Zoomed-In Sampling Points") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  ) +
  coord_quickmap(xlim = long_range, ylim = lat_range)  # Zoom into the area where the points are located


# Create basemap
baltic_sea_map <- basemap(
  limits = c(min(unique_points$LONGI) - 1, max(unique_points$LONGI) + 1, min(unique_points$LATIT) - 1, max(unique_points$LATIT) + 1),
  land.col = "#eeeac4",
  land.border.col = "black",
  rotate = TRUE,
  bathymetry = FALSE
)

map <- baltic_sea_map + 
  geom_spatial_point(aes(x=unique_points$LONGI, y =unique_points$LATIT), color = "red") +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.ontop = FALSE) 

ggsave("map.png", map, device = "png")
