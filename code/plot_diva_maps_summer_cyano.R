library(tidyverse)
library(JuliaCall)
library(rnaturalearth)
library(sf)
library(oce)

# Set up Julia
if (Sys.info()["sysname"] == "Windows") {
  julia_setup(JULIA_HOME = "C:/Julia-1.10.4/bin/")
} else {
  julia_setup(JULIA_HOME = "/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
}
julia_install_package("DIVAnd")
julia_command("using DIVAnd")

# Define paths
data_path <- file.path("multiyear", "SHARK_IFCB_2022_2024_Baltic_SMHI")
ifcb_path <- Sys.getenv("ifcb_path")

# Define taxa
selected_taxa <- c("Nodularia spumigena", "Aphanizomenon", "Dolichospermum")

# Read data
data <- read_tsv(file.path(ifcb_path, "shark", data_path, "processed_data", "data.txt"))

# Filter relevant data
data_filtered <- data %>%
  filter(LATNM %in% selected_taxa) %>%
  filter(IMAGE_VERIFICATION == "PredictedByMachine") %>%
  mutate(YEAR_CRUISE = paste(MYEAR, CRUISE_NO, sep = "_")) %>%
  filter(!is.na(CRUISE_NO))

# Only use SMHI cruises from June, July and August 2023-2024
data_filtered <- data_filtered %>%
  filter(YEAR_CRUISE %in% c("2023_011", "2023_012", "2023_013", "2024_013", "2024_016", "2024_017"))

# Get land polygons from Natural Earth
land_polygons <- ne_countries(scale = "medium", returnclass = "sf")

# To be populated in loop
interpolated_list <- list()
sampling_points <- list()
months <- list()

for (i in seq_along(unique(data_filtered$YEAR_CRUISE))) {
  
  data_cruise <- data_filtered %>%
    filter(YEAR_CRUISE == unique(data_filtered$YEAR_CRUISE)[i])
  
  coordinates <- data_cruise %>%
    select(LATIT, LONGI) %>%
    distinct()
  
  month <- unique(month(data_cruise$SDATE))
  
  sampling_points[[unique(data_filtered$YEAR_CRUISE)[i]]] <- coordinates
  months[[unique(data_filtered$YEAR_CRUISE)[i]]] <- month
  
  for (j in seq_along(unique(data_cruise$LATNM))) {
    data_cruise_taxa <- data_cruise %>%
      filter(LATNM == unique(data_cruise$LATNM)[j])
    
    # Assign the dimensions
    x <- data_cruise_taxa$LONGI
    y <- data_cruise_taxa$LATIT
    f <- data_cruise_taxa$C_CONC
    
    # Set the interpolation grid
    NX <- 100
    NY <- 100

    xx <- seq(11.8, 22, length.out=NX)
    yy <- seq(53.5, 60, length.out=NY)
    
    # Assign values in julia
    julia_assign("xx", xx)
    julia_assign("yy", yy)
    
    julia_command("xi, yi = ndgrid(xx, yy)");
    xi = julia_eval("xi")
    yi = julia_eval("yi")

    # Create a grid of spatial points from xx and yy
    grid_points <- expand.grid(x = xx, y = yy)
    
    # Convert grid to sf object with the correct CRS
    grid_points_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = st_crs(land_polygons))
    
    # Check which grid points intersect the land polygons
    land_mask <- st_intersects(grid_points_sf, land_polygons, sparse = FALSE)
    
    # 'land_mask' is a logical vector indicating which points are on land (TRUE means it's on land)
    land_mask <- apply(land_mask, 1, any)
    
    # Apply the additional condition: Latitude > 60 and Longitude < 15
    additional_mask <- with(grid_points, y > 55.6 & x < 13)

    # Combine the two masks: land and the additional conditions
    combined_mask <- land_mask | additional_mask  # TRUE if it's either on land or meets the condition
    
    # Reshape the combined_mask vector into a matrix with dimensions NX by NY
    mask <- matrix(!combined_mask, nrow = NX, ncol = NY, byrow = FALSE)
    
    # Calculate the resolution (spacing) for longitude, latitude, and time
    dx <- xi[2,1] - xi[1,1]  # Longitude resolution
    dy <- yi[1,2] - yi[1,1]  # Latitude resolution

    # Metrics: pm and pn for spatial dimensions
    pm <- matrix(1, NX, NY) / dx  # Longitude metric
    pn <- matrix(1, NX, NY) / dy  # Latitude metric
    
    # Assign metrics to Julia
    julia_assign("pm", pm)
    julia_assign("pn", pn)

    # Analysis parameters
    # Correlation length
    len <- .75
    
    # Obs. error variance normalized by the background error variance
    epsilon2 <- 1.0
    
    # From R variable to Julia variable
    julia_assign("xi", xi)
    julia_assign("yi", yi)
    julia_assign("x", x)
    julia_assign("y", y)
    julia_assign("f", f)
    julia_assign("len", len)
    julia_assign("mask", mask)
    julia_assign("epsilon2", epsilon2)
    
    tryCatch(
      {
        # Run the DIVAnd interpolation in 2D (longitude, latitude)
        julia_command("fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2);")
        
        # Call DIVAnd_errormap to get the error map
        julia_command("e, errormap = DIVAnd_errormap(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2, s; method = :cheap, Bscale = false);")
      
        fi <- julia_eval("fi")
        e <- julia_eval("e")
        
        interpolated_list[[paste(unique(data_cruise_taxa$YEAR_CRUISE), unique(data_cruise_taxa$LATNM), sep = ";")]] <- fi
        },
        error = function(e) {
        # Print the error message (optional)
        message("Error occurred: ", e$message)
        
        # Skip to the next iteration
        next
      }
    )
  }
}

# Extract cruise numbers
cruise_numbers <- sub(";.*", "", names(interpolated_list))

# Group by cruise number
grouped_list <- split(interpolated_list, cruise_numbers)

# Loop through each cruise number and its associated taxa
plots <- list()

# Extract year from names
years <- sub("_.*", "", names(interpolated_list))

# Group interpolated_list by year
grouped_by_year <- split(interpolated_list, years)

# Pre-calculate global zlim
global_zlim <- range(unlist(lapply(grouped_by_year, as.vector)), na.rm = TRUE)

# Pre-calculate zlim for each year
year_zlim <- lapply(grouped_by_year, function(year_data) {
  # Flatten all matrices for the year and compute range
  all_values <- unlist(lapply(year_data, as.vector))
  range(all_values, na.rm = TRUE)
})

for (cruise_number in names(grouped_list)) {
  # Extract year for the cruise number
  year <- sub("_.*", "", cruise_number)
  
  # Get the zlim for this year
  zlim <- year_zlim[[year]]
  
  # Extract taxa list for the current cruise
  taxa_list <- grouped_list[[cruise_number]]
  
  # Extract sampling points for the current cruise
  sampling_data <- sampling_points[[cruise_number]]
  
  # Extract month information
  month <- months[[cruise_number]]
  if (is.null(sampling_data)) next  # Skip if no sampling points for this cruise
  
  # Create a list to store data frames for each taxon
  taxa_data <- list()
  
  for (taxon_name in names(taxa_list)) {
    # Get the 2D matrix for the current taxon
    taxon_matrix <- taxa_list[[taxon_name]]
    
    # Skip if the matrix is empty or full of NaN
    if (all(is.nan(taxon_matrix))) next
    
    # Convert the matrix to a data frame for plotting
    interp_grid <- expand.grid(
      lon = xx,
      lat = yy
    )
    interp_grid$fi <- as.vector(taxon_matrix)
    
    # Extract cruise numbers
    taxa <- sub("^[^;]*;", "", taxon_name)
    
    # Add the taxon name for faceting
    interp_grid$taxon <- taxa
    
    # Append to the taxa data list
    taxa_data[[taxon_name]] <- interp_grid
  }
  
  # Combine all taxa data into a single data frame
  if (length(taxa_data) == 0) next  # Skip if no valid data
  cruise_data <- bind_rows(taxa_data)
  
  # Define color palette and limits
  pal <- oce.colorsTemperature(100)

  # Generate the faceted plot
  faceted_plot <- ggplot() +
    # Add raster layer for the interpolated grid
    geom_raster(data = cruise_data, aes(x = lon, y = lat, fill = fi)) +
    # Add sampling points as a point layer
    geom_point(data = sampling_data, aes(x = LONGI, y = LATIT), 
               size = 1.5, shape = 21, fill = alpha("white", .5), stroke = 0) +
    # Add land polygons from rnaturalearth
    geom_sf(data = land_polygons, fill = "#eeeac4", color = "black") +
    # Set color scale for interpolation
    scale_fill_gradientn(colors = pal, limits = global_zlim, na.value = "transparent", name = "Biomass (µg C/L)") +
    # Adjust the plot limits based on your longitude and latitude ranges
    coord_sf(xlim = c(min(cruise_data$lon, na.rm = TRUE) - 1, max(cruise_data$lon, na.rm = TRUE) + 1), 
             ylim = c(min(cruise_data$lat, na.rm = TRUE) - 1, max(cruise_data$lat, na.rm = TRUE) + 1), 
             expand = FALSE) +
    # Add labels and titles
    labs(
      title = paste0("Cyanobacterial biomass, cruise: ", cruise_number, ", month(s): ", paste(month.name[month], collapse = ",")),
      x = "Longitude", 
      y = "Latitude"
    ) +
    # Create facets for each taxon
    facet_wrap(~ taxon) +
    # Use a minimal theme
    theme_minimal()
  
  # Store the plot
  plots[[cruise_number]] <- faceted_plot
}

# Create a directory to save the plots (if it doesn't already exist)
output_dir <- "plots/summer"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through the plots and save each one
for (cruise_number in names(plots)) {
  # Define the output file name
  output_file <- file.path(output_dir, paste0("plot_", cruise_number, ".png"))
  
  # Save the plot as a PNG
  ggsave(filename = output_file, 
         plot = plots[[cruise_number]], 
         width = 10, height = 8, dpi = 300,
         bg = "white")
}

# Read bucket sample data from planktologist onboard Svea July 2024
bucket_samples <- read_tsv("data/bucket_samples_july_2024.txt")

# Prepare the data
bucket_samples_long <- bucket_samples %>%
  pivot_longer(
    cols = c(Aphanizomenon, `Nodularia spumigena`, Dolichospermum),
    names_to = "Taxon",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))  # Remove rows with NA values

# Convert `Value` into a categorical variable
bucket_samples_long <- bucket_samples_long %>%
  mutate(
    Abundance = case_when(
      Value == 1 ~ "Present",
      Value == 2 ~ "Present (2)",
      Value == 3 ~ "Common",
      Value == 4 ~ "Very Common",
      TRUE ~ NA_character_
    ),
    Abundance = factor(Abundance, levels = c("Present", "Present (2)", "Common", "Very Common"))
  )

# Plot
bucket_plot <- ggplot(bucket_samples_long, aes(x = Lon, y = Lat)) +
  # Add points with size representing Abundance
  geom_point(aes(size = Abundance, color = Taxon), shape = 21, alpha = 0.8) +
  # Add land polygons from rnaturalearth
  geom_sf(data = land_polygons, fill = "#eeeac4", color = "black", inherit.aes = FALSE) +
  # Customize scale for size and color
  scale_size_manual(
    name = "Abundance class",
    values = c("Present" = 1, "Present (2)" = 2, "Common" = 4, "Very Common" = 6)
  ) +
  scale_color_manual(
    name = "Taxon",
    values = c("Aphanizomenon" = "blue", 
               "Nodularia spumigena" = "red", 
               "Dolichospermum" = "green")
  ) +
  # Add labels and titles
  labs(
    title = "Cyanobacterial presence (microscopy), cruise: 2024_016, month(s): July",
    x = "Longitude",
    y = "Latitude"
  ) +
  # Adjust the coordinate system
  coord_sf(xlim = c(min(cruise_data$lon, na.rm = TRUE) - 1, max(cruise_data$lon, na.rm = TRUE) + 1), 
           ylim = c(min(cruise_data$lat, na.rm = TRUE) - 1, max(cruise_data$lat, na.rm = TRUE) + 1), 
           expand = FALSE) +
  # Use faceting for separate plots by Taxon
  facet_wrap(~ Taxon) +
  # Use a minimal theme
  theme_minimal() +
  theme(legend.position = "right")

# Print the plot
bucket_plot

# Save the plot as a PNG
ggsave(filename = "plots/summer/plot_2024_016_bucket.png", 
       plot = bucket_plot, 
       width = 10, height = 8, dpi = 300,
       bg = "white")



### Create aggregated cyanobacteria maps

# Aggregate data
data_aggregated <- data_filtered %>%
  group_by(YEAR_CRUISE, SDATE, LONGI, LATIT, SMPNO) %>%
  summarise(C_CONC = sum(C_CONC, na.rm = TRUE)) %>%
  ungroup()

# To be populated in loop
interpolated_list <- list()
sampling_points <- list()
months <- list()
years <- list()

for (i in seq_along(unique(data_aggregated$YEAR_CRUISE))) {
  
  data_cruise <- data_aggregated %>%
    filter(YEAR_CRUISE == unique(data_aggregated$YEAR_CRUISE)[i])
  
  coordinates <- data_cruise %>%
    select(LATIT, LONGI) %>%
    distinct()
  
  month <- unique(month(data_cruise$SDATE))
  year <- unique(year(data_cruise$SDATE))
  
  sampling_points[[unique(data_aggregated$YEAR_CRUISE)[i]]] <- coordinates
  months[[unique(data_aggregated$YEAR_CRUISE)[i]]] <- month
  years[[unique(data_aggregated$YEAR_CRUISE)[i]]] <- year
  
  # Assign the dimensions
  x <- data_cruise$LONGI
  y <- data_cruise$LATIT
  f <- data_cruise$C_CONC
  
  # Set the interpolation grid
  NX <- 100
  NY <- 100
  
  xx <- seq(11.8, 22, length.out=NX)
  yy <- seq(53.5, 60, length.out=NY)
  
  # Assign values in julia
  julia_assign("xx", xx)
  julia_assign("yy", yy)
  
  julia_command("xi, yi = ndgrid(xx, yy)");
  xi = julia_eval("xi")
  yi = julia_eval("yi")
  
  # Create a grid of spatial points from xx and yy
  grid_points <- expand.grid(x = xx, y = yy)
  
  # Convert grid to sf object with the correct CRS
  grid_points_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = st_crs(land_polygons))
  
  # Check which grid points intersect the land polygons
  land_mask <- st_intersects(grid_points_sf, land_polygons, sparse = FALSE)
  
  # 'land_mask' is a logical vector indicating which points are on land (TRUE means it's on land)
  land_mask <- apply(land_mask, 1, any)
  
  # Apply the additional condition: Latitude > 60 and Longitude < 15
  additional_mask <- with(grid_points, y > 55.6 & x < 13)
  
  # Combine the two masks: land and the additional conditions
  combined_mask <- land_mask | additional_mask  # TRUE if it's either on land or meets the condition
  
  # Reshape the combined_mask vector into a matrix with dimensions NX by NY
  mask <- matrix(!combined_mask, nrow = NX, ncol = NY, byrow = FALSE)
  
  # Calculate the resolution (spacing) for longitude, latitude, and time
  dx <- xi[2,1] - xi[1,1]  # Longitude resolution
  dy <- yi[1,2] - yi[1,1]  # Latitude resolution
  
  # Metrics: pm and pn for spatial dimensions
  pm <- matrix(1, NX, NY) / dx  # Longitude metric
  pn <- matrix(1, NX, NY) / dy  # Latitude metric
  
  # Assign metrics to Julia
  julia_assign("pm", pm)
  julia_assign("pn", pn)
  
  # Analysis parameters
  # Correlation length
  len <- .75
  
  # Obs. error variance normalized by the background error variance
  epsilon2 <- 1.0
  
  # From R variable to Julia variable
  julia_assign("xi", xi)
  julia_assign("yi", yi)
  julia_assign("x", x)
  julia_assign("y", y)
  julia_assign("f", f)
  julia_assign("len", len)
  julia_assign("mask", mask)
  julia_assign("epsilon2", epsilon2)
  
  tryCatch(
    {
      # Run the DIVAnd interpolation in 2D (longitude, latitude)
      julia_command("fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2);")
      
      # Call DIVAnd_errormap to get the error map
      julia_command("e, errormap = DIVAnd_errormap(mask, (pm, pn), (xi, yi), (x, y), f, len, epsilon2, s; method = :cheap, Bscale = false);")
      
      fi <- julia_eval("fi")
      e <- julia_eval("e")
      
      interpolated_list[[unique(data_cruise$YEAR_CRUISE)]] <- fi
    },
    error = function(e) {
      # Print the error message (optional)
      message("Error occurred: ", e$message)
      
      # Skip to the next iteration
      next
    }
  )
}

# Extract year from names
years <- sub("_.*", "", names(interpolated_list))

# Group interpolated_list by year
grouped_by_year <- split(interpolated_list, years)

# Pre-calculate global zlim
global_zlim <- range(unlist(lapply(grouped_by_year, as.vector)), na.rm = TRUE)

# Pre-calculate zlim for each year
year_zlim <- lapply(grouped_by_year, function(year_data) {
  # Flatten all matrices for the year and compute range
  all_values <- unlist(lapply(year_data, as.vector))
  range(all_values, na.rm = TRUE)
})

# Initialize a list for storing plots
yearly_plots <- list()

for (year in names(grouped_by_year)) {
  # Extract cruises for the year
  cruises <- grouped_by_year[[year]]
  
  # Get the zlim for this year
  zlim <- year_zlim[[year]]
  
  # Prepare data for plotting
  year_data <- list()
  sampling_points_year <- list()
  
  for (cruise_number in names(cruises)) {
    # Extract the matrix for the current cruise
    cruise_matrix <- cruises[[cruise_number]]
    
    # Get the month name for the cruise
    month <- month.name[months[[cruise_number]]]
    
    # Skip if the matrix is empty or full of NaN
    if (all(is.nan(cruise_matrix))) next
    
    # Convert the matrix to a data frame for plotting
    interp_grid <- expand.grid(
      lon = xx,
      lat = yy
    )
    interp_grid$fi <- as.vector(cruise_matrix)
    
    # Add cruise identifier for faceting
    interp_grid$cruise <- paste0(cruise_number, ": ", month)
    
    # Append to the year's data
    year_data[[cruise_number]] <- interp_grid
    
    # Collect sampling points for the cruise
    if (!is.null(sampling_points[[cruise_number]])) {
      sampling_points_year[[cruise_number]] <- sampling_points[[cruise_number]] %>%
        mutate(
          cruise = paste0(cruise_number, ": ", month),  # Add cruise identifier
          month = month            # Add month name
        )
    }
  }
  
  # Combine all cruise data into a single data frame
  if (length(year_data) == 0) next  # Skip if no valid data
  combined_data <- bind_rows(year_data)
  
  # Combine sampling points for all cruises
  combined_sampling <- bind_rows(
    sampling_points_year
  )
  
  # Define color palette
  pal <- oce.colorsTemperature(100)
  
  # Generate the faceted plot for the year
  yearly_plot <- ggplot() +
    # Add raster layer for the interpolated grid
    geom_raster(data = combined_data, aes(x = lon, y = lat, fill = fi)) +
    # Add sampling points as a point layer
    geom_point(data = combined_sampling, aes(x = LONGI, y = LATIT, group = cruise), 
               size = 1.5, shape = 21, fill = alpha("white", .5), stroke = 0) +
    # Add land polygons from rnaturalearth
    geom_sf(data = land_polygons, fill = "#eeeac4", color = "black") +
    # Set color scale for interpolation
    scale_fill_gradientn(colors = pal, limits = global_zlim, na.value = "transparent", name = "Biomass (µg C/L)") +
    # Adjust the plot limits based on your longitude and latitude ranges
    coord_sf(xlim = c(min(combined_data$lon, na.rm = TRUE) - 1, max(combined_data$lon, na.rm = TRUE) + 1), 
             ylim = c(min(combined_data$lat, na.rm = TRUE) - 1, max(combined_data$lat, na.rm = TRUE) + 1), 
             expand = FALSE) +
    # Add labels and titles
    labs(
      title = paste0("Surface cyanobacterial biomass for year: ", year),
      x = "Longitude", 
      y = "Latitude"
    ) +
    # Create facets for each cruise
    facet_wrap(~ cruise) +
    # Use a minimal theme
    theme_minimal()
  
  # Store the plot
  yearly_plots[[year]] <- yearly_plot
}

# Create a directory to save the plots (if it doesn't already exist)
output_dir <- "plots/yearly"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save the yearly plots
for (year in names(yearly_plots)) {
  # Define the output file name
  output_file <- file.path(output_dir, paste0("plot_", year, ".png"))
  
  # Save the plot as a PNG
  ggsave(filename = output_file, 
         plot = yearly_plots[[year]], 
         width = 12, height = 10, dpi = 300,
         bg = "white")
}

