library(tidyverse)
library(SHARK4R)

ifcb_path <- Sys.getenv("ifcb_path")

baltic_2024 <- read_tsv(file.path(ifcb_path, "shark", "shark2024_Baltic_v9/SHARK_IFCB_2024_Baltic_SMHI/processed_data/data.txt"))
baltic_2025 <- read_tsv(file.path(ifcb_path, "shark", "shark2025_Baltic_v9/SHARK_IFCB_2025_Baltic_SMHI/processed_data/data.txt"))
westcoast_2024 <- read_tsv(file.path(ifcb_path, "shark", "shark2024_Skagerrak-Kattegat_v7/SHARK_IFCB_2024_Skagerrak-Kattegat_SMHI/processed_data/data.txt"))
westcoast_2025 <- read_tsv(file.path(ifcb_path, "shark", "shark2025_Skagerrak-Kattegat_v7/SHARK_IFCB_2025_Skagerrak-Kattegat_SMHI/processed_data/data.txt"))

all_data <- bind_rows(baltic_2024, westcoast_2024, baltic_2025, westcoast_2025)

all_data <- all_data %>%
  filter(!is.na(LATNM)) %>%
  mutate(plankton_group = assign_phytoplankton_group(LATNM, APHIA_ID)$plankton_group)

all_data_grouped <- all_data %>%
  group_by(MYEAR, CRUISE_NO, SDATE, STIME, SMPNO, plankton_group) %>%
  summarise(biomass = sum(C_CONC, na.rm = TRUE),
            biovol = sum(BIOVOL, na.rm = TRUE))

# Plot for 2024
plot_2024 <- plot_data %>%
  filter(MYEAR == 2024) %>%
  ggplot(aes(x = datetime, y = biomass, fill = plankton_group)) +
  geom_area(alpha = 0.8, position = "stack") +
  facet_wrap(~ CRUISE_NO, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 70)) +
  labs(
    x = "Date-Time",
    y = "Biomass",
    fill = "Plankton Group",
    title = "Plankton Group Biomass Time Series - 2024"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Plot for 2025
plot_2025 <- plot_data %>%
  filter(MYEAR == 2025) %>%
  ggplot(aes(x = datetime, y = biomass, fill = plankton_group)) +
  geom_area(alpha = 0.8, position = "stack") +
  facet_wrap(~ CRUISE_NO, scales = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(ylim = c(0, 300)) +
  labs(
    x = "Date-Time",
    y = "Biomass",
    fill = "Plankton Group",
    title = "Plankton Group Biomass Time Series - 2025"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Display plots
plot_2024
plot_2025

# Save the 2024 plot
ggsave("plankton_biomass_2024.png", plot = plot_2024,
       width = 12, height = 6, dpi = 300)

# Save the 2025 plot
ggsave("plankton_biomass_2025.png", plot = plot_2025,
       width = 12, height = 6, dpi = 300)
