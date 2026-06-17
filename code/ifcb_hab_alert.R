# -------------------------------
# Load libraries
# -------------------------------

suppressPackageStartupMessages({
  library(SHARK4R, quietly = TRUE)
  library(iRfcb, quietly = TRUE)
  library(curl, quietly = TRUE)
  library(hdf5r, quietly = TRUE)
})

source("code/utils/clean_taxa_fn.R")

# -------------------------------
# Configuration
# -------------------------------

cnn_model <- "SMHI-NIVA-SYKE-SAMS-SZN-ResNet50-V6" # niva_smhi_baltic or NA for MATLAB

year <- as.numeric(format(Sys.Date(), "%Y"))
date <- format(Sys.Date()-1, "%Y%m%d") # Yesterday

ifcb_path <- Sys.getenv("ifcb_path")
ifcb_base_path <- Sys.getenv("ifcb_base_path")
havgem_path <- Sys.getenv("havgem_path")
repo_dir <- Sys.getenv("repo_path")
smtp_server <- Sys.getenv("SMTP_SERVER")
at_email <- Sys.getenv("at_email")
bk_email <- Sys.getenv("bk_email")
alg1_email <- Sys.getenv("alg1_email")
ifcb_email <- Sys.getenv("ifcb_email")

feature_folder <- file.path(ifcb_path, "features", "v4", year)
raw_folder <- file.path(ifcb_path, "data", year, paste0("D", date))

class_folder <- file.path(ifcb_path, "classified", cnn_model, paste0("class", year, "_v3"))

# Taxa that should be forced to be regarded as diatoms, as they are not according to WoRMS
diatom_include = c("Actinocyclus_spp", "Actinocyclus", "Navicula-like", "Navicula", "cf_Proboscia_rhizosolenia")

# -------------------------------
# Get data
# -------------------------------

# List already merged class files
class_files <- list.files(class_folder, ".h5", full.names = TRUE, recursive = TRUE)
combined_class_files <- class_files[grepl(paste0("D", date), basename(class_files))]
class_bins <- gsub("_class.h5", "", basename(combined_class_files))

if (length(class_bins) > 0) {

  cat("Calculating biovolume and biomass data...\n")
  
  biovolumes <- ifcb_summarize_biovolumes(feature_folder = feature_folder,
                                          class_files = combined_class_files,
                                          hdr_folder = raw_folder,
                                          micron_factor = 1 / 2.77,
                                          diatom_class = "Bacillariophyceae",
                                          diatom_include = diatom_include,
                                          use_python = FALSE,
                                          feature_version = 4,
                                          verbose = TRUE,
                                          drop_zero_volume = TRUE)
  
} else {
  message("No bins to process, stopping script.")
  quit(save = "no", status = 0)
}

# -------------------------------
# Correct taxa names
# -------------------------------

hab_list <- get_hab_list(species_only = FALSE) %>%
  filter(taxonRank %in% c("Species", "Genus"))

warning_levels <- read_csv(file.path(repo_dir, "data", "taxa_lookup.csv"),
                        progress = FALSE, col_types = cols()) %>%
  select(AphiaID, warning_level, HAB) %>%
  filter(HAB)

hab_aphia_ids <- unique(c(hab_list$AphiaID, warning_levels$AphiaID))

# Extract unique taxa names from biovolumes and manual_biovolumes
taxa_names <- unique(biovolumes$class)

class_names <- clean_taxa(taxa_names) %>%
  mutate(trophic_type = ifcb_get_trophic_type(cleaned_name),
         hab_taxa = aphia_id %in% hab_aphia_ids) %>%
  rename(class_clean = cleaned_name,
         scientificname = reported_name,
         class = class_name,
         sflag = flag) %>%
  distinct() %>%
  left_join(warning_levels, by = c("aphia_id" = "AphiaID")) %>%
  distinct()

# -------------------------------
# Join data
# -------------------------------

hab_data <- biovolumes %>%
  left_join(class_names, by = "class") %>%
  filter(hab_taxa)

max_counts <- hab_data %>%
  group_by(aphia_id) %>%
  slice_max(counts_per_liter, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    aphia_id,
    class_clean,
    sample,
    counts,
    counts_per_liter,
    warning_level
  ) %>%
  mutate(counts_per_liter = round(counts_per_liter)) %>%
  arrange(class_clean)

biomass_data <- biovolumes %>%
  group_by(sample) %>%
  summarize(total_carbon = sum(carbon_ug_per_liter, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    bloom_category = case_when(
      total_carbon > 120  ~ "Major bloom",
      total_carbon >= 70 ~ "Bloom",
      total_carbon >= 40 ~ "Minor bloom",
      TRUE               ~ "No bloom"
    ),
    bloom_category = factor(
      bloom_category,
      levels = c("No bloom", "Minor bloom", "Bloom", "Major bloom")
    )
  )

cyano_counts <- biovolumes %>%
  left_join(class_names, by = "class") %>%
  filter(worms_phylum == "Cyanobacteria") %>%
  group_by(aphia_id) %>%
  slice_max(counts_per_liter, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    aphia_id,
    class_clean,
    sample,
    counts,
    counts_per_liter,
    warning_level
  ) %>%
  mutate(counts_per_liter = round(counts_per_liter)) %>%
  arrange(class_clean)

# -------------------------------
# Bloom notification
# -------------------------------
bloom_samples  <- biomass_data %>% filter(bloom_category != "No bloom")
bloom_detected <- nrow(bloom_samples) > 0

if (bloom_detected) {
  top_sample <- biomass_data %>%
    slice_max(total_carbon, n = 1, with_ties = FALSE)
  
  bloom_link <- paste0(
    '<a href="https://ifcb-dashboard-utv.smhi.se/timeline?dataset=RV_Svea&bin=',
    top_sample$sample, '">', top_sample$sample, '</a>'
  )
  
  bloom_notice <- paste0(
    '<p style="background-color:#ffcccc;font-weight:bold;padding:8px;">',
    nrow(bloom_samples), ' sample(s) reached bloom level. ',
    'Highest particle carbon: ',
    format(round(top_sample$total_carbon, 1), big.mark = " ", trim = TRUE),
    ' &micro;g C/L (', as.character(top_sample$bloom_category),
    ') in sample ', bloom_link, '.',
    '</p>'
  )
} else {
  bloom_notice <- ""
}

# -------------------------------
# Send E-mail notification
# -------------------------------

# Dynamic emails
# names <- c("Anders Torstensson", "Bengt Karlson")
# emails <- c(at_email, bk_email)
names <- c("Anders Torstensson")
emails <- c(at_email)
to_header <- paste0('"', names, '" <', emails, '>')  # no spaces added
to_header <- paste(to_header, collapse = ",")        # no spaces after comma

date_string <- format(as.Date(date, "%Y%m%d"), "%Y-%m-%d")

warning_detected <- any(
  max_counts$counts_per_liter > max_counts$warning_level,
  na.rm = TRUE
)

subject <- if (warning_detected) {
  paste0("IFCB HAB ALERT - ", date_string)
} else {
  paste0("IFCB HAB Summary - ", date_string)
}

email_table <- max_counts %>%
  mutate(
    sample = paste0(
      '<a href="https://ifcb-dashboard-utv.smhi.se/timeline?dataset=RV_Svea&bin=',
      sample,
      '">',
      sample,
      '</a>'
    ),
    alert = !is.na(warning_level) & counts_per_liter > warning_level,
    row_style = ifelse(
      alert,
      ' style="background-color:#ffcccc;font-weight:bold;"',
      ''
    ),
    counts_per_liter = format(counts_per_liter, big.mark = " ", trim = TRUE),
    warning_level = ifelse(
      is.na(warning_level),
      "",
      format(warning_level, big.mark = " ", trim = TRUE)
    )
  )

table_rows <- paste0(
  apply(email_table, 1, function(x) {
    paste0(
      "<tr",
      x["row_style"],
      ">",
      "<td>", x["aphia_id"], "</td>",
      "<td>", x["class_clean"], "</td>",
      "<td>", x["sample"], "</td>",
      "<td>", x["counts"], "</td>",
      "<td>", x["counts_per_liter"], "</td>",
      "<td>", ifelse(is.na(x["warning_level"]), "", x["warning_level"]), "</td>",
      "</tr>"
    )
  }),
  collapse = "\n"
)

html_table <- paste0(
  '<table border="1" cellpadding="5" cellspacing="0">',
  '<tr>',
  '<th>Aphia ID</th>',
  '<th>Scientific name</th>',
  '<th>Sample</th>',
  '<th>Images</th>',
  '<th>Counts/L</th>',
  '<th>Warning level</th>',
  '</tr>',
  table_rows,
  '</table>'
)

message <- paste0(
  'From: "IFCB" <ifcb.u@smhi.se>
To: ', to_header, '
Subject: ', subject, '
MIME-Version: 1.0
Content-Type: text/html; charset=UTF-8
<html>
<body>
<p>Hej,</p>
', bloom_notice, '
<p>
The table below summarizes the maximum observed abundance for each HAB taxon
during ', date_string, '. For each taxon, the sample shown is the sample in
which the highest concentration (counts per liter) was observed during that day.
</p>
', html_table, '
<p>
Rows highlighted in red indicate that the observed concentration exceeded the
configured warning level.
</p>
<p>
Mvh<br>
ifcb.u
</p>
</body>
</html>'
)

# Deliver mail to recipients
send_mail(
  mail_from = ifcb_email,
  mail_rcpt = emails,
  message = message,
  smtp_server = smtp_server,
  use_ssl = "force",
  verbose = FALSE
)


# -------------------------------
# Send BAWS E-mail notification
# -------------------------------

# Dynamic emails
names <- c("Alg.gbg1", "Anders Torstensson")
emails <- c(alg1_email, at_email)
to_header <- paste0('"', names, '" <', emails, '>')  # no spaces added
to_header <- paste(to_header, collapse = ",")        # no spaces after comma

date_string <- format(as.Date(date, "%Y%m%d"), "%Y-%m-%d")

d <- as.Date(date_string)

month_day <- format(d, "%m-%d")

warning_detected <- any(
  cyano_counts$counts_per_liter > cyano_counts$warning_level,
  na.rm = TRUE
)

subject <- if (warning_detected) {
  paste0("IFCB BAWS ALERT - ", date_string)
} else {
  paste0("IFCB BAWS Summary - ", date_string)
}

priority_species <- c("Aphanizomenon flosaquae", "Nodularia spumigena", "Dolichospermum")

email_table <- cyano_counts %>%
  mutate(
    sample = paste0(
      '<a href="https://ifcb-dashboard-utv.smhi.se/timeline?dataset=RV_Svea&bin=',
      sample,
      '">',
      sample,
      '</a>'
    ),
    alert = !is.na(warning_level) & counts_per_liter > warning_level,
    is_priority = class_clean %in% priority_species,
    row_style = case_when(
      alert        ~ ' style="background-color:#ffcccc;font-weight:bold;"',
      is_priority  ~ ' style="background-color:#d4edda;color:#155724;font-weight:bold;"',
      TRUE         ~ ''
    ),
    counts_per_liter = format(counts_per_liter, big.mark = " ", trim = TRUE),
  ) %>%
  arrange(desc(is_priority), class_clean)  # priority species first, then alphabetical

table_rows <- paste0(
  apply(email_table, 1, function(x) {
    paste0(
      "<tr",
      x["row_style"],
      ">",
      "<td>", x["aphia_id"], "</td>",
      "<td>", x["class_clean"], "</td>",
      "<td>", x["sample"], "</td>",
      "<td>", x["counts"], "</td>",
      "<td>", x["counts_per_liter"], "</td>",
      "</tr>"
    )
  }),
  collapse = "\n"
)

html_table <- paste0(
  '<table border="1" cellpadding="5" cellspacing="0">',
  '<tr>',
  '<th>Aphia ID</th>',
  '<th>Scientific name</th>',
  '<th>Sample</th>',
  '<th>Images</th>',
  '<th>Counts/L</th>',
  '</tr>',
  table_rows,
  '</table>'
)

message <- paste0(
  'From: "IFCB" <ifcb.u@smhi.se>
To: ', to_header, '
Subject: ', subject, '
MIME-Version: 1.0
Content-Type: text/html; charset=UTF-8
<html>
<body>
<p>Hej,</p>
', bloom_notice, '
<p>
The table below summarizes the maximum observed abundance of each cyanobacterial 
taxon on ', date_string, '. For each taxon, the sample shown is the sample in
which the highest concentration (counts per liter) was observed during that day.
</p>
', html_table, '
<p style="font-size:0.9em;color:#155724;">
Rows highlighted in green indicate filamentous cyanobacteria
(<em>Aphanizomenon flosaquae</em>, <em>Nodularia spumigena</em>, <em>Dolichospermum</em>)
that are of particular concern for bloom formation and toxin production.
</p>
<p>
Mvh<br>
ifcb.u
</p>
</body>
</html>'
)

if (month_day >= "05-01" && month_day <= "10-31") {
  # between 1 May and 31 October
  # Deliver mail to recipients
  send_mail(
    mail_from = ifcb_email,
    mail_rcpt = emails,
    message = message,
    smtp_server = smtp_server,
    use_ssl = "force",
    verbose = FALSE
  )
}
