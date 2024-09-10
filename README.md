# SMHI IFCB Data Pipelines

## Overview

This repository contains the SMHI IFCB Data Pipelines, an R-based tool for processing and exporting manually labelled and pre-classified data from the Imaging FlowCytobot (IFCB) to [SHARK](https://sharkweb.smhi.se/) and [Figshare](https://figshare.scilifelab.se/). The data has already been classified using a Random Forest algorithm with the MATLAB package [`ifcb-analysis`](https://github.com/hsosik/ifcb-analysis) (Sosik and Olson, 2007). The pipelines focuses on subsequent steps like data cleaning, biovolume calculation, quality control, image extraction and export to the SHARK database and to the [`SMHI IFCB plankton image reference library`](https://doi.org/10.17044/scilifelab.25883455) (Torstensson et al. 2024). The pipelines relies heavily on functions from the [`iRfcb`](https://github.com/EuropeanIFCBGroup/iRfcb) R package (> v0.3.9 is required).

### Key Features
- **Read Data**: Imports and loads manually labelled and pre-classified IFCB data, including metadata such as position, timestamp, and class scores.
- **Biovolume Calculation**: Computes biovolume and carbon content based on the size, taxonomic group and shape of detected particles.
- **Quality Control**: Implements automated checks to ensure data accuracy before export.
- **SHARK Format Export**: Prepares and exports the processed data in the required format for SHARK.
- **Image Export**: Extracts and archives annotated IFCB images for inclusion in the SMHI IFCB plankton image reference library.

## Repository Structure
- **`ifcb-data-pipeline.Rmd`**: The main RMarkdown file that runs the data pipeline to a standardized SHARK format. Currently placed in the root directory.
- **`ifcb-image-export-pipeline.Rmd`**: The main RMarkdown file that runs the image export pipeline to the `SMHI IFCB plankton image reference library`. Currently placed in the root directory.
- **`data/`**: Folder to store input data files.
- **`output/`**: Folder to store out data files, including HTML reports.
- **`templates/`**: Folder to store template README files for use in the `ifcb-image-export-pipeline.Rmd`.

### IFCB Data Folder Structure
The IFCB data folder, containing results processed by the `ifcb-analysis` MATLAB package, is linked within the R environment (see [Step 4](#step-4) of [Installation and Setup](#installation-and-setup)). This recommended folder structure is fully compatible with the data pipeline.

```
├── classified/
│   ├── classifier_name/
│   │   ├── classYYYY_vX/
│   │   │   └── *.mat        # Classified data files for a specific classifier and version
│
├── config/
│   └── class2use_classifier_name.mat   # Configuration file mapping classes to classifier names
│
├── data/
│   ├── YYYY/               # Yearly data directory
│   │   ├── DYYYYMMDDD/     # Subdirectory for each day of the year
│   │   │   ├── *.hdr       # Header files containing metadata for raw data
│   │   │   ├── *.roi       # Region of interest files
│   │   │   └── *.adc       # ADC files with ancillary data for each sample
│
├── features/
│   └── YYYY/               # Feature files directory for each year
│       └── *.csv           # CSV files containing extracted features
│
├── manual/
│   ├── classifier_name/    # Directory for manually classified data
│   │   ├── *.mat           # Manual classification files
│   │   └── summary/        # Summary files for manual classification
│   │       ├── classifier_name.mat     # MATLAB Random Forest results file
│   │       ├── classifier_name.csv     # CSV file with class scores from OOB analysis (no threshold)
│   │       └── classifier_name_opt.csv # CSV file with class scores from OOB analysis (opt threshold)
```

## Installation and Setup

1. **Clone the repository**:  
   ```
   git clone https://github.com/nodc-sweden/ifcb-data-pipeline.git
   ```

2. **Install Python**:  
   The pipeline uses Python for some data processing tasks using `reticulate`. To prepare Python:
   - Install [Python](https://www.python.org/downloads/) (the pipeline has been tested using versions > 3.10)
   - A virtual environment will be created when the pipeline is initialized using `iRfcb::ifcb_py_install()`

3. **Install necessary packages**:  
   In R: Ensure you have the required R packages installed by running:
   ```r
   install.packages(c("tidyverse", "worrms", "knitr", "rmarkdown", "leaflet", 
                      "htmltools", "patchwork", "reticulate"))
   
   # install.packages("devtools")
   devtools::install_github("EuropeanIFCBGroup/iRfcb", dependencies = TRUE)
   ```
   
4. <a name="step-4"></a> **Edit environmental variables**:  
   Open the project and set your paths and e-mail address as .Renviron variables:
   ```r
   # Copy template
   file.copy(".Renviron-template", ".Renviron")
   
   # Edit file
   usethis::edit_r_environ("project")
   ```

5. **Define parameters**:  
   Update the parameters inside the `ifcb-data-pipeline.Rmd` and `ifcb-image-export-pipeline.Rmd` files to suit your data processing needs.

## Usage

### IFCB data pipeline

To run the data pipeline, open `ifcb-data-pipeline.Rmd` in RStudio and knit it to execute the steps. The `ifcb-data-pipeline.Rmd` will proceed through the following stages:

1. **Data Import**: Reads and loads pre-classified IFCB data.
2. **Data Processing**: Conducts biovolume calculations, applies quality control checks, and prepares the data.
3. **Export**: Outputs the final processed data in the format required by SHARK.

To run the image export pipeline, open `ifcb-image-export-pipeline.Rmd` in RStudio and knit it to execute the steps.

### IFCB image export pipeline

The `ifcb-image-export-pipeline.Rmd` will proceed through the following stages:

1. **Extract PNG Images**: Extracts annotated `.png` images organized into subfolders based on the specified classifiers.
2. **Zip Files**: Creates zip archives containing the `.png` images, MATLAB files, and raw data.
3. **README and Manifest Creation**: Updates README files with image counts and generates a `MANIFEST.txt` for the final output.

This process currently supports multiple classifiers, such as: "Baltic", "Skagerrak-Kattegat", "Tångesund", and "iRfcb" and is essential for updating the [`SMHI IFCB plankton image reference library`](https://doi.org/10.17044/scilifelab.25883455).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References
- Sosik, H.M., & Olson, R.J. (2007). Automated submersible flow cytometer for imaging phytoplankton. *Limnology and Oceanography: Methods*.
- Torstensson (2024). I 'R' FlowCytobot (iRfcb): Tools for Analyzing and Processing Data from the IFCB. R package version 0.3.9. https://doi.org/10.5281/zenodo.12533225
- Torstensson, Anders; Skjevik, Ann-Turi; Mohlin, Malin; Karlberg, Maria; Karlson, Bengt (2024). SMHI IFCB plankton image reference library. SciLifeLab. Dataset. doi:10.17044/scilifelab.25883455
