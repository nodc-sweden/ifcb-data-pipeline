# SMHI IFCB Data Pipeline

## Overview

This repository contains the SMHI IFCB Data Pipeline, an R-based tool for processing and exporting pre-classified data from the Imaging FlowCytobot (IFCB) to [SHARK](https://sharkweb.smhi.se/). The data has already been classified using a Random Forest algorithm with the MATLAB package [`ifcb-analysis`](https://github.com/hsosik/ifcb-analysis) (Sosik and Olson, 2007). This pipeline focuses on subsequent steps like data cleaning, biovolume calculation, quality control, and export to the SHARK database. The pipeline relies heavily on functions from the [`iRfcb`](https://github.com/EuropeanIFCBGroup/iRfcb) R package.

### Key Features
- **Biovolume Calculation**: Computes biovolume based on the size and shape of detected particles.
- **Quality Control**: Implements checks to ensure data accuracy before export.
- **SHARK Format Export**: Prepares and exports the processed data in the required format for SHARK.
- **Annual and Multi-year Processing**: Allows for the processing of data for single years or across multiple years.

## Repository Structure
- **`ifcb-data-pipeline.Rmd`**: The main RMarkdown file that runs the entire pipeline. Currently placed in the root directory.
- **`data/`**: Folder to store input data files.
- **`output/`**: Folder to store out data files, including HTML reports.

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
   
4. **Edit environmental variables**:  
   Open the project and set your paths as .Renviron variables:
   ```r
   # Copy template
   file.copy(".Renviron-template", ".Renviron")
   
   # Edit file
   usethis::edit_r_environ("project")
   ```

5. **Define parameters**:  
   Update the parameters inside the `ifcb-data-pipeline.Rmd` file to suit your data processing needs, including:
   - Year to process
   - Classifier thresholds
   - Data export options

## Usage

To run the data pipeline, open `ifcb-data-pipeline.Rmd` in RStudio and knit it to execute the steps. The pipeline will proceed through the following stages:
1. **Data Import**: Reads and loads pre-classified IFCB data.
2. **Data Processing**: Conducts biovolume calculations, applies quality control checks, and prepares the data.
3. **Export**: Outputs the final processed data in the format required by SHARK.

### Customizing the Pipeline

You can modify the following parameters in the RMarkdown file:
- **Year and Classifier Thresholds**: Specify the year and thresholds for processing data.
- **Data Export Options**: Configure options for annual or multi-year data exports.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References
- Sosik, H.M., & Olson, R.J. (2007). Automated submersible flow cytometer for imaging phytoplankton. *Limnology and Oceanography: Methods*.
- Torstensson (2024). I 'R' FlowCytobot (iRfcb): Tools for Analyzing and Processing Data from the IFCB. R package version 0.3.7. https://doi.org/10.5281/zenodo.12533225
