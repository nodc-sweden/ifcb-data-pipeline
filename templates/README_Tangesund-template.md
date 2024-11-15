## General information

- Author: Anders Torstensson, Ann-Turi Skjevik, Malin Mohlin, Maria Karlberg, Bengt Karlson
- Contact e-mail: <E-MAIL>
- DOI: 10.17044/scilifelab.25883455
- License: CC BY 4.0
- Version: <VERSION>
- This readme file was last updated: <DATE>

Please cite as: Torstensson, Anders; Skjevik, Ann-Turi; Mohlin, Malin; Karlberg, Maria; Karlson, Bengt (<YEAR>). SMHI IFCB plankton image reference library. SciLifeLab. Dataset. https://doi.org/10.17044/scilifelab.25883455.v<VERSION>

## Dataset description

This dataset includes manually annotated plankton images by phytoplankton experts at the Swedish Meteorological and Hydrological Institute (SMHI). These images can be used for training classifiers to identify various plankton species. The images were captured using an Imaging FlowCytobot (IFCB, McLane Research Laboratories). In 2016, an IFCB was deployed in situ at depths between 3 and 18 meters, near a mussel farm in Tångesund, Mollösund (Skagerrak). This dataset contains <N_IMAGES> annotated images from <CLASSES> different classes.

## Available data

There are two zip-packages available from this dataset:

- <ZIP_NAME>_annotated_images.zip - contains .png images that are manually annotated and organized into subfolders for each class
- <ZIP_NAME>_matlab_files.zip - includes raw data files (.roi, .hdr, .adc) and MATLAB files for creating a random forest image classifier using the code available at https://github.com/hsosik/ifcb-analysis
