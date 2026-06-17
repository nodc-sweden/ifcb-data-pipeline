## General information

- Author: Wenche Eikrem
- Contact e-mail: <E-MAIL>
- DOI: 
- License: CC BY 4.0
- Version: <VERSION>
- This readme file was last updated: <DATE>

Please cite as: 

## Dataset description

This repository includes a dataset of manually annotated plankton images by phytoplankton experts at the Norwegian Institute for Water Research (NIVA). These images can be used for training classifiers to identify various plankton species. The images were captured using an Imaging FlowCytobot (IFCB, McLane Research Laboratories) from different locations and seasons in the Skagerrak and Kattegat. Images were gathered during regular cruises from <YEAR_START> to <YEAR_END>, utilizing an IFCB mounted as part of the underway FerryBox system on the M/S Color Fantasy. This collection consists of <N_IMAGES> annotated images across <CLASSES> different classes.

## Available data

There are two zip-packages available from this dataset:

- <ZIP_NAME>_annotated_images.zip - contains .png images that are manually annotated and organized into subfolders for each class
- <ZIP_NAME>_matlab_files.zip - includes raw data files (.roi, .hdr, .adc) and MATLAB files for creating a random forest image classifier using the code available at https://github.com/hsosik/ifcb-analysis
