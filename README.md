# Cavity-nesting-Hymenoptera-in-deadwood
This is the R script for statistical analyses pertaining to the manuscript: "Moisture reduces saproxylic bee, wasp, and parasitoid diversity in lying and standing deadwood" (Martini et al. submitted).
Files include the main R script, and a secondary script with helper functions sourced into the main analysis

# README for data "log_data.csv", "deadwd_insect_data.csv" and associated R scripts

Access this dataset on Figshare: https://doi.org/10.6084/m9.figshare.31979481

Data and analysis of host-parasitoid diversity in deadwood traps from the 'BEF-China' platform.
Collected in 2024 using deadwood traps for cavity-nesting Hymenoptera and their parasitoids.

## Description of the data and file structure

"log_data.csv":
Individual deadwood-level data on the abundance and species richness of cavity-nesting Hymenoptera and their natural enemies (parasites, parasitoids, kleptoparasites).
Includes deadwood-lelvel variables such as ID, treatment, position, ant exclusion (binary), presence (binary) and abundance of ants on the deadwood, and the gravimetric moisture content (GMC) of the wood. Also includes plot-level forest stand variables including tree species richness, canopy cover (as % pr proportion), coarse woody debris (CWD) volume, ant occurrence from the leaf litter, and geomorphological categorization (also numerical) of each plot. Note that some dataopoints for GMC (10%), CWD (6%), and leaf litter ant occurrence (5%) were imputed prior to this analysis using GLMM model-based predictions fitted to complete cases following Nakagawa and Freckleton, 2008; and Van Buuren, S. 2018. 

"deadwd_insect_data.csv":
Raw host and parasitoid data from each deadwood trap. Used for community composition and other analyses.

Main processed data file is "log_data.csv"
Raw data file is "deadwd_insect_data.csv"
Metadata file is "log_data_metadata.csv"
Main data analysis R script is "deadwd_script.R"
Secondary R script with helper functions is "deadwd_functions.R". This secondary script is sourced into the main script for data analysis.

Data is stored in Figshare, scripts are stored on Github + Zenodo

For questions please contact:
Massimo Martini - massimo.martini@nature.uni-freiburg.de

## Code/Software

All code was written and runs in R. 
R version 4.5.1 was used for the data analysis.
