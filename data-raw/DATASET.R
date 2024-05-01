## code to prepare `DATASET` dataset goes here

# Load raw data from .csv file
hnsc_data<- readRDS("C:/Users/farza/Documents/RESEARCH/SSTP21/Proposal defense/final data analysis/R package for marginal analysis/DSFDRC/data-raw/hnsc.rds")
# Apply preprocessing...
# Save the cleaned data in the required R package location
usethis::use_data(hnsc_data,overwrite = TRUE)
