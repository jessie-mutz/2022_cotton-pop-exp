## USDA cotton experiment
## Chemical data setup
## Code written and annotated by JM (last updated 2024-10-16)

# The code in this script file:
# - reformat dataframe to one row per plant;
# - standardize abundances to sample mass to get relative abundances;
# - save cleaned data as .csv.

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("./USDA-cotton-exp")

# Load required packages.
require(reshape2)

# Read in chemical data
chem_data <- read.csv("raw-data/chem-data-formatted.csv")

# Convert to long format with one row per plant.
chem_data_melt <- reshape2::melt(data = chem_data, id.vars = "compound")
chem_data <- reshape2::dcast(variable ~ compound, value.var = "value", 
                             data = chem_data_melt)

# Extract plant ID from variable column
chem_data$variable <- with(chem_data, as.character(variable))
chem_data$Plant_ID <- with(chem_data, substr(variable,
                                              start = nchar(variable) - 3,
                                              stop = nchar(variable)))
chem_data[,c(2:5)] <- sapply(chem_data[,c(2:5)], FUN = as.numeric)

# Add abundances for H20 blanks
chem_data$H2OBlank_3.6 <- 150526763
chem_data$H2OBlank_10.5 <- 600869354

# Indicate whether to include an observation based on whether the relative
# abundance is at least 2x greater than the H20 blank.
chem_data$H_3.6_include <- with(chem_data,
                                ifelse(Hemigossypolone_3.6 < 2*H2OBlank_3.6, 
                                       "N", "Y"))
chem_data$H_10.5_include <- with(chem_data,
                                 ifelse(Hemigossypolone_10.5 < 2*H2OBlank_10.5, 
                                        "N", "Y"))

# For observations to include, standardize abundances by sample mass. For
# observations not to include, set relative abundance to NA.
chem_data$H_3.6_relabun <- with(chem_data,
                                ifelse(H_3.6_include == "Y", 
                                       Hemigossypolone_3.6/sample_mass, NA))
chem_data$H_10.5_relabun <- with(chem_data,
                                 ifelse(H_10.5_include == "Y", 
                                        Hemigossypolone_10.5/sample_mass, NA))

# Save cleaned data as .csv.
write.csv(chem_data, "cleaned-data/chem-data.csv", row.names = FALSE)

rm(chem_data, chem_data_melt)
