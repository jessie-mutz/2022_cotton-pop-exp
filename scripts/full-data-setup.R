## USDA cotton experiment
## Full data setup
## Code written and annotated by JM (last updated 2024-10-21)

# The code in this script file:
# - reads in data from raw and cleaned data files;
# - creates one dataframe with all defense, resistance, and damage measurements 
#   for each experimental plant, columns for population ID, herbivore density
#   treatment, and timepoint of sampling;
# - saves cleaned data as .csv.

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("./USDA-cotton-exp")

# Read in cleaned data files
full_data <- read.csv("raw-data/plant-guide.csv")
phys_traits <- read.csv ("cleaned-data/phys-trait-std.csv")
chem_traits <- read.csv("cleaned-data/chem-data.csv")
bioassay_data <- read.csv("cleaned-data/bioassay-data.csv")
field_obs <- read.csv("cleaned-data/field-data.csv")

# Set up full_data with days after the start of herbivory and herbivore density.
times <- data.frame(cbind(
  Timepoint = seq(1, 5),
  Days = c(2, 4, 7, 12, 21)
))

herbs <- data.frame(cbind(
  Herb_trt = c("control", "low", "mid", "high"),
  Herb_dens = c(0, 6, 12, 25)
))

full_data$Herb_add <- herbs[match(full_data$Herb_dens, herbs$Herb_trt), 'Herb_dens']
full_data$Herb_add <- as.numeric(full_data$Herb_add)
full_data$Day <- times[match(full_data$Timepoint, times$Timepoint), 'Days']
full_data$Timepoint <- as.factor(full_data$Timepoint)

## Add measurements from each dataframe (phys_traits, chem_traits, bioassay_data,
## field_obs) to full_data

# Specify which columns to add to full_data.
phys_traits_cols <- colnames(phys_traits)[-which(colnames(phys_traits) == 'Plant_ID')]
chem_cols <- c('H_3.6_relabun', 'H_10.5_relabun')
bioassay_cols <- c('Dead_1', 'Dead_2', 'Damage_1', 'Damage_2',
                   'RGR_1', 'RGR_2', 'Initial_leaf_area',
                   'Area_consumed_1', 'Area_consumed_2', 'Mean_area_consumed'
                   )
field_cols <- c('Initial_Percent_damage', 'Sampling_Percent_damage',
                'Sampling_Aphids', 'Harvest_Percent_damage',
                'Initial_biomass_log', 'Harvest_biomass_log', 'Plant_RGR')

# Add physical traits to full_data
full_data[,phys_traits_cols] <- phys_traits[
  match(full_data$Plant_ID, phys_traits$Plant_ID), phys_traits_cols]
# Add chemical traits to full_data
full_data[,chem_cols] <- chem_traits[
  match(full_data$Plant_ID, chem_traits$Plant_ID), chem_cols]
# Add bioassay responses to full_data
full_data[,bioassay_cols] <- bioassay_data[
  match(full_data$Plant_ID, bioassay_data$Plant_ID), bioassay_cols]
full_data$RGR_mean <- rowMeans(full_data[,c('RGR_1', 'RGR_2')], na.rm = T)

# For bioassay data:
# Add column for whether any leaf tissue was consumed during the bioassay period.
full_data$Any_consumed_1 <- with(full_data,
                                 ifelse(Area_consumed_1 == 0, 0, 1))
full_data$Any_consumed_2 <- with(full_data,
                                 ifelse(Area_consumed_2 == 0, 0, 1))
# Add columns for the total number of leaves with any sign of damage during the
# bioassay period (up to 2) and the number of observations of leaf consumption 
# (up to 2).
full_data$Total_any_consumed <- rowSums(full_data[,c('Any_consumed_1',
                                                     'Any_consumed_2')],
                                        na.rm = TRUE)
full_data$Obs_any_consumed <- 2 - rowSums(is.na( full_data[,c('Any_consumed_1', 
                                                              'Any_consumed_2')]))
# Number of bioassay caterpillars that died during the bioassay period.  
full_data$Total_dead <- rowSums(full_data[,c('Dead_1', 'Dead_2')], na.rm = TRUE)


# Add damage and plant biomass data to full_data
full_data[,field_cols] <- field_obs[
  match(full_data$Plant_ID, field_obs$Plant_ID), field_cols]


# Save compiled data as .csv.
write.csv(full_data, "cleaned-data/full-data.csv", row.names = FALSE)
