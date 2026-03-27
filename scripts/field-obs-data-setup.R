## USDA cotton experiment
## Field data setup
## Code written and annotated by JM (last updated 2026-03-25)

# The code in this script file:
# - determines the relationship between non-destructive measurements of stem length
#   and dry biomass;
# - uses this relationship to express initial and final plant sizes in terms of
#   dry biomass and estimate RGR for each experimental plant across the season;
# - converts data from long form (separate entries for initial, sampling, and
#   harvest measurements) into wide form (one row per plant);
# - adds information about herbivore density and sampling timepoint to the ant
#   attack data;
# - saves cleaned data as .csv.

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("./USDA-cotton-exp")

# Load required packages.
require(doBy)
require(ggplot2)

## Dataframe for plant damage and biomass.======================================

## First, find relationship between summed stem lengths and dry biomass, measured
##  for a subset of experimental plants.

# Read in data.
stem.to.biomass <- read.csv("raw-data/stem-lengths-to-biomass.csv")

# Model log dry mass as a function of summed stem lengths.
stem.to.biomass.mod <- lm(log(Dry_mass) ~ Stem_lengths, data = stem.to.biomass)
summary(stem.to.biomass.mod)

# Visualize model best-fit.
ggplot(data = stem.to.biomass, aes(x = Stem_lengths, y = log(Dry_mass))) +
  geom_point() +
  geom_abline(intercept = 4.23, slope = 0.0013)

## Then, prep data for dry biomass at harvest.

# Read in data on dry biomass at harvest.
final_biomass <- read.csv("raw-data/final-biomass.csv")
final_biomass <- final_biomass[-which(final_biomass$Exclude == "Y"),]

# Some plants had biomass contained in multiple bags -- sum these per plant.
final_biomass_sum <- doBy::summary_by(Dry_mass ~ Plot_ID + Plant_ID,
                                      data = final_biomass,
                                      FUN = c(sum, length))
hist(log(final_biomass_sum$Dry_mass.sum))

## Put it all together and convert from long to wide form.

# Read in damage measurements.
field_obs <- read.csv("raw-data/field-obs-data.csv")
# Specify columns to keep
cols_want <- c("Plot_ID", "Treatment", "Timepoint", "Herb_dens", 
               "Plant_ID", "Geno", "Date", "Stem_lengths", "Percent_damage")
# Subset to just initial plant measurements (one per plant)...
field_obs2 <- field_obs[which(field_obs$Census_ID == "initial"),
                        cols_want]
# ...then rename date, stem lengths, and damage columns to indicate "initial" measurements
colnames(field_obs2) <- c(cols_want[1:6], "Initial_Date", "Initial_Stem_lengths",
                          "Initial_Percent_damage")

# Add columns for date and damage estimate at time of sampling
field_obs2[,c("Sampling_Plant_ID", "Sampling_Date", "Sampling_Percent_damage", "Sampling_Aphids")] <-
  field_obs[which(field_obs$Census_ID == "sampling"),][match(field_obs2$Plant_ID, field_obs[which(field_obs$Census_ID == "sampling"), 'Plant_ID']),
                                                       c("Plant_ID", "Date", "Percent_damage", "Aphids")]
# Add columns for date and damage estimate at harvest
field_obs2[,c("Harvest_Plant_ID", "Harvest_Date", "Harvest_Percent_damage")] <-
  field_obs[which(field_obs$Census_ID == "final"),][match(field_obs2$Plant_ID, field_obs[which(field_obs$Census_ID == "final"), 'Plant_ID']),
                                                    c("Plant_ID", "Date", "Percent_damage")]
# Using dataframe with summed dry biomass, add column for dry biomass at harvest
field_obs2[,c("Harvest2_Plant_ID", "Harvest_Dry_mass")] <-
  final_biomass_sum[match(field_obs2$Plant_ID, final_biomass_sum$Plant_ID), c("Plant_ID", "Dry_mass.sum")]
# Convert initial measurements of stem length into log dry biomass
field_obs2$Initial_biomass_log <- with(field_obs2,
                                       4.23 + 0.0013*Initial_Stem_lengths)
# Get log dry biomass at harvest
field_obs2$Harvest_biomass_log <- log(field_obs2$Harvest_Dry_mass)
# Estimate plant RGR across the season using these dry biomass measurements
field_obs2$Plant_RGR <- with(field_obs2, Harvest_biomass_log - Initial_biomass_log)

# Check that all plant IDs match (data associated with correct plant during merge)
field_obs2[which(field_obs2$Plant_ID != field_obs2$Sampling_Plant_ID),]
field_obs2[which(field_obs2$Plant_ID != field_obs2$Harvest_Plant_ID),]
field_obs2[which(field_obs2$Plant_ID != field_obs2$Harvest2_Plant_ID),]
# All good, so get rid of the extra plant ID columns
field_obs <- field_obs2[,!(names(field_obs2) %in% c("Sampling_Plant_ID", "Harvest_Plant_ID", "Harvest2_Plant_ID"))]

# Save cleaned file
write.csv(field_obs, "cleaned-data/field-data.csv", row.names = FALSE)

## Ant attack dataframe.========================================================

## Read in ant attack data and full data (for metadata on herbivore density and
##  timepoint associated with each experimental population)
ants <- read.csv("raw-data/plants-attacked.csv")
full_data <- read.csv("cleaned-data/full-data.csv")

# Make a guide to the sampling timepoint and herbivore density associated with
#  each experimental population.
pop_guide <- doBy::summaryBy(Day ~ Cage_pos,
                             id = c("Timepoint", "Herb_add"),
                             data = full_data,
                             keep.names = TRUE)

# Add these metadata to the ants dataframe.
ants[,c('Herb_add', 'Timepoint', 'Day')] <- pop_guide[match(ants$Cage_pos,
                                                            pop_guide$Cage_pos),
                                                      c('Herb_add', 'Timepoint', 'Day')]

# Save cleaned ant attack data.
write.csv(ants, "cleaned-data/ant-attack.csv")
