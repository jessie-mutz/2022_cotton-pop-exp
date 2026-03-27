## USDA cotton experiment
## Bioassay data setup
## Code written and annotated by JM (last updated 2024-10-16)

# The code in this script file:
# - changes RGR data into wide-format so that there is one line of data per plant;
# - compiles data on initial and final leaf area measurements and leaf area 
#   consumed by each bioassay caterpillar;
# - combines these into a single dataframe and saves to 'cleaned-data' as .csv.

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("./USDA-cotton-exp")

# Load required packages.
require(reshape2)

## Change RGR dataframe into wide format so that there is one line per plant.
# Read in RGR data.
bioassay_RGR <- read.csv("raw-data/bioassay-RGR-data.csv")

# Make bioassay_RGR wide-format.
bioassay_RGR_melt <- reshape2::melt(bioassay_RGR, id.vars = 1:9)

bioassay_RGR_cast <- reshape2::dcast(bioassay_RGR_melt, 
                                     Plant_ID ~ variable + Spod_ID,
                                     value.var = "value")

## Initial and final leaf area data.
# Read in leaf area data.
bioassay_initial <- read.csv("raw-data/bioassay-leaves-initial.csv")
bioassay_final <- read.csv("raw-data/bioassay-leaves-final.csv")

# Remove extraneous columns
bioassay_initial <- bioassay_initial[, -c(1:6)]
bioassay_final <- bioassay_final[,-c(1:6)]
bioassay_final <- na.omit(bioassay_final)

# For final leaf area data:
# Get plant ID from "Notes" column
bioassay_final$Plant_ID <- substr(bioassay_final$Notes, start = 1, stop = 4)
# Get leaf rep (half 1 or half 2) from "Notes" column
bioassay_final$Leaf_rep <- substr(bioassay_final$Notes, start = 6, stop = 6)

# Change bioassay_final dataframe to wide format
bioassay_final_melt <- reshape2::melt(bioassay_final, 
                                      measure.vars = c("Consumed_Leaf_Area_.cm2.",
                                                       "Total_Leaf_Area_.cm2."))

bioassay_final_cast <- reshape2::dcast(bioassay_final_melt, 
                                       Plant_ID ~ variable + Leaf_rep,
                                       value.var = "value")

## Now, we want to add leaf area details to the wide-format RGR dataframe.
# Add initial leaf area column to bioassay_RGR_cast dataframe.
bioassay_RGR_cast$Initial_leaf_area <- bioassay_initial[match(
  bioassay_RGR_cast$Plant_ID, bioassay_initial$Notes), 'Total_Leaf_Area_.cm2.']

# Add columns for final leaf area, leaf area consumed (for each leaf rep),
# and total area consumed to bioassay_RGR_cast dataframe.
bioassay_RGR_cast[,c('Area_consumed_1', 'Area_consumed_2')] <-
  bioassay_final_cast[match(bioassay_RGR_cast$Plant_ID, bioassay_final_cast$Plant_ID),
                 c('Consumed_Leaf_Area_.cm2._1', 'Consumed_Leaf_Area_.cm2._2')]

# Average leaf area consumed across the two leaf area reps for each plant.
bioassay_RGR_cast$Mean_area_consumed <- rowMeans(bioassay_RGR_cast[,c('Area_consumed_1', 'Area_consumed_2')],
                                                 na.rm = TRUE)
# Converted bioassay caterpillar mortality to 0 (survived) and 1 (died)
bioassay_RGR_cast$Dead_1 <- with(bioassay_RGR_cast,
                                 as.numeric(as.factor(Dead_1)) - 1)
bioassay_RGR_cast$Dead_2 <- with(bioassay_RGR_cast,
                                 as.numeric(as.factor(Dead_2)) - 1)

# Save cleaned dataframe as csv.
write.csv(bioassay_RGR_cast, "cleaned-data/bioassay-data.csv",
          row.names = FALSE)

rm(bioassay_RGR, bioassay_RGR_melt, bioassay_RGR_cast,
   bioassay_initial, bioassay_final, bioassay_final_melt, bioassay_final_cast)
