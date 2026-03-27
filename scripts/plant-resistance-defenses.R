## USDA cotton experiment
## Data analysis: Resistance ~ defenses
## Code written and annotated by JM (last updated 2026-03-26)

# The code in this script file:
# - runs multiple regression models for plant-scale observations of resistance 
#   (as bioassay caterpillar RGR and non-zero leaf area consumed) as a function
#   of plant-scale defenses (trichome density, pigment gland proportional leaf area, 
#   HGQ 3.6 relative abundance, and HGQ 10.5 relative abundance);
# - saves as output:
#    *.txt file of data used to run the model;
#    *.csv summary of model posterior for all monitored parameters;
#    *.csv model estimate of population mean response variable for all experimental 
#      populations;
#    *.pdf of trace plots for a subset of monitored parameters to verify convergence.
# **Note that this code uses functions in 'scripts/helper-funcs.R' to compile
#   data for models and summarize model output.
# **This code takes a long time to run. JM ran the code using the Dartmouth
#   Discovery cluster.

# Load required packages.
require(rjags)
require(jagsUI)
require(doBy)
require(bayestestR)

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("USDA-cotton-exp")

# Read in data.
full_data <- read.csv("cleaned-data/full-data.csv")

# Add a column for population ID (Cage_pos) as integers 1-70.
full_data$Pop_ID <- with(full_data, substr(Plant_ID - 1, start = 2, stop = 3))
full_data$Pop_ID <- with(full_data, as.numeric(Pop_ID) + 1)

# Scale leaf size measurements (initial leaf area or midvein length) for use as
#  covariates
full_data$Initial_leaf_area_sc <- scale(full_data$Initial_leaf_area)[,1]
full_data$Midvein_length_sc <- scale(full_data$Midvein_length)[,1]

# Load models
source('scripts/resistance-defenses-models.R')

# Load helper functions for compiling data lists ('data_vecs_resp') and summarizing
# model output ('save_out')
source("scripts/helper-funcs.R")

## =============================================================================
## Plant scale
## =============================================================================

# Specify model parameters.
nc <- 12 #number of chains
ni <- 75000 #number of iterations
nb <- 10000 #burn-in
nt <- 2 #thinning rate

# Parameters to estimate.
params.est <- c( "int", "slope_resp", "slope_tri", "slope_gnd", 
                 "slope_H10", "slope_H36" )

## RGR model ===================================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'RGR_mean', 'Trichomes_std_cm2', 'Gland_prop',
                       'H_10.5_relabun', 'H_3.6_relabun', 'Initial_leaf_area_sc')]
mod_df <- na.omit(mod_df)

# Get vectors of sequential integers, one designating plant genotype and one
#  designating population ID.
mod_df$Geno_ordered <- with(mod_df, as.numeric(as.factor(Geno)))
mod_df$Pop_ordered <- with(mod_df, as.numeric(as.factor(Pop_ID)))

# Assign data vectors
nGeno <- with(mod_df, length(unique(Geno)))
nTimept <- with(mod_df, length(unique(Timepoint)))
nObs <- dim(mod_df)[1] #number of observations

resp <- mod_df$RGR_mean #resistance response variable at plant scale

tri <- with(mod_df, sqrt(Trichomes_std_cm2)) #trichome density
gnd <- with(mod_df, sqrt(Gland_prop))
H10 <- with(mod_df, log(H_10.5_relabun)) #hemigossypolone 10.5
H36 <- with(mod_df, log(H_3.6_relabun)) #hemigossypolone 3.6
leaf_size <- with(mod_df, Initial_leaf_area_sc)

timept <- mod_df$Timepoint
geno <- mod_df$Geno_ordered

# List of data vectors for models.
dat.list <- list( nObs = nObs, nTimept = nTimept, nGeno = nGeno,
                  geno = geno, timept = timept, leaf_size = leaf_size,
                  resp = resp, tri = tri, gnd = gnd, H10 = H10, H36 = H36
)

## (2) Run model and check convergence
plant_rgr_defen <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                                model.file = "scripts/plant-resistance-defense.txt", 
                                n.thin = nt, n.chains = nc, 
                                n.burnin = nb, n.iter = ni, DIC = TRUE)

print("RGR ~ defenses (plant scale)")
plant_rgr_defen

# Check traceplots and density plots to visually confirm convergence
pdf(file = "output/plant-rgr-defen.pdf")
plot(plant_rgr_defen)

## (3) Save output
sink("output/plant-rgr-defen-dat.txt")
print(dat.list)
sink()

plant_rgr_defen_out <- save_out(model = plant_rgr_defen, dataframe = mod_df, pop_means = FALSE)
write.csv(plant_rgr_defen_out[['samples']], "output/plant-rgr-defen-samples.csv")
write.csv(plant_rgr_defen_out[['summary']], "output/plant-rgr-defen-summary.csv")
rm(plant_rgr_defen_out)


## Leaf area consumed model ====================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[which(full_data$Mean_area_consumed > 0),]
mod_df <- mod_df[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                    'Mean_area_consumed', 'Trichomes_std_cm2', 'Gland_prop',
                    'H_10.5_relabun', 'H_3.6_relabun')]
mod_df <- na.omit(mod_df)

# Get vectors of sequential integers, one designating plant genotype and one
#  designating population ID.
mod_df$Pop_ordered <- with(mod_df, as.numeric(as.factor(Pop_ID)))

# Assign data vectors
nTimept <- with(mod_df, length(unique(Timepoint)))
nObs <- dim(mod_df)[1] #number of observations

resp <- with(mod_df, log(Mean_area_consumed)) #resistance response variable at plant scale

tri <- with(mod_df, sqrt(Trichomes_std_cm2)) #trichome density
gnd <- with(mod_df, sqrt(Gland_prop))
H10 <- with(mod_df, log(H_10.5_relabun)) #hemigossypolone 10.5
H36 <- with(mod_df, log(H_3.6_relabun)) #hemigossypolone 3.6

timept <- mod_df$Timepoint

# List of data vectors for models.
dat.list <- list( nObs = nObs, nTimept = nTimept, nGeno = nGeno,
                  geno = geno, timept = timept, leaf_size = leaf_size,
                  resp = resp, tri = tri, gnd = gnd, H10 = H10, H36 = H36
)

## (2) Run model and check convergence
plant_con_defen <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                                model.file = "scripts/plant-resistance-defense.txt", 
                                n.thin = nt, n.chains = nc, 
                                n.burnin = nb, n.iter = ni, DIC = TRUE)

print("Leaf area consumed (cont) ~ defenses (plant scale)")
plant_con_defen

# Check traceplots and density plots to visually confirm convergence
pdf(file = "output/plant-con-defen.pdf")
plot(plant_con_defen)

## (3) Save output
sink("output/plant-con-defen-dat.txt")
print(dat.list)
sink()

plant_con_defen_out <- save_out(model = plant_con_defen, dataframe = mod_df, pop_means = FALSE)
write.csv(plant_con_defen_out[['samples']], "output/plant-con-defen-samples.csv")
write.csv(plant_con_defen_out[['summary']], "output/plant-con-defen-summary.csv")
rm(plant_con_defen_out)
