## USDA cotton experiment
## Data analysis: Resistance ~ defenses
## Code written and annotated by JM (last updated 2026-03-26)

# The code in this script file:
# - runs hierarchical models specified in Eq. 2 and in 'scripts/resistance-defense-models'
#   for relationship between mean population resistance (as bioassay caterpillar 
#   RGR and non-zero leaf area consumed) and mean population defenses (trichome 
#   density, pigment gland proportional leaf area, HGQ 3.6 relative abundance, 
#   and HGQ 10.5 relative abundance);
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
setwd("USDA-cotton-exp")

# Read in data.
full_data <- read.csv("cleaned-data/full-data.csv")

# Add a column for population ID (Cage_pos) as integers 1-70.
full_data$Pop_ID <- with(full_data, substr(Plant_ID - 1, start = 2, stop = 3))
full_data$Pop_ID <- with(full_data, as.numeric(Pop_ID) + 1)

# Make a guide to the sampling timepoint and herbivore density associated with
#  each experimental population.
pop_guide <- doBy::summaryBy(Day ~ Pop_ID,
                             id = c("Timepoint", "Herb_add", "Cage_pos"),
                             data = full_data,
                             keep.names = TRUE)

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
## Population scale: full path
## =============================================================================

# Specify model fitting parameters.
nc <- 15 #number of chains
ni <- 150000 #number of iterations
nb <- 15000 #burn-in
nt <- 5 #thinning length

# Parameters to estimate
params.est <- c( "intPop", "slope_densPop", 
                 "slope_triPop", "slope_gndPop", 
                 "slope_H10Pop", "slope_H36Pop",
                 "slope_densPop_tri", "slope_densPop_gnd",
                 "slope_densPop_H10", "slope_densPop_H36",
                 "respPop", "tri_respPop", "gnd_respPop", "H10_respPop", "H36_respPop")

params.fit <- c( "intPop", "slope_densPop", 
                 "slope_triPop", "slope_gndPop", 
                 "slope_H10Pop", "slope_H36Pop",
                 "slope_densPop_tri", "slope_densPop_gnd",
                 "slope_densPop_H10", "slope_densPop_H36")


# Scale defenses to entire dataset (before subsetting)
pop_path_data <- full_data
pop_path_data$Trichomes_sc <- with(pop_path_data, scale(sqrt(Trichomes_std_cm2))[,1])
pop_path_data$Gland_prop_sc <- with(pop_path_data, scale(sqrt(Gland_prop))[,1])
pop_path_data$H_10.5_sc <- with(pop_path_data, scale(log(H_10.5_relabun))[,1])
pop_path_data$H_3.6_sc <- with(pop_path_data, scale(log(H_3.6_relabun))[,1])

## Bioassay leaf area consumed model ===========================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- pop_path_data[which(pop_path_data$Mean_area_consumed > 0),]
mod_df <- mod_df[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                    'Mean_area_consumed', 'Trichomes_sc', 'Gland_prop_sc',
                    'H_10.5_sc', 'H_3.6_sc',
                    'Midvein_length_sc', 'Initial_leaf_area_sc')]
mod_df <- na.omit(mod_df)

# View table of how many plants per population have associated data.
table(mod_df$Pop_ID)

# Run models and save output, varying the threshold for including a population from 4-7 plants.
thresholds <- 4:7

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
  mod_df <- mod_df[which(mod_df$Pop_ID %in% names(which(table(mod_df$Pop_ID) >= thresholds[i]))),]

# Get vectors of sequential integers, one designating plant genotype and one
#  designating population ID.
mod_df$Geno_ordered <- with(mod_df, as.numeric(as.factor(Geno)))
mod_df$Pop_ordered <- with(mod_df, as.numeric(as.factor(Pop_ID)))

# Assign data vectors
nGeno <- with(mod_df, length(unique(Geno))) #number of genotypes
nTimept <- with(mod_df, length(unique(Timepoint))) #number of timepoints
nPop <- with(mod_df, length(unique(Pop_ID))) #number of populations
nObs <- dim(mod_df)[1] #number of observations

resp <- log(mod_df$Mean_area_consumed) #resistance response variable at plant scale
leaf_size <- mod_df$Midvein_length_sc #leaf size as covariate for defenses
leaf_size_resp <- mod_df$Initial_leaf_area_sc #leaf size as covariate for resistance

tri <- mod_df$Trichomes_sc #scaled trichome density
gnd <- mod_df$Gland_prop_sc #scaled terpenoid gland area
H10 <- mod_df$H_10.5_sc #scaled hemigossypolone 10.5
H36 <- mod_df$H_3.6_sc #scaled hemigossypolone 3.6

pop <- mod_df$Pop_ordered #population ID associated with each observation
geno <- mod_df$Geno_ordered #genotype associated with each observation

# Get dataframe with one row for each population, with columns for associated 
#  sampling timepoint and herbivore density.
pop_sum <- doBy::summaryBy( Timepoint + Herb_add ~ Pop_ID, data = mod_df,
                            keep.names = TRUE )

densPop <- pop_sum$Herb_add #herbivore density associated with each population
timeptPop <- pop_sum$Timepoint #sampling timepoint for each population

# List of data vectors
dat.list <- list( nGeno = nGeno, nTimept = nTimept, nPop = nPop, nObs = nObs,
                  resp = resp, leaf_size = leaf_size, pop = pop, geno = geno, 
                  leaf_size_resp = leaf_size_resp, tri = tri, gnd = gnd,
                  H10 = H10, H36 = H36,
                  densPop = densPop, timeptPop = timeptPop )

pop_con_defen4 <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                               model.file = "scripts/pop-resis-NoTimexDens-1.txt", 
                               n.thin = nt, n.chains = nc, 
                               n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Leaf area consumed (cont) ~ defenses (pop scale", 
           paste0("Min. plants = ", thresholds[i]), 
           pop_con_defen4$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/pop-con-defen-plots-min", thresholds[i], ".pdf"))
plot(pop_con_defen4, parameters = params.fit)

## (3) Save output
pop_con_defen_out <- save_out(model = pop_con_defen4, dataframe = mod_df,
                              pop_means = TRUE)
write.csv(pop_con_defen_out[['samples']], paste0("output/pop-con-defen-samples-min",
                                                 thresholds[i], ".csv"), row.names = FALSE)
write.csv(pop_con_defen_out[['summary']], paste0("output/pop-con-defen-summary-min",
                                                 thresholds[i], ".csv"))
write.csv(pop_con_defen_out[['pop_means']], paste0("output/pop-con-defen-pop-means-min",
                                                   thresholds[i], ".csv"))
rm(pop_con_defen_out)

}

## Bioassay caterpillar RGR model ==============================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- pop_path_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                           'RGR_mean', 'Trichomes_sc', 'Gland_prop_sc',
                           'H_10.5_sc', 'H_3.6_sc',
                           'Midvein_length_sc', 'Initial_leaf_area_sc')]
mod_df <- na.omit(mod_df)
# View table of how many plants per population have associated data.
table(mod_df$Pop_ID)

# Run models and save output, varying the threshold for including a population from 4-7 plants.
thresholds <- 4:7

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
  mod_df <- mod_df[which(mod_df$Pop_ID %in% names(which(table(mod_df$Pop_ID) >= thresholds[i]))),]
  
# Get vectors of sequential integers, one designating plant genotype and one
#  designating population ID.
mod_df$Geno_ordered <- with(mod_df, as.numeric(as.factor(Geno)))
mod_df$Pop_ordered <- with(mod_df, as.numeric(as.factor(Pop_ID)))

# Assign data vectors
nGeno <- with(mod_df, length(unique(Geno))) #number of genotypes
nTimept <- with(mod_df, length(unique(Timepoint))) #number of timepoints
nPop <- with(mod_df, length(unique(Pop_ID))) #number of populations
nObs <- dim(mod_df)[1] #number of observations

resp <- mod_df$RGR_mean #resistance response variable on plant scale
leaf_size <- mod_df$Midvein_length_sc #leaf size as covariate for defenses
leaf_size_resp <- mod_df$Initial_leaf_area_sc #leaf size as covariate for resistance

tri <- mod_df$Trichomes_sc #scaled trichome density
gnd <- mod_df$Gland_prop_sc #scaled terpenoid gland area
H10 <- mod_df$H_10.5_sc #scaled hemigossypolone 10.5
H36 <- mod_df$H_3.6_sc #scaled hemigossypolone 3.6

pop <- mod_df$Pop_ordered #population ID associated with each observation
geno <- mod_df$Geno_ordered #genotype associated with each observation

# Get dataframe with one row for each population, with columns for associated 
#  sampling timepoint and herbivore density.
pop_sum <- doBy::summaryBy( Timepoint + Herb_add ~ Pop_ID, data = mod_df,
                            keep.names = TRUE )

densPop <- pop_sum$Herb_add #herbivore density associated with each population
timeptPop <- pop_sum$Timepoint #sampling timepoint for each population

# List of data vectors
dat.list <- list( nGeno = nGeno, nTimept = nTimept, nPop = nPop, nObs = nObs,
                  resp = resp, leaf_size = leaf_size, pop = pop, geno = geno, 
                  leaf_size_resp = leaf_size_resp, tri = tri, gnd = gnd,
                  H10 = H10, H36 = H36,
                  densPop = densPop, timeptPop = timeptPop )

pop_rgr_defen4 <- jags(data = dat.list, parameters.to.save = params.est,
                       model.file = "scripts/pop-resis-NoTimexDens-1.txt", 
                       n.thin = nt, n.chains = nc, 
                       n.burnin = nb, n.iter = ni, DIC = TRUE)

pop_rgr_defen4

## (3) Save output
pop_rgr_defen_out <- save_out(model = pop_rgr_defen4, dataframe = mod_df,
                              pop_means = TRUE)
write.csv(pop_rgr_defen_out[['samples']], paste0("output/pop-rgr-defen-samples-min",
                                                 thresholds[i], ".csv"), row.names = FALSE)
write.csv(pop_rgr_defen_out[['summary']], paste0("output/pop-rgr-defen-summary-min",
                                                 thresholds[i], ".csv"))
write.csv(pop_rgr_defen_out[['pop_means']], paste0("output/pop-rgr-defen-pop-means-min",
                                                   thresholds[i], ".csv"))
rm(pop_rgr_defen_out)

}