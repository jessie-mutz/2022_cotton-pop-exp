## USDA cotton experiment
## Data analysis: KCV of resistance and defenses ~ herbivore density
## Code written and annotated by JM (last updated 2026-03-26)

# The code in this script file:
# - standardizes plant defense and resistance measurements to mean leaf size and 
#   uses these to calculate KCV for each experimental population;
# - runs the model in 'scripts/kcv-models.R' for experimental population KCV of
#   resistance (as bioassay caterpillar RGR and leaf area consumed Y/N and non-zero),
#   defenses (trichome density, pigment gland proportional leaf area, HGQ 3.6 
#   relative abundance, and HGQ 10.5 relative abundance) and damage (at sampling
#   and at harvest) using data sets determined by n_threshold = 4, 5, 6, or 7
#   (where n_threshold is the minimum number of plants with data for that response
#   variable required for 'population' inclusion in analysis);
# - saves as output:
#    *.txt file of data used to run the model;
#    *.csv summary of model posterior for all monitored parameters;
#    *.csv model estimate of population mean response variable for all experimental 
#      populations;
#    *.pdf of trace plots for a subset of monitored parameters to verify convergence.
# **Note that this code uses functions in 'scripts/helper-funcs.R' to compile
#   data for models and summarize model output.

# Load required packages.
require(rjags)
require(jagsUI)
require(doBy)
require(bayestestR)

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("USDA-cotton-exp")

# The models used in the following analyses are specified in 'scripts/induced-
#  response-Bayesian-models.R':
source("scripts/kcv-models.R")

# Load helper functions for compiling data lists ('data_vecs_resp') and summarizing
# model output ('save_out')
source("scripts/helper-funcs.R")


## =============================================================================
## Data setup for KCV analyses
## =============================================================================

# Read in data for defenses, resistance, and damage.
full_data <- read.csv("cleaned-data/full-data.csv")

# Add a column for population ID (Cage_pos) as integers 1-70.
full_data$Pop_ID <- with(full_data, substr(Plant_ID - 1, start = 2, stop = 3))
full_data$Pop_ID <- with(full_data, as.numeric(Pop_ID) + 1)

# Scale leaf size measurements (initial leaf area or midvein length) for use as
#  covariates
full_data$Initial_leaf_area_sc <- scale(full_data$Initial_leaf_area)[,1]
full_data$Midvein_length_sc <- scale(full_data$Midvein_length)[,1]

# Read in summaries from single trait defense or resistance models. Estimates for
#  the effect of leaf size will be used to standardize plant-scale defense and 
#  resistance measurements.
rgr_summary <- read.csv("output/full-rgr-summary-min5.csv")
con_summary <- read.csv("output/full-con-summary-min5.csv")
tri_summary <- read.csv("output/full-tri-summary-min5.csv")
H10_summary <- read.csv("output/full-H10-summary-min5.csv")
H36_summary <- read.csv("output/full-H36-summary-min5.csv")
gna_summary <- read.csv("output/full-gland-area-summary-min5.csv")


## =============================================================================
## Standardize observations to mean leaf size
## =============================================================================

# For each response, determine whether there is an effect of leaf size
rgr_summary[which(rgr_summary$X == "slope_resp"), ]
con_summary[which(con_summary$X == "slope_resp"), ] #no effect here
tri_summary[which(tri_summary$X == "slope_resp"), ]
H10_summary[which(H10_summary$X == "slope_resp"), ] #no effect here
H36_summary[which(H36_summary$X == "slope_resp"), ]
gna_summary[which(gna_summary$X == "slope_resp"), ]

# For response variables with significant effect of leaf size, standardize trait 
#  values to mean leaf size
full_data$rgr_sizestd <- with(full_data, RGR_mean - 
                                rgr_summary[which(rgr_summary$X == "slope_resp"), 'mean']*Initial_leaf_area_sc)
full_data$con_sizestd <- with(full_data, ifelse(Mean_area_consumed == 0, NA,
                                                log(Mean_area_consumed)))
full_data$Trichomes_sizestd <- with(full_data, sqrt(Trichomes_std_cm2) - 
                                      tri_summary[which(tri_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$H_10.5_sizestd <- with(full_data, log(H_10.5_relabun))
full_data$H_3.6_sizestd <- with(full_data, log(H_3.6_relabun) - 
                                  H36_summary[which(H36_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)
full_data$Gland_prop_sizestd <- with(full_data, sqrt(Gland_prop) - 
                                       gna_summary[which(gna_summary$X == "slope_resp"), 'mean']*Midvein_length_sc)

# Get leaf size-standardized, untransformed values of each defenses or resistance variable.
full_dat2 <- full_data
full_dat2$Gland_prop_sizestd <- with(full_dat2, Gland_prop_sizestd^2)
full_dat2$H_10.5_sizestd <- with(full_dat2, exp(H_10.5_sizestd))
full_dat2$H_3.6_sizestd <- with(full_dat2, exp(H_3.6_sizestd))
full_dat2$con_sizestd <- with(full_dat2, exp(con_sizestd))
full_dat2$Trichomes_sizestd <- with(full_dat2, Trichomes_sizestd^2)


## =============================================================================
## Estimate KCVs for experimental populations
## =============================================================================

# Functions for estimating Kvalseth's CV from Lobry et al. 2023

sdn <- function(x, i, ...){
  n <- sum(!is.na(x)) # number of non missing values
  return(sqrt((n - 1)/n)*sd(x[i], ...))
}

PCV <- function(x, i, ...){
  barx <- mean(x[i], ...)
  if(isTRUE(all.equal(barx, 0))) return(0)
  return(sdn(x[i], ...)/barx)
}

KCV <- function(x, i, ...){
  PCV2 <- PCV(x, i, ...)^2
  return(sqrt(PCV2/(1 + PCV2)))
}

# Estimate KCV 
out_sum <- doBy::summaryBy(Gland_prop_sizestd + H_10.5_sizestd + H_3.6_sizestd +
                             rgr_sizestd + con_sizestd + Trichomes_sizestd +
                             Sampling_Percent_damage + Harvest_Percent_damage +
                             Plant_RGR ~ Cage_pos, 
                           id =~ Herb_add + Day + Timepoint, data = full_dat2, 
                           FUN = c(function(x) KCV(x, na.rm = TRUE), 
                                   function(x) sum(!is.na(x))),
                           fun.names = c("KCV", "notNA"),
                           var.names = c("gna", "H10", "H36", "rgr", "con", "tri",
                                         "samp_dam", "harv_dam", "plant_rgr"))

# Save KCV estimates.
write.csv(out_sum, 'output/KCV-estimates.csv', row.names = FALSE)

## =============================================================================
## Analysis for KCV ~ herbivore density
## =============================================================================

# Parameters for model fitting (set globally):
nc <- 8 #number of chains
ni <- 15000 #number of iterations
nb <- 5000 #burn-in
nt <- 2 #thinning rate

# Parameters to estimate from the model.
params.est <- c( "intPop", "slope_densPop", "pred")

thresholds <- 4:7

## Bioassay caterpillar RGR ====================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'rgr.KCV', 'rgr.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
  mod_df <- mod_df[which(mod_df$rgr.notNA >= thresholds[i]),]
  
# Assign data vectors.
nTimept <- length(unique(mod_df$Timepoint)) #number of timepoints
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
timeptPop <- mod_df$Timepoint
respPop <- mod_df$rgr.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
rgr_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("RGR", paste0("Min. plants = ", thresholds[i]), rgr_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-rgr-plots-min", thresholds[i], ".pdf"))
plot(rgr_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
rgr_out <- save_out(model = rgr_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(rgr_out[['samples']], paste0("output/kcv-rgr-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(rgr_out[['summary']], paste0("output/kcv-rgr-summary-min",
                                       thresholds[i], ".csv"))
rm(rgr_out)

}

## Bioassay leaf area consumed =================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'con.KCV', 'con.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$con.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$con.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
con_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Leaf area consumed (cont)", paste0("Min. plants = ", 
                                               thresholds[i]), con_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-con-plots-min", thresholds[i], ".pdf"))
plot(con_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
con_out <- save_out(model = con_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(con_out[['samples']], paste0("output/kcv-con-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(con_out[['summary']], paste0("output/kcv-con-summary-min",
                                       thresholds[i], ".csv"))
rm(con_out)

}

## Hemigossypolone 10.5 ========================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'H10.KCV', 'H10.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$H10.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$H10.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
H10_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("HGQ 10.5", paste0("Min. plants = ", thresholds[i]), H10_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-H10-plots-min", thresholds[i], ".pdf"))
plot(H10_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
H10_out <- save_out(model = H10_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(H10_out[['samples']], paste0("output/kcv-H10-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(H10_out[['summary']], paste0("output/kcv-H10-summary-min",
                                       thresholds[i], ".csv"))
rm(H10_out)

}

## Hemigossypolone 3.6 =========================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'H36.KCV', 'H36.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$H36.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$H36.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
H36_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("HGQ 3.6", paste0("Min. plants = ", thresholds[i]), H36_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-H36-plots-min", thresholds[i], ".pdf"))
plot(H36_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
H36_out <- save_out(model = H36_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(H36_out[['samples']], paste0("output/kcv-H36-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(H36_out[['summary']], paste0("output/kcv-H36-summary-min",
                                       thresholds[i], ".csv"))
rm(H36_out)

}

## Terpenoid gland density =====================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'gna.KCV', 'gna.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$gna.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$gna.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
gna_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Gland prop. area", paste0("Min. plants = ", thresholds[i]), gna_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-gna-plots-min", thresholds[i], ".pdf"))
plot(gna_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
gna_out <- save_out(model = gna_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(gna_out[['samples']], paste0("output/kcv-gland-area-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(gna_out[['summary']], paste0("output/kcv-gland-area-summary-min",
                                       thresholds[i], ".csv"))
rm(gna_out)

}

## Trichome density ============================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'tri.KCV', 'tri.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$tri.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$tri.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
tri_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Trichome density", paste0("Min. plants = ", thresholds[i]), tri_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-tri-plots-min", thresholds[i], ".pdf"))
plot(tri_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
tri_out <- save_out(model = tri_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(tri_out[['samples']], paste0("output/kcv-tri-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(tri_out[['summary']], paste0("output/kcv-tri-summary-min",
                                       thresholds[i], ".csv"))
rm(tri_out)

}

## Percent damage at sampling ==================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'samp_dam.KCV', 'samp_dam.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$samp_dam.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$samp_dam.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
dam_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Percent damage at sampling", 
           paste0("Min. plants = ", thresholds[i]), dam_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-dam-plots-min", thresholds[i], ".pdf"))
plot(dam_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
dam_out <- save_out(model = dam_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(dam_out[['samples']], paste0("output/kcv-damage-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(dam_out[['summary']], paste0("output/kcv-damage-summary-min",
                                       thresholds[i], ".csv"))
rm(dam_out)

}

## Percent damage at harvest ===================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'harv_dam.KCV', 'harv_dam.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$harv_dam.notNA >= thresholds[i]),]

# Assign data vectors.
nPop <- length(mod_df[,1]) #number of populations
densPop <- mod_df$Herb_add
respPop <- mod_df$harv_dam.KCV

dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                  nDens = 4, densPop_newdat = c(0, 6, 12, 25) )

## (2) Run model and check convergence
har_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Percent damage at harvest", 
           paste0("Min. plants = ", thresholds[i]), har_kcv$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/kcv-har-plots-min", thresholds[i], ".pdf"))
plot(har_kcv, parameters = params.est[1:2])

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
har_out <- save_out(model = har_kcv, dataframe = mod_df, pop_means = FALSE)

# Save each list element as .csv in 'output' folder.
write.csv(har_out[['samples']], paste0("output/kcv-harvest-damage-samples-min",
                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(har_out[['summary']], paste0("output/kcv-harvest-damage-summary-min",
                                       thresholds[i], ".csv"))

rm(har_out)

}

## Plant RGR ===================================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- out_sum[,c('Timepoint', 'Herb_add', 'plant_rgr.KCV', 'plant_rgr.notNA')]
mod_df <- na.omit(mod_df)

# Run models and save output, varying the threshold for including a population from 4-7 plants.

for (i in 1:length(thresholds)) {
  
  # Include only populations with at least thresholds[i] plant observations.
  mod_df <- mod_df[which(mod_df$plant_rgr.notNA >= thresholds[i]),]
  
  # Assign data vectors.
  nPop <- length(mod_df[,1]) #number of populations
  densPop <- mod_df$Herb_add
  respPop <- mod_df$plant_rgr.KCV
  
  dat.list <- list( nPop = nPop, respPop = respPop, densPop = densPop, 
                    nDens = 4, densPop_newdat = c(0, 6, 12, 25) )
  
  ## (2) Run model and check convergence
  plg_kcv <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/NoTimexDens-KCV-beta.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)
  
  print(list("Plant RGR", 
             paste0("Min. plants = ", thresholds[i]), plg_kcv$summary))
  
  # Check traceplots and density plots to visually confirm convergence
  pdf(file = paste0("output/kcv-plg-plots-min", thresholds[i], ".pdf"))
  plot(plg_kcv, parameters = params.est[1:2])
  
  ## (3) Save output
  # Use 'save_out' function to get posterior samples, summary statistics of 
  #  monitored variables, and population mean estimates.
  plg_out <- save_out(model = plg_kcv, dataframe = mod_df, pop_means = FALSE)
  
  # Save each list element as .csv in 'output' folder.
  write.csv(plg_out[['samples']], paste0("output/kcv-plant-rgr-damage-samples-min",
                                         thresholds[i], ".csv"), row.names = FALSE)
  write.csv(plg_out[['summary']], paste0("output/kcv-plant-rgr-damage-summary-min",
                                         thresholds[i], ".csv"))
  
  rm(plg_out)
  
}
