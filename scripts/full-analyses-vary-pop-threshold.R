## USDA cotton experiment
## Data analysis: Resistance and defenses ~ herbivore density
## Code written and annotated by JM (last updated 2026-03-25)

# The code in this script file:
# - runs the models specified by Eq. 1 and in 'scripts/induced-response-Bayesian-models.R'
#   for experimental G. hirsutum population resistance (as bioassay caterpillar
#   RGR and leaf area consumed Y/N and non-zero) and defenses (trichome density, 
#   pigment gland proportional leaf area, HGQ 3.6 relative abundance, and HGQ 10.5 
#   relative abundance) using data sets determined by n_threshold = 4, 5, 6, or 7
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

attr(scale(full_data$Initial_leaf_area), "scaled:center")
attr(scale(full_data$Midvein_length), "scaled:center")

# The models used in the following analyses are specified in 'scripts/induced-
#  response-Bayesian-models.R':
source("scripts/induced-response-Bayesian-models.R")

# Load helper functions for compiling data lists ('data_vecs_resp') and summarizing
# model output ('save_out')
source("scripts/helper-funcs.R")

# Parameters for model fitting (set globally):
nc <- 8 #number of chains
ni <- 15000 #number of iterations
nb <- 5000 #burn-in
nt <- 2 #thinning rate

# Parameters to estimate from the model.
params.fit <- c( "intPop", "slope_densPop", "slope_resp" ) #for checking convergence
#params.est <- c( "intPop", "slope_densPop", "slope_resp", "respPop" )
params.est <- c( "intPop", "slope_densPop", "slope_resp", "respPop",
                 "densT1T2", "densT1T3", "densT1T4", "densT1T5", 
                 "densT2T3", "densT2T4", "densT2T5", "densT3T4", 
                 "densT3T5", "densT4T5", "pred")


## =============================================================================
## Effects of herbivore density on resistance
## =============================================================================

## Bioassay caterpillar RGR ====================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'RGR_mean', 'Initial_leaf_area_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- mod_df$RGR_mean 
dat.list[['leaf_size']] <- mod_df$Initial_leaf_area_sc

## (2) Run model and check convergence
rgr_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("RGR", paste0("Min. plants = ", thresholds[i]), rgr_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-rgr-plots-min", thresholds[i], ".pdf"))
plot(rgr_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
rgr_out <- save_out(model = rgr_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-rgr-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(rgr_out[['samples']], paste0("output/full-rgr-samples-min", 
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(rgr_out[['summary']], paste0("output/full-rgr-summary-min",
                                       thresholds[i], ".csv"))
write.csv(rgr_out[['pop_means']], paste0("output/full-rgr-pop-means-min",
                                         thresholds[i], ".csv"))
rm(rgr_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(mod_df$RGR_mean)[,1]

rgr_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

rgr_out <- save_out(model = rgr_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(rgr_out[['summary']], paste0("output/scaled-rgr-summary-min",
                                       thresholds[i], ".csv"))
write.csv(rgr_out[['pop_means']], paste0("output/scaled-rgr-pop-means-min",
                                         thresholds[i], ".csv"))
rm(rgr_out)

}

## Bioassay leaf area consumed =================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[which(full_data$Mean_area_consumed > 0),]
mod_df <- mod_df[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                    'Mean_area_consumed', 'Initial_leaf_area_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- log(mod_df$Mean_area_consumed)
dat.list[['leaf_size']] <- mod_df$Initial_leaf_area_sc

## (2) Run model and check convergence
con_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Leaf area consumed (cont)", paste0("Min. plants = ", thresholds[i]), con_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-con-plots-min", thresholds[i], ".pdf"))
plot(con_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
con_out <- save_out(model = con_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-con-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(con_out[['samples']], paste0("output/full-con-samples-min",
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(con_out[['summary']], paste0("output/full-con-summary-min",
                                       thresholds[i], ".csv"))
write.csv(con_out[['pop_means']], paste0("output/full-con-pop-means-min",
                                         thresholds[i], ".csv"))
rm(con_out)

## (4) Rerun model with scaled values of leaf area consumed & save output
dat.list[['resp']] <- scale(log(mod_df$Mean_area_consumed))[,1]

con_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

con_out <- save_out(model = con_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(con_out[['summary']], paste0("output/scaled-con-summary-min",
                                       thresholds[i], ".csv"))
write.csv(con_out[['pop_means']], paste0("output/scaled-con-pop-means-min",
                                         thresholds[i], ".csv"))
rm(con_out)
}

## Bioassay leaf consumed Y/N ==================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[-which(full_data$Obs_any_consumed == 0),]
mod_df <- mod_df[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                    'Total_any_consumed', 'Obs_any_consumed', 'Initial_leaf_area_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- mod_df$Total_any_consumed
dat.list[['leaf_size']] <- mod_df$Initial_leaf_area_sc
dat.list[['num_trials']] <- mod_df$Obs_any_consumed

## (2) Run model and check convergence
conBin_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                             model.file = "scripts/TimexDens-LinThr-Geno-bin.txt", 
                             n.thin = nt, n.chains = nc, 
                             n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Leaf area consumed (bin)", paste0("Min. plants = ", thresholds[i]), conBin_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-conBin-plots-min", thresholds[i], ".pdf"))
plot(conBin_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
conBin_out <- save_out(model = conBin_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
#write.csv(conBin_out[['samples']], paste0("output/full-conBin-samples-min",
#                                          thresholds[i], ".csv"), row.names = FALSE)
write.csv(conBin_out[['summary']], paste0("output/full-conBin-summary-min",
                                          thresholds[i], ".csv"))
write.csv(conBin_out[['pop_means']], paste0("output/full-conBin-pop-means-min",
                                            thresholds[i], ".csv"))
rm(conBin_out)

}

## =============================================================================
## Effects of herbivore density on defense traits
## =============================================================================

## Trichome density ============================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'Trichomes_std_cm2', 'Midvein_length_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- sqrt(mod_df$Trichomes_std_cm2) #response variable at plant scale
dat.list[['leaf_size']] <- mod_df$Midvein_length_sc #leaf size as covariate

## (2) Run model and check convergence
tri_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Trichome density", paste0("Min. plants = ", thresholds[i]), 
           tri_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-tri-plots-min", thresholds[i], ".pdf"))
plot(tri_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
tri_out <- save_out(model = tri_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-tri-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(tri_out[['samples']], paste0("output/full-tri-samples-min",
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(tri_out[['summary']], paste0("output/full-tri-summary-min", 
                                       thresholds[i], ".csv"))
write.csv(tri_out[['pop_means']], paste0("output/full-tri-pop-means-min",
                                         thresholds[i], ".csv"))
rm(tri_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(sqrt(mod_df$Trichomes_std_cm2))[,1]

tri_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

tri_out <- save_out(model = tri_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(tri_out[['summary']], paste0("output/scaled-tri-summary-min",
                                       thresholds[i], ".csv"))
write.csv(tri_out[['pop_means']], paste0("output/scaled-tri-pop-means-min",
                                         thresholds[i], ".csv"))
rm(tri_out)

}

## Hemigossypolone 10.5 ========================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'H_10.5_relabun', 'Midvein_length_sc')]
mod_df <- na.omit(mod_df)

# View table of how many plants per population have associated data.
table(mod_df$Pop_ID)

# Run models and save output, varying the threshold for including a population from 4-7 plants.
thresholds <- 4:6 #using 6 instead of 7 because data are too sparse with 7

for (i in 1:length(thresholds)) {

# Include only populations with at least thresholds[i] plant observations.
mod_df <- mod_df[which(mod_df$Pop_ID %in% names(which(table(mod_df$Pop_ID) >= thresholds[i]))),]

# Get vectors of sequential integers, one designating plant genotype and one
#  designating population ID.
mod_df$Geno_ordered <- with(mod_df, as.numeric(as.factor(Geno)))
mod_df$Pop_ordered <- with(mod_df, as.numeric(as.factor(Pop_ID)))

# Assign data vectors
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- log(mod_df$H_10.5_relabun) #response variable at plant scale
dat.list[['leaf_size']] <- mod_df$Midvein_length_sc #leaf size as covariate

## (2) Run model and check convergence
H10_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("HGQ10.5", paste0("Min. plants = ", thresholds[i]), H10_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-H10-plots-min", thresholds[i], ".pdf"))
plot(H10_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
H10_out <- save_out(model = H10_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-H10-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(H10_out[['samples']], paste0("output/full-H10-samples-min",
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(H10_out[['summary']], paste0("output/full-H10-summary-min",
                                       thresholds[i], ".csv"))
write.csv(H10_out[['pop_means']], paste0("output/full-H10-pop-means-min",
                                         thresholds[i], ".csv"))
rm(H10_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(log(mod_df$H_10.5_relabun))[,1]

H10_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

H10_out <- save_out(model = H10_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(H10_out[['summary']], paste0("output/scaled-H10-summary-min",
                                       thresholds[i], ".csv"))
write.csv(H10_out[['pop_means']], paste0("output/scaled-H10-pop-means-min",
                                         thresholds[i], ".csv"))
rm(H10_out)

}


## Hemigossypolone 3.6 =========================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'H_3.6_relabun', 'Midvein_length_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- log(mod_df$H_3.6_relabun) #response variable at plant scale
dat.list[['leaf_size']] <- mod_df$Midvein_length_sc #leaf size as covariate

## (2) Run model and check convergence
H36_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("HGQ3.6", paste0("Min. plants = ", thresholds[i]), H36_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-H36-plots-min", thresholds[i], ".pdf"))
plot(H36_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
H36_out <- save_out(model = H36_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-H36-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(H36_out[['samples']], paste0("output/full-H36-samples-min", 
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(H36_out[['summary']], paste0("output/full-H36-summary-min", 
                                       thresholds[i], ".csv"))
write.csv(H36_out[['pop_means']], paste0("output/full-H36-pop-means-min",
                                         thresholds[i], ".csv"))
rm(H36_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(log(mod_df$H_3.6_relabun))[,1]
H36_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

H36_out <- save_out(model = H36_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(H36_out[['summary']], paste0("output/scaled-H36-summary-min",
                                       thresholds[i], ".csv"))
write.csv(H36_out[['pop_means']], paste0("output/scaled-H36-pop-means-min",
                                         thresholds[i], ".csv"))
rm(H36_out)

}

## Terpenoid gland density =====================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'Glands_std_cm2', 'Midvein_length_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- sqrt(mod_df$Glands_std_cm2) #response variable at plant scale
dat.list[['leaf_size']] <- mod_df$Midvein_length_sc #leaf size as covariate

## (2) Run model and check convergence
gnd_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Gland density", paste0("Min. plants = ", thresholds[i]), gnd_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-gnd-plots-min", thresholds[i], ".pdf"))
plot(gnd_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
gnd_out <- save_out(model = gnd_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-gnd-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(gnd_out[['samples']], paste0("output/full-gland-dens-samples-min",
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(gnd_out[['summary']], paste0("output/full-gland-dens-summary-min",
                                       thresholds[i], ".csv"))
write.csv(gnd_out[['pop_means']], paste0("output/full-gland-dens-pop-means-min",
                                         thresholds[i], ".csv"))
rm(gnd_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(sqrt(mod_df$Glands_std_cm2))[,1]

gnd_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

gnd_out <- save_out(model = gnd_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(gnd_out[['summary']], paste0("output/scaled-gland-dens-summary-min",
                                       thresholds[i], ".csv"))
write.csv(gnd_out[['pop_means']], paste0("output/scaled-gland-dens-pop-means-min",
                                         thresholds[i], ".csv"))
rm(gnd_out)

}

## Terpenoid gland area ========================================================

## (1) Data prep
# Subset the full data frame to columns that will be included in the model, 
#  then remove any rows with NAs.
mod_df <- full_data[,c('Pop_ID', 'Plant_ID', 'Geno', 'Timepoint', 'Herb_add', 
                       'Gland_prop', 'Midvein_length_sc')]
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
dat.list <- data_vecs_resp(df = mod_df)
dat.list[['resp']] <- sqrt(mod_df$Gland_prop) #response variable at plant scale
dat.list[['leaf_size']] <- mod_df$Midvein_length_sc #leaf size as covariate

## (2) Run model and check convergence
gna_1_lin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

print(list("Gland prop. area", paste0("Min. plants = ", thresholds[i]), 
           gna_1_lin$summary))

# Check traceplots and density plots to visually confirm convergence
pdf(file = paste0("output/full-gna-plots-min", thresholds[i], ".pdf"))
plot(gna_1_lin, parameters = params.fit)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
gna_out <- save_out(model = gna_1_lin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink(paste0("output/full-gna-dat-min", thresholds[i], ".txt"))
print(dat.list)
sink()

#write.csv(gna_out[['samples']], paste0("output/full-gland-area-samples-min",
#                                       thresholds[i], ".csv"), row.names = FALSE)
write.csv(gna_out[['summary']], paste0("output/full-gland-area-summary-min",
                                       thresholds[i], ".csv"))
write.csv(gna_out[['pop_means']], paste0("output/full-gland-area-pop-means-min",
                                         thresholds[i], ".csv"))
rm(gna_out)

## (4) Rerun model with scaled values of RGR & save output
dat.list[['resp']] <- scale(sqrt(mod_df$Gland_prop))[,1]

gna_sc <- jagsUI:: jags(data = dat.list, parameters.to.save = params.est,
                        model.file = "scripts/TimexDens-LinThr-Geno.txt", 
                        n.thin = nt, n.chains = nc, 
                        n.burnin = nb, n.iter = ni, DIC = TRUE)

gna_out <- save_out(model = gna_sc, dataframe = mod_df, pop_means = TRUE)
write.csv(gna_out[['summary']], paste0("output/scaled-gland-area-summary-min",
                                       thresholds[i], ".csv"))
write.csv(gna_out[['pop_means']], paste0("output/scaled-gland-area-pop-means-min",
                                         thresholds[i], ".csv"))
rm(gna_out)

}

