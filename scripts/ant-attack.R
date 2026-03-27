## USDA cotton experiment
## Data analysis: Ant tending ~ herbivore density
## Code written and annotated by JM (last updated 2026-03-25)

# Load required packages.
require(rjags)
require(jagsUI)
require(doBy)
require(bayestestR)

# Set working directory to 'USDA-cotton-exp' repository:
#setwd("USDA-cotton-exp")

# Load helper functions for compiling data lists ('data_vecs_resp') and summarizing
# model output ('save_out')
source("scripts/helper-funcs.R")

## Specify model.===============================================================

sink("scripts/ant-attack-bin.txt")
cat("
model {
  
  for (i in 1:nObs) {
    resp[i] ~ dbin(prob_resp[i], num_trials[i])
    logit(prob_resp[i]) <- int + slope_dens*dens[i]
  }
  
  slope_dens ~ dnorm(0, 0.001)
  int ~ dnorm(0, 0.001)

  for (i in 1:nDens) {
    logit(pred[i]) <- int + slope_dens*dens_newdat[i]
  }

}
", fill = TRUE)
sink()

## Specify model.===============================================================

ants <- read.csv("cleaned-data/ant-attack.csv")

resp <- ants$No_plants_attacked
num_trials <- ants$No_plants_obs
dens <- as.numeric(ants$Herb_add)
nDens <- 3
nObs <- dim(ants)[1]
dens_newdat <- c(6, 12, 25)

dat.list <- list( resp = resp, num_trials = num_trials, dens = dens,
                  nDens = nDens, nObs = nObs, dens_newdat = dens_newdat )

params.est <- c( 'slope_dens', 'pred', 'int' )

# Parameters for model fitting (set globally):
nc <- 8 #number of chains
ni <- 15000 #number of iterations
nb <- 5000 #burn-in
nt <- 2 #thinning rate

## (2) Run model and check convergence
ant_1_bin <- jagsUI::jags(data = dat.list, parameters.to.save = params.est,
                          model.file = "scripts/ant-attack-bin.txt", 
                          n.thin = nt, n.chains = nc, 
                          n.burnin = nb, n.iter = ni, DIC = TRUE)

ant_1_bin

# Check traceplots and density plots to visually confirm convergence
pdf(file = "output/ant-attack.pdf")
plot(ant_1_bin, parameters = params.est)

## (3) Save output
# Use 'save_out' function to get posterior samples, summary statistics of 
#  monitored variables, and population mean estimates.
ant_out <- save_out(model = ant_1_bin, dataframe = mod_df, pop_means = TRUE)

# Save each list element as .csv in 'output' folder.
sink("output/ant-attack.txt")
print(dat.list)
sink()

write.csv(ant_out[['summary']], "output/ant-attack.csv")
write.csv(ant_out[['pop_means']], "output/ant-attack.csv")
rm(ant_out)
