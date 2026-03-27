## USDA cotton experiment
## Data analysis: Resistance ~ defenses Bayesian models
## Code written and annotated by JM (last updated 2026-03-26)

# The code in this script file:
# - Specifies models for relationship between defenses and resistance at the 
#   plant and population scales;
# - Saves models as .txt files.

## =============================================================================
## Plant-scale resistance ~ defenses multiple regression
## =============================================================================

sink("scripts/plant-resistance-defense.txt")
cat("
model {

  for (i in 1:nTimept) {
    int[i] ~ dnorm(mu_int, tau_int)
  }
  
#  for (i in 1:nGeno) {
#    int_geno[i] ~ dnorm(mu_int_geno, tau_int_geno)
#  }
  
  for (i in 1:nObs) {
    resp[i] ~ dnorm(mu_resp[i], tau_resp)
    mu_resp[i] <- int[timept[i]] + 
                  #int_geno[geno[i]] +
                  slope_tri*tri[i] + slope_gnd*gnd[i] + 
                  slope_H10*H10[i] + slope_H36*H36[i] + 
                  slope_resp*leaf_size[i]
  }
  
  slope_tri ~ dnorm(0, 0.001)
  slope_gnd ~ dnorm(0, 0.001)
  slope_H10 ~ dnorm(0, 0.001)
  slope_H36 ~ dnorm(0, 0.001)
  slope_resp ~ dnorm(0, 0.001)
  
  sigma_resp ~ dunif(0, 10)
  tau_resp <- 1/(sigma_resp*sigma_resp)
  
  mu_int ~ dnorm(0, 0.001)
  sigma_int ~ dunif(0,10)
  tau_int <- 1/(sigma_int*sigma_int)
  
#  mu_int_geno ~ dnorm(0, 0.001)
#  sigma_int_geno ~ dunif(0,10)
#  tau_int_geno <- 1/(sigma_int*sigma_int)

}
", fill = TRUE)
sink()


## =============================================================================
## Population-scale resistance ~ defenses
## =============================================================================

sink("scripts/pop-resis-NoTimexDens-1.txt")
cat("
model {

  for (i in 1:nGeno) {
    int_tri[i] ~ dnorm(0, tau_int_tri)
    int_gnd[i] ~ dnorm(0, tau_int_gnd)
    int_H10[i] ~ dnorm(0, tau_int_H10)
    int_H36[i] ~ dnorm(0, tau_int_H36)
    int_resp[i] ~ dnorm(0, tau_int_resp)
  }
  
  for (i in 1:nTimept) {
    intPop[i] ~ dnorm(0, 0.001)
    
    intPop_tri[i] ~ dnorm(0, 0.001)
    slope_densPop_tri[i] ~ dnorm(0, 0.001)
    
    intPop_gnd[i] ~ dnorm(0, 0.001)
    slope_densPop_gnd[i] ~ dnorm(0, 0.001)
    
    intPop_H10[i] ~ dnorm(0, 0.001)
    slope_densPop_H10[i] ~ dnorm(0, 0.001)
    
    intPop_H36[i] ~ dnorm(0, 0.001)
    slope_densPop_H36[i] ~ dnorm(0, 0.001)
  }
  
  for (i in 1:nPop) {
    respPop[i] ~ dnorm(mu_respPop[i], tau_respPop)
    mu_respPop[i] <- intPop[timeptPop[i]] + 
                     slope_triPop*tri_respPop[i] + 
                     slope_gndPop*gnd_respPop[i] +
                     slope_H10Pop*H10_respPop[i] + 
                     slope_H36Pop*H36_respPop[i]
                     
    tri_respPop[i] ~ dnorm(mu_triPop[i], tau_triPop)
    mu_triPop[i] <- intPop_tri[timeptPop[i]] + slope_densPop_tri[timeptPop[i]]*densPop[i]
    
    gnd_respPop[i] ~ dnorm(mu_gndPop[i], tau_gndPop)
    mu_gndPop[i] <- intPop_gnd[timeptPop[i]] + slope_densPop_gnd[timeptPop[i]]*densPop[i]
    
    H10_respPop[i] ~ dnorm(mu_H10Pop[i], tau_H10Pop)
    mu_H10Pop[i] <- intPop_H10[timeptPop[i]] + slope_densPop_H10[timeptPop[i]]*densPop[i]
    
    H36_respPop[i] ~ dnorm(mu_H36Pop[i], tau_H36Pop)
    mu_H36Pop[i] <- intPop_H36[timeptPop[i]] + slope_densPop_H36[timeptPop[i]]*densPop[i]
  }
  
  for (i in 1:nObs) {
    resp[i] ~ dnorm(mu_resp[i], tau_resp)
    mu_resp[i] <- respPop[pop[i]] + slope_resp*leaf_size_resp[i] + int_resp[geno[i]]
  
    tri[i] ~ dnorm(mu_tri[i], tau_tri)
    mu_tri[i] <- tri_respPop[pop[i]] + slope_tri*leaf_size[i] + int_tri[geno[i]]
  
    gnd[i] ~ dnorm(mu_gnd[i], tau_gnd)
    mu_gnd[i] <- gnd_respPop[pop[i]] + slope_gnd*leaf_size[i] + int_gnd[geno[i]]
    
    H10[i] ~ dnorm(mu_H10[i], tau_H10)
    mu_H10[i] <- H10_respPop[pop[i]] + slope_H10*leaf_size[i] + int_H10[geno[i]]
    
    H36[i] ~ dnorm(mu_H36[i], tau_H36)
    mu_H36[i] <- H36_respPop[pop[i]] + slope_H36*leaf_size[i] + int_H36[geno[i]]
  }
  
  # Priors for variance in estimate of intercept varying randomly among genotypes
  sigma_int_resp ~ dunif(0, 10)
  tau_int_resp <- 1/(sigma_int_resp*sigma_int_resp)
  sigma_int_tri ~ dunif(0, 10)
  tau_int_tri <- 1/(sigma_int_tri*sigma_int_tri)
  sigma_int_gnd ~ dunif(0, 10)
  tau_int_gnd <- 1/(sigma_int_gnd*sigma_int_gnd)
  sigma_int_H10 ~ dunif(0, 10)
  tau_int_H10 <- 1/(sigma_int_H10*sigma_int_H10)
  sigma_int_H36 ~ dunif(0, 10)
  tau_int_H36 <- 1/(sigma_int_H36*sigma_int_H36)
  
  # Priors for variance in plant-scale observations
  sigma_resp ~ dunif(0, 10)
  tau_resp <- 1/(sigma_resp*sigma_resp)
  sigma_tri ~ dunif(0, 10)
  tau_tri <- 1/(sigma_tri*sigma_tri)
  sigma_gnd ~ dunif(0, 10)
  tau_gnd <- 1/(sigma_gnd*sigma_gnd)
  sigma_H10 ~ dunif(0, 10)
  tau_H10 <- 1/(sigma_H10*sigma_H10)
  sigma_H36 ~ dunif(0, 10)
  tau_H36 <- 1/(sigma_H36*sigma_H36)
  
  # Priors for variance in population mean
  sigma_respPop ~ dunif(0, 10)
  tau_respPop <- 1/(sigma_respPop*sigma_respPop)
  sigma_triPop ~ dunif(0, 10)
  tau_triPop <- 1/(sigma_triPop*sigma_triPop)
  sigma_gndPop ~ dunif(0, 10)
  tau_gndPop <- 1/(sigma_gndPop*sigma_gndPop)
  sigma_H10Pop ~ dunif(0, 10)
  tau_H10Pop <- 1/(sigma_H10Pop*sigma_H10Pop)
  sigma_H36Pop ~ dunif(0, 10)
  tau_H36Pop <- 1/(sigma_H36Pop*sigma_H36Pop)
  
  # Priors for effects of leaf size.
  slope_resp ~ dnorm(0, 0.001)
  slope_tri ~ dnorm(0, 0.001)
  slope_gnd ~ dnorm(0, 0.001)
  slope_H10 ~ dnorm(0, 0.001)
  slope_H36 ~ dnorm(0, 0.001)

  # Priors for effects of herbivore density on defense traits.
  slope_triPop ~ dnorm(0, 0.001)
  slope_gndPop ~ dnorm(0, 0.001)
  slope_H10Pop ~ dnorm(0, 0.001)
  slope_H36Pop ~ dnorm(0, 0.001)
}
", fill = TRUE)
sink()
