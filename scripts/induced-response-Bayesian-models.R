## USDA cotton experiment
## Data analysis: Resistance and defenses ~ herbivore density Bayesian models
## Code written and annotated by JM (last updated 2024-10-18)

# The code in this script file:
# - Specifies models as described in Eq. 1, for normally distributed or binomially
#   distributed response variable;
# - Saves models as .txt files.

## Normally-distributed response variable ======================================

sink("scripts/TimexDens-LinThr-Geno.txt")
cat("
model {

  for (i in 1:nGeno) {
    int[i] ~ dnorm(0, tau_int)
  }
  
  for (i in 1:nTimept) {
    intPop[i] ~ dnorm(0, 0.001)
    slope_densPop[i] ~ dnorm(0, 0.001)
  }
  
  for (i in 1:nPop) {
    respPop[i] ~ dnorm(mu_respPop[i], tau_respPop)
    mu_respPop[i] <- intPop[timeptPop[i]] + slope_densPop[timeptPop[i]]*densPop[i]
  }
  
  for (i in 1:nObs) {
    resp[i] ~ dnorm(mu_resp[i], tau_resp)
    mu_resp[i] <- respPop[pop[i]] + slope_resp*leaf_size[i] + int[geno[i]]
  }
  
  sigma_int ~ dunif(0, 10)
  tau_int <- 1/(sigma_int*sigma_int)
  
  slope_resp ~ dnorm(0, 0.001)
  
  sigma_respPop ~ dunif(0, 10)
  tau_respPop <- 1/(sigma_respPop*sigma_respPop)
  
  sigma_resp ~ dunif(0, 10)
  tau_resp <- 1/(sigma_resp*sigma_resp)
  
  for (i in 1:nTimept) {
  for (j in 1:nDens) {
    pred[i,j] <- slope_densPop[i]*densPop_newdat[j]
  }}

densT1T2 <- slope_densPop[1] - slope_densPop[2]
densT1T3 <- slope_densPop[1] - slope_densPop[3]
densT1T4 <- slope_densPop[1] - slope_densPop[4]
densT1T5 <- slope_densPop[1] - slope_densPop[5]
densT2T3 <- slope_densPop[2] - slope_densPop[3]
densT2T4 <- slope_densPop[2] - slope_densPop[4]
densT2T5 <- slope_densPop[2] - slope_densPop[5]
densT3T4 <- slope_densPop[3] - slope_densPop[4]
densT3T5 <- slope_densPop[3] - slope_densPop[5]
densT4T5 <- slope_densPop[4] - slope_densPop[5]
  
}
", fill = TRUE)
sink()


## Binomially-distributed response variable (logit link) =======================

sink("scripts/TimexDens-LinThr-Geno-bin.txt")
cat("
model {

  for (i in 1:nGeno) {
    int[i] ~ dnorm(0, tau_int)
  }
  
  for (i in 1:nTimept) {
    intPop[i] ~ dnorm(0, 0.001)
    slope_densPop[i] ~ dnorm(0, 0.001)
  }
  
  for (i in 1:nPop) {
    respPop[i] ~ dnorm(mu_respPop[i], tau_respPop)
    mu_respPop[i] <- intPop[timeptPop[i]] + slope_densPop[timeptPop[i]]*densPop[i]
  }
  
  for (i in 1:nObs) {
    resp[i] ~ dbin(prob_resp[i], num_trials[i])
    logit(prob_resp[i]) <- respPop[pop[i]] + slope_resp*leaf_size[i] + int[geno[i]]
  }
  
  sigma_int ~ dunif(0, 10)
  tau_int <- 1/(sigma_int*sigma_int)
  
  slope_resp ~ dnorm(0, 0.001)
  
  sigma_respPop ~ dunif(0, 10)
  tau_respPop <- 1/(sigma_respPop*sigma_respPop)
  
  for (i in 1:nTimept) {
  for (j in 1:nDens) {
    pred[i,j] <- slope_densPop[i]*densPop_newdat[j]
  }}

densT1T2 <- slope_densPop[1] - slope_densPop[2]
densT1T3 <- slope_densPop[1] - slope_densPop[3]
densT1T4 <- slope_densPop[1] - slope_densPop[4]
densT1T5 <- slope_densPop[1] - slope_densPop[5]
densT2T3 <- slope_densPop[2] - slope_densPop[3]
densT2T4 <- slope_densPop[2] - slope_densPop[4]
densT2T5 <- slope_densPop[2] - slope_densPop[5]
densT3T4 <- slope_densPop[3] - slope_densPop[4]
densT3T5 <- slope_densPop[3] - slope_densPop[5]
densT4T5 <- slope_densPop[4] - slope_densPop[5]
  
}
", fill = TRUE)
sink()
