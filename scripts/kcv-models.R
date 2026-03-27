## =============================================================================
## Models for KCV ~ herbivore density
## =============================================================================

# ## For normally-distributed response variable ==================================
# # (note that 'pred' parameter to estimate is prediction with intercept, i.e.,
# # not the strength of induction)
# sink("scripts/TimexDens-KCV.txt")
# cat("
# model {
#   
#   for (i in 1:nTimept) {
#     intPop[i] ~ dnorm(0, 0.001)
#     slope_densPop[i] ~ dnorm(0, 0.001)
#   }
#   
#   for (i in 1:nPop) {
#     respPop[i] ~ dnorm(mu_respPop[i], tau_respPop)
#     mu_respPop[i] <- intPop[timeptPop[i]] + slope_densPop[timeptPop[i]]*densPop[i]
#   }
#   
#   sigma_respPop ~ dunif(0, 10)
#   tau_respPop <- 1/(sigma_respPop*sigma_respPop)
#   
#   for (i in 1:nTimept) {
#   for (j in 1:nDens) {
#     pred[i,j] <- intPop[i] + slope_densPop[i]*densPop_newdat[j]
#   }}
#   
# }
# ", fill = TRUE)
# sink()
# 
# ## For beta-distributed response variable ======================================
# # (note that 'pred' parameter to estimate is prediction with intercept, i.e.,
# # not the strength of induction)
# sink("scripts/TimexDens-KCV-beta.txt")
# cat("
# model {
#   
#   for (i in 1:nTimept) {
#     intPop[i] ~ dnorm(0, 0.001)
#     slope_densPop[i] ~ dnorm(0, 0.001)
#   }
#   
#   for (i in 1:nPop) {
#     respPop[i] ~ dbeta(alpha[i], beta[i])
#     alpha[i] <- mu_respPop[i] * phi
#     beta[i] <- ( 1 - mu_respPop[i] ) * phi
#     logit(mu_respPop[i]) <- intPop[timeptPop[i]] + slope_densPop[timeptPop[i]]*densPop[i]
#   }
#   
#   for (i in 1:nTimept) {
#   for (j in 1:nDens) {
#     logit(pred[i,j]) <- intPop[i] + slope_densPop[i]*densPop_newdat[j]
#   }}
# 
# phi ~ dgamma(0.1,0.1)
#   
# }
# ", fill = TRUE)
# sink()

## For beta-distributed response variable, without effect of timepoint =========
# (note that 'pred' parameter to estimate is prediction with intercept, i.e.,
# not the strength of induction)
sink("scripts/NoTimexDens-KCV-beta.txt")
cat("
model {
  
  intPop ~ dnorm(0, 0.001)
  slope_densPop ~ dnorm(0, 0.001)
  
  for (i in 1:nPop) {
    respPop[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu_respPop[i] * phi
    beta[i] <- ( 1 - mu_respPop[i] ) * phi
    logit(mu_respPop[i]) <- intPop + slope_densPop*densPop[i]
  }


  for (j in 1:nDens) {
    logit(pred[j]) <- intPop + slope_densPop*densPop_newdat[j]
  }

phi ~ dgamma(0.1,0.1)
  
}
", fill = TRUE)
sink()

