## USDA cotton experiment
## Data analysis: Helper functions for Bayesian data analysis
## Code written and annotated by JM (last updated 2024-10-21)

require(bayestestR)

## Function for assigning data vectors, given dataframe 'df' subsetted from 
##  'full_data' and 'pop_summary' dataframe that summarizes population-level
##  experimental treatments (i.e., sampling timepoint and herbivore density).
data_vecs_resp <- function(df, pop_summary = pop_guide) {
  # Assign data vectors.
  nGeno <- length(unique(df$Geno)) #number of genotypes
  nTimept <- length(unique(df$Timepoint)) #number of timepoints
  nPop <- length(unique(df$Pop_ID)) #number of populations
  nObs <- dim(df)[1] #number of observations
  
  pop <- df$Pop_ordered #population ID associated with each observation
  geno <- df$Geno_ordered #genotype associated with each observation
  
  # Use 'pop_guide' to get herbivore density and timepoint associated with each
  #  population in the subsetted dataframe ('pop_guide' is already ordered by Pop_ID).
  densPop <- pop_summary[which(pop_summary$Pop_ID %in% df$Pop_ID), 'Herb_add']
  timeptPop <- pop_summary[which(pop_summary$Pop_ID %in% df$Pop_ID), 'Timepoint']  
  
  # Return compiled list of data vectors.
  return(list( nGeno = nGeno, nTimept = nTimept, nPop = nPop, nObs = nObs,
               pop = pop, geno = geno, 
               densPop = densPop, timeptPop = timeptPop,
               nDens = 3, densPop_newdat = c(6, 12, 25) ))
}


## Function for saving output
save_out <- function ( model, dataframe, pop_summary = pop_guide, pop_means ) {
  ## Samples from posterior distribution.
  # Add chain ID column to each list element in samples
  samps <- mapply(cbind, model$samples, "Chain_ID" = seq(1,nc), SIMPLIFY = F)
  # Compile samples from all chains into single dataframe
  samps <- do.call(rbind, samps)
  
  ## Summary statistics for monitored variables, including median, 89% CRI, 95% CRI
  post <- sapply(colnames(model$samples[[1]]), 
                 function(x) { bayestestR::describe_posterior(unlist(model$samples[,x]), 
                                                              ci = c(0.89, 0.95)) })
  
  post_sum <- data.frame(median = do.call(rbind, post['Median',])[,1],
                         CI_low_89 = do.call(rbind, post['CI_low',])[,1],
                         CI_high_89 = do.call(rbind, post['CI_high',])[,1],
                         CI_low_95 = do.call(rbind, post['CI_low',])[,2],
                         CI_high_95 = do.call(rbind, post['CI_high',])[,2])
  
  post_sum <- cbind(summary(model), post_sum)
  
  if (pop_means == TRUE) {
    ## Estimates of population means.
    # Get the estimates of population means from the model summary (any rows containing
    #  'respPop')
    respPop <- post_sum[which(sapply(dimnames(post_sum)[1], 
                                     FUN = function(x){grepl("respPop", x)})),]
    # Make dataframe with actual Pop_ID (from mod_df) and associated estimates
    pop <- data.frame(cbind(
      Pop_ordered = unique(dataframe$Pop_ordered),
      Pop_actual = unique(dataframe$Pop_ID),
      respPop[unique(dataframe$Pop_ordered),]
    ))
    # Add information about herbivore density and timestep associated with each 
    #  population
    pop <- cbind(pop, pop_guide[match(pop$Pop_actual, pop_guide$Pop_ID),])
  }
  
  if (pop_means == FALSE) { pop <- NA }
  
  return(list(samples = samps, summary = post_sum, pop_means = pop))
}
