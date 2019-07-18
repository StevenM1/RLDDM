### Model 
# a ~ condition
# DDM as choice model
# SARSA as update rule
# st0 estimated

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('a' = 1,
                                     'm' = 'condition',
                                     't0' = 1,
                                     'eta1' = 1,
                                     'st0' = 1),
                 'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='SARSA',
                 'choiceFunction'='DDM')

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition) {
  # "Declare"
  eta1 <- eta2 <- rep(pars[['eta1']], length(condition))
  return(list(eta1=eta1, eta2=eta2))
}

### the following function gets trial-by-trial DDM pars
transformChoicePars <- function(pars, condition, ev) {
  ### Gets trial-by-trial DDM parameters ###
  delta_ev = ev[,1]-ev[,2]
  nTrials = length(condition)
  a <- v <- t0 <- z <- sv <- sz <- s <- rep(NA, nTrials)
  
  # all current models have no variability in sz, sv, s, t0
  t0 = rep(pars[['t0']], nTrials)
  z = rep(pars[['z']], nTrials)
  sv <- rep(pars[['sv']], nTrials)
  sz <- rep(pars[['sz']], nTrials)
  s <- rep(pars[['s']], nTrials)
  a <- rep(pars[['a']], nTrials)
  st0 <- rep(pars[['st0']], nTrials)
  
  # all models assume a linear relation between delta_ev and v
  scalingFactor <- rep(NA, nTrials)
  scalingFactor[condition == 'SPD'] <- pars[['m.SPD']]
  scalingFactor[condition == 'ACC'] <- pars[['m.ACC']]
  v = delta_ev*scalingFactor

  # rescale z from [0, 1] to [0, a]
  z = z*a
  return(list(t0=t0, a=a, v=v, z=z, sz=sz, sv=sv, s=s, st0=st0))
}

defaultBounds <- function() {
  list('a'=c(1e-3, 5),
       'v'=c(1e-3, 5),
       'm'=c(1e-3, 50),
       'z'=c(.2, .8),
       't0'=c(0, .5),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'sz'=c(0, 5),
       'sv'=c(0, 2),
       's'=c(0, 5),
       'vmax'=c(0, 100),
       'k'=c(0, 50),
       'st0'=c(0, .5),
       'beta'=c(0, 50))
}
