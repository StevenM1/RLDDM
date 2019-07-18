### Model 1
# a ~ Cue
# 

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('B' = 'condition',
                                     'm' = 1,
#                                     'k' = 1,
                                     't0' = 1,
                                     'eta1' = 1),
                 'constants'=c('st0'=0, 'A'=0, 'eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='SARSA')

obj <- objRLRWMultiCond

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
  nTrials = length(condition)
  B <- t0 <- A <- m <- rep(NA, nTrials)
  
  # all current models have no variability in sz, sv, s, t0
  A = rep(pars[['A']], nTrials)
  t0 = rep(pars[['t0']], nTrials)
  st0 <- rep(pars[['st0']], nTrials)
  
  # all models assume a linear relation between delta_ev and v
  v1 = exp(ev[,1])
  v2 = exp(ev[,2])
  v1 = v1/(v1+v2)*pars[['m']]
  v2 = v2/(v1+v2)*pars[['m']]
  
  # B differs by condition
  B[condition=='SPD'] <- pars[['B.SPD']]
  B[condition=='ACC'] <- pars[['B.ACC']]
  
  # rescale z from [0, 1] to [0, a]
  return(list(t0=t0, B=B, A=A, v1=v1, v2=v2, st0=st0))
}


defaultBounds <- function() {
  list('a'=c(1e-3, 5),
       'A'=c(0, 5),
       'B'=c(1e-3, 5),
       'v1'=c(1e-3, 5),
       'v2'=c(1e-3, 5),
       'v'=c(1e-3, 5),
       'm'=c(-100, 100),
       'z'=c(.2, .8),
       't0'=c(0, .5),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'sz'=c(0, 5),
       'sv'=c(0, 5),
       's'=c(0, 5),
       'vmax'=c(0, 100),
       'k'=c(0, 250),
       'st0'=c(0, .5),
       'beta'=c(0, 50),
       'startValue'=c(0, 2))
}
