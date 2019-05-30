### Model 2
# eta ~ cue

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('a' = 1,
                                     'm' = 1,
                                     't0' = 1,
                                     'eta1' = 'condition'),
                 'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1, 'eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'))

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition) {
  # "Declare"
  eta1 <- eta2 <- rep(NA, length(condition))
  eta1[condition=='SPD'] = pars[['eta1.SPD']]
  eta1[condition=='ACC'] = pars[['eta1.ACC']]
  eta2 <- eta1
  
  return(list(eta1=eta1, eta2=eta2))
}

### the following function gets trial-by-trial DDM pars
transformDDMPars <- function(pars, condition, delta_ev) {
  ### Gets trial-by-trial DDM parameters ###
  
  nTrials = length(condition)
  a <- v <- t0 <- z <- sv <- sz <- s <- rep(NA, nTrials)
  
  # all current models have no variability in sz, sv, s, t0
  t0 = rep(pars[['t0']], nTrials)
  z = rep(pars[['z']], nTrials)
  sv <- rep(pars[['sv']], nTrials)
  sz <- rep(pars[['sz']], nTrials)
  s <- rep(pars[['s']], nTrials)
  
  # all models assume a linear relation between delta_ev and v
  v = delta_ev*pars[['m']]
  
  # Model 2: constant threshold across cues
  a = rep(pars[['a']], nTrials)
  
  return(list(t0=t0, a=a, v=v, z=z, sz=sz, sv=sv, s=s))
}
