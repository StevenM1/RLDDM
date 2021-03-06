### Model 1qbnl
# a ~ Cue, Q-learning, 2 learning rates, v ~ nonlin(ev)

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('a' = 'condition',
                                     'm' = 1,
                                     't0' = 1,
                                     'eta1' = 1,
                                     'eta2' = 1,
                                     'vmax' = 1),
                 'constants'=c('sv'=0, 'sz'=0, 'z'=0.5, 's'=1),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='Qlearning')

obj <- objRLDDMMultiCond

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition) {
  eta1 <- rep(pars[['eta1']], length(condition))
  eta2 <- rep(pars[['eta2']], length(condition))
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
  v = (2*pars[['vmax']])/(1+exp(-delta_ev*pars[['m']]))-pars[['vmax']]
  
  # a differs by condition
  a[condition=='SPD'] <- pars[['a.SPD']]
  a[condition=='ACC'] <- pars[['a.ACC']]
  
  # rescale z from [0, 1] to [0, a]
  z = z*a
  return(list(t0=t0, a=a, v=v, z=z, sz=sz, sv=sv, s=s))
}

defaultBounds <- function() {
  list('a'=c(1e-3, 10),
       'v'=c(1e-3, 10),
       'm'=c(1e-3, 1000),
       'z'=c(.2, .8),
       't0'=c(0, .5),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'sz'=c(0, 5),
       'sv'=c(0, 5),
       's'=c(0, 10),
       'vmax'=c(0, 100),
       'k'=c(-10, 10))
}