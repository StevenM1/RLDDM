### Model LBA1
# a depends on cue
# mean_v_1 ~ 1
# mean_v_2 ~ 1


### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('B' = 'condition',
                                     'A' = 1,
                                     'mean_v1' = 1,
                                     'mean_v2' = 1,
                                     'sd_v2' = 1,
                                     't0' = 1),
                 'constants'=c('sd_v1'=1, 'eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'))

obj <- objLBAMultiCond

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition, ...) {
  # "Declare"
  eta1 <- eta2 <- rep(pars[['eta1']], length(condition))
  eta2 <- eta1
  return(list(eta1=eta1, eta2=eta2))
}

### the following function gets trial-by-trial LBA pars
transformLBAPars <- function(pars, condition, ...) {
  ### Gets trial-by-trial DDM parameters ###
  
  nTrials = length(condition)
  t0 <- A <- b <- sd_v2 <- rep(NA, nTrials)

  # no trial-by-trial variability in these params
  t0 = rep(pars[['t0']], nTrials)
  A = rep(pars[['A']], nTrials)
  mean_v1 = rep(pars[['mean_v1']], nTrials)
  mean_v2 = rep(pars[['mean_v2']], nTrials)
  sd_v1 = rep(pars[['sd_v1']], nTrials)
  sd_v2 = rep(pars[['sd_v2']], nTrials)
  
  # a differs
  b[condition=='SPD'] <- pars[['A']] + pars[['B.SPD']]
  b[condition=='ACC'] <- pars[['A']] + pars[['B.ACC']]
  
  # to matrix
#  mean_v = matrix((mean_v1, mean_v2), ncol=2)
#  sd_v = matrix(c(sd_v1, sd_v2), ncol=2)
  
  return(list(t0=t0, A=A, mean_v=list(mean_v1, mean_v2), sd_v=list(sd_v1, sd_v2), b=b))
}

defaultBounds <- function() {
  list('A'=c(1e-3, 10),
       'B'=c(1e-3, 10),
       'mean_v1'=c(-3, 10),
       'mean_v2'=c(-3, 10),
       'sd_v1'=c(.1, 5),
       'sd_v2'=c(.1, 5),
       'm'=c(1e-3, 100),
       't0'=c(0, .6),
       'eta1'=c(0, 1),
       'eta2'=c(0, 1),
       'b'=c(0, 5),
       'sv'=c(0, 5),
       's'=c(0, 10),
       'vmax'=c(0, 100),
       'k'=c(-10, 10))
}