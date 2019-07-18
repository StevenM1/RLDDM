### Model 2
# eta1 ~ Cue
# 

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('beta' = 1,
                                     'eta1' = 'condition'),
                 'constants'=c('eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='Qlearning')

obj <- objSoftmaxMultiCond

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
transformChoicePars <- function(pars, condition) {
  ### Gets trial-by-trial beta parameters ###
  nTrials = length(condition)
  beta <- rep(pars[['beta']], nTrials)
#  beta[condition=='SPD'] <- pars[['beta.SPD']]
#  beta[condition=='ACC'] <- pars[['beta.ACC']]
  return(list(beta=beta))
}
