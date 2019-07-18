### Model 
# 1 ~ condition
# Softmax (SM) as choice model
# Q-learning as update rule

### modelSpec is a list containing:
# 1. The parameters to fit, and the factors they depend on
# 2. constants in the model
# 3. The factors from (1), and their levels
modelSpec = list('variablePars'=list('beta' = 1,
                                     'eta1' = 1),
                 'constants'=c('eta2'=-Inf),
                 'condition'=c('SPD', 'ACC'),
                 'learningRule'='Qlearning',
                 'choiceFunction'='SM')

### transformLearningRate is a function transforming
### "global" parameters to trial-by-trial values, dependent 
### on the condition
transformLearningRate <- function(pars, condition) {
  # "Declare"
  eta1 <- eta2 <- rep(pars[['eta1']], length(condition))
  return(list(eta1=eta1, eta2=eta2))
}

### the following function gets trial-by-trial choice pars
transformChoicePars <- function(pars, condition, ...) {
  ### Gets trial-by-trial beta parameters ###
  nTrials = length(condition)
  beta <- rep(pars[['beta']], nTrials)
  return(list(beta=beta))
}
