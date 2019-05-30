rm(list=ls())  # fresh start

# Libraries, directories ---------------------------------------------------------------
.libPaths(c(.libPaths(), '/home/stevenm/rpackages')) # for tux server
library(rtdists)
library(DEoptim)
library(RLDDM)
fitModel <- '~/surfdrive/data/learningTask/Barbara_preprocessed'
workDir <- '/Users/steven/Sync/PhDprojects/RLDDM'

# Model to fit -------------------------------------------------------------------
modelType = 'RLDDM'
source(file.path(workDir, 'analysis/fittingFunctions.R'))
exp <- 'exp3'
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
dir.create(resDir, showWarnings = FALSE, recursive=TRUE)

# Optimization options ---------------------------------------------------------
nCores <- 2
# Load data
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
# very minor changes to data
dat$cue <- as.character(dat$cue)
dat[is.na(dat$cue), 'cue'] <- 'NEU'
dat$cue <- as.factor(dat$cue)
ppsToFit <- unique(dat$pp)
modelsToFit <- 1:4


# Collect all parameters --------------------------------------------------
## Get all parameters, mean RT and mean accuracy into single dataframe
fit <- data.frame(expand.grid(pp=ppsToFit, modelN=modelsToFit, cue=c('SPD', 'ACC')),
                  eta1=NA, eta2=NA, a=NA, t0=NA, m=NA, mrt=NA, accObjective=NA, accSubjective=NA, 
                  negLL=NA, nPars=NA, nObs=NA, BIC=NA, minBIC=NA, wBIC=NA, AIC=NA, minAIC=NA, wAIC=NA)

for(modelN in 1:4) {
  source(file.path(workDir, 'analysis', 'models', paste0('model', modelN, '.R')))
  modelSetup <- model(modelSpec)
  
  for(pp in ppsToFit) {
    # load data
    d = prepareForFitting(dat[dat$pp==pp&as.character(dat$Block)=='Miniblocks',]) 
    df=d$df
    outcomes=d$outcomes
    VVchoiceIdx=d$VVchoiceIdx
    values=d$values
    choice=d$choice

    # load model
    modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_Miniblocks.rdat'))
    load(modelFn)
    names(bestPars) <- modelSetup$p.vector
    fitModel <- objRLDDMMultiCond(bestPars,
                               rt=df$rt,
                               parNames=modelSetup$p.vector,
                               constants=modelSetup$constants,
                               choice=choice,
                               condition=df$cue,
                               values=values,
                               outcomes=outcomes,
                               VVchoiceIdx=VVchoiceIdx,
                               returnType='full')
    df$EVchoice <- fitModel$VV[,1]-fitModel$VV[,2]
    df$modelChoice <- ifelse(df$EVchoice > 0, 2, 1)
    # "subjective" accuracy: choice in line with highest EV-stimulus
    df$accuracySubjective <- df$modelChoice == df$choice
    # "objective" accuracy: Choice is optimal
    df$accuracyObjective <- df$choice == 2
    df$PE <- apply(fitModel$PE, 1, sum, na.rm=TRUE)
    
    # some metrics
    mrt <- aggregate(rt ~ cue, df, mean)
    accObjective <- aggregate(accuracyObjective ~ cue, df, mean)
    accSubjective <- aggregate(accuracySubjective ~ cue, df, mean)
    idx1 <- fit$pp==pp & fit$modelN==modelN
    fit[idx1, 'nObs'] <- nrow(df)
    for(cue in mrt$cue) {
      fit[idx1&fit$cue==cue, 'mrt'] <- mrt[mrt$cue==cue, 'rt']
      fit[idx1&fit$cue==cue, 'accObjective'] <- accObjective[accObjective$cue==cue, 'accuracyObjective']
      fit[idx1&fit$cue==cue, 'accSubjective'] <- accSubjective[accSubjective$cue==cue, 'accuracySubjective']
    }
    
    fit[idx1, 'negLL'] <- -fitModel$LL
    names(bestPars) <- modelSetup$p.vector
    fit[idx1, 't0'] <- bestPars[['t0']]
    fit[idx1, 'm'] <- bestPars[['m']]
    fit[idx1, 'nPars'] <- length(bestPars)
    if(modelN == 1) {
      fit[idx1, 'eta1'] <- bestPars[['eta1']]
      fit[idx1&fit$cue=='SPD', 'a'] <- bestPars[['a.SPD']]
      fit[idx1&fit$cue=='ACC', 'a'] <- bestPars[['a.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'a'] <- bestPars[['a.NEU']]
    } else if(modelN == 2) {
      fit[idx1, 'a'] <- bestPars[['a']]
      fit[idx1&fit$cue=='SPD', 'eta1'] <- bestPars[['eta1.SPD']]
      fit[idx1&fit$cue=='ACC', 'eta1'] <- bestPars[['eta1.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'eta1'] <- bestPars[['eta1.NEU']]
    } else if(modelN == 3) {
      fit[idx1, 'a'] <- bestPars[['a']]
      fit[idx1, 'eta1'] <- bestPars[['eta1']]
    } else if(modelN == 4) {
      fit[idx1&fit$cue=='SPD', 'a'] <- bestPars[['a.SPD']]
      fit[idx1&fit$cue=='ACC', 'a'] <- bestPars[['a.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'a'] <- bestPars[['a.NEU']]
      fit[idx1&fit$cue=='SPD', 'eta1'] <- bestPars[['eta1.SPD']]
      fit[idx1&fit$cue=='ACC', 'eta1'] <- bestPars[['eta1.ACC']]
      #      fit[idx1&fit$cue=='NEU', 'eta1'] <- bestPars[['eta1.NEU']]
    }
    fit[idx1, 'BIC'] <- 2*fit[idx1, 'negLL'] + log(fit[idx1, 'nObs'])*fit[idx1, 'nPars']
    fit[idx1, 'AIC'] <- 2*fit[idx1, 'negLL'] + 2*fit[idx1, 'nPars']
    # fit[idx1, 'minBIC'] <- fit[idx1, 'BIC']-min(fit[idx1, 'BIC'])
  }
}





###
### only if fit on single block
#fit[fit$modelN %in% c(1,2), 'nPars'] <- fit[fit$modelN %in% c(1,2), 'nPars']-1
#fit[fit$modelN==4, 'nPars'] <- fit[fit$modelN==4, 'nPars']-2

fitWide = reshape(data=fit,#[fit$modelN!=3,],
                  v.names=c('eta1', 'eta2', 'a', 't0', 'm', 'mrt', 'accObjective', 'accSubjective'), 
                  direction='wide', idvar = c('pp', 'modelN'), timevar=c('cue'))
for(pp in fitWide$pp) {
  idx <- fitWide$pp==pp
  fitWide[idx, 'minBIC'] = fitWide[idx, 'BIC'] - min(fitWide[idx, 'BIC'])
  fitWide[idx, 'minAIC'] = fitWide[idx, 'AIC'] - min(fitWide[idx, 'AIC'])
}

fitWide[fitWide$minBIC==0, 'modelN']
fitWide[fitWide$minAIC==0, 'modelN']

winningModels <- data.frame(pp=unique(fitWide$pp), BIC=NA, AIC=NA, nObs=NA)
for(pp in fitWide$pp) {
  winningModels[winningModels$pp==pp, 'BIC'] = fitWide[fitWide$pp==pp & fitWide$minBIC==0, 'modelN']
  winningModels[winningModels$pp==pp, 'AIC'] = fitWide[fitWide$pp==pp & fitWide$minAIC==0, 'modelN']
  winningModels[winningModels$pp==pp, 'nObs'] = fitWide[fitWide$pp==pp & fitWide$minBIC==0, 'nObs']
}
table(winningModels$BIC)
table(winningModels$BIC)/sum(table(winningModels$BIC))


#table(winningModels$AIC)



####
fitWide$accObjDiff <- fitWide$accObjective.ACC-fitWide$accObjective.SPD
fitWide$accSubjDiff <- fitWide$accSubjective.ACC-fitWide$accSubjective.SPD
fitWide$mrtDiff <- fitWide$mrt.ACC-fitWide$mrt.SPD

mean(fitWide[fitWide$modelN==1, 'accObjDiff'])
mean(fitWide[fitWide$modelN==1, 'accSubjDiff'])