# Compare parameters
rm(list=ls())  # fresh start

# Libraries, directories ---------------------------------------------------------------
source('./setup_locale.R') # gets dataDir, workDir, dependent on server
library(rtdists)
library(DEoptim)
library(RLDDM)
library(plyr)

# Model to fit -------------------------------------------------------------------
modelType = 'RLRW'
source(file.path(workDir, 'analysis/src/fittingFunctions.R'))
source(file.path(workDir, 'analysis/src/racingWaldDists.R'))
exp <- 'exp2'
resDir <- file.path(workDir, 'fits', 'Barbara', exp, paste0('model-', modelType))
dir.create(resDir, showWarnings = FALSE, recursive=TRUE)

# Load data
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
block <- 'Miniblocks'
ppsToFit <- unique(dat$pp)

# Create dataframe with all parameters
pardf <- NULL
for(block in c('Calibration', 'Miniblocks', 'Trialwise')) {
  if(block == 'Calibration') {
    modelN <- 'racingWaldAccAdvantages'
  } else {
    modelN <- '2racingWaldAccAdvantages'
  }
  source(file.path(workDir, 'analysis', 'models', paste0('model', modelN, '.R')))
  modelSetup <- model(modelSpec)
  pardfc <- NULL
  
  for(pp in ppsToFit) {
    modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat'))
    load(modelFn)
    names(bestPars) <- modelSetup$p.vector
    pardfc <- rbind(pardfc, c('sub'=pp, 'model'=modelN, 'block'=block, bestPars))
  }
  pardfc <- data.frame(pardfc)
  if(!is.null(pardf)) {
    pardf <- rbind.fill(pardf, pardfc)
  } else {
    pardf <- pardfc
  }
}

# cast to right types
for(col in colnames(pardf)) {
  if(!col %in% c('sub', 'model', 'block')) {
    pardf[,col] <- as.numeric(as.character(pardf[,col]))
  }
}


# Plot --------------------------------------------------------------------
par(mfrow=c(3,2))
#includeRowIdx <- rep(TRUE, nrow(pardf)) #
includeRowIdx <- pardf$eta1 > 0.001 & pardf$b1 < 20 & pardf$eta1 < .5
tmp <- aggregate(includeRowIdx ~ sub, pardf, sum)
includeSubs <- tmp$sub[tmp$includeRowIdx==3]
pardfIncl <- pardf[pardf$sub %in% includeSubs,]
xName = 'Miniblocks'
yName = 'Trialwise'
idx1 <- pardfIncl$block == xName
idx2 <- pardfIncl$block == yName
parnames <- c('b0.SPD', 'b0.ACC', 'b1', 'eta1', 't0', 'B')#, 'bU')
for(parname in parnames) {
  x = pardfIncl[idx1, parname]
  if(parname %in% c('b0.SPD', 'b0.ACC') & (yName == 'Calibration')) {
    y = pardfIncl[idx2, 'b0']
  } else {
    y = pardfIncl[idx2, parname]
  }
  
  if(parname %in% c('b1')) {
    xlim <- ylim <- c(0, 15)
  } else {
    xlim <- ylim <- NULL
  }
  plot(x, y, main=parname, ylab=yName, xlab=xName, xlim=xlim, ylim=ylim)
  abline(reg=lm(y~x), lty=2)
  abline(a=0, b=1, col='grey')
  loc <- ifelse(parname %in% c('eta1', 'b1'), 'right', 'bottomright')
  legend(loc, c(paste0('r=', round(cor(x, y), 3))))
}


###### Between-subject correlations?
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(pardf[pardf$model=='racingWaldAccAdvantages', c('B', 't0', 'eta1', 'b0', 'b1')], panel = panel.smooth,
      cex = 1.5, horOdd=TRUE,
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)


###### Approximation of within-subject correlations
# eta1 ~ b1
pp <- 32
modelN <- 'racingWaldAccAdvantages'
block <- 'Calibration'
load(file.path(dataDir, paste0('data_', exp, '.Rdata')))
dat$cue <- as.character(dat$cue)
dat[is.na(dat$cue), 'cue'] <- 'NEU'
dat$cue <- as.factor(dat$cue)
thisPpDat <- dat[dat$pp==pp&as.character(dat$Block)==block,]
d = prepareForFitting(thisPpDat)
df=d$df
outcomes=d$outcomes
VVchoiceIdx=d$VVchoiceIdx
values=d$values
choice=d$choice
choice = ifelse(choice==2, 1, 2)

fitPars <- pardf[pardf$sub==pp&pardf$model==modelN, c('B', 't0', 'eta1', 'b0', 'b1')]
modelFn = file.path(resDir, paste0('sub-', pp, '_model-', modelN, '_', block, '.rdat'))
load(modelFn)
names(bestPars) <- modelSetup$p.vector

eta1s <- seq(0, 1, .01)
b1s <- seq(0, 100, .25)
LLdf <- expand.grid(eta1=eta1s, b1=b1s, LL=NA)
for(i in 1:nrow(LLdf)) {
  thisPars <- fitPars
  thisPars[['b1']] <- LLdf[i,'b1']
  thisPars[['eta1']] <- LLdf[i,'eta1']
  LLdf[i,'LL'] <- obj(thisPars, 
                  rt=df$rt,
                  parNames=modelSetup$p.vector,
                  constants=modelSetup$constants,
                  choice=choice,
                  condition=df$cue,
                  values=values,
                  outcomes=outcomes, 
                  VVchoiceIdx=VVchoiceIdx,
                  returnType='negLL',
                  learningRule = modelSetup$learningRule)
}
#LLdf$LL[LLdf$LL > 100] <- 100
LLdf$posLL <- -LLdf$LL  # flip sign
LLdf$posLL[LLdf$posLL < -100] <- -100
bestLL <- max(LLdf$posLL)
breaks <- bestLL - c(100, 50, 10, 5, 1)

ggplot(LLdf, aes(x = b1, y = eta1, z = posLL, fill=posLL)) +
  geom_tile() +
#  geom_contour(breaks=breaks, colour=c('red')) +
#  scale_fill_gradient2(low="green", mid="lightblue", high="red",  # colors in the scale
#                        midpoint=90,                              # same midpoint for plots (mean of the range)
#                        breaks=seq(80,110,5),                     # breaks in the scale bar
#                        limits=c(90, 110), 
#                        guide = "colorbar") + 
  annotate("point", x=bestPars[['b1']], y=bestPars[['eta1']],
          alpha = 1, color='purple')

