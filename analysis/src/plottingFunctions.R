plotCDF <- function(qps, add=FALSE, ...) {
  
  if(!'def.prob' %in% names(qps)) {
    def.props <- qps$prop*as.numeric(gsub("[^0-9\\.]", "", names(qps$qps)))/100
  } else {
    def.props <- qps$def.prob
  }
  if(add) {
    lines(qps$qps, def.props, ...)
  } else {
    plot(qps$qps, def.props, ...)
  }
}

getQuantiles <- function(responses, rts, goodResponse='upper', badResponse='lower', quants=seq(.1, .9, .2)) {
  # Get cum. def. quantiles for a dataset, as well as some other variables of interest
  #
  qps1 <- quantile(rts[responses==goodResponse], quants)
  qps2 <- quantile(rts[responses==badResponse], quants)
  df <- data.frame(rt=rts, responses=responses, trialN=1:length(rts))
  df$bin <- cut(df$trialN, 10, labels=FALSE)
  
  list('r2'=list(prop=mean(responses==goodResponse),
                 qps=qps1,
                 def.probability=quants*mean(responses==goodResponse)),
       'r1'=list(prop=mean(responses==badResponse),
                 qps=qps2,
                 def.probability=quants*mean(responses==badResponse)),
       'rtByBin'=aggregate(rt~bin, df, mean)[,2], 
       'accByBin'=aggregate((responses==goodResponse) ~ bin, df, mean)[,2],
       'skewByBin'=aggregate(rt~bin, df, skewness)[,2],
       'iqrByBin'=aggregate(rt~bin, df, IQR)[,2])
}

plotSingleSub <- function(nSim=50, byColumn=NULL) {
  ## for now, assume all model/parameter pars are in globals
  
  #  if(outDEoptim$optim$bestval > 1e3) { print(pp); dev.off(); next }
  #  modelSetup <- model(modelSpec)
  
  names(bestPars) <- modelSetup$p.vector
  bestparsNoConstants <- bestPars
  fitModel <- obj(bestPars, 
                  rt=df$rt,
                  parNames=modelSetup$p.vector,
                  constants=modelSetup$constants,
                  choice=choice,
                  condition=df$cue,
                  values=values,
                  outcomes=outcomes, 
                  VVchoiceIdx=VVchoiceIdx,
                  returnType='full',
                  learningRule = modelSetup$learningRule,
                  choiceFunction = modelSetup$choiceFunction)
  
  for(i in 1:length(modelSetup$constants)) bestPars[[names(modelSetup$constants[i])]] <- modelSetup$constants[[i]]
  choicePars = transformChoicePars(condition = df$cue, pars=bestPars, fitModel$VV)
  
#  cat('simulating', sep='')
  quants <- c(.1, .3, .5, .7, .9)
  
#  sfExportAll()
#  sims.qps <- sfLapply(1:nSim, function(x) {
  sims.qps <- lapply(1:nSim, function(x) {
#    cat('.', sep='')
    if(modelSetup$choiceFunction=='RW') {
      simdat <- rWaldRaceSM(n=nrow(df), 
                            A=choicePars[['A']], 
                            B=choicePars[['B']], 
                            v=list(choicePars[['v1']], choicePars[['v2']]), 
                            t0=choicePars[['t0']])
      getQuantiles(rts=simdat$RT, responses=simdat$R, goodResponse=1, badResponse=2)
    } else if(modelSetup$choiceFunction == 'DDM') {
      simdat <- rdiffusion(n=nrow(df), 
                           a=choicePars[['a']], 
                           v=choicePars[['v']], 
                           z=choicePars[['z']], 
                           s=choicePars[['s']], 
                           t0=choicePars[['t0']],
                           sv=choicePars[['sv']],
                           sz=choicePars[['sz']],
                           st0=choicePars[['st0']])
      getQuantiles(rts=simdat$rt, responses=simdat$response, goodResponse='upper', badResponse='lower')
    } else if(modelSetup$choiceFunction == 'LBA') {
      stop('LBA plotting not yet implemented')
    }
  })
  
  dat.qps <- getQuantiles(rts=df$rt, responses=df$choiceIsHighP, goodResponse=1, badResponse=0)
  
  # determine x/y limit
  xlim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$qps)))),
                  max(sapply(dat.qps, function(x) sapply(dat.qps, function(x) if('qps' %in% names(x)) {max(x$qps)} else {0}))))
  ylim.max <- max(max(sapply(sims.qps, function(x) sapply(c('r1', 'r2'), function(y) max(x[[y]]$prop)))),
                  max(sapply(dat.qps, function(x) sapply(dat.qps, function(x) if('qps' %in% names(x)) {max(x$qps)} else {0}))))
  
  if(is.na(xlim.max)) {
    warning('xlim.max could not be determined')
    xlim.max = 1
  }
  
  # Start plot
  plot(0, 0, type='n', xlim=c(0, xlim.max), ylim=c(0, ylim.max), xlab='RT', ylab='Def. prob.')
  abline(v=seq(0, 10, .25), lty='dotted', col='grey')
  abline(h=seq(0, 1, .1), lty='dotted', col='grey')
  abline(v=choicePars[['t0']][1], lty=2, col='black')
  abline(v=choicePars[['t0']][1]+choicePars[['st0']][1], lty=2, col='black')
  
  # Plot individual simulations (grey lines)
  tmp = lapply(sims.qps, function(x) lapply(1:2, function(y) plotCDF(qps=x[[y]], add=TRUE, type='l', lty=1, col='grey')))
  
  # mean over simulations
  sims.mean <- list('r1'=list(qps=apply(do.call(rbind, lapply(sims.qps, function(x) x$r1$qps)), 2, mean),
                              def.prob=apply(do.call(rbind, lapply(sims.qps, function(x) x$r1$def.probability)), 2, mean),
                              prop=mean(sapply(sims.qps, function(x) x$r1$prop))),
                    'r2'=list(qps=apply(do.call(rbind, lapply(sims.qps, function(x) x$r2$qps)), 2, mean),
                              def.prob=apply(do.call(rbind, lapply(sims.qps, function(x) x$r2$def.probability)), 2, mean),
                              prop=mean(sapply(sims.qps, function(x) x$r2$prop))),
                    'rtByBin'=apply(sapply(sims.qps, function(x) x$rtByBin), 1, mean),
                    'accByBin'=apply(sapply(sims.qps, function(x) x$accByBin), 1, mean),
                    'skewByBin'=apply(sapply(sims.qps, function(x) x$skewByBin), 1, mean),
                    'iqrByBin'=apply(sapply(sims.qps, function(x) x$iqrByBin), 1, mean))
  
  # plot data, mean simulation
  plotCDF(qps=dat.qps$r1, add=TRUE, type='b', lty=3, lwd=3, col=1)
  plotCDF(qps=dat.qps$r2, add=TRUE, type='b', lty=3, lwd=3, col=2)
  plotCDF(sims.mean$r1, add=TRUE, type='b', pch=4, lwd=2, col=1)
  plotCDF(sims.mean$r2, add=TRUE, type='b', pch=4, lwd=2, col=2)
  legend('bottomright', c('Data', 'Model'), lty=c(3, 1), lwd=c(3,2), pch=c(1, 4), bty='n')
  
  # round parameters for title
  if('outDEoptim' %in% ls(envir = .GlobalEnv)) {
    bestparsNoConstants <- c(bestparsNoConstants, 'nLL'=outDEoptim$optim$bestval)
  }
  pars.rounded <- c()
  for(p in 1:length(bestparsNoConstants)) {
    if(names(bestparsNoConstants)[p] %in% c('bU', 'eta1', 'eta2', 'bB')) {
      pars.rounded <- c(pars.rounded, round(bestparsNoConstants[p], 4))
    } else {
      pars.rounded <- c(pars.rounded, round(bestparsNoConstants[p], 2))
    }
  }
  title(main=paste('Sub', pp), sub=paste(names(bestparsNoConstants), '=', pars.rounded, collapse='  '))
  
  # plot mRT, acc, skew, IQR over time
  plot(x=1:10, y=sims.mean$rtByBin, type='l', ylim=c(0.3, .9), xlab='Trial (binned)', ylab='Mean RT')
  lines(x=1:10, y=dat.qps$rtByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  
  plot(x=1:10, y=sims.mean$accByBin, type='l', ylim=c(0.3, 1), xlab='Trial (binned)', ylab='Accuracy')
  lines(x=1:10, y=dat.qps$accByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  
  plot(x=1:10, y=sims.mean$skewByBin, type='l', ylim=c(0, 3), xlab='Trial (binned)', ylab='Skewness')
  lines(x=1:10, y=dat.qps$skewByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  
  plot(x=1:10, y=sims.mean$iqrByBin, type='l', ylim=c(0, 1), xlab='Trial (binned)', ylab='IQR')
  lines(x=1:10, y=dat.qps$iqrByBin, type='b')
  legend('topright', c('Data', 'Model'), lty=c(2, 1), pch=c(1, NA),  bty='n')
  
  # return info on quantiles
  return(list(all.dat.qps=dat.qps, all.sims.qps=sims.mean))
}


getMeanSimQps <- function(all.qps) {
  list('r1' = list(qps=apply(t(sapply(all.qps, function(x) x$r1$qps)), 2, mean, na.rm=TRUE),
                   def.prob=apply(t(sapply(all.qps, function(x) x$r1$def.prob)), 2, mean, na.rm=TRUE)),
       'r2' = list(qps=apply(t(sapply(all.qps, function(x) x$r2$qps)), 2, mean, na.rm=TRUE),
                   def.prob=apply(t(sapply(all.qps, function(x) x$r2$def.prob)), 2, mean, na.rm=TRUE)),
       'rtByBin' = apply(sapply(all.qps, function(x) x$rtByBin), 1, mean, na.rm=TRUE),
       'accByBin' = apply(sapply(all.qps, function(x) x$accByBin), 1, mean, na.rm=TRUE),
       'skewByBin' = apply(sapply(all.qps, function(x) x$skewByBin), 1, mean, na.rm=TRUE),
       'iqrByBin' = apply(sapply(all.qps, function(x) x$iqrByBin), 1, mean, na.rm=TRUE))
}

getMeanDatQps <- function(all.aps) {
  list('r1' = list(qps=apply(t(sapply(all.aps, function(x) x$r1$qps)), 2, mean, na.rm=TRUE),
                   def.prob=apply(t(sapply(all.aps, function(x) x$r1$def.prob)), 2, mean, na.rm=TRUE)),
       'r2' = list(qps=apply(t(sapply(all.aps, function(x) x$r2$qps)), 2, mean, na.rm=TRUE),
                   def.prob=apply(t(sapply(all.aps, function(x) x$r2$def.prob)), 2, mean, na.rm=TRUE)),
       'meanRtByBin' = apply(sapply(all.aps, function(x) x$rtByBin), 1, mean, na.rm=TRUE),
       'sdRtByBin' = apply(sapply(all.aps, function(x) x$rtByBin), 1, sd, na.rm=TRUE),
       'meanAccByBin' = apply(sapply(all.aps, function(x) x$accByBin), 1, mean, na.rm=TRUE),
       'sdAccByBin' = apply(sapply(all.aps, function(x) x$accByBin), 1, sd, na.rm=TRUE),
       'meanSkewByBin' = apply(sapply(all.aps, function(x) x$skewByBin), 1, mean, na.rm=TRUE),
       'sdSkewByBin' = apply(sapply(all.aps, function(x) x$skewByBin), 1, sd, na.rm=TRUE),
       'meanIqrByBin' = apply(sapply(all.aps, function(x) x$iqrByBin), 1, mean, na.rm=TRUE),
       'sdIqrByBin' = apply(sapply(all.aps, function(x) x$iqrByBin), 1, sd, na.rm=TRUE))
}