# this script contains helper analysis functions 
library(coin)
# check the distribution of scheduled delays
# ...as measured in number of key presses (for the instrumental version of the task)
scheduledDelays <- function(blockData,label) {
  cat(sprintf('Scheduled delays for %s\n',blockLabel))
  bkDelays = blockData$scheduledWait
  print(summary(bkDelays))
  # empirical cumulative distribution of scheduled delays
  fn <- ecdf(blockData$scheduledWait)
  plot(fn, main = sprintf('Scheduled delays: %s',blockLabel), xlab='Scheduled delay (s)',
       ylab='Cumulative proportion', xlim=c(0,30))
  # autocorrelation function
  # acfOutput <- acf(bkDelays, lag.max=20, main = sprintf('Scheduled delays: %s',blockLabel))
}


# plot trialwise responses in detail
trialPlots <- function(thisTrialData,label = " ") {
  # change the trialNum to accumulated trialNum if needed
  if(length(unique(thisTrialData$blockNum)) > 1){
    nTrial1 = sum(thisTrialData$blockNum == 1)
    nTrial1n2 = sum(thisTrialData$blockNum == 1 | thisTrialData$blockNum == 2)
  }
  # values to be plotted
  rwdIdx = thisTrialData$trialEarnings > loseValue
  quitIdx = thisTrialData$trialEarnings <= loseValue
  rwdTrialNo = thisTrialData$trialNum[rwdIdx]
  quitTrialNo = thisTrialData$trialNum[quitIdx]
  rwdSchedDelay = thisTrialData$scheduledWait[rwdIdx]
  quitSchedDelay = thisTrialData$scheduledWait[quitIdx]
  waitDuration = thisTrialData$timeWaited
  quitTime = waitDuration[quitIdx]
  # other parameters
  nTrials = length(thisTrialData$trialEarnings)
  # prepare plotData
  plotData = data.frame("trialNum" = c(rwdTrialNo, quitTrialNo, quitTrialNo),
                        "trialDuration" = c(rwdSchedDelay, quitTime, quitSchedDelay),
                        "condition" = rep(c('reward', 'quit', 'quitSchedule'), time = 
                                        c(length(rwdTrialNo), length(quitTrialNo),
                                          length(quitTrialNo))))
  plotData$condition = factor( plotData$condition, levels = c('reward', 'quit', 'quitSchedule'))
  # plot the main figure
  p = ggplot(plotData, aes(trialNum, trialDuration, color = condition)) + geom_point() +
  geom_line(data = plotData[plotData$condition != 'quitSchedule',],
            aes(trialNum, trialDuration, color = condition)) +
    scale_color_manual(values = c('blue', 'red', 'black')) + 
    xlab('Trial') + ylab('Waiting duration (s)') + ggtitle(label) + myTheme
  # add block lines if we have multiple blocks 
  if(length(unique(thisTrialData$blockNum)) > 1){
    p = p + geom_vline(xintercept = c(nTrial1,nTrial1n2),  linetype='dashed',
                   color = "grey", size = 1)
  }
  print(p)
  return(p)
}


# calculate kaplan-meier and area under the curve
kmsc <- function(thisTrialData,tMax,label='',plotKMSC=FALSE,grid=0) {
  library(survival)
  waitDuration = thisTrialData$timeWaited
  quitIdx = (thisTrialData$trialEarnings == 0)
  # for rewarded trials, base the duration on the reward delivery time (not the subsequent response)
  waitDuration[!quitIdx] <- thisTrialData$scheduledWait[!quitIdx]
  # fit the survival function
  kmfit <- survfit(Surv(waitDuration, quitIdx, type='right') ~ 1, 
                 type='kaplan-meier', conf.type='none', start.time=0, se.fit=FALSE)
  # extract elements of the survival curve object (?survfit.object)
  kmT = kmfit$time
  kmF = kmfit$surv
  # add a point at zero, since "kaplan-meier" starts from the first event
  kmT = c(0, kmT)
  kmF = c(1, kmF)
  # keep only points up through tMax 
  # if you use the same tMax for both condition, otherwise you can skip this
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # extend the last value to exactly tMax
  # notice that kmT is not evenly spaced
  kmT <- c(kmT, tMax)
  kmF <- c(kmF, tail(kmF,1))
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  kmFmy = kmF
  kmFmy[length(kmFmy)] = 0
  # sum(tail(kmT, -1) * diff((1 - kmFmy)))
  stdWd = sqrt(sum((tail(kmT, -1) - auc)^2 * diff((1 - kmFmy))))
  
  # calculate variance
  # plot if requested
  if (plotKMSC) {
   plotData = data.frame(kmT = kmT, kmF = kmF)
   p = ggplot(plotData, aes(kmT, kmF)) + geom_line() + xlab('Delay (s)') +
      ylab('Survival rate') + ylim(c(0,1)) + xlim(c(0,tMax)) +
        ggtitle(sprintf('KMSC: %s (AUC = %1.1f)',label,auc)) + 
        displayTheme
   print(p)
  }
  # put the survival curve on a standard grid
  kmOnGrid = vector()
  for (gIdx in 1:length(grid)) {
    g = grid[gIdx]
    # use the last point where t is less than or equal to the current grid value
    kmOnGrid[gIdx] = kmF[max(which(kmT<=g))]
  }
  return(list(kmT=kmT, kmF=kmF, auc=auc, kmOnGrid=kmOnGrid, stdWd = stdWd))
}

# this function can truncate trials in the simualtion object
# which enables us to zoom in and look and specific trials
truncateTrials = function(data, startTidx, endTidx){
  nVar = length(data)
  varNames = names(data)
  outputs = vector(mode = "list", length = nVar)
  for(i in 1 : nVar){
    junk = data[[varNames[i]]]
    if(is.matrix(junk)) outputs[[i]] = junk[, startTidx:endTidx]
    else outputs[[i]] = junk[startTidx:endTidx]
  }
  names(outputs) = varNames
  return(outputs)
}

# willingness to wait time-series
wtwTS <- function(thisTrialData, tGrid, wtwCeiling, label = "", plotWTW = F) {
  trialWTW = numeric(length = length(thisTrialData$trialEarnings)) # initialize the per-trial estimate of WTW
  quitIdx = thisTrialData$trialEarnings == 0
  # use either the rewardTime (for reward trials) or time waited (for quit trials)
  #   (not using time waited for reward trials because this includes the post-reward RT)
  timeWaited = thisTrialData$scheduledWait # use rewardtime make more sense but sometime nan
  timeWaited[quitIdx] = thisTrialData$timeWaited[quitIdx]
  ### find the longest time waited up through the first quit trial
  #   (or, if there were no quit trials, the longest time waited at all)
  #   that will be the WTW estimate for all trials prior to the first quit
  firstQuit = which(quitIdx)[1]
  if (is.na(firstQuit)) {firstQuit = length(thisTrialData$trialEarnings)} # if no quit, set to the last trial
  currentWTW = max(timeWaited[1:firstQuit])
  # start from the trial before the current quit
  thisTrialIdx = firstQuit - 1
  trialWTW[1:thisTrialIdx] = currentWTW
  ### iterate through the remaining trials, updating currentWTW
  while (thisTrialIdx < length(thisTrialData$trialEarnings)) {
    thisTrialIdx = thisTrialIdx + 1
    if (quitIdx[thisTrialIdx]) {currentWTW = timeWaited[thisTrialIdx]}
    else {currentWTW = max(currentWTW, timeWaited[thisTrialIdx])}
    trialWTW[thisTrialIdx] = currentWTW
  }
  ### impose a ceiling value, since trial durations exceeding some value may be infrequent
  trialWTW = pmin(trialWTW, wtwCeiling)
  ### convert from per-trial to per-second over the course of the block
  endTimeTrial = thisTrialData$sellTime
  timeWTW = trial2sec(trialWTW, endTimeTrial, tGrid)
  ### for testing
  # for testing: plot trialWTW on top of an individual's trialwise plot
  if(plotWTW){
    # for testing: plot timeWTW
    p = ggplot(data.frame(tGrid, timeWTW), aes(tGrid, timeWTW)) + geom_line() +
      xlab("Time in block (s)") + ylab("WTW (s)") + ggtitle(sprintf('WTW : %s', label)) +
      displayTheme
    print(p)
  }
  outputs = list(timeWTW = timeWTW,  trialWTW = trialWTW)
  return(outputs)
}
# this function maps trial-wise data to a continous time scale
trial2sec = function(dataTrial, endTimeTrial, tGrid){
  dataTime = numeric(length = length(tGrid))
  nTrial = length(dataTrial)
  binStartIdx = 1
  for(tIdx in 1 : nTrial){
    binEndTime = endTimeTrial[tIdx] 
    binEndIdx = max(which(tGrid < binEndTime)) # last grid point that falls within this trial
    dataTime[binStartIdx:binEndIdx] = dataTrial[tIdx]
    binStartIdx = binEndIdx + 1
  }
  # extend the final value to the end of the vector
  dataTime[binStartIdx:length(dataTime)] = dataTrial[nTrial]  
  return(dataTime)
}
# correlation plot
# the first col of plotData is x, the second col is y, the third col is the group
plotCorrelation = function(data, dotColor = "black",isRank){
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "cond")
  
  # calculate correlations
  corTests = lapply(1:2, function(i) cor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              method = "spearman"))
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~ data[data$cond == conditions[i], "y"]))
  rhos = sapply(1:2, function(i) round(as.numeric(corTests[[i]]$estimate), 3))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(as.numeric(pvalue(corTestsPerm[[i]])), 3))
  
  textColors = ifelse(ps < 0.05, "red", "blue")
  textData = data.frame(label = paste(rhos, "(p =", ps, ")"),
                        cond= c("HP", "LP"), color = textColors)
  # prepare rank 
  if(isRank){
    plotData = data %>% group_by(cond) %>% mutate(xRank = rank(x), yRank = rank(y))
  }
  
  # plot
  if(isRank){
    p0 = ggplot(plotData, aes(xRank, yRank)) + geom_point(size = 3, color = dotColor, fill = dotColor)
  }else{
    p0 = ggplot(plotData, aes(x, y)) + geom_point(size = 3, color = dotColor, fill = dotColor)
  }
  p = p0  + geom_text(data = textData,aes(x = -Inf,y = -Inf, label = label),
              hjust   = -0.2,vjust = -1,color = textColors, size = 5, fontface = 2, color = "#252525") +
    facet_grid(~cond)
 return(p)
} 


getCorrelation = function(data){
  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "cond")

  # calculate correlations
  # since we can't get rho from the later
  corTests = lapply(1:2, function(i) cor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              method = "kendall") 
                    )
  # supposedly, kendall can deal with data with a lot of ties
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~
  #                                                    data[data$cond == conditions[i], "y"]))
  rhos = sapply(1:2, function(i) as.numeric(corTests[[i]]$estimate))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(pvalue(corTestsPerm [[i]]), 3))
  return(list(rhos = rhos, ps = ps))
}

getPartCorrelation = function(data){
  library("ppcor")

  conditions = c("HP", "LP")
  colnames(data) = c("x", "y", "z", "cond")
  
  # calculate correlations
  # since we can't get rho from the later
  corTests = lapply(1:2, function(i) pcor.test(data[data$cond == conditions[i], "x"],
                                              data[data$cond == conditions[i], "y"],
                                              data[data$cond == conditions[i], "z"],
                                              method = "kendall") 
  )
  # supposedly, kendall can deal with data with a lot of ties
  # corTestsPerm = lapply(1:2, function(i) spearman_test(data[data$cond == conditions[i], "x"] ~
  #                                                    data[data$cond == conditions[i], "y"]))
  rhos = sapply(1:2, function(i) as.numeric(corTests[[i]]$estimate))
  ps = sapply(1:2, function(i) round(corTests[[i]]$p.value, 3))
  # ps = sapply(1:2, function(i) round(pvalue(corTestsPerm [[i]]), 3))
  return(list(rhos = rhos, ps = ps))
}

# convert data of multiple blocks into one session
block2session = function(tempt){
  nBlock = length(unique(tempt$blockNum))
  nTrials = sapply(1:nBlock, function(i) sum(tempt$blockNum == i))
  thisTrialData = within(tempt, {trialNum = trialNum + rep(c(0, cumsum(nTrials)[1:nBlock - 1]), time = nTrials);
  sellTime = sellTime + rep((1:nBlock-1) * blockSecs, time = nTrials);
  trialStartTime = trialStartTime + rep((1:nBlock-1) * blockSecs, time = nTrials);
  totalEarnings = totalEarnings +  rep(c(0, totalEarnings[cumsum(nTrials)[1:nBlock - 1]]),
                                       time = nTrials)
  })
  return(thisTrialData)
}


# kmscMoving 
kmscMoving = function(thisTrialData, tMax, label, plotKMSC, tGrid, window, by){
  library(survival)
  # initialize
  nTrial = length(thisTrialData$scheduledWait)
  # (nTrial - window + 1) / by is the number of gap, +1 to get nWindow, +2 to make one more
  nWindow = ceiling(((nTrial - window + 1) / by)) +1  
  endTimes = numeric(length = nWindow)
  survFits = vector(mode = "list", length = nWindow)
  # loop over windows
  nCore = parallel::detectCores() -1 # only for the local computer
  library("doMC")
  library("foreach")
  registerDoMC(nCore)
  survFits = foreach(i = 1 : nWindow) %dopar% {
    startTrial = 1 + by * (i-1)
    endTrial = min(startTrial + window - 1, nTrial)
    # survival analysis
    trialData = truncateTrials(thisTrialData, startTrial, endTrial)
    waitDuration = trialData$timeWaited
    quitIdx = (trialData$trialEarnings == 0)
    # for rewarded trials, base the duration on the reward delivery time (not the subsequent response)
    waitDuration[!quitIdx] <- trialData$scheduledWait[!quitIdx]
    tempt = fitSurv(waitDuration, quitIdx)
  }
  winEndTimes = sapply(1 : nWindow, function(i) thisTrialData$sellTime[min(i*by+window-by,
                                                                           nTrial)])
  winAUCs = sapply(1 : nWindow, function(i) survFits[[i]]$auc)
  winStdWds = sapply(1 : nWindow, function(i) survFits[[i]]$stdWd)
  # map to the continous time scale 
  timeAUCs = trial2sec(winAUCs, winEndTimes, tGrid)
  timeStdWds =  trial2sec(winStdWds, winEndTimes, tGrid)
  return(list(timeAUCs = timeAUCs, timeStdWds = timeStdWds,
              winAUCs = winAUCs, winStdWds = winStdWds))
}

fitSurv = function(waitDuration, quitIdx){
  # fit the survival function
  kmfit <- survfit(Surv(waitDuration, quitIdx, type='right') ~ 1, 
                   type='kaplan-meier', conf.type='none', start.time=0, se.fit=FALSE)
  # extract elements of the survival curve object (?survfit.object)
  kmT = kmfit$time
  kmF = kmfit$surv
  # add a point at zero, since "kaplan-meier" starts from the first event
  kmT = c(0, kmT)
  kmF = c(1, kmF)
  # keep only points up through tMax 
  # if you use the same tMax for both condition, otherwise you can skip this
  keepIdx = kmT<=tMax
  kmT <- kmT[keepIdx]
  kmF <- kmF[keepIdx]
  # extend the last value to exactly tMax
  # notice that kmT is not evenly spaced
  kmT <- c(kmT, tMax)
  kmF <- c(kmF, tail(kmF,1))
  # calculate auc
  auc <- sum(diff(kmT) * head(kmF,-1))
  kmFmy = kmF
  kmFmy[length(kmFmy)] = 0
  # sum(tail(kmT, -1) * diff((1 - kmFmy)))
  stdWd = sqrt(sum((tail(kmT, -1) - auc)^2 * diff((1 - kmFmy))))
  # return 
  return(list(kmT=kmT, kmF=kmF, auc=auc, stdWd = stdWd))
}

deMean = function(input){
  output = input - mean(input)
  return(output)
}

stand = function(input){
  output = (input - mean(input)) / sd(input)
  return(output)
}

lastTrunc = function(thisTrialData){
  excludedTrials = lapply(1 : nBlock, function(i)
    which(thisTrialData$trialStartTime > (blockSecs - tMaxs[i]) &
            (thisTrialData$blockNum == i)))
  includeStart = which(thisTrialData$trialNum == 1)
  includeEnd = sapply(1 : nBlock, function(i){
    if(length(excludedTrials[[i]] > 0)){
      min(excludedTrials[[i]])-1
    }else{
      max(which(thisTrialData$blockNum ==i))  
    }
  })
  tempt = lapply(1 : nBlock, function(i)
    truncateTrials(thisTrialData, includeStart[i], includeEnd[i]))
  thisTrialData = do.call("rbind", tempt)
  
}
