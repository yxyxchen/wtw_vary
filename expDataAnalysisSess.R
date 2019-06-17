# in this dataset, only trials within the 7 mins will be kept. Therefore, we don't need to delete any data
# determine whether to truncate data
isTrun = T
# load libraries
source('subFxs/loadFxs.R') # for loading data 
source('subFxs/analysisFxs.R') # for analysis 
source("subFxs/plotThemes.R")
library("ggplot2")
library('dplyr')
dir.create("genData")
dir.create("genData/expDataAnalysis")

# load setting parameters 
load("wtwSettings.RData")
nBlock = 4
if(isTrun){
  tGrid = seq(0, nBlock * blockSecs, by = 1)
}

if(isTrun){
  tGrid = seq(0, blockSecs * nBlock, by = 1) # here I use a truncated tGrid, according to max(sellTime) 
}else{
  tGrid = seq(0, blockSecs * nBlock, by = 1)
}
# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for',n,'subjects.\n')

# define nBlock
nBlock = 4

# control which individual-level plots to generate
plotTrialwiseData = F
plotKMSC = F
plotWTW = F

# initialize outputs, organised by block
AUC = numeric(length =n * nBlock)
conditions = vector(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
nAction = numeric(length =n * nBlock)
wtwEarly = numeric(length =n * nBlock)
timeWTW_ = vector(mode = "list", length = n * nBlock)
trialWTW_ = vector(mode = "list", length = n * nBlock)
kmOnGrid_ = vector(mode = "list", length = n * nBlock)
stdQuitTime = numeric(length =n * nBlock)
cvQuitTime = numeric(length =n * nBlock)
muQuitTime = numeric(length =n * nBlock)
nQuit = numeric(length =n * nBlock)
nTrial = numeric(length =n * nBlock)
stdWd = numeric(length =n * nBlock) # standard deviation from the survival curve for the whole block
cvWd =  numeric(length =n * nBlock)
# descriptive statistics for individual subjects and blocks
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  # select data 
  thisTrialData = trialData[[thisID]]
  # generate arguments for later analysis 
  label = sprintf('Subject %s, Cond %s',thisID, unique(thisTrialData$condition))
  tMax = min(tMaxs)
  # trunc
  if(isTrun){
    excluedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2] * 0.5))
    nExclude[[sIdx]] = length(excluedTrials)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
  }
  thisTrialData = block2session(thisTrialData)
  
  # calcualte totalEarsIdxisIdxgs
  conditions[sIdx] = unique(thisTrialData$condition)
  totalEarnings[sIdx] =  sum(thisTrialData$trialEarnings)
  timeWaited = thisTrialData$timeWaited
  trialEarnings = thisTrialData$trialEarnings
  scheduledWait = thisTrialData$scheduledWait
  timeWaited[trialEarnings > loseValue] = scheduledWait[trialEarnings > loseValue]
  nAction[sIdx] = sum(round(ifelse(trialEarnings > loseValue, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
  nTrial[sIdx] = length(timeWaited)
  # calculate varQuitTime
  stdQuitTime[sIdx] = ifelse(totalEarnings[sIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]))
  cvQuitTime[sIdx] = ifelse(totalEarnings[sIdx] == 0, NA, sd(timeWaited[trialEarnings == 0]) / mean(timeWaited[trialEarnings == 0]))
  muQuitTime[sIdx] = mean(timeWaited[trialEarnings == 0])
  nQuit[sIdx] = sum(trialEarnings == 0)
      
  # plot trial-by-trial data
  if (plotTrialwiseData) {
    trialPlots(thisTrialData,label)
    readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
    graphics.off()
  }
    
  # survival analysis
  kmscResults = kmsc(thisTrialData,min(tMaxs),label,plotKMSC, kmGrid)
  AUC[noIdx] = kmscResults[['auc']]
  kmOnGrid_[[noIdx]] = kmscResults$kmOnGrid
  stdWd[noIdx] = kmscResults$stdWd
  cvWd[noIdx] = kmscResults$stdWd / kmscResults$auc
    
  if (plotKMSC) {
    readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
    graphics.off()
  }

  # WTW time series
  wtwCeiling = min(tMaxs)
  wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
  timeWTW_[[noIdx]] = wtwtsResults$timeWTW
  trialWTW_[[noIdx]] = wtwtsResults$trialWTW
  wtwEarly[noIdx] =   wtwtsResults$trialWTW[1]
  # wait for input before continuing, if individual plots were requested
  if (plotWTW) {
    readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
    graphics.off()
  }
}

# save data
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridBlock.RData')
sessionData = data.frame(id = rep(allIDs, each = nBlock), blockNum = rep( t(1 : nBlock), n),
                       condition = factor(conditions, levels = c("Rising", "Falling")),
                       AUC = AUC, wtwEarly = wtwEarly,totalEarnings = totalEarnings,
                       nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                       muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial, stdWd = stdWd, cvWd = cvWd)
save(sessionData, file = 'genData/expDataAnalysis/sessionData.RData')

# descriptive statistics for individual subjects and blocks
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  thisTrialData = trialData[[thisID]]
  label = sprintf('Subject %s',thisID)
  trialPlots(block2session(thisTrialData),label)
  readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
  graphics.off()
}

# plot AUC in two conditions
library("ggpubr")
load("wtwSettings.RData")
blockData %>% ggplot(aes(condition, AUC)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) +
  scale_fill_manual(values = conditionColors) + 
  xlab("") + ylab("Average WTW(s)") + myTheme +
  stat_compare_means(comparisons = list(c("Rising", "Falling")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6) + ylim(c(0, 36))
dir.create("figures")
dir.create("figures/expDataAnalysis")
ggsave(sprintf("figures/expDataAnalysis/AUC.png"), width = 4, height = 3)

# plot wtw 
plotData = data.frame(wtw = unlist(timeWTW_), time = rep(tGrid, n),
           cbal = rep(as.factor(hdrData$cbal), each = length(tGrid))) %>% group_by(cbal, time) %>%
  summarise(mean = mean(wtw), se = sd(wtw) / sqrt(length(wtw)), min = mean - se, max = mean + se) 

policy = data.frame(cbal = c("1", "2"), wt = unlist(optimWaitTimes))
plotData %>% ggplot(aes(time, mean, color = cbal, fill = cbal)) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha = 0.3,  colour=NA) +
  geom_line(size = 1) + facet_wrap(~cbal, scales = "free") +
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors) +
  xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60 * 10),
                     labels = seq(0, max(tGrid), by = 60 * 10) / 60) + 
  ylab("Willingness to wait (s)") +
  myTheme + ylim(c(0, 32))
ggsave("figures/expDataAnalysisSess/wtw_timecourse.png", width = 6, height = 3)

# plot survival curve
condition =  rep(blockData$condition)
data.frame(kmsc = unlist(kmOnGrid_), time = rep(kmGrid, n * nBlock),
                      condition = factor(rep(condition, each = length(kmGrid))), levels = conditions) %>%
  group_by(condition, time) %>%
  summarise(mean = mean(kmsc), se = sd(kmsc) / sqrt(length(kmsc)), min = mean - se, max = mean + se) %>% 
  ggplot(aes(time, mean, color = condition, fill = condition)) + 
  geom_ribbon(aes(ymin=min, ymax=max),alpha = 0.3, colour=NA)+
  geom_line(size = 1.5) + myTheme + scale_fill_manual(values = conditionColors) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + scale_color_manual(values = conditionColors)

#   geom_line(data = ideal, aes(time, kmsc, color = condition), linetype = "dotted", size = 1)
ggsave("figures/expDataAnalysis/kmsc_timecourse.png", width = 5, height = 4) 



###
blockData$cbal = rep(hdrData$cbal, each = 4)
R2F = blockData %>% filter(cbal == 1) %>% with(AUC[blockNum == 3] - AUC[blockNum == 4])
F2R = blockData %>% filter(cbal == 2) %>% with(AUC[blockNum == 3] - AUC[blockNum == 4])

hist(R2F)
hist(F2R)

plotData = data.frame(value = c(abs(R2F), abs(F2R)), condition = c("R2F", "F2R"), time = c(length(R2F), length(F2R)))
library(ggpubr)
compare_means(value ~ condition, plotData)

library("ggplot2")
ggplot(plotData, aes(condition, value)) + geom_boxplot() + stat_compare_means() + myTheme + 
  ylab("AUC abs adaption (s)")
