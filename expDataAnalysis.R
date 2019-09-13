# in this dataset, only trials within the 7 mins will be kept. Therefore, we don't need to delete any data
datasetColors = c('#c53932', '#529D3E', '#3976AF')
# determine whether to truncate data
isTrun = F

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

# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for',n,'subjects.\n')

# define nBlock
nBlock = 2

# control which individual-level plots to generate
plotTrialwiseData = F
plotKMSC = F
plotWTW = F

# initialize outputs, organised by block
AUC = numeric(length =n * nBlock)
conditions = vector(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
nExclude =  numeric(length =n * nBlock)
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
  #if(blockData[blockData$id == thisID, "AUC"] > 20 & blockData$condition[blockData$id == thisID] == "LP"){
  for (bkIdx in 1: nBlock){
    noIdx = (sIdx - 1) * nBlock + bkIdx 
    # select data 
    thisTrialData = trialData[[thisID]]
    thisBlockIdx = (thisTrialData$blockNum == bkIdx)
    thisTrialData = thisTrialData[thisBlockIdx,]
    # truncate the last min(tMaxs) seconds
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "Rising", 1, 2)
    if(isTrun){
      excludedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
      thisTrialData = thisTrialData[! (1 : nrow(thisTrialData) %in% excludedTrials),]
      nExclude[[noIdx]] = length(excludedTrials)
    }
    # generate arguments for later analysis 
    label = sprintf('Subject %s, Cond %s',thisID, unique(thisTrialData$condition))
    tMax = min(tMaxs)
    
    # calcualte totalEarnings
    conditions[noIdx] = unique(thisTrialData$condition)
    totalEarnings[noIdx] =  sum(thisTrialData$trialEarnings)
    timeWaited = thisTrialData$timeWaited
    trialEarnings = thisTrialData$trialEarnings
    scheduledWait = thisTrialData$scheduledWait
    timeWaited[trialEarnings == loseValue] = scheduledWait[trialEarnings == loseValue]
    nAction[noIdx] = sum(round( ceiling(timeWaited / stepDuration)))
    nTrial[noIdx] = length(timeWaited)
    # calculate varQuitTime
    stdQuitTime[noIdx] = ifelse(totalEarnings[noIdx] == loseValue, NA, sd(timeWaited[trialEarnings == loseValue]))
    cvQuitTime[noIdx] = ifelse(totalEarnings[noIdx] == loseValue, NA, sd(timeWaited[trialEarnings == loseValue]) / mean(timeWaited[trialEarnings == loseValue]))
    muQuitTime[noIdx] = mean(timeWaited[trialEarnings == loseValue])
    nQuit[noIdx] = sum(trialEarnings == loseValue)
      
    # plot trial-by-trial data
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
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
  } # loop over blocks
}

# save data
blockData = data.frame(id = rep(allIDs, each = nBlock), blockNum = rep( t(1 : nBlock), n),
                       cbal = rep(hdrData$cbal,each = 2),
                       condition = factor(conditions, levels = c("Rising", "Falling")),
                       AUC = AUC, wtwEarly = wtwEarly,totalEarnings = totalEarnings,
                       nAction = nAction, stdQuitTime = stdQuitTime, cvQuitTime = cvQuitTime,
                       muQuitTime = muQuitTime, nQuit = nQuit, nTrial = nTrial, stdWd = stdWd, cvWd = cvWd,
                       nExclude = nExclude)
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridBlock.RData')
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')

# descriptive statistics for individual subjects and blocks
# for (sIdx in 1 : n) {
#   thisID = allIDs[sIdx]
#   thisTrialData = trialData[[thisID]]
#   label = sprintf('Subject %s',thisID)
#   trialPlots(block2session(thisTrialData),label)
#   readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
#   graphics.off()
# }

# demonstrate the optimism bias
data.frame(diff = blockData$AUC[blockData$condition == "Falling"] -
             blockData$AUC[blockData$condition == "Rising"],
           cbal = factor(blockData$cbal[blockData$condition == "Falling"],
                         labels = c("Rise-Fall", "Fall-Rise"))) %>% 
  ggplot(aes(x = cbal, y = diff)) + geom_boxplot() +
  stat_compare_means(comparisons = list(c("Rise-Fall", "Fall-Rise")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6) + ylim(c(-30, 30)) +
  xlab("") + ylab("Fall - Rise (s)") + myTheme
ggsave(sprintf("figures/expDataAnalysis/optimism.png"), width = 4, height = 3)
  

blockData %>% 
  mutate(cbal = factor(blockData$cbal,
                       labels = c("Rise-Fall", "Fall-Rise"))) %>%
  ggplot(aes(condition, AUC)) + geom_boxplot() +
  facet_grid(~cbal) + myTheme
ggsave(sprintf("figures/expDataAnalysis/AUC_cond.png"), width = 4, height = 3)

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
ggsave(sprintf("figures/expDataAnalysis/zTruc_AUC.png"), width = 4, height = 3)

# plot for sumData
wTest = wilcox.test( blockData[blockData$condition == "Rising", "AUC"],
                     blockData[blockData$condition == "Falling", "AUC"], paired = T)
datasetColors = c('#c53932', '#529D3E', '#3976AF')
data.frame(HPAUC = blockData$AUC[blockData$condition == 'Rising'],
           LPAUC = blockData$AUC[blockData$condition == 'Falling']) %>% 
  ggplot(aes(LPAUC, HPAUC))  + geom_point(shape = 21, size = 4, colour = "black", fill = datasetColors[3]) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0, min(tMaxs)))+ ylim(c(0, min(tMaxs))) +
  xlab("LP AUC (s)") + ylab("HP AUC (s)") +
  annotate("text", x = 12, y = 3, label = sprintf('p < 0.001***', wTest$p.value))+
  sumTheme
ggsave("figures/expDataAnalysis/zTruc_AUC_Cmp_sum.eps", width = 4, height = 3) 
  
# plot wtw 
# select = (blockData$blockNum == 2)
select = rep(T, n * nBlock)
plotData = data.frame(wtw = unlist(timeWTW_[select]), time = rep(tGrid, sum(select)),
           condition = rep(blockData$condition[select], each = length(tGrid))) %>% group_by(condition, time) %>%
  summarise(mean = mean(wtw), se = sd(wtw) / sqrt(length(wtw)), min = mean - se, max = mean + se) 

policy = data.frame(condition = c("Rising", "Falling"), wt = unlist(optimWaitTimes))
plotData %>% ggplot(aes(time, mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max),alpha = 0.3,  colour=NA) +
  geom_line(size = 1) + facet_wrap(~condition, scales = "free") +
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors) +
  xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60 * 5),
                     labels = paste(seq(0, 10, by = 5))) + 
  ylab("Willingness to wait (s)") +
  myTheme + ylim(c(7, 18))
ggsave("figures/expDataAnalysis/zTruc_wtw_timecourse.png", width = 6, height = 3)

# sum WTW
tempt = vector(length = nrow(plotData))
tempt[plotData$condition == "Rising"] = "HP"
tempt[plotData$condition == "Falling"] = "LP"
policy = data.frame(cond= c("HP", "LP"), wt = unlist(optimWaitTimes))
plotData %>% mutate(cond = factor(condition, levels = conditions, labels = c('HP', 'LP'))) %>%
  ggplot(aes(time, mean)) +
  geom_ribbon(aes(ymin=min, ymax=max),fill = '#a6bddb', colour= '#a6bddb') +
  geom_line(size = 1, color = datasetColors[3]) + facet_wrap(~cond, scales = "free")  +
  xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60 * 5),
                     labels = paste(seq(0, 10, by = 5))) + 
  ylab("WTW (s)") +
  sumTheme + ylim(c(7, 18))
ggsave("figures/expDataAnalysis/zTruc_wtw_timecourse.eps", width = 6, height = 2)
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
ggsave("figures/expDataAnalysis/zTruc_kmsc_timecourse.png", width = 5, height = 4) 


