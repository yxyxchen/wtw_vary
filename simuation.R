# 
dir.create("genData")
dir.create("genData/simulation")
library('ggplot2')
library('plyr')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # called by simulate 
source("subFxs/helpFxs.R") # getParas
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load scheduledWait from empirical data
source("subFxs/analysisFxs.R") 

# modelName 
modelName = "RL2"
repFun = getRepFun(modelName)

# load expData
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$ID          
nSub = length(ids)      

# load expPara
paraNames = getParaNames(modelName)
parentDir ="genData/expModelFitting"; dirName = sprintf("%s/%sdb",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)

set.seed(123)
simTrialData = list()
for(sIdx in 1 : nSub){
  id = ids[sIdx]
  paras = as.double(expPara[expPara$id == id, 1 : length(paraNames)])
  # prepare input
  thisTrialData = trialData[[id]] # here we id instead of sIdx
  # excluded some trials
  excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                           thisTrialData$condition == conditions[1])
  excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                           thisTrialData$condition == conditions[2])
  excluedTrials = c(excluedTrials1, excluedTrials2)
  thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials & 
                                  thisTrialData$blockNum <= 2,]
  condition = thisTrialData$condition
  scheduledWait = thisTrialData$scheduledWait
  blockNum = thisTrialData$blockNum
  scheduledReward = thisTrialData$trialEarnings
  scheduledReward[scheduledReward == 0] = ifelse(scheduledWait[scheduledReward == 0] < 6.6195,
                                                 ifelse(condition[scheduledReward == 0] == "Rising", min(tokenValue), max(tokenValue)),
                                                 ifelse(condition[scheduledReward == 0]  == "Rising", max(tokenValue), min(tokenValue)))
  scheduledReward = sapply(scheduledReward, function(x) which(tokenValue == x))
  id = factor(ids[sIdx], levels = hdrData$ID)
  simTrialData[[id]] = repFun(paras, condition, scheduledWait, scheduledReward, blockNum)
}
hdrData$ID = factor(hdrData$ID, levels = hdrData$ID)
save(simTrialData, hdrData, file = "genData/simulation/simTrialData.RData")

