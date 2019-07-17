# libraries and scripts
library("stringr")
library("ggplot2")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/modelComparisonFxs.R")
source("subFxs/plotThemes.R")
library("dplyr")
library("tidyr")
load("wtwSettings.RData")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs) 
load("genData/expDataAnalysis/blockData.RData")
# select common useID
idList = hdrData$ID
modelNames = factor(c("PRbs", "PRbsNC", "Rlearn", "RlearnL"),
                    levels = c("PRbs", "PRbsNC", "Rlearn", "RlearnL"))
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paras = getParas(modelName)
  expPara = loadExpPara(paras, sprintf("genData/expModelFitting/%sdb", modelName))
  useID_[[i]] = factor(getUseID(expPara, paras), levels = levels(hdrData$ID))
}
useID = idList[apply(sapply(1 : nModel, function(i )idList %in% useID_[[i]]), MARGIN = 1,
              all)]
nUse = length(useID)

# extract logEvidece_ from loo 
logEvidence_ = matrix(NA, nUse, nModel)
logLik_ = matrix(NA, nUse, nModel)
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    fileName = sprintf("genData/expModelFitting/%sdb/s%s_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic # here is like loglikelyhood, larger the better 
    logLik_[sIdx, m] = WAIC$elpd_waic  + WAIC$p_waic / 2
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
# save output for modelComparision 
output = data.frame(logEvidence_)
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

# participants best desribed by 
library("ggpubr")
bestNums = sapply(1 : nModel, function(i) sum(apply(logEvidence_[,1:nModel], MARGIN = 1, FUN = function(x) which.max(x) == i)))
data.frame(model = modelNames, bestNums = bestNums) %>%  ggplot(aes(x="", y=bestNums, fill=model)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + ylab("") + xlab("") + ggtitle(sprintf("Best described (n = %d)", nUse))+ 
  myTheme
dir.create("figures/expModelComparison")
ggsave("figures/expModelComparison/loo_nBest.png", width = 5, height =3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% ggplot(aes(model, pwaic)) + geom_boxplot() + myTheme
ggsave("figures/expModelComparison/loo_pwaic.png", width = 5, height =3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% 
  summarise(muData = mean(pwaic), seData = sd(pwaic) / sqrt(length(pwaic)),
            minData = muData - seData, maxData = muData + seData) 

# extract logEvidence, cross validation
modelNames = c("PRbs", "PRbsNC", "Rlearn", "RlearnL")
nModel = length(modelNames)
ids = hdrData$ID
nSub = length(ids)
nFold = 10
logEvidence = matrix(nrow = length(ids), ncol= nModel) 
logEvidenceTrain = list(length = nModel)
for(mIdx in 1 : nModel){
  modelName = modelNames[mIdx]
  paras = getParas(modelName)
  nPara = length(paras)
  logLikFun = getLogLikFun(modelName)
  thisLogEvidenceTrain = matrix(nrow = nFold, ncol = nSub)
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    load(sprintf("genData/expModelFittingCV/split/s%d.RData", id))
    thisTrialData = trialData[[id]]
    excluedTrialsRise = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                                thisTrialData$condition == conditions[1])
    excluedTrialsFall = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                                thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrialsRise, excluedTrialsFall)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    thisTrialData = thisTrialData[thisTrialData$blockNum <= 2,]
    # prepare
    cond = thisTrialData$cond
    nTrial = length(thisTrialData$trialEarnings)
    trialEarnings = thisTrialData$trialEarnings
    timeWaited = pmin(thisTrialData$timeWaited, ifelse(cond == conditions[1], tMaxs[1], tMaxs[2]))
    Ts = round(ceiling(timeWaited / stepDuration) + 1)
    scheduledWait = thisTrialData$scheduledWait
    # 
    cvPara = loadCVPara(paras,
                      sprintf("genData/expModelFittingCV/%sdb",modelName),
                      pattern = sprintf("s%d_f[0-9]{1,2}_summary.txt", id))
    # initialize 
    LL_ = vector(length = nFold)
    if(length(getUseID(cvPara, paras)) == 10){
      for(f in 1 : nFold){
        # determine training end testing trials
        trials = partTable[f,]
        trials = trials[trials < nTrial]
        junk = 1 : nTrial
        trialsTrain = junk[!junk %in% trials]
        
        thisParas = as.double(cvPara[f,1:nPara])
        lik_ = logLikFun(thisParas, cond, trialEarnings, timeWaited)$lik_
        LL_[f] = sum(sapply(1 : length(trials), function(i){
          trial = trials[i]
          if(trialEarnings[trial] > 0){
            junk = log(lik_[1 : max(Ts[trial]-1, 1), trial])
            junk[is.infinite(junk)] = -10000
            sum(junk)
          }else{
            junk = c(log(lik_[1:max(Ts[trial] - 2,1), trial]), log(1-lik_[Ts[trial] - 1, trial]))
            junk[is.infinite(junk)] = -10000
            sum(junk)
          }
        }))
        thisLogEvidenceTrain[f, sIdx] = sum(sapply(1 : length(trialsTrain), function(i){
          trial = trialsTrain[i]
          if(trialEarnings[trial] > 0){
            junk = log(lik_[1 : max(Ts[trial]-1, 1), trial])
            junk[is.infinite(junk)] = -10000
            sum(junk)
          }else{
            junk = c(log(lik_[1:max(Ts[trial] - 2,1), trial]), log(1-lik_[Ts[trial] - 1, trial]))
            junk[is.infinite(junk)] = -10000
            sum(junk)         
          }
        }))
          
      }
      logEvidence[sIdx, mIdx] = sum(LL_)
      logEvidenceTrain[[mIdx]] = thisLogEvidenceTrain
    }
  }
}
select = apply(sapply(1 : nModel, function(i) !is.na(logEvidence[,i])), MARGIN = 1, FUN = all)
useID = ids[select]

output = data.frame(cvLik = logEvidence[select,])
f= "genData/expModelFitting/logEvidenceListCV.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

bestNums = sapply(1 : nModel, function(i) sum(apply(logEvidence[,1:nModel], MARGIN = 1, FUN = function(x) which.max(x) == i)))
data.frame(model = modelNames, bestNums = bestNums) %>%  ggplot(aes(x="", y=bestNums, fill=model)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + ylab("") + xlab("") + ggtitle(sprintf("Best described (n = %d)", nUse))+ 
  myTheme
dir.create("figures/expModelComparison")
ggsave("figures/expModelComparison/CV_nBest.png", width = 5, height =3.5)

