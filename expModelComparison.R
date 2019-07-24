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
modelNames =  factor(c("QL1", "QL2", "RL1", "RL2", "BL"),
                     levels = c("QL1", "QL2", "RL1", "RL2", "BL"))
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paraNames = getParaNames(modelName)
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFitting/%sdb", modelName))
  useID_[[i]] = factor(getUseID(expPara, paraNames), levels = levels(hdrData$ID))
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
# extract logEvidence, cross validation
modelNames = factor(c("QL1", "QL2", "RL1", "RL2", "BL"),
                    levels = c("QL1", "QL2", "RL1", "RL2", "BL"))
nModel = length(modelNames)
ids = hdrData$ID
nSub = length(ids)
nFold = 10
logEvidence = matrix(nrow = length(ids), ncol= nModel) 
for(mIdx in 1 : nModel){
  modelName = modelNames[mIdx]
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  likFun = getLikFun(modelName)
  for(sIdx in 1 : nSub){
    id = ids[sIdx]
    load(sprintf("genData/expModelFittingCV/split/s%s.RData", id))
    thisTrialData = trialData[[id]]
    # excluded some trials
    excluedTrialsHP = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                              thisTrialData$condition == conditions[1])
    excluedTrialsLP = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                              thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrialsHP, excluedTrialsLP)
    thisTrialData = thisTrialData[(!(1 : nrow(thisTrialData)) %in% excluedTrials) &
                                    thisTrialData$blockNum <= 2,]
    # prepare the data
    nTrial = length(thisTrialData$trialEarnings)
    cond = thisTrialData$cond
    trialEarnings = thisTrialData$trialEarnings
    timeWaited = thisTrialData$timeWaited
    timeWaited[trialEarnings != 0] = thisTrialData$scheduledWait[trialEarnings != 0]
    Ts = round(ceiling(timeWaited / stepDuration) + 1)
    cvPara = loadCVPara(paraNames,
                        sprintf("genData/expModelFittingCV/%sdb",modelName),
                        pattern = sprintf("s%s_f[0-9]{1,2}_summary.txt", id))
    # test
    # tempt = read.csv(sprintf("genData/expModelFitting/%s/s%s.txt", modelName, id), header = F)
    # 
    # paras = as.double(tempt[1,1:nPara]) # use the para estmates from the training set
    # lik_ = likFun(paras, cond, trialEarnings, timeWaited)$lik_
    # LL_all = sum(sapply(1 : length(trialEarnings), function(i){
    #   if(trialEarnings[i] != 0){
    #     junk = log(lik_[1 : max(Ts[i]-1, 1), i])
    #     junk[is.infinite(junk)] = -10000
    #     sum(junk)
    #   }else{
    #     if(Ts[i] > 2){
    #       junk = c(log(lik_[1:max(Ts[i] - 2,1), i]), log(1-lik_[Ts[i] - 1, i]))
    #     }else{
    #       junk = log(1-lik_[Ts[i] - 1, i])
    #     }
    #     junk[is.infinite(junk)] = -10000
    #     sum(junk)
    #   }
    # }))
    # 
    # initialize 
    LL_ = vector(length = nFold)
    if(length(getUseID(cvPara, paraNames)) == 10){
      for(f in 1 : nFold){
        # determine training end testing trials
        trials = partTable[f,]
        trials = trials[trials < nTrial]
        junk = 1 : nTrial
        
        paras = as.double(cvPara[f,1:nPara]) # use the para estmates from the training set
        lik_ = likFun(paras, cond, trialEarnings, timeWaited)$lik_ # the likelyhood of waitting on the whole data sets
        # sum first all steps, then overl trials within a fold
        LL_[f] = sum(sapply(1 : length(trials), function(i){
          trial = trials[i]
          # change for this version
          if(trialEarnings[trial] != 0){
            junk = log(lik_[1 : max(Ts[trial]-1, 1), trial])
            junk[is.infinite(junk)] = -10000
            sum(junk)
          }else{
            if(Ts[trial] > 2){
              junk = c(log(lik_[1:max(Ts[trial] - 2,1), trial]), log(1-lik_[Ts[trial] - 1, trial]))
              junk[is.infinite(junk)] = -10000
              sum(junk)
            }else{
              junk = log(1-lik_[Ts[trial] - 1, trial])
              junk
            }
          }
        }))
      }
      logEvidence[sIdx, mIdx] = sum(LL_) # sum over folds 
    }
  }
}
select = apply(sapply(1 : nModel, function(i) !is.na(logEvidence[,i])), MARGIN = 1, FUN = all)
useID = ids[select]
nUse = length(useID)
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



