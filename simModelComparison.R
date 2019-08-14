# libraries and scripts
library("stringr")
library("ggplot2")
library("dplyr")
library("tidyr")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/modelComparisonFxs.R")
source("subFxs/plotThemes.R")
load("wtwSettings.RData")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs) 
load("genData/expDataAnalysis/blockData.RData")
# select common useID
encode = "RL2"
idList = hdrData$ID
modelNames = factor(c("QL1", "QL2", "RL1", "RL2", "BL"),
                    levels = c("QL1", "QL2", "RL1", "RL2", "BL"))
nModel = length(modelNames)
useID_ = vector(mode = "list", length = nModel)
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  paraNames = getParaNames(modelName)
  expPara = loadExpPara(paraNames, sprintf("genData/simModelFitting/%s/%sdb", encode, modelName))
  useID_[[i]] = getUseID(expPara, paraNames)
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
    fileName = sprintf("genData/simModelFitting/%s/%sdb/s%s_waic.RData", encode, modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic # here is like loglikelyhood, larger the better 
    logLik_[sIdx, m] = WAIC$elpd_waic  + WAIC$p_waic / 2
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
# save output for modelComparision 
output = data.frame(logEvidence_)
f= "genData/simModelFitting/logEvidenceList.csv"
write.table(file = f, output, sep = ",", col.names = F, row.names = F)

# participants best desribed by 
library("ggpubr")
bestNums = sapply(1 : nModel, function(i) sum(apply(logEvidence_[,1:nModel], MARGIN = 1, FUN = function(x) which.max(x) == i)))
data.frame(model = modelNames, bestNums = bestNums) %>%  ggplot(aes(x="", y=bestNums, fill=model)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + ylab("") + xlab("") + ggtitle(sprintf("Best described (n = %d)", nUse))+ 
  myTheme
dir.create("figures/simModelComparison")
dir.create(sprintf("figures/simModelComparison/%s", encode))
ggsave(sprintf("figures/simModelComparison/%s/loo_nBest.png", encode), width = 5, height =3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% ggplot(aes(model, pwaic)) + geom_boxplot() + myTheme
ggsave(sprintf("figures/simModelComparison/%s/loo_pwaic.png", encode), width = 5, height =3.5)

data.frame(pwaic = as.vector(pWaic_), model = rep(modelNames, each = nUse)) %>%
  group_by(model) %>% 
  summarise(muData = mean(pwaic), seData = sd(pwaic) / sqrt(length(pwaic)),
              minData = muData - seData, maxData = muData + seData) 

