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
idList = hdrData$ID

# test useID
decodes = c("RL1", "RL2", "BL", "QL1", "QL2")
encodes = c("RL2", "QL2")
for(e in 1 : length(encodes)){
  encode = encodes[e]
  # load simData
  load(sprintf("genData/simulation/%s.RData", encode))
  ids = hdrData$ID
  nSub = length(ids)    
  for(d in 1 : length(decodes)){
    decode = decodes[d]
    # determine paraNames
    paraNames = getParaNames(decode)
    nPara = length(paraNames)
    if(paraNames == "wrong model name"){
      print(paraNames)
      break
    }
    # determine excID
    expPara = loadExpPara(paraNames,
                          sprintf("genData/simModelFitting/%s/%sdb", encode, decode))
    useID = getUseID(expPara, paraNames)
    excID = ids[!ids %in% useID]
    txt = sprintf("%s_%s:%s", encode, decode, length(useID))
    print(txt)
    readline("continue")
  }
}

# 
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, sprintf("genData/expModelFitting/%sdb", modelName))
simPara = loadExpPara(paraNames, sprintf("genData/simModelFitting/%s/%sdb", modelName, modelName))
expID = getUseID(expPara, paraNames)
simID = getUseID(simPara, paraNames)
useID = idList[idList %in% expID & idList %in% simID]

tempt  = cbind(exp = as.vector(as.matrix(expPara[,1:nPara])), sim = as.vector(as.matrix(simPara[,1 : nPara])))
data.frame(tempt, para = factor(rep(paraNames, each = nrow(simPara)), levels = paraNames),
           id = rep(simPara$id,  nPara)) %>% filter(id %in% useID) %>% 
  ggplot(aes(exp, sim)) + geom_point() + facet_wrap(.~para, scales = "free") +
  geom_abline(aes(slope = 1, intercept = 0), color = "red", linetype = "dashed") + 
  myTheme
dir.create("figures/paraRecovery")
ggsave(sprintf("figures/paraRecovery/%s.png", modelName) ,width = 6, height = 4)

