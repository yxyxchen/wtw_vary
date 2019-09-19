# libraries and scripts
library("ggplot2")
library("dplyr")
library("tidyr")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/plotThemes.R")

# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$id[hdrData$stress == "no_stress"]                 # column of subject IDs
nSub = length(ids) 

# check fit
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s/%s", modelName, modelName))
simPara = loadExpPara(paraNames, sprintf("genData/simModelFit/%s/%s", modelName, modelName)) 
passCheckIds = simPara$id[checkFit(paraNames, simPara)]
nPassCheck = length(passCheckIds)

data.frame(rbind(expPara[,c(1:nPara,25)], simPara[,c(1:nPara,25)])) %>%
  filter(id %in% passCheckIds) %>% 
  mutate(source = rep(c("exp", "sim"), each = nPassCheck )) %>%
  group_by(source) %>% 
  gather(key = "paraName", value = "paraValue", -source) %>%
  ggplot(aes())
  
  
  
  
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

