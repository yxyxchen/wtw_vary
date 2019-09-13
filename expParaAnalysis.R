library("ggplot2")
library("dplyr")
library("tidyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")
# load trialData since we need scheduledWait 
allData = loadAllData()
hdrData = allData$hdrData 
trialData = allData$trialData       
allIDs = hdrData$ID 

modelName = "RL2"

# create output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "genData/expModelFitting"
dirName = sprintf("%s/%sdb",parentDir, modelName)
expPara= loadExpPara(paraNames, dirName)
useID = getUseID(expPara, paraNames)

# plot hist
expPara %>% filter(id %in% useID) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 3)

expPara %>% filter(id %in% useID) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>% 
  group_by(para) %>% summarise(mu = mean(value), median = median(value))

load('genData/expDataAnalysis/blockData.RData')
expPara$nQuit = blockData$nQuit[1 : 80 %% 2 == 0] + blockData$nQuit[1 : 80 %% 2 == 1]
wilcox.test(expPara$nega[expPara$nQuit >= 10] - 1)
expPara %>% filter(nQuit >= 10) %>%ggplot(aes(nega)) + geom_histogram(bins = 10) +
  geom_vline(xintercept = 1) + myTheme 
ggsave(sprintf('figures/expParaAnalysis/optimism_%s.png', modelName),
       width = 4, height = 4)

median(expPara$nega[expPara$nQuit >= 10])
