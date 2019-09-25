load("expParas.RData")
library("ggplot2"); library("Hmisc"); library("coin")
library("dplyr"); library("tidyr")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getparaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
source('MFAnalysis.R')

# model Name
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)

# output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# 
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]
# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "genData/expModelFit"
dirName = sprintf("%s/%s",parentDir, modelName)
expPara = loadExpPara(paraNames, dirName)
passCheck = checkFit(paraNames, expPara)


# plot hist 
# paraNames = c("LR", "LP", expression(tau), expression(gamma), "P")
# paraNames = c("LR", "LP", expression(tau), "P")
paraNames = paraNames
expPara$condition = sumStats$condition[1 : nrow(expPara)]
expPara %>% filter(passCheck ) %>% select(c(paraNames, "condition")) %>%
  gather(-c("condition"), key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(condition ~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 4)

# summary stats for expPara
expPara %>% filter(passCheck) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>% 
  group_by(para) %>% summarise(mu = mean(value), median = median(value))


# optimism bias
wilcoxResults = wilcox.test(expPara$phi_pos[passCheck] - expPara$phi_neg[passCheck])
# we use the 1.5 x the IQR for the larger of the two values 
# (phi_pos, phi_neg) as the criteria 
# for exclude outliers
junk = expPara %>% filter(passCheck) %>%
  select(c("phi_pos", "phi_neg")) %>%
  gather("paraName", "paraValue") %>%
  group_by(paraName) %>% 
  summarise(qLower = quantile(paraValue)[1],
            qUpper = quantile(paraValue)[3],
            IQR = qUpper - qLower,
            limLower = qLower - IQR * 1.5,
            limUpper = qUpper + IQR * 1.5)
limLower = max(0, min(junk$limLower))
limUpper = max(junk$limUpper)
# optimism bias
expPara %>% filter(passCheck) %>% 
  ggplot(aes(phi_neg, phi_pos)) +
  geom_point(color = themeColor,  fill = "#9ecae1", shape= 21, stroke = 1, size = 5) + 
  xlim(c(limLower, limUpper)) + ylim(c(limLower, limUpper)) + 
  geom_abline(slope = 1, intercept = 0) + 
  annotate("text", x = 0.015, y = 0.015, label = sprintf("p = %.3f", wilcoxResults$p.value)) +
  myTheme
ggsave("figures/expParaAnalysis/optimism.eps", width = 6, height = 6)

