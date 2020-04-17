set.seed(123)
# create output directories
dir.create("figures/expModelRep/")
dir.create(sprintf("figures/expModelRep/%s",modelName))

# load experiment parameters
load('expParas.RData')

# load packages and sub functions 
library("tidyverse")
library("latex2exp")
source("expSchematics.R")
source("subFxs/plotThemes.R")
source("subFxs/helpFxs.R") 
source("subFxs/loadFxs.R") # 
source("subFxs/analysisFxs.R") 
source("subFxs/modelRep.R")

# load empirical data 
allData = loadAllData()
hdrData = allData$hdrData 
trialData = allData$trialData       
ids = hdrData$id
nSub = length(ids)

## WTW from empirical data 
source("MFAnalysis.R")
MFResults = MFAnalysis(isTrct = T)
blockStats = MFResults[['blockStats']]
muWTWEmp = blockStats$muWTW
stdWTWEmp = blockStats$stdWTW

# modelNames
models = c("optim_noise", "optim_noise_bias", 
           "QL1", "QL2", "RL1", "RL2")
modelLabels =  c("optim_noise", "optim_noise_bias", 
                 "Q1", "Q2", "R1", "R2")
nModel = length(models)

# num of repetitions 
nRep = 10

# loop over models
outputs = list()
passCheck_ = matrix(NA, nrow = nSub, ncol = nModel)
for(m in 1 : nModel){
  modelName = models[m]
  paraNames = getParaNames(modelName)
  nPara = length(paraNames)
  
  ## check fit
  expPara = loadExpPara(paraNames, sprintf("genData/expModelFit/%s", modelName))
  passCheck = checkFit(paraNames, expPara)
  passCheck_[,m] = passCheck
  # compare willingness to wait (WTW) from empirical and replicated data
  
  ## rep
  load(sprintf("genData/expModelRep/%s_trct.RData", modelName))
  output =  data.frame(
    muAUC = repOutputs$muWTWRep_mu,
    muCIP = repOutputs$stdWTWRep_mu
  )
  outputs[[m]] = output
  rm(repOutputs)
  rm(output)
}
allPass = apply(passCheck_, 1, sum) == nModel
sum(allPass)
######################### predict AUC #######################
plotData = bind_rows(outputs) %>%
  mutate(model = rep(modelLabels, each = nSub * nBlock),
         empAUC = rep(muWTWEmp, nModel),
         empCIP = rep(stdWTWEmp, nModel),
         condition = rep(blockStats$condition, nModel)) %>% 
  filter(rep(allPass, nModel * nBlock)) %>%
  arrange(model, condition, empAUC)  %>%
  mutate(rank = rep(c(1 : sum(blockStats$condition == "HP" & allPass), 1 : sum(blockStats$condition == "LP" & allPass)), nModel),
         type = ifelse(model %in% c("optim_noise", "optim_noise_bias"),
                        "baseline", "RL")) 
# main text version
plotData %>% filter(model %in% c("optim_noise", "optim_noise_bias", "Q2")) %>%
  ggplot(aes(rank, muAUC)) +
  geom_bar(data = plotData[plotData$model == "Q1",],
           aes(rank, empAUC), inherit.aes = F, stat = "identity",
           fill = NA, color = "grey") +
  geom_line(aes(color = model)) + 
  geom_point(aes(color = model)) + 
  myTheme + 
  ylab("AUC (s)") + xlab("") +
  facet_grid(~condition) +
  scale_x_continuous(breaks = c()) +
  xlab('') +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e7298a")) +
  theme(legend.title = element_blank())
fileName = "figures/expModelRep/AUC_all.eps"
ggsave(filename = fileName,  width = 6, height = 3)
fileName = "figures/expModelRep/AUC_all.png"
ggsave(filename = fileName,  width = 6, height = 3)

# full version
plotData %>%
  ggplot(aes(rank, muAUC)) +
  geom_bar(data = plotData[plotData$model == "Q1",],
           aes(rank, empAUC), inherit.aes = F, stat = "identity",
           fill = NA, color = "grey") +
  geom_line(aes(color = model, linetype = type)) + 
  geom_point(aes(color = model, shape = type)) + 
  myTheme + 
  ylab("AUC (s)") + xlab("") +
  facet_grid(~condition) +
  scale_x_continuous(breaks = c()) +
  xlab('') +
  scale_color_manual(values = c("#1f78b4", "#33a02c",
                                "#e7298a", "#e41a1c", "#d95f02", "#984ea3")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(NA, 16)) + 
  theme(legend.title = element_blank())
fileName = "figures/expModelRep/AUC_all_full.eps"
ggsave(filename = fileName,  width = 6, height = 3)
fileName = "figures/expModelRep/AUC_all_full.png"
ggsave(filename = fileName,  width = 6, height = 3)

######################### predict CIP #######################
plotData = bind_rows(outputs) %>%
  mutate(model = rep(modelLabels, each = nSub * nBlock),
         empAUC = rep(muWTWEmp, nModel),
         empCIP = rep(stdWTWEmp, nModel),
         condition = rep(blockStats$condition, nModel)) %>% 
  filter(rep(allPass, nModel * nBlock)) %>%
  arrange(model, condition, empCIP)  %>%
  mutate(rank = rep(c(1 : sum(blockStats$condition == "HP" & allPass), 1 : sum(blockStats$condition == "LP" & allPass)), nModel),
         type = ifelse(model %in% c("optim_noise", "optim_noise_bias"),
                       "baseline", "RL")) 
# full version
plotData %>%
  ggplot(aes(rank, muCIP)) +
  geom_bar(data = plotData[plotData$model == "Q1",],
           aes(rank, empCIP), inherit.aes = F, stat = "identity",
           fill = NA, color = "grey") +
  geom_line(aes(color = model, linetype = type)) + 
  geom_point(aes(color = model, shape = type)) + 
  myTheme + 
  ylab(TeX("CIP ($s^2$)")) + xlab("") +
  facet_grid(~condition) +
  scale_x_continuous(breaks = c()) +
  xlab('') +
  scale_color_manual(values = c("#1f78b4", "#33a02c",
                                "#e7298a", "#e41a1c", "#d95f02", "#984ea3")) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_shape_manual(values = c(NA, 16)) + 
  theme(legend.title = element_blank())
fileName = "figures/expModelRep/CIP_all_full.eps"
ggsave(filename = fileName,  width = 6, height = 3)
fileName = "figures/expModelRep/CIP_all_full.png"
ggsave(filename = fileName,  width = 6, height = 3)



