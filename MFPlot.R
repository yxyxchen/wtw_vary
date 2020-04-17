MFPlot = function(){
  library('dplyr')
  library("tidyr")
  library("ggplot2")
  library("latex2exp")
  library("ggpubr")
  source("subFxs/plotThemes.R")
  source('MFAnalysis.R')
  source("expSchematics.R")
  
  # normative analysis 
  iti = 2
  normResults = expSchematics(0, iti, F)
  optimWaitThresholds = normResults$optimWaitThresholds
  
  # output dir
  dir.create('figures/MFplot')
  
  # load experiment parameters
  load("expParas.RData")
  
  # plot WTW timecourses in two environments
  nBlock = 4
  MFResults = MFAnalysis(isTrct = F)
  blockStats = MFResults[['blockStats']]
  timeWTW_ = MFResults[['timeWTW_']]
  nSub = length(unique(blockStats$id))
  cbal = blockStats$cbal; cbal[cbal == 1] = "HPLP"; cbal[cbal == 2] = "LPHP"
  # define background colors
  greenData = data.frame(
    xmin = c(0, blockSec * 2, blockSec, blockSec * 3),
    xmax = c(blockSec, blockSec* nBlock, blockSec * 2, blockSec * 4) - max(delayMaxs),
    cbal = c("HPLP", "HPLP", "LPHP", "LPHP")
  )
  purpleData = data.frame(
    xmin = c(blockSec,  blockSec * 3, 0, blockSec * 2),
    xmax = c(blockSec* 2, blockSec * 4, blockSec, blockSec* 3) - max(delayMaxs),
    cbal = c("HPLP", "HPLP", "LPHP", "LPHP")
  )
  greyData = data.frame(
    xmin = 1:nBlock * blockSec - max(delayMaxs), xmax = 1:nBlock * blockSec
  )
  plotData = data.frame(wtw = unlist(timeWTW_),
                        time = rep(seq(0, blockSec * nBlock - 1, by = 1), nSub),
                        condition = rep(blockStats$condition, each = length(tGrid)),
                        blockIdx = rep(rep(1 : nBlock, length(tGrid)), nSub),
                        cbal = rep(cbal, each = length(tGrid))) %>%
    group_by(time, cbal) %>% 
    dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),min = mu- se, max = mu + se) 
  plotData$mu[mod(plotData$time, blockSec) == 0] = NA
  plotData %>%
    ggplot(aes(time, mu)) +
    geom_rect(greenData, mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20),
              fill = conditionColorBacks[1], inherit.aes = F) +
    geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20),
              fill = "#d9d9d9", inherit.aes = F) +
    geom_rect(purpleData, mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 20),
              fill = conditionColorBacks[2], inherit.aes = F) +
    geom_ribbon(aes(ymin=min, ymax=max), fill = "#737373") +
    geom_line(size = 1) +
    xlab("Task time (min)") + ylab("WTW (s)") + 
    myTheme +
    facet_grid(.~cbal) +
    scale_x_continuous(breaks = seq(0, blockSec * nBlock, by = blockSec),
                       labels = seq(0, blockSec * nBlock, by = blockSec) / 60) +
    scale_y_continuous(breaks = c(0,  10, 20), limits = c(0, 20)) 
  ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 3)
  ggsave("figures/MFPlot/wtw_timecourse.png", width = 5, height = 3) 
  
  ####################### separate analysis for the first and second two blocks ############
  MFResults = MFAnalysis(isTrct = T)
  # plot the first two blocks and the second two blocks separately 
  for(k in 1 : 2){
    blockStats = MFResults[['blockStats']]
    survCurve_ = MFResults$survCurve_
    
    # filter 
    if(k == 2){
      survCurve_= survCurve_[blockStats$blockNum > 2]
      blockStats = blockStats[blockStats$blockNum > 2,]
    }else{
      survCurve_= survCurve_[blockStats$blockNum <= 2]
      blockStats = blockStats[blockStats$blockNum <= 2,]
    }
    ############## plot AUCs in two environments #############
    # test
    wTest = wilcox.test( blockStats[blockStats$condition == "HP", "muWTW"],
                         blockStats[blockStats$condition == "LP", "muWTW"],paired = T)
    blockStats %>% group_by(condition) %>% summarise(median = median(muWTW))
    # plot
    data.frame(muWTWHP = blockStats$muWTW[blockStats$condition == 'HP'],
               muWTWLP = blockStats$muWTW[blockStats$condition == 'LP'],
               cbal = blockStats$cbal[blockStats$condition == "HP"]) %>%
      ggplot(aes(muWTWLP, muWTWHP)) +
      geom_point(size = 5, shape = 21, stroke =1) +
      geom_abline(slope = 1, intercept = 0)  +
      xlab("LP muAUC / (s)") + ylab("HP muAUC / (s)") + 
      myTheme + xlim(c(-1,31)) + ylim(c(-1,31)) + 
      xlab("LP AUC (s)") + ylab("HP AUC (s)")
    ggsave(sprintf("figures/MFPlot/muWTW_comparison_%d.eps", k), width = 4, height = 4)
    ggsave(sprintf("figures/MFPlot/muWTW_comparison_%d.png", k), width = 4, height = 4)
    
    ################### plot CIPs in two environments ###################
    # test
    wTest = wilcox.test( blockStats[blockStats$condition == "HP", "stdWTW"],
                         blockStats[blockStats$condition == "LP", "stdWTW"],paired = T)
    blockStats %>% group_by(condition) %>% summarise(median = median(stdWTW))
    # plot
    data.frame(stdWTWHP = blockStats$stdWTW[blockStats$condition == 'HP'],
               stdWTWLP = blockStats$stdWTW[blockStats$condition == 'LP'],
               cbal = blockStats$cbal[blockStats$condition == "HP"]) %>%
      ggplot(aes(stdWTWLP, stdWTWHP)) +
      geom_point(size = 5, shape = 21, stroke =1) +
      geom_abline(slope = 1, intercept = 0)  +
      xlab(TeX("LP CIP (s^2)")) + ylab(TeX("HP CIP (s^2)")) + 
      myTheme + xlim(c(-1,16)) + ylim(c(-1,16)) 
    ggsave(sprintf("figures/MFPlot/stdWTW_comparison_%d.eps", k), width = 4, height = 4)
    ggsave(sprintf("figures/MFPlot/stdWTW_comparison_%d.png", k), width = 4, height = 4)
    
    ################### plot CIP and AUC correlations ###################
    cor.test(blockStats$muWTW, blockStats$stdWTW, method = "spearman")
    cor.test(blockStats$muWTW[blockStats$condition == "HP"], blockStats$stdWTW[blockStats$condition == "HP"], method = "spearman")
    cor.test(blockStats$muWTW[blockStats$condition == "LP"], blockStats$stdWTW[blockStats$condition == "LP"], method = "spearman")
    blockStats %>% ggplot(aes(muWTW, stdWTW)) +
      geom_point(aes(color = condition), size = 3) +
      facet_grid(~condition) + myTheme +
      xlab("AUC (s)") + ylab(TeX("CIP ($s^2$)")) +
      scale_color_manual(values = conditionColors) +
      theme(legend.position = "none")
    ggsave(sprintf("figures/MFPlot/stdWTW_muWTW_%d.eps", k), width = 8, height = 4)
    ggsave(sprintf("figures/MFPlot/stdWTW_muWTW_%d.png", k), width = 8, height = 4)
    
    ################################ plot survival curves #####################
    # optimal strategy
    optim = data.frame(
      t = rep(kmGrid,  2),
      surv = rep(1, length(kmGrid) * 2),
      condition = rep(conditions, each = length(kmGrid)),
      select = rep(1:2, each = length(kmGrid))
    ) 
    optim$surv[optim$condition == "LP" & kmGrid> optimWaitThresholds$LP] = 0 # quit after 2.2 s
    optim$surv[optim$condition == "HP" & kmGrid> optimWaitThresholds$LP] = NA # don't plot after 2.2 s
    optim$select[optim$condition == "HP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors 
    optim$select[optim$condition == "LP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors
    # stats test
    plotData = data.frame(survCurve = unlist(survCurve_),
                          time = rep(kmGrid, nSub * nCondition),
                          condition = rep(blockStats$condition, each = length(kmGrid)))
    isSig = sapply(1 : length(kmGrid) , function(i)
    {
      t = kmGrid[i]
      HP = plotData$survCurve[plotData$condition == "HP" & plotData$time == t]
      LP = plotData$survCurve[plotData$condition == "LP" & plotData$time == t]
      tempt = wilcox.test(HP, LP, paired = T)
      ifelse(tempt$p.value < 0.05, 1.01, NA)
    }
    )
    sigData = data.frame(
      isSig = isSig,
      t = kmGrid
    )
    ## plot
    plotData %>%
      group_by(condition, time) %>%
      dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
                       min = mu- se, max = mu + se) %>%
      ggplot(aes(time, mu, color = condition, fill = condition)) + geom_line() +
      geom_ribbon(aes(time, ymin = min, ymax = max), alpha = 0.5, color = NA) +
      geom_line(data = optim, aes(t, surv, color = condition, linetype = condition, alpha = condition), size = 1.2) +
      geom_line(data = data.frame(t = kmGrid[kmGrid > 2],surv = 1),
                aes(t, surv), color = conditionColors[1], size = 1.2, inherit.aes = F, alpha = 0.8) + 
      geom_point(data = sigData, aes(t, isSig), inherit.aes = F, color = "black", shape = 4, size = 0.8) + 
      scale_fill_manual(values = conditionColors) +
      scale_color_manual(values = conditionColors) +
      scale_linetype_manual(values = c("solid", "dotted")) +
      scale_alpha_manual(values = c(0.8, 1))+
      xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
      theme(legend.position = "none") 
    ggsave(sprintf("figures/MFPlot/survival_curve_%d.eps", k), width = 4, height = 4)
    ggsave(sprintf("figures/MFPlot/survival_curve_%d.png", k), width = 4, height = 4) 
  }
  
}

# mixed effect anova, I should be more interested in non-stressed. 
library("lme4")
library(lmerTest)
blockStats$blockNum = factor(rep(c(-1, 1), 40))
blockStats$cIdx = ifelse(c)
blockStats %>% group_by(blockNum, condition) %>% summarise(mu = mean(muWTW))
options(contrasts = rep ("contr.sum", 2))
fit = lmer(muWTW ~  condition * blockNum + (1|id),
           data = blockStats)
summary(fit)
library("car")
Anova(fit, type="III")


