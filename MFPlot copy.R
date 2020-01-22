library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
source("subFxs/plotThemes.R")
source('MFAnalysis.R')

# output dir
dir.create('figures/MFplot')

# load experiment parameters
load("expParas.RData")

# plot WTW timecourses in two environments
MFResults = MFAnalysis(isTrct = F)
sumStats = MFResults[['sumStats']]
timeWTW_ = MFResults[['timeWTW_']]
nSub = nrow(sumStats)
cbal = sumStats$cbal; cbal[cbal == 1] = "HPLP"; cbal[cbal == 2] = "LPHP"
yellowData = data.frame(
  xmin = 0:3 * blockSec, xmax = 0:3 * blockSec + blockSec - max(tMaxs)
)
greyData = data.frame(
  xmin = 1:4 * blockSec - max(tMaxs), xmax = 1:4 * blockSec
)

data.frame(wtw = unlist(timeWTW_),
           time = rep(seq(0, blockSec * 4 - 1, by = 1), nSub),
           condition = rep(sumStats$condition, each = length(tGrid)),
           cbal = rep(cbal, each = length(tGrid))) %>%
  group_by(time, cbal) %>% 
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu)) +
  geom_rect(yellowData, mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 16), 
            fill = "#ffffcc",inherit.aes = F) +
  geom_rect(data = greyData, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 16),
            fill = "#d9d9d9", inherit.aes = F) +
  geom_ribbon(aes(ymin=min, ymax=max), fill = '#9ecae1') +
  geom_line(color = themeColor, size = 1) +
  xlab("Task time (s)") + ylab("WTW (s)") + 
  myTheme + ggtitle(expName) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  facet_grid(cbal~.) + scale_x_continuous(breaks = c(0, 600, 1200)) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 16)) +
  theme(
    strip.text.y = element_blank()
  ) + geom_text(data = data.frame(cbal = unique(cbal)),
                aes(label = cbal), x = 300, y = 5, inherit.aes = F, size = 6)
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 6)
ggsave("figures/MFPlot/wtw_timecourse.png", width = 5, height = 6) 


# plot average WTWs in two environments
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]
wTest1 = wilcox.test( sumStats[sumStats$condition == "HP" & sumStats$cbal == 1, "muWTW"],
             sumStats[sumStats$condition == "LP" & sumStats$cbal == 1, "muWTW"], paired = T)
wTest2 = wilcox.test( sumStats[sumStats$condition == "HP" & sumStats$cbal == 2, "muWTW"],
             sumStats[sumStats$condition == "LP" & sumStats$cbal == 2, "muWTW"], paired = T)
pData = data.frame(
  label = sprintf('p = %.3f', c(wTest1$p.value, wTest2$p.value)),
  cbal = c(1, 2)
)
data.frame(muWTWHP = sumStats$muWTW[sumStats$condition == 'HP'],
           muWTWLP = sumStats$muWTW[sumStats$condition == 'LP'],
           cbal = sumStats$cbal[sumStats$condition == "HP"]) %>%
  ggplot(aes(muWTWLP, muWTWHP)) +
  geom_point(color = themeColor, size = 5, shape = 21, fill = '#9ecae1', stroke =1) +
  geom_abline(slope = 1, intercept = 0) + facet_grid(~cbal) + 
  geom_text(data = pData, aes(label = label), x = 15, y = 3, inherit.aes = F) +
  xlab("LP muAUC / (s)") + ylab("HP muAUC / (s)") + 
    myTheme + xlim(c(-1,31)) + ylim(c(-1,31)) 
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 4, height = 3)
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 4, height = 3)


