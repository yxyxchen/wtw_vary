
library('ggplot2')
library('tidyr'); library('dplyr')
source('subFxs/plotThemes.R')
dir.create("figures/expSchematics")

load("expParas.RData")
# plot exp
# cdf from 1 
x = seq(0.5, tMaxs[2], by = 0.5);
fastCDF = pgamma(x,  2, scale = 2, log = FALSE) / pgamma(tMaxs[2],  2, scale = 2)
slowCDF = pgamma(x,  6, scale = 2, log = FALSE)/ pgamma(tMaxs[2],  6, scale = 2)
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

data.frame(cdf = c(0, fastCDF, 0, slowCDF),
           time = rep(c(0, x), 2),
           cond = rep(c('Early', 'Late'), each = length(fastCDF) + 1)) %>%
  ggplot(aes(time, cdf)) + geom_line(color = themeColor, size = 3) + facet_grid(~cond) + 
  myTheme + xlab('Delay duration (s)') + ylab('CDF') + ggtitle(expName) + 
  theme(plot.title = element_text(hjust = 0.5, color = themeColor)) +
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)),limits = c(0, max(tMaxs) * 1.1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  geom_text(data = data.frame(label = c("HP = -1¢", "HP = +8¢"), cond = c("Early", "Late"), x = c(16, 22)),
            aes(label = label, x = x), y = 0.4, size = 5) + 
  geom_text(data = data.frame(label = c("LP = +8¢", "LP = -1¢"), cond = c("Early", "Late"), x = c(16, 22)),
            aes(label = label, x = x), y = 0.2, size = 5)   
ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
ggsave('figures/expSchematics/CDF.png', width =4, height = 3)


# plot reward rates
slowMeanDelay = cumsum((x - 0.5 * x[1]) * slowPDF) / cumsum(slowPDF)
fastMeanDelay =  cumsum((x - 0.5 * x[1]) * fastPDF) / cumsum(fastPDF)

policy = data.frame(cond = c("HP", "LP"), rewardRate = c(x[which.max(HPRate)], x[which.max(LPRate)]))
HPRate = (8 *  slowCDF + (-1) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + x * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + x * (1 - fastCDF)) / 2 + iti)

LPRate = ((-1) *  slowCDF + (8) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + x * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + x * (1 - fastCDF)) / 2 + iti)

data.frame(rate = c(0, HPRate, 0, LPRate), 
           time = rep(c(0,x), 2),
           cond = rep(c('HP', 'LP'), each = length(x) + 1)) %>%
  ggplot(aes(time, rate)) + geom_line(size = 2) + facet_grid(~cond) +
  myTheme + 
  ylab(expression(bold("Reward rate (cent s"^"-1"*")"))) + xlab("Waiting policy (s)")  +
  geom_vline(data = policy, aes(xintercept = rewardRate),
             linetype = "dashed", size = 1.5) 
ggsave("figures/exp/reward_rate.png", width = 6, height = 3)