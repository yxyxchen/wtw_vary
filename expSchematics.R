
library('ggplot2')
library('tidyr'); library('dplyr')
source('subFxs/plotThemes.R')

load("expParas.RData")
# plot exp
# cdf from 1 
x = seq(0.5, tMaxs[2], by = 0.5);
fastCDF = pgamma(x,  2, scale = 2, log = FALSE) / pgamma(tMaxs[2],  2, scale = 2)
slowCDF = pgamma(x,  6, scale = 2, log = FALSE)/ pgamma(tMaxs[2],  6, scale = 2)
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

data.frame(pdf = c(0, fastPDF, 0, slowPDF, 0, slowPDF, 0, fastPDF), time = rep(c(0, x), 4),
           cond = rep(c('HP', 'LP'), each = 2 * length(x) + 2),
           reward = factor(rep(c(tokenValue, tokenValue), each = length(x) + 1))) %>%
ggplot(aes(time, pdf, fill = reward)) + facet_grid(~cond) + geom_line(color = NA) + 
  geom_ribbon(aes(ymin=0, ymax=pdf)) + 
  scale_fill_manual(values = c("#565656", "#9ecae1")) + myTheme + 
  xlab('Delay duration(s)') + ylab('PDF') 
dir.create('figures/expSchematics')





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