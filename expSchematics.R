# this script plots delay distributions and reward rates in two environments

# load experiment parameters 
load("expParas.RData")

# for display purposes, all variables on the continous time scale
# are discretized into 0.1 second time bins
bin = 0.1 # width of a time bin
time = list(
  HP = seq(bin, tMaxs[1], by = bin),
  LP = seq(bin, tMaxs[2], by = bin)
) 

# delay CDFS
## early delay distribution: gamma(k = 2, theta = 2), truncated at 30
## late delay distribution: gamma(k = 6, theta = 2), truncated at 30
fastCDF = pgamma(time$HP,  2, scale = 2, log = FALSE) / pgamma(tMaxs[2],  2, scale = 2)
slowCDF = pgamma(time$LP,  6, scale = 2, log = FALSE)/ pgamma(tMaxs[2],  6, scale = 2)

# delay PDFs
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

# average waiting durations given different policies
# Here we assume rewards occur at the middle of each time bin
slowMeanDelay = cumsum((time$HP - 0.5 * bin) * slowPDF) / cumsum(slowPDF)
fastMeanDelay =  cumsum((time$LP - 0.5 * bin) * fastPDF) / cumsum(fastPDF)


# rewardRates given different policies
## might be different from the values used in expParas.R, 
## which are calcuated with a higher temporal resoluation
HPRate = (8 *  slowCDF + (-1) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + time$HP * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + time$HP * (1 - fastCDF)) / 2 + iti)
LPRate = ((-1) *  slowCDF + (8) * fastCDF) / 2 / 
  ((slowMeanDelay *  slowCDF + time$LP * (1 - slowCDF)) / 2 +
     (fastMeanDelay *  fastCDF + time$LP * (1 - fastCDF)) / 2 + iti)

# optimal raward rates and optimal policies
optimWaitThresholds = list()
optimWaitThresholds$HP = time$HP[which.max(HPRate)]
optimWaitThresholds$LP = time$LP[which.max(LPRate)]
optimRewardRates = list()
optimRewardRates$HP = max(HPRate)
optimRewardRates$LP = max(LPRate)

# plot cdf 
library('ggplot2')
library('tidyr'); library('dplyr')
source('subFxs/plotThemes.R')
dir.create("figures/expSchematics")
data.frame(cdf = c(0, fastCDF, 0, slowCDF),
           time = c(0, time$HP, 0, time$LP),
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
data.frame(rate = c(0, HPRate, 0, LPRate), 
           time = c(0, time$HP, 0, time$LP),
           cond = rep(c('HP', 'LP'), each = length(time$HP) + 1)) %>%
  ggplot(aes(time, rate)) + geom_line(size = 2, color = themeColor) + facet_grid(~cond) +
  myTheme + 
  ylab(expression(bold("Reward rate ( s"^"-1"*")"))) + xlab("Waiting policy (s)") +
  ggtitle(expName) + 
  theme(plot.title = element_text(hjust = 0.5, color = themeColor)) +
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)),
                     limits = c(0, max(tMaxs) * 1.1))
ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)


