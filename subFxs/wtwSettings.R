# 
######## condition varibles #########
conditions = c("Rising", "Falling")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(30, 30) # trial durations
blockMins = 10 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 1)
kmGrid = seq(0, min(tMaxs), 0.1)
stepDuration = 1

######### reward variable ########
tokenValue = c(-1, 8) #value of the token
loseValue = 0

########## waiting time and reward rate vairbales ########
# in this task, two truncated gamma distributions were equally sampled
# the gamma document in Matlab is a little bit strange, where a is the shape, namely k
# while b is not the rate but the scale, which is usually presented by theta.
# (fast:k = 2, theta = 2; slow: a= 6, theta = 2). btw, mu_fast = 4, mu_slow = 12
# To sample the distribution evenly, they first divided each distribution by its quartiles into 4 components 
# and sampled equally from 4 components.

optimWaitTimes = list(HP = 30, LP = 5) # roughly estimations from Joe's grant. 
HP = mean(tokenValue) / (30 + iti)
LP = mean(tokenValue) / (5 + iti)
optimRewardRates = list(HP = HP, LP = LP) 
save("conditions", "conditionColors", "tMaxs", "blockMins", "blockSecs", "iti", "tGrid", 
     "tokenValue", "stepDuration", "optimRewardRates", 
     "optimWaitTimes", "loseValue", "kmGrid", file = "wtwSettings.RData")

# plot exp
# cdf from 1 
x = seq(0.5, tMaxs[2], by = 0.5);
fastCDF = pgamma(x,  2, scale = 2, log = FALSE) / pgamma(tMaxs[2],  2, scale = 2)
slowCDF = pgamma(x,  6, scale = 2, log = FALSE)/ pgamma(tMaxs[2],  6, scale = 2)
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

library('ggplot2')
library('tidyr'); library('dplyr')
source('subFxs/plotThemes.R')
data.frame(pdf = c(0, fastPDF, 0, slowPDF, 0, slowPDF, 0, fastPDF), time = rep(c(0, x), 4),
           cond = rep(c('HP', 'LP'), each = 2 * length(x) + 2),
           reward = factor(rep(c(tokenValue, tokenValue), each = length(x) + 1))) %>%
ggplot(aes(time, pdf, color = reward, fill = reward)) + facet_grid(~cond) + geom_line() + 
  geom_ribbon(aes(ymin=0, ymax=pdf), alpha = 0.5) + 
  scale_fill_manual(values = c("#969696", "#fed976")) +
  scale_color_manual(values = c("#969696", "#fed976")) + myTheme + 
  xlab('Delay duration(s)') + ylab('PDF') 
dir.create('figures/exp')
ggsave('figures/exp/cdp.png', width =6, height = 3)


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