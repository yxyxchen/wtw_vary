# this script plots delay distributions, reward rates and subjective values in two environments
expSchematics = function(smallReward, iti, isPlot){
# small reward 
smallReward = 0 # get 0 when the agent quits
  
# load experiment parameters 
load("expParas.RData")

# for display purposes, all variables on the continous time scale
# are discretized into 0.1 second time bins
bin = 0.1 # width of a time bin
time = list(
  HP = seq(bin, delayMaxs[1], by = bin),
  LP = seq(bin, delayMaxs[2], by = bin)
) 

# delay CDFS
## early delay distribution: gamma(k = 2, theta = 2), truncated at 30
## late delay distribution: gamma(k = 6, theta = 2), truncated at 30
fastCDF = pgamma(time$HP,  2, scale = 2, log = FALSE) 
fastCDF[length(fastCDF)] = 1 
slowCDF = pgamma(time$LP,  6, scale = 2, log = FALSE)
slowCDF[length(slowCDF)] = 1 


# delay PDFs
fastPDF = diff(c(0,fastCDF))
slowPDF = diff(c(0,slowCDF))

# average waiting durations given different policies
# Here we assume rewards occur at the middle of each time bin
slowMeanDelay = cumsum((time$HP - 0.5 * bin) * slowPDF) / cumsum(slowPDF)
fastMeanDelay =  cumsum((time$LP - 0.5 * bin) * fastPDF) / cumsum(fastPDF)


# rewardRates given different policies
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

# calculate the value of waiting as a function of elapsed time 
subjectValues = list()
cIdx = 1
i = 1
for(cIdx in 1 : 2){
  condition = conditions[cIdx]
  delayMax = delayMaxs[cIdx]
  thisTime = time[[cIdx]]
  Rstar = optimRewardRates[[cIdx]]
  ts = seq(0, delayMax, by = 0.1) # elapsed time
  
  # initialize 
  thisSubjectValues = vector(length = length(ts)) # value of waiting given the elapsed time value
  Tstars = vector(length = length(ts)) # optimal waiting policy given the elapsed time value
  
  # loop over different elapsed time values
  for(i in 1 : length(ts)){
    t = ts[i] # this elapsed time
    trctTime = thisTime[thisTime > t]
    
    # loop over different waiting policies to find Tstar and gt_max
    if(t == delayMax){
      Tstar = t
      gt_max = max((ifelse(condition == "HP", -1, 8) * tail(fastPDF, 1) +
                  ifelse(condition == "HP", 8, -1)  * tail(slowPDF, 1)) / (tail(fastPDF, 1) + tail(slowPDF, 1)), 0)
    }else{
      Tstar = t
      gt_max = -100
      for(T in seq(t, delayMax, by = 0.1)){
        # normalize the truncated distributions
        trctFastPDF = fastPDF[thisTime > t] / sum(fastPDF[thisTime > t] + slowPDF[thisTime > t])
        trctSlowPDF = slowPDF[thisTime > t] / sum(fastPDF[thisTime > t] + slowPDF[thisTime > t])
        # expected reward
        at = ifelse(condition == "HP", -1, 8) * sum(trctFastPDF[trctTime <= T]) + 
          ifelse(condition == "HP", 8, -1) * sum(trctSlowPDF[trctTime <= T])
        # remaining delay 
        bt = sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctFastPDF[trctTime <= T]) + 
          sum((trctTime[trctTime <= T] - 0.5 * bin - t) * trctSlowPDF[trctTime <= T])+
          (T - t) * sum(trctFastPDF[trctTime > T] + trctSlowPDF[trctTime > T]) 
        # net value
        gt = at - bt * Rstar
        if(gt > gt_max){
          gt_max = gt
          Tstar = T
        }
      }
    }
    thisSubjectValues[i] = gt_max 
    Tstars[i] = Tstar
  }
  subjectValues[[condition]] = thisSubjectValues
}

if(isPlot){
  # plot cdf 
  library("tidyverse")
  library("latex2exp")
  source('subFxs/plotThemes.R')
  dir.create("figures/expSchematics")
  data.frame(cdf = c(0, fastCDF, 0, slowCDF),
             time = c(0, time$HP, 0, time$LP),
             cond = rep(c('Early', 'Late'), each = length(fastCDF) + 1)) %>%
    ggplot(aes(time, cdf)) + geom_line(size = 3) + facet_grid(~cond) + 
    myTheme + xlab('Delay duration (s)') + ylab('CDF') + 
    theme(plot.title = element_text(hjust = 0.5, color = themeColor)) +
    scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),limits = c(0, max(delayMaxs) * 1.1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) + 
    geom_text(data = data.frame(label = c("HP = -1¢", "HP = +8¢"), cond = c("Early", "Late"), x = c(16, 22)),
              aes(label = label, x = x), y = 0.4, size = 5, color = conditionColors[1]) + 
    geom_text(data = data.frame(label = c("LP = +8¢", "LP = -1¢"), cond = c("Early", "Late"), x = c(16, 22)),
              aes(label = label, x = x), y = 0.2, size = 5, color = conditionColors[2])   
  ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
  ggsave('figures/expSchematics/CDF.png', width =4, height = 3)
  
  # plot reward rates
  data.frame(rate = c(0, HPRate, 0, LPRate), 
             time = c(0, time$HP, 0, time$LP),
             condition = rep(c('HP', 'LP'), each = length(time$HP) + 1)) %>%
    ggplot(aes(time, rate)) + geom_line(size = 2, aes(color = condition)) + facet_grid(~cond) +
    myTheme + 
    ylab(TeX('Reward rate (¢ $s^{-1}$)')) + xlab("Waiting policy (s)") + 
    theme(plot.title = element_text(hjust = 0.5, color = themeColor)) +
    scale_x_continuous(breaks = c(0, max(delayMaxs)/ 2, max(delayMaxs)),
                       limits = c(0, max(delayMaxs) * 1.1)) +
    scale_color_manual(values = conditionColors) +
    theme(legend.position = "none") + facet_grid(~condition)
  ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
  ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)
  
  # plot subjective value of waiting 
  data.frame(
    value =  c(subjectValues$HP, subjectValues$LP),
    t = c(seq(0, delayMaxs[1], by = 0.1), seq(0, delayMaxs[2], by = 0.1)),
    condition = rep(conditions, c(length(subjectValues$HP), length(subjectValues$LP))))%>%
    ggplot(aes(t, value)) +
    geom_line(aes(color = condition), size = 2) +
    myTheme+
    scale_color_manual(values = conditionColors) +
    scale_linetype_manual(values = c(1, 2)) +
    xlab("Elapsed time (s)") + ylab("Subjective value (¢)")  + 
    theme(legend.position = "none") + facet_grid(~condition)
  ggsave("figures/expSchematics/subjective.eps", width = 4, height = 3)
  ggsave("figures/expSchematics/subjective.png", width = 4, height = 3)      
}    


# return outputs 
outputs = list(
  "optimWaitThresholds" = optimWaitThresholds,
  "optimRewardRates" = optimRewardRates,
  "time" = time,
  "subjectValues" = subjectValues,
  "fastPDF" = fastPDF,
  "slowPDF" = slowPDF)
}