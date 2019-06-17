# in this script, we try to calculate the optimWaitingTimes and the optimRewardRates
# to be as close as to the normal analysis (like integrate the prob density)
# we assume rewards happen at the middle of the gap(therefore, the meanRewardDelay would be unbiased)
# yet in wtwSettingsEnd.R, to unfiy different algorithms, we assume rewards happen at the end of the gap
# however, results for LP still change with the stepDuration
# we do all the calculation by stepDuration = 0.5, and the optimWaitTime, you know is not that...

# we use this script to get stepDuration
# we don't use the reward rate here, it is close to the normal analysis, but not that good.
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

######### reward variable ########
tokenValue = 10 #value of the token
loseValue = 0
stepDuration = 1
########## supporting vairbales ########
optimWaitTimes = list(HP = 30, LP = 5)
HP = tokenValue / ((2 + 16) / 2 + iti)
LP = tokenValue / ((2 + 16) / 2 + iti)
optimRewardRates = list(HP = HP, LP = LP) 

save("conditions", "conditionColors", "tMaxs", "blockMins", "blockSecs", "iti", "tGrid", 
     "tokenValue", "stepDuration", "optimRewardRates", 
     "optimWaitTimes", "loseValue", "kmGrid", file = "wtwSettings.RData")

