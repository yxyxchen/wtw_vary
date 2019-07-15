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

