conditions = c("HP", "LP")
tMaxs = c(30, 30) # max trial durations in secs
nBlock = 2
blockMin = 10 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
iti = 2 # iti duration in secs
tokenValue = c(-1, 8)
# optimal waiting thresholds, unit = sec
optimWaitThresholds = list(HP = 30, LP = 5) # roughly estimations from Joe's grant. 
# optimal reward rates
HP = mean(tokenValue) / (30 + iti) # this is obviously wrong. We know what 
LP = mean(tokenValue) / (5 + iti)
optimRewardRates = list(HP = HP, LP = LP) 
# analyses parameters
tGrid = seq(0, blockSec, by = 2) # time grid for wtw time courses
kmGrid = seq(0, min(tMaxs), by = 0.1) # time grid for Kaplan-Meier survival curves
save("conditions" = conditions,
     "tMaxs" = tMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "iti" = iti,
     "tokenValue" = tokenValue,
     "optimRewardRates" = optimRewardRates,
     "optimWaitThresholds" = optimWaitThresholds,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")