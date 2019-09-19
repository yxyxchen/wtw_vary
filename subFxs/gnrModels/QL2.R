# our reinfocement learning generative models simulate persistence behavior as wait-or-quit choices 
# on discrete time steps:
# QL1: Q-learning model with a single learning rate
# QL2: Q-learning model with separate learning rates for rewards and punishments
# RL1: R-learning model with a single learning rate 
# RL2: R-learning model with separate learning rates for rewards and punishments

# inputs:
# paras: parameter values 
# condition_: [nTrialx1 factor] HP or LP
# scheduledWait_ : [nTrialx1 num] trial-wise reward delays

# outputs
# trialNum : [nTrialx1 int] 1 : nTrial
# condition : [nTrialx1 factor] from inputs
# scheduledWait : [nTrialx1 num] from inputs 
# trialEarnings : [nTrialx1 int] payment for each trial, either 10 or 0
# timeWaited : [nTrialx1 num] waiting duration for each trial 
# sellTime : [nTrialx1 num]  when the agent sells then token on each trial 
# Qwaits_ : [nStepMax x nTrial num] action value of waiting for each time step at each trial
# Viti_ : [nTrialx1 num] state value of the iti stage at each trial 

QL2 = function(paras, condition_, scheduledWait_, scheduledReward_){
  load('expParas.RData')
  
  # extract parameter values
  phi_pos = paras[1]; phi_neg = paras[2]; tau = paras[3]; gamma = paras[4]; prior = paras[5]
  
  # num of trials
  nTrial = length(scheduledWait_) # num of trials 
  
  # parameters for discrete temproal states 
  stepSec = 1 # duration of one time step 
  nStepMax =  max(tMaxs) / stepSec # maximal number of steps in a trial
  
  # initialize action values 
  Viti = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
  Qwaits = (prior  - 1 : nStepMax) * 0.1 + Viti
  
  # initialize output variables
  Qwaits_ = matrix(NA, nStepMax, nTrial); Qwaits_[,1] = Qwaits 
  Viti_ = vector(length = nTrial); Viti_[1] = Viti
  trialEarnings_ = rep(0, nTrial)
  timeWaited_ = rep(0, nTrial)
  sellTime_ = rep(0, nTrial)
  
  # track elpased time from the beginning of the task 
  elapsedTime = 0 
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # current scheduledWait 
    scheduledWait = scheduledWait_[tIdx]
    scheduledReward = scheduledReward_[tIdx]
    
    # loop over steps until a trial ends 
    t = 1
    while(t <= nStepMax){
      # decide whether to wait or quit
      pWait =  1 / sum(1  + exp((Viti - Qwaits[t])* tau))
      action = ifelse(runif(1) < pWait, 'wait', 'quit')
      
      # if a reward occurs and the agent is still waiting, the agent gets the reward
      rewardOccur = scheduledWait <= (t * stepSec) && scheduledWait > ((t-1) * stepSec)
      getReward = (action == 'wait' && rewardOccur)
      
      # a trial ends if the agent gets a reward or quits. if the trial ends, 
      # proceed to the next trial.Otherwise, proceed to the next step 
      isTerminal = (getReward || action == "quit")
      if(isTerminal){
        # update trial-wise variables 
        T = t + 1 # terminal state
        trialEarnings = ifelse(getReward, scheduledReward, 0) # payment 
        timeWaited =  ifelse(getReward, scheduledWait, t * stepSec) # waiting duration 
        sellTime = elapsedTime + timeWaited # elapsed task time when the agent sells the token
        elapsedTime = elapsedTime + timeWaited + iti # elapsed task time before the next trial
        # record trial-wise variables
        trialEarnings_[tIdx] = trialEarnings
        timeWaited_[tIdx] = timeWaited
        sellTime_[tIdx] = sellTime
        break
      }else{
        t = t + 1
      }
    }
    
    # update action values at the end of each trial
    if(tIdx < nTrial){
      # calculate the reward signal for updating action value which equals 
      # trialEarnings + discounted value of the successor state gamma. Noticably,
      # the successor state at the end of trial is always the iti state before the next trial
      rwdSignal = trialEarnings + Viti * gamma
      # discounted reward signals for step 1 - (T-1)
      stepRwdSignals = sapply(1 : (T-1), function(t) gamma^(T-t-1) * rwdSignal)
      # discounted reward signals for the iti state
      itiRwdSignal = rwdSignal * gamma^(T-2 + iti / stepSec)
      
      # update Qwaits 
      if(getReward){
        if(trialEarnings > 0){
          Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi_pos *(stepRwdSignals[1 : (T-1)] - Qwaits[1 : (T-1)])
        }else{
          Qwaits[1 : (T-1)] = Qwaits[1 : (T-1)] + phi_neg *(stepRwdSignals[1 : (T-1)] - Qwaits[1 : (T-1)])
        }
 
      }else{
        if(T > 2){
          # in non-rewarded trials, Qwait in the last step will not be updated
          # since the agent chooses to quit on that step 
          Qwaits[1 : (T-2)] = Qwaits[1 : (T-2)] + phi_neg * (stepRwdSignals[1 : (T-2)] - Qwaits[1 : (T-2)])
        }
      }
      
      # update Viti
      itiDelta = itiRwdSignal - Viti
      Viti = ifelse(getReward, Viti + phi_pos * itiDelta, Viti + phi_neg * itiDelta)
      
      # record updated values
      Qwaits_[,tIdx + 1] = Qwaits
      Viti_[tIdx + 1] = Viti
    }# end of the value update section
  } # end of the loop over trials
  
  # return outputs
  outputs = list( 
    "trialNum" = 1 : nTrial, 
    "condition" = condition_,
    "trialEarnings" = trialEarnings_, 
    "timeWaited" = timeWaited_,
    "sellTime" = sellTime_,
    "scheduledWait" = scheduledWait_,
    "Qwaits_" = Qwaits_, 
    "Viti_" = Viti_
  )
  return(outputs)
}