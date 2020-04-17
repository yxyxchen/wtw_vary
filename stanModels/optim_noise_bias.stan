data {
  // experiment parameters
  real stepSec;// duration between two decision points
  real iti;// iti duration
  int  nDecPoint; // number of possible decision points
  real tWaits[nDecPoint]; // time for each decision point 
  
  // empirical data
  int N; // number of trials
  int Rs[N]; // payoff in each trial
  real Ts[N]; // a trial ends at t == T
  int lastDecPoints[N];// at which decision point each trial ends
  int cIdxs[N];
  
  // theoretical values
  real HPvalues[nDecPoint];
  real LPvalues[nDecPoint];
}
transformed data {
  // total number of decision points in all trials
  int nDecPointTotal = sum(lastDecPoints);
  
    ## the theoretic present value of the awaited reward sampled at 1 hz
  HPtheoreticValues = c(presentValues$HP[seq(0, delayMaxs[1], by = 0.1) %in% seq(0, delayMax, by = stepSec)],
                        rep(0, (nDecPoint -delayMaxs[1]) / stepSec - 1))
  LPtheoreticValues = presentValues$LP[seq(0, delayMaxs[2], by = 0.1) %in% seq(0, delayMax, by = stepSec)]
}
parameters {
  // parameters:
  // tau: action consistency 
  // theta : decision bias
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_theta;

}
transformed parameters{
  // transfer paras
  real tau = (raw_tau + 0.5) * 21.9 + 0.1 ; // tau ~ unif(0.1, 22)
  real theta = (raw_theta + 0.5) * 40 - 20; // tau ~ unif(-20, 20)
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues;
  real values[nDecPoint];
  actionValues[2] = 0; // decision threshold 
  
  
  // sample
  raw_tau ~ uniform(-0.5, 0.5);
  raw_theta ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // determine decision values 
    if(cIdxs[tIdx] == 1){
      values = HPvalues;
    }else{
      values = LPvalues;
    }
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      // the agent wait in every decision point in rewarded trials
      // and wait except for the last decision point in non-rewarded trials
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood 
      actionValues[1] = values[i] * tau + theta;
      target +=  categorical_logit_lpmf(action | actionValues);
    } 
  }
}
generated quantities {
 // initialize variables
  int action;
  vector[2] actionValues;
  real values[nDecPoint];
  vector[nDecPointTotal] log_lik = rep_vector(0, nDecPointTotal); // trial-wise log likelihood 
  real totalLL; // total log likelihood
  int no = 1; // action index
  
  actionValues[2] = 0; // decision threshold 
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial

    // determine decision values 
    if(cIdxs[tIdx] == 1){
      values = HPvalues;
    }else{
      values = LPvalues;
    }
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood 
      actionValues[1] = values[i] * tau + theta;
      log_lik[no] = categorical_logit_lpmf(action | actionValues);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  totalLL =sum(log_lik);
}



