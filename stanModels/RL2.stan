data {
  // experiment parameters
  real stepSec;// duration between two decision points
  real iti;// iti duration
  int  nDecPoint; // number of possible decision points
  real tWaits[nDecPoint]; // time for each decision point 
  
  // initial value for reward rate
  real rewardRate_ini; 
  
  // empirical data
  int N; // number of trials
  int Rs[N]; // payoff in each trial
  real Ts[N]; // a trial ends at t == T
  int lastDecPoints[N];// at which decision point each trial ends
}
transformed data {
  // total number of decision points in all trials
  int nDecPointTotal = sum(lastDecPoints);
}
parameters {
  // parameters:
  // alpha : learning rate
  // tau : action consistency, namely the soft-max temperature parameter
  // prior: prior belief parameter
  // beta : learning rate for the task reward rate
  
  // for computational efficiency,we sample raw parameters from unif(-0.5, 0.5)
  // which are later transformed into actual parameters
  real<lower = -0.5, upper = 0.5> raw_alphaR;
  real<lower = -0.5, upper = 0.5> raw_alphaU_alphaR;
  real<lower = -0.5, upper = 0.5> raw_tau;
  real<lower = -0.5, upper = 0.5> raw_prior;
  real<lower = -0.5, upper = 0.5> raw_beta_alphaR;
}
transformed parameters{
  // scale raw parameters into real parameters
  real alphaR = (raw_alphaR + 0.5) * 0.3; // alpha ~ unif(0, 0.3)
  real alphaU_alphaR = (raw_alphaU_alphaR + 0.5) * 5; // the ratio between alphaU and alphaR
  real alphaU =  min([alphaR * alphaU_alphaR, 1]');// alphaU
  real tau = (raw_tau + 0.5) * 21.9 + 0.1; // tau ~ unif(0.1, 22)
  real prior = (raw_prior + 0.5) * 6.5; // prior ~ unif(0, 6.5)
  real beta_alphaR = (raw_beta_alphaR + 0.5) * 1; // the ratio between beta and alphaR
  real beta = beta_alphaR * alphaR; // beta ~ unif(0, 0.3)
  
  // declare variables 
  // // state value of t = 0
  real V0; 
  // // action value of waiting in each decision points 
  vector[nDecPoint] Qwaits; 
  // // variables to record action values 
  matrix[nDecPoint, N] Qwaits_ = rep_matrix(0, nDecPoint, N);
  vector[N] V0_ = rep_vector(0, N);
  // // expected return 
  real G0;
  // // prediction error
  real delta0;
  // reward rate 
  real rewardRate;
  
  
  // initialize action values 
  //// the initial value of t = 0 
  V0 = 0; 
  rewardRate = rewardRate_ini;
  // the initial waiting value delines with elapsed time 
  // and the prior parameter determines at which step it falls below V0
  for(i in 1 : nDecPoint){
    Qwaits[i] = - i * 0.1 + prior + V0;
  }
  
  // record initial action values
  Qwaits_[,1] = Qwaits;
  V0_[1] = V0;
 
  //loop over trials
  for(tIdx in 1 : (N - 1)){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // determine alpha
    real alpha;
    if(R > 0){
      alpha = alphaR;
    }else{
      alpha = alphaU;
    }
    
    // update Qwaits towards the discounted returns
    for(i in 1 : lastDecPoint){
      real t = tWaits[i]; // time for this decision points 
      real Gt = R - rewardRate * (T - t) + V0;
      Qwaits[i] = Qwaits[i] + alpha * (Gt - Qwaits[i]);
    }
    
    // update V0 towards the discounted returns 
    G0 = R - rewardRate * T + V0;
    delta0 = G0 - V0;
    V0 = V0 + alpha * delta0;
    
    // update reward rate 
    rewardRate = rewardRate + beta * delta0;
    
    // save action values
    Qwaits_[,tIdx+1] = Qwaits;
    V0_[tIdx+1] = V0;
  }
}
model {
  // delcare variables 
  int action; 
  vector[2] actionValues; 
  // distributions for raw parameters
  raw_alphaR ~ uniform(-0.5, 0.5);
  raw_alphaU_alphaR ~ uniform(-0.5, 0.5);
  raw_tau ~ uniform(-0.5, 0.5);
  raw_prior ~ uniform(-0.5, 0.5);
  raw_beta_alphaR ~ uniform(-0.5, 0.5);
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      // the agent wait in every decision point in rewarded trials
      // and wait except for the last decision point in non-rewarded trials
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[i, tIdx] * tau;
      actionValues[2] = V0_[tIdx] * tau;
      target += categorical_logit_lpmf(action | actionValues);
    } 
  }
}
generated quantities {
 // initialize variables
  vector[2] actionValues;
  int action;
  vector[nDecPointTotal] log_lik = rep_vector(0, nDecPointTotal); // trial-wise log likelihood 
  real totalLL; // total log likelihood
  int no = 1; // action index
  
  // loop over trials
  for(tIdx in 1 : N){
    real T = Ts[tIdx]; // this trial ends on t = T
    int R = Rs[tIdx]; // payoff in this trial
    int lastDecPoint = lastDecPoints[tIdx]; // last decision point in this trial
    
    // loop over decision points
    for(i in 1 : lastDecPoint){
      if(R == 0 && i == lastDecPoint){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      // calculate the likelihood using the soft-max function
      actionValues[1] = Qwaits_[i, tIdx] * tau;
      actionValues[2] = V0_[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | actionValues);
      no = no + 1;
    }
  }
  // calculate total log likelihood
  totalLL =sum(log_lik);
}
