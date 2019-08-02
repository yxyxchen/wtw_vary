data {
  // depending on the task
  real stepDuration;
  real iti;
  
  // depending on the condition
  real wIni;
  int nTimeSteps; // nTimeSteps = tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] trialEarnings;
  int Ts[N]; // terminal time step index 
  
  // debug
  int nPara;
  vector[nPara] low;
  vector[nPara] up;
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = low[1], upper = up[1]> phi;
  real<lower = low[2], upper = up[2]> nega; 
  real<lower = low[3], upper = up[3]> tau;
  real<lower = low[4], upper = up[4]> gamma;
  real<lower = low[5], upper = up[5]> prior; 
}
transformed parameters{
  // initialize action values 
  real Viti = wIni;
  vector[nTimeSteps] Qwait;
  
    // initialize variables to record action values 
  matrix[nTimeSteps, N] Qwaits = rep_matrix(0, nTimeSteps, N);
  vector[N] Vitis = rep_vector(0, N);

  // initialize caching variables
  real G1;
  
  // enter the initial values
  for(i in 1 : nTimeSteps){
    Qwait[i] = prior*0.1 - 0.1*(i - 1) + Viti;
  }
  Qwaits[,1] = Qwait;
  Vitis[1] = Viti;
 
  //loop over trials
  for(tIdx in 1 : (N -1)){
    // determine the termial timestep T 
    int T = Ts[tIdx];
    real RT = trialEarnings[tIdx];
    
    // update action values for rewarded trials
    if(RT != 0){
      for(t in 1 : (T - 1)){
        real G = RT * gamma^(T - t -1) + Viti * gamma^(T - t);
        if(RT > 0){
          Qwait[t] = Qwait[t] + phi * (G - Qwait[t]);
        }else{
          Qwait[t] = Qwait[t] + nega * phi * (G - Qwait[t]);
        }
      }
    }else{
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  RT  * gamma^(T - t -1) + Viti * gamma^(T - t);
          Qwait[t] = Qwait[t] + nega * phi * (G - Qwait[t]);    
        }
      }
    }
    // update Qquit by counterfactual thiking
    G1 =  RT  * gamma^(T - 2) + Viti * gamma^(T - 1);
    // update Viti
    if(RT > 0){
       Viti = Viti + phi * (G1 * gamma^(iti / stepDuration) - Viti);
    }else{
       Viti = Viti + nega * phi * (G1 * gamma^(iti / stepDuration) - Viti);
    }
   
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Vitis[tIdx+1] = Viti;
  }// end of the loop
}
model {
  phi ~ uniform(low[1], up[1]);
  nega ~ uniform(low[2], up[2]);
  tau ~ uniform(low[3], up[3]);
  gamma ~ uniform(low[4], up[4]);
  prior ~ uniform(low[5], up[5]);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Vitis[tIdx] * tau;
      //action ~ categorical_logit(values);
      target += categorical_logit_lpmf(action | values);
    } 
  }
}
generated quantities {
// initialize log_lik
  vector[totalSteps] log_lik = rep_vector(0, totalSteps);
  vector[2] values;
  real LL_all;
  int no = 1;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    int T = Ts[tIdx];
    for(i in 1 : (T - 1)){
      if(trialEarnings[tIdx] == 0 && i == (T-1)){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Vitis[tIdx] * tau;
      log_lik[no] =categorical_logit_lpmf(action | values);
      no = no + 1;
    }
  }// end of the loop
  LL_all =sum(log_lik);
}
