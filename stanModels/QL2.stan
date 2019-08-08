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
}
transformed data {
  int totalSteps = sum(Ts) - N;
}
parameters {
  real<lower = 0, upper = 0.3> phi;
  real<lower = 0, upper = min([5, 1 / phi]')> nega; 
  real<lower = 0.1, upper = 22> tau;
  real<lower = 0.5, upper = 1> gamma;
  real<lower = 0, upper = 65> prior; 
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
          Qwait[t] = Qwait[t] + phi * nega * (G - Qwait[t]);
        }
      }
    }else{
      if(T > 2){
        for(t in 1 : (T-2)){
          real G =  RT  * gamma^(T - t -1) + Viti * gamma^(T - t);
          Qwait[t] = Qwait[t] + phi * nega * (G - Qwait[t]);    
        }
      }
    }
    // update Qquit by counterfactual thiking
    G1 =  RT  * gamma^(T - 2) + Viti * gamma^(T - 1);
    // update Viti
    if(RT > 0){
       Viti = Viti + phi * (G1 * gamma^(iti / stepDuration) - Viti);
    }else{
       Viti = Viti + phi * nega * (G1 * gamma^(iti / stepDuration) - Viti);
    }
   
    // save action values
    Qwaits[,tIdx+1] = Qwait;
    Vitis[tIdx+1] = Viti;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  nega ~ uniform(0, min([5, 1 / phi]'));
  tau ~ uniform(0.1, 22);
  gamma ~ uniform(0.5, 1);
  prior ~ uniform(0, 65);
  
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
