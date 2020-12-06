// Project Name: YUKON RIVER CHINOOK IMEG WORKING GROUP - STAN version of age-structured state-space spawner-recruitment model with AR-1 process variation (Fleischman et al. CJFAS. 2013)
// Creator: Brendan Connors, Fisheries and Oceans Canada and Curry Cunningham, College of Fisheries and Ocean Sciences, UAF
// Date: 3.11.2020
// Purpose: State-space spawner recruitment model

data{
  int nyrs;  // number of calender years
  int a_min;  // minimum age class
  int a_max;  // maximum age class
  int A; // number of age classes
  int nRyrs; // number of recruitment years
  int H_comps[nyrs, A]; // catch observed age composition in counts by age class
  int S_comps[nyrs, A]; // spawners observed age composition in counts by age class
  vector[nyrs] S_obs; // observed spawners
  vector[nyrs] H_obs; // observed harvest
  vector[nyrs] S_cv; // spawner observation error CV
  vector[nyrs] H_cv; // harvest observation error CV
}

parameters{
  vector<lower=0>[nRyrs] lnR; // log recruitment states
  real<lower=0> lnalpha; // log Ricker a
  real<lower=0> beta; // Ricker b
  real<lower=0> sigma_R; // process error
  real<lower=0> sigma_R0; // process error for first a.max years with no spawner link
  real<lower=-1,upper=1> phi; // lag-1 correlation in process error
  real lnresid_0;  // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0; // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[nyrs] U;  // harvest rate
  vector<lower=0,upper=1>[3] prob;  // maturity schedule probs
  real<lower=0,upper=1> D_scale; // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[nRyrs, A] g; // individual year/age class gamma variates for generating age at maturity proportions
  
  //BC new
  vector<lower=0>[A] Sel; // fishery selectivity by age
}

transformed parameters{
  vector<lower=0>[nyrs] N;  // run size states
  vector<lower=0>[nyrs] S;  // spawner states
  vector[nyrs] lnS;  // log spawner states
  vector<lower=0>[nyrs] C;  // catch states
  vector[nyrs] lnC;  // log catch states
  vector<lower=0>[nRyrs] R;  // recruitment states
  vector[nRyrs] lnresid;  // log recruitment residuals
  vector[nRyrs] lnRm_1;  // log recruitment states in absence of lag-one correlation
  vector[nRyrs] lnRm_2;  // log recruitment states after accounting for lag-one correlation
  matrix<lower=0>[nyrs, A] N_ta; // returns by age matrix
  matrix<lower=0, upper=1>[nRyrs, A] p; // age at maturity proportions
  vector<lower=0,upper=1>[4] pi;  // maturity schedule probs
  real<lower=0> D_sum; // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[A] Dir_alpha; // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions

  // BC new
  matrix<lower=0, upper=1>[nyrs, A] C_pa; // age composition in catch by year/age class 
  matrix<lower=0, upper=1>[nyrs, A] S_pa; // age composition in spawners by year/age class 
  matrix<lower=0>[nyrs, A] C_ta; // catch by age matrix
  matrix<lower=0>[nyrs, A] S_ta; // spawners by age matrix
  
  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = prob[3] * (1 - pi[1] - pi[2]);
  pi[4] = 1 - pi[1] - pi[2] - pi[3];
  D_sum = 1/D_scale^2;

  for (a in 1:A) {
    Dir_alpha[a] = D_sum * pi[a];
    for (y in 1:nRyrs) {
      p[y,a] = g[y,a]/sum(g[y,]);
    }
  }

  R = exp(lnR);

  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  // Calculate catch and spawners by age and return year 
  for (t in 1:nyrs) {
    for(a in 1:A){
      N_ta[t,a] = R[t+A-a] * p[t+A-a,a];
      C_ta[t,a] = Sel[a] * N_ta[t,a] * U[t];
      S_ta[t,a] = N_ta[t,a] - C_ta[t,a];
    }
  }

  // Calculate returns, spawners and catch by return year
  for(t in 1:nyrs) {
    N[t] = sum(N_ta[t,1:A]);
    S[t] = sum(S_ta[t,1:A]);
    lnS[t] = log(S[t]);
    C[t] = sum(C_ta[t,1:A]);
    lnC[t] = log(C[t]);
  }

  // Calculate age proportions in harvest and spawners by return year
  for (t in 1:nyrs) {
    for(a in 1:A){
      C_pa[t,a] = C_ta[t,a]/C[t];
      S_pa[t,a] = S_ta[t,a]/S[t];
    }
  }

  // Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
  for (i in 1:nRyrs) {
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
    lnRm_2[i] = 0.0;
  }

  for (y in (A+a_min):nRyrs) {
    lnRm_1[y] = lnS[y-a_max] + lnalpha - beta * S[y-a_max];
    lnresid[y] = lnR[y] - lnRm_1[y];
  }

  lnRm_2[A+a_min] =  lnRm_1[A+a_min] + phi * lnresid_0;

  for (y in (A+a_min+1):nRyrs) {
    lnRm_2[y] =  lnRm_1[y] + phi * lnresid[y-1];
  }
}

model{
  // Priors
  lnalpha ~ normal(0,3);
  beta ~ normal(0.1638,0.025); // PR based prior, inverse of PR based Smax in units of 100,000s, so if Smax = 6.10354 then mean of prior is 1/6.10354 = 0.1638 hist(1/rnorm(1000,0.1638,0.02))
  sigma_R ~ normal(0,2);
  lnresid_0 ~ normal(0,20);
  mean_ln_R0 ~ normal(0,20);
  sigma_R0 ~ inv_gamma(1.25,0.4); // made this an informative prior based on metanalysis of other AK chinook stocks (Fleischman et al. 2013), less informative priors resulted in divergent tranistions
  prob[1] ~ beta(1,1);
  prob[2] ~ beta(1,1);
  prob[3] ~ beta(1,1);
  D_scale ~ beta(1,1);
  Sel[1] ~ normal(0,2);
  Sel[2] ~ normal(0,2);
  Sel[3] ~ normal(0,2);
  Sel[4] ~ normal(0,2);


  // Gamma variates for each year and age class which are used to determine age at maturity proportions
  for (y in 1:nRyrs) {
    for (a in 1:A) {
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y,a]|Dir_alpha[a],1);
    }
  }

  // First `a.max` years of recruits, for which there is no spawner link
  lnR[1:a_max] ~ normal(mean_ln_R0, sigma_R0);

  // State model
  lnR[(A+a_min):nRyrs] ~ normal(lnRm_2[(A+a_min):nRyrs], sigma_R);

  // Observation model
  for(t in 1:nyrs){
  //C_comps[t,1:A]) ~ multinomial(C_pa[t,1:A]);
  //S_comps[t,1:A]) ~ multinomial(S_pa[t,1:A]);
    target += multinomial_lpmf(H_comps[t,1:A]|to_vector(C_pa[t,1:A]));
    target += multinomial_lpmf(S_comps[t,1:A]|to_vector(S_pa[t,1:A]));
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}

generated quantities {
  // biological reference points
  real<lower=0> S_max; // spawner abundance that produces maximum recruitment
  real<lower=0> S_eq; // equilibrium spawner abundance
  real<lower=0> S_msy; // spawner abundance that produces maximum sustainable yield
  real<lower=0, upper=1> U_msy; // harvest rate that maximizes yield
  vector[nRyrs] lnalpha_y;
  real<lower=0> lnalpha_c; // bias corrected log alpha
  real<lower=0> S_eq_c; // bias corrected equilibrium spawner abundance
  real S_msy_c; // bias corrected spawner abundance that produces maximum sustainable yield
  real U_msy_c; // bias corrected harvest rate that maximizes yield

  S_max = 1/beta;
  S_eq = lnalpha * S_max;
  S_msy = S_eq * (0.5 - 0.07 * lnalpha);
  U_msy = lnalpha * (0.5 - 0.07 * lnalpha);
  lnalpha_y = lnalpha +lnresid;
  lnalpha_c = lnalpha + (sigma_R * sigma_R)/2/(1-phi * phi);
  S_eq_c = lnalpha_c * S_max;
  S_msy_c = S_eq * (0.5 - 0.07 * lnalpha_c);
  U_msy_c = lnalpha_c * (0.5 - 0.07 * lnalpha_c);

}



