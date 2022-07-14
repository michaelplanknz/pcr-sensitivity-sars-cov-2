data {
  int<lower=0> N; // number of tests
  int<lower=0> P; // number of individuals
  int patient_ID[N]; // case IDs
  vector[N] day_of_test; // day of test (integer)
  int test_result[N]; // boolean test result
  vector[P] time_first_symptom; // day of symptom onset
  vector[P] te_upper_bound; // maximum time at which infection must have occured before
}

transformed data{
  vector<lower = 0>[P] inc_lower_bound = time_first_symptom - te_upper_bound; // lower bounds for truncated incubation periods, i.e. symptom onset time - min(symptom onset time, time of first positive test)
}

parameters {
  real beta1;
  real<upper = 0> beta2;
  real<upper = 0> beta3;
  //vector<lower = inc_lower_bound>[P] Inc; // requires stan version 2.24 to use vector as lower bound
  vector<lower = 0>[P] Inc_raw; // see https://mc-stan.org/docs/2_28/stan-users-guide/vectors-with-varying-bounds.html
  real<lower = 0> cutpoint;
}

transformed parameters {
  vector[P] Inc = inc_lower_bound + Inc_raw;
  vector <lower = 0> [P] T_e =  time_first_symptom - Inc; // infection times for each individual
  vector[N] diff = day_of_test - T_e[patient_ID];
  vector[N] lp;
   
  // PCR positivity log-likelihoods for each test result
  for(i in 1:N) {
  // piecewise quadratic: b0 + b1*(t-c)^2 + b2*H(t-c)*(t-c)^2 
  lp[i] = bernoulli_logit_lpmf(test_result[i] | beta1 + beta2 * pow((diff[i] - cutpoint), 2) + (-beta2 + beta3) * step(diff[i] - cutpoint) * pow((diff[i] - cutpoint), 2));
  }
}

model {

  // Prior log-density of incubation period for each individual j:
  for(j in 1:P) {
    // incubation period params: point estimate of logmean and logsd from Lauer et al. (2020)
    target += lognormal_lpdf(Inc_raw[j] + inc_lower_bound[j]  | 1.621,0.418) - lognormal_lccdf(inc_lower_bound[j] | 1.621,0.418); 
    //target += lognormal_lpdf(Inc_raw[j] + inc_lower_bound[j]  | 1.685,0.4871) - lognormal_lccdf(inc_lower_bound[j] | 1.685,0.4871); // sensitivity analysis lmean+1SD, lsd+1SD
    //target += lognormal_lpdf(Inc_raw[j] + inc_lower_bound[j]  | 1.557,0.4871) - lognormal_lccdf(inc_lower_bound[j] | 1.557,0.4871); // sensitivity analysis lmean-1SD, lsd+1SD
    //target += lognormal_lpdf(Inc_raw[j] + inc_lower_bound[j]  | 1.557,0.3489) - lognormal_lccdf(inc_lower_bound[j] | 1.557,0.3489); // sensitivity analysis lmean-1SD, lsd-1SD
    //target += lognormal_lpdf(Inc_raw[j] + inc_lower_bound[j]  | 1.685,0.3489) - lognormal_lccdf(inc_lower_bound[j] | 1.685,0.3489); // sensitivity analysis lmean+1SD, lsd-1SD
  }
 
  // Priors
  beta1 ~ normal(0, 5); // mean 0, sd 5
  beta2 ~ normal(0, 1);  
  beta3 ~ normal(0, 1); 
  cutpoint ~ normal(5, 5) T[0, ]; // mean 5, sd 5
  
  // Joint PCR positivity log-likelihood
  target += sum(lp);
}

generated quantities {
  vector[501] p;
  real k;
  for(j in 1:501) {
    k = (j * 0.1) - 0.1;
    p[j] = inv_logit(beta1 + beta2 * pow((k - cutpoint), 2) + (-beta2 + beta3) * step(k - cutpoint) * pow((k - cutpoint), 2));
  }
}
