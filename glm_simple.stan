data {
  int<lower=1> N; // number of observations
  //int<lower=1> S; //number of species
  int<lower=0,upper=1> presence[N]; //obs of presence absence
  vector[N] hsm;
  vector[N] fsm;
  //vector<upper=S> [N] species; //species of each obs
}

parameters {
  real intercept;
  real<lower=0> b_hsm;
  real<lower=0> b_fsm;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  presence ~ bernoulli_logit(intercept + b_hsm * hsm + b_fsm * fsm);
}
