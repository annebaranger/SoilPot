// input data
data {
  int<lower=1> N; // number of observation
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  real prior_K;
}

// the parameters to be estimated from the data
parameters {
  real <lower=0,upper=1> K_int;
}

model {
  //priors
  K_int~normal(prior_K,0.1);

  // How the data are expected to have been generated:
  presence ~ bernoulli(K_int);
}
