// input data
data {
  int<lower=1> N; // number of observation
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  vector[N] fsm; // predictor : frost safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=0,upper=1> K_int;
  real <lower=0> r_fsm;
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  proba = K_int + r_fsm * fsm;

}

model {
  // How the data are expected to have been generated:
  presence ~ bernoulli(proba);
}
