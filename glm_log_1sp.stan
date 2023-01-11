// input data
data {
  int<lower=1> N; // number of observation
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  real prior_t_fsm;
  real prior_t_hsm;
  vector[N] fsm; // predictor : frost safety margins
  vector[N] hsm; // predictor : hydraulic safety margins
}

// the parameters to be estimated from the data
parameters {
  real <lower=0> K_int;
  real <lower=0.1> r_fsm;
  real<lower=0.1> r_hsm;
  real <lower=min(fsm),upper=max(fsm)> t_fsm; //point where regression changes
  real <lower=min(hsm),upper=max(hsm)> t_hsm; //point where regression changes
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  vector [N] K_vect;
  K_vect=rep_vector(K_int,N);
  proba = K_vect ./
                (
                  (1 + exp(-r_fsm * (fsm - t_fsm))).*
                  (1 + exp(-r_hsm * (hsm - t_hsm)))
                  );
  // for (i in 1:N) {
  //   proba[i]= K_int / 
  //               (
  //                 (1 + exp(-r_fsm * (fsm[i] - t_fsm)))*
  //                 (1 + exp(-r_hsm * (hsm[i] - t_hsm)))
  //                 );
  // }
}
// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  K_int~uniform(0,1);
  t_fsm~normal(prior_t_fsm,5);
  t_hsm~normal(prior_t_hsm,5);

  // How the data are expected to have been generated:
  presence ~ bernoulli(proba);
}

// // predictions
// generated quantities {
//   int presence_pred [N];
//   for (i in 1:N){
//     presence_pred[i] = bernoulli_rng(K_sp[species[i]] /(
//                           (1 + exp(-r_fsm * (fsm[i] - t_fsm)))*
//                           (1+exp(-r_hsm * (hsm[i] - t_hsm)))
//                           ));
//   }
// }
