// input data
data {
  int<lower=1> N; // number of observations
  int<lower=1> S; //number of species
  int <lower=0,upper=1> presence [N]; //obs of presence absence
  int <lower=1,upper=S> species [N];
  real prior [4];
  vector[N] mat; // predictor : frost safety margins
  vector[N] wai; // predictor : hydraulic safety margins
  // real mean_mat;
  // real mean_wai;
}

// the parameters to be estimated from the data
parameters {
  real <lower=0> K_int;
  vector <lower=0,upper=1> [S] K_sp; // plateau of the second segment, or threshold value
  real <lower=0.1> r_mat;
  real<lower=0.1> r_wai;
  real <lower=min(mat),upper=max(mat)> t_mat; //point where regression changes
  real <lower=min(wai),upper=max(wai)> t_wai; //point where regression changes
}

transformed parameters {
  vector <lower=0,upper=1> [N] proba;
  proba = K_sp[species] ./
                (
                  (1 + exp(-r_mat * (mat - t_mat))).*
                  (1 + exp(-r_wai * (wai - t_wai)))
                  );
}
// The model to be estimated.
// presence are distributed according to Bernoulli law
// logit used as link-function
model {
  //priors
  K_int~normal(0.5,1);
  K_sp~normal(K_int,2); 
  r_mat~normal(prior[1],1);
  r_wai~normal(prior[2],1);
  t_mat~normal(prior[3],1);
  t_wai~normal(prior[4],1);
  
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
