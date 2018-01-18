data { 
  int W; // number of records
  int K; // number of islands  
  int O_t[W]; // number of reported cases
  real Ostar[W]; // exposition (past incidence weighted by the discretized serial interval)
  int sumO_t[W]; // cumulative number of reported cases
  int island[W]; // island index (1 to K)
  int pop[K]; // island population
}

parameters {
  // declaring hyperparameters (normally distributed random effects)
  real<lower=0> mu_R0_CHIKV; // non-centered
  real<lower=0> sigma_R0_CHIKV;
  real IL_mu_rho_CHIKV; // non-centered, with logit transformation
  real<lower=0> IL_sigma_rho_CHIKV;

  // declaring random parameters by island (non-centered)
  real R0_CHIKV_TILDE[K]; // island-level base transmission for CHIKV
  real IL_rho_CHIKV_TILDE[K]; // reporting rate by island (inverse-logit transformed)
  
  // declaring dispersion parameter
  real<lower=0> phi;
}

transformed parameters {
  // declaring rescaled random parameters
  real<lower=0> R0_CHIKV[K]; // island-level base transmission for CHIKV
  real<lower=0,upper=1> rho_CHIKV[K]; // reporting rate by island
  
  // declaring model intermediates
  real<lower=0> lp[W]; // mean prediction
  real<lower=0> sampledisp[W]; // dispersion = mean/phi

  // rescaling random parameters
  for(i in 1:K) {
    R0_CHIKV[i] = mu_R0_CHIKV + sigma_R0_CHIKV * R0_CHIKV_TILDE[i];
    rho_CHIKV[i] = inv_logit(IL_mu_rho_CHIKV + IL_sigma_rho_CHIKV * IL_rho_CHIKV_TILDE[i]);
  }
  
  // building negative binomial model
  for(i in 1:W) {
    lp[i] = R0_CHIKV[island[i]] * Ostar[i] * ( 1 - sumO_t[i] / ( rho_CHIKV[island[i]] * pop[island[i]] ) );
    if(lp[i]==0) lp[i] = 0.0001;
    sampledisp[i] = lp[i]/phi;
    if(sampledisp[i]==0) sampledisp[i] = 0.0001;
  }
}

model {
  // setting priors 
  mu_R0_CHIKV ~ gamma(1,0.2); // implies: R0 strictly positive, away from 0, probably between 1 and 5
  sigma_R0_CHIKV ~ cauchy(0,2.5); 
  IL_mu_rho_CHIKV ~ normal(0,1.5); // implies: mu_rho_CHIKV uniform on [0,1]
  IL_sigma_rho_CHIKV ~ cauchy(0,2.5);

  for(i in 1:K) {
    R0_CHIKV_TILDE[i] ~ normal(0,1); // implies R0_CHIKV ~ normal( mu_R0_CHIKV,sigma_R0_CHIKV )
    IL_rho_CHIKV_TILDE[i] ~ normal(0,1); // implies rho_CHIKV ~ inverse_logit( normal( IL_mu_rho_CHIKV,IL_sigma_rho_CHIKV ) )
  }
  
  phi ~ cauchy(0,2.5);

  // likelihood
  target += neg_binomial_2_lpmf(O_t|lp,sampledisp);
}

generated quantities {
  real log_lik[W]; // save the likelihood
  real pred_lp[W]; // predicted values around lp
  real resid_lp[W]; // residuals
  real mu_rho_CHIKV;

  // diagnostics
  for (i in 1:W) {
    log_lik[i] = neg_binomial_2_lpmf(O_t[i]|lp[i],sampledisp[i]);
    pred_lp[i] = neg_binomial_2_rng(lp[i],sampledisp[i]);
    resid_lp[i] = O_t[i] - lp[i];
  }
  
  // rescaling hyperparameter
  mu_rho_CHIKV = inv_logit(IL_mu_rho_CHIKV);
}
