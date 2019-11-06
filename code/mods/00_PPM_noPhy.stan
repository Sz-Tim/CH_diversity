data {
  int<lower=0> nCell; // number of grid cells
  int<lower=0> nSpp;  // number of species
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> W[nCell,nSpp];  // cell-scale citizen science counts
  int<lower=0> Y[nCell,nSpp];  // cell-scale citizen science counts
  matrix[nCell,R] X;  // cell-scale citizen science covariates
  real effort_prSoil;
}

parameters {
  vector[R] beta;
  vector[R] eta;
}

transformed parameters {
  matrix<lower=0>[nCell,nSpp] LAMBDA;
  vector<lower=0, upper=1>[nCell] E = inv_logit(X * eta);
  for(s in 1:nSpp) {
    LAMBDA[,s] = exp(X * beta);
  }
}

model {
  beta[1] ~ normal(12, 1);
  eta[1] ~ normal(-10, 1);
  for(r in 2:R) {
    beta[r] ~ normal(0, 1);
    eta[r] ~ normal(0, 1);
  }
  for(s in 1:nSpp) {
    Y[,s] ~ poisson(LAMBDA[,s]*effort_prSoil);
    W[,s] ~ poisson(LAMBDA[,s].*E);
  }
}

