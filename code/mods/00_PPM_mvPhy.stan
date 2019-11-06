data {
  int<lower=0> nCell; // number of grid cells
  int<lower=0> nSpp;  // number of species
  int<lower=0> nGen;  // number of species
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> tax_i[nSpp,2];  // species-genus lookup
  int<lower=0> W[nCell,nSpp];  // cell-scale citizen science counts
  int<lower=0> Y[nCell,nSpp];  // cell-scale citizen science counts
  matrix[nCell,R] X;  // cell-scale citizen science covariates
  real effort_prSoil;
}

parameters {
  matrix[R,nSpp] b;  // environmental covariate species slopes
  real<lower=0> sigma_b;
  matrix[R,nGen] B;  // environmental covariate genus slopes
  vector[R] beta;
  cholesky_factor_corr[nGen] L_Omega;
  vector[R] eta;
}

transformed parameters {
  matrix<lower=0>[nCell,nSpp] LAMBDA = exp(X * b);
  vector<lower=0, upper=1>[nCell] E = inv_logit(X * eta);
}

model {
  L_Omega ~ lkj_corr_cholesky(2);
  beta[1] ~ normal(12, 1);
  eta[1] ~ normal(-10, 1);
  for(r in 2:R) {
    beta[r] ~ normal(0, 1);
    eta[r] ~ normal(0, 1);
  }
  for(r in 1:R) {
    b[r,] ~ normal(B[r,tax_i[,2]], sigma_b);
    B[r,] ~ multi_normal_cholesky(rep_vector(beta[r], nGen), L_Omega);
  }
  for(s in 1:nSpp) {
    Y[,s] ~ poisson(LAMBDA[,s]*effort_prSoil);
    W[,s] ~ poisson(LAMBDA[,s].*E);
  }
}

