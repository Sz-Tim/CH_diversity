data {
  // dimensions and indices
  int<lower=0> nCell_W; // grid cells for W (train)
  int<lower=0> nCell_W_; // grid cells for W (test)
  int<lower=0> nCell_Y; // grid cells for Y (train)
  int<lower=0> nCell_Y_; // grid cells for Y (test)
  int<lower=0> nSpp;  // number of species
  int<lower=0> nGen;  // number of species
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  // taxonomy
  int<lower=0> tax_i[nSpp,2];  // species-genus lookup
  // observed data
  int<lower=0> W[nCell_W,nSpp];  // W counts (train)
  int<lower=0> Y[nCell_Y,nSpp];  // Y counts (train)
  // covariates
  matrix[nCell_W,R] X_W;  // W covariates (train)
  matrix[nCell_W_,R] X_W_;  // W covariates (test)
  matrix[nCell_Y,R] X_Y;  // Y covariates (train)
  matrix[nCell_Y_,R] X_Y_;  // Y covariates (test)
  matrix[nCell_W,Q] V_W;  // W effort covariates (train)
  matrix[nCell_W_,Q] V_W_;  // W effort covariates (test)
  real effort_prSoil;  // proportion of cell sampled for Y
}

parameters {
  matrix[R,nSpp] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b;  // conspecific sd 
  matrix[R,nGen] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[nGen] L_Omega;  // genus-level correlation matrix 
  vector[Q] eta;  // W sampling effort slopes 
  vector<lower=0>[nSpp] spp_effort;  // species bias in W
}

transformed parameters {
  matrix[R,nSpp] b;  // species slopes
  matrix[R,nGen] B;  // genus slopes
  matrix<lower=0>[nCell_W,nSpp] LAMBDA_W;  // latent lambda W
  matrix<lower=0>[nCell_Y,nSpp] LAMBDA_Y;  // latent lambda Y
  vector<lower=0, upper=1>[nCell_W] E = inv_logit(V_W * eta);  // sampling effort
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega;  // B ~ mvNorm(beta, L_Omega)
  }
  b = B[,tax_i[,2]] + b_std * sigma_b; // b ~ Norm(B, sigma_b)
  LAMBDA_W = exp(X_W * b);
  LAMBDA_Y = exp(X_Y * b);
}

model {
  spp_effort ~ normal(1, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  sigma_b ~ normal(0, 1);
  beta[1] ~ normal(10, 1);
  eta[1] ~ normal(-10, 1);
  for(r in 2:R) {
    beta[r] ~ normal(0, 1);
  }
  for(q in 2:Q) {
    eta[q] ~ normal(0, 1);
  }
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    B_std[r,] ~ normal(0, 1);
  }
  for(s in 1:nSpp) {
    W[,s] ~ poisson(LAMBDA_W[,s].*E*spp_effort[s]);
    Y[,s] ~ poisson(LAMBDA_Y[,s]*effort_prSoil);
  }
}

generated quantities {
  // int<lower=0> Wpp[nCell_W,nSpp];  // W counts (train)
  // int<lower=0> Ypp[nCell_Y,nSpp];  // Y counts (train)
  matrix<lower=0>[nCell_W,nSpp] What;  // latent lambda W
  matrix<lower=0>[nCell_Y,nSpp] Yhat;  // latent lambda Y
  // int<lower=0> W_[nCell_W_,nSpp];  // W counts (train)
  // int<lower=0> Y_[nCell_Y_,nSpp];  // Y counts (train)
  matrix<lower=0>[nCell_W_,nSpp] LAMBDA_W_ = exp(X_W_ * b);
  matrix<lower=0>[nCell_Y_,nSpp] LAMBDA_Y_ = exp(X_Y_ * b);
  vector<lower=0, upper=1>[nCell_W_] E_ = inv_logit(V_W_ * eta);
  matrix<lower=0>[nCell_W_,nSpp] W_hat;  // latent lambda W
  matrix<lower=0>[nCell_Y_,nSpp] Y_hat;  // latent lambda Y

  for(s in 1:nSpp) {
    What[,s] = LAMBDA_W[,s].*E*spp_effort[s];
    Yhat[,s] = LAMBDA_Y[,s]*effort_prSoil;
    // Wpp[,s] = poisson_rng(What[,s]);
    // Ypp[,s] = poisson_rng(Yhat[,s]);
    W_hat[,s] = LAMBDA_W_[,s].*E_*spp_effort[s];
    Y_hat[,s] = LAMBDA_Y_[,s]*effort_prSoil;
    // W_[,s] = poisson_rng(W_hat[,s]);
    // Y_[,s] = poisson_rng(Y_hat[,s]);
  }
}



