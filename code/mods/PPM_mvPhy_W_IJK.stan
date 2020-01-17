data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> K_; // grid cells for W (test)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> J_; // grid cells for Y (test)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  real<lower=0.5, upper=1.5> D_prior[S];
  
  // observed data
  int<lower=0> W[K,S];  // W counts (train)
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[K_+J_,R] X_;  // W grid cell covariates (test)
  matrix[K,Q] U;  // W effort covariates (train)
  matrix[K_,Q] U_;  // W effort covariates (test)
  
}



parameters {
  
  // slopes: cell level
  matrix[R,S] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b;  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega_B;  // genus-level correlation matrix 
  
  // slopes: W effort
  vector[Q] eta;  // W sampling effort slopes 
  
  // bias: W species random effects
  vector<lower=0>[S] D;  // species bias in W

}



transformed parameters {
  
  // slopes
  matrix[R,S] b;  // species cell level
  matrix[R,G] B;  // genus cell level
  
  // Lambda/lambda
  matrix<lower=0>[K+J,S] LAMBDA;  // cell level
  vector<lower=0, upper=1>[K] E = inv_logit(U * eta);  // sampling effort
  
  // cell level
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega_B;  // B ~ mvNorm(beta, L_Omega_B)
  }
  b = B[,tax_i[,2]] + b_std * sigma_b; // b ~ Norm(B, sigma_b)
  LAMBDA = exp(X * b);
  
}

model {
  
  // effort and species bias priors
  D ~ normal(D_prior, 1);
  eta[1] ~ normal(-10, 1);
  eta[2:Q] ~ normal(0, 1);
  
  // cell level priors
  beta[1] ~ normal(10, 1);
  beta[2:R] ~ normal(0, 1);
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    B_std[r,] ~ normal(0, 1);
  }
  sigma_b ~ normal(0, 1);
  L_Omega_B ~ lkj_corr_cholesky(2);
  
  // likelihood
  for(s in 1:S) {
    W[,s] ~ poisson(LAMBDA[1:K,s].*E*D[s]);
  }
  
}

generated quantities {
  // matrix<lower=0>[K_W,S] What;  // latent lambda W
  // matrix<lower=0>[I_Y,S] Yhat;  // latent lambda Y
  matrix<lower=0>[K_+J_,S] LAMBDA_ = exp(X_ * b);
  vector<lower=0, upper=1>[K_] E_ = inv_logit(U_ * eta);
  // matrix<lower=0>[K_W_,S] W_hat;  // latent lambda W
  // matrix<lower=0>[I_Y_,S] Y_hat;  // latent lambda Y
  
  // for(s in 1:S) {
  //   What[,s] = LAMBDA_W[,s].*E*D[s];
  //   Yhat[,s] = lambda_Y[,s];
  //   W_hat[,s] = LAMBDA_W_[,s].*E_*D[s];
  //   Y_hat[,s] = lambda_Y_[,s];
  // }
}



