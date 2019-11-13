data {
  // dimensions and indices
  int<lower=0> K_Y; // grid cells for Y (train)
  int<lower=0> K_Y_; // grid cells for Y (test)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  // taxonomy
  int<lower=0> tax_i[S,3];  // species-genus lookup
  // observed data
  int<lower=0> Y[K_Y,S];  // Y counts (train)
  // covariates
  matrix[K_Y,R] X_Y;  // Y covariates (train)
  matrix[K_Y_,R] X_Y_;  // Y covariates (test)
  real h;  // proportion of cell sampled for Y
}

parameters {
  matrix[R,S] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b;  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega;  // genus-level correlation matrix 
}

transformed parameters {
  matrix[R,S] b;  // species slopes
  matrix[R,G] B;  // genus slopes
  matrix<lower=0>[K_Y,S] LAMBDA_Y;  // latent lambda Y
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega;  // B ~ mvNorm(beta, L_Omega)
  }
  b = B[,tax_i[,2]] + b_std * sigma_b; // b ~ Norm(B, sigma_b)
  LAMBDA_Y = exp(X_Y * b);
}

model {
  L_Omega ~ lkj_corr_cholesky(2);
  sigma_b ~ normal(0, 1);
  beta[1] ~ normal(10, 1);
  for(r in 2:R) {
    beta[r] ~ normal(0, 1);
  }
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    B_std[r,] ~ normal(0, 1);
  }
  for(s in 1:S) {
    Y[,s] ~ poisson(LAMBDA_Y[,s]*h);
  }
}

generated quantities {
  // int<lower=0> Ypp[K_Y,S];  // Y counts (train)
  matrix<lower=0>[K_Y,S] Yhat;  // latent lambda Y
  // int<lower=0> Y_[K_Y_,S];  // Y counts (train)
  matrix<lower=0>[K_Y_,S] LAMBDA_Y_ = exp(X_Y_ * b);
  matrix<lower=0>[K_Y_,S] Y_hat;  // latent lambda Y
  matrix<lower=0>[K_Y,S] prPres_Y = 1 - exp(-LAMBDA_Y);
  matrix<lower=0>[K_Y_,S] prPres_Y_ = 1 - exp(-LAMBDA_Y_);

  for(s in 1:S) {
    Yhat[,s] = LAMBDA_Y[,s]*h;
    // Ypp[,s] = poisson_rng(Yhat[,s]);
    Y_hat[,s] = LAMBDA_Y_[,s]*h;
    // Y_[,s] = poisson_rng(Y_hat[,s]);
  }
}



