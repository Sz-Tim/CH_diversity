data {
  // dimensions and indices
  int<lower=0> K_W; // grid cells for W (train)
  int<lower=0> K_W_; // grid cells for W (test)
  int<lower=0> K_Y; // grid cells for Y (train)
  int<lower=0> K_Y_; // grid cells for Y (test)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  // taxonomy
  int<lower=0> tax_i[S,3];  // species-genus lookup
  // observed data
  int<lower=0> W[K_W,S];  // W counts (train)
  int<lower=0> Y[K_Y,S];  // Y counts (train)
  // covariates
  matrix[K_W,R] X_W;  // W covariates (train)
  matrix[K_W_,R] X_W_;  // W covariates (test)
  matrix[K_Y,R] X_Y;  // Y covariates (train)
  matrix[K_Y_,R] X_Y_;  // Y covariates (test)
  matrix[K_W,Q] U_W;  // W effort covariates (train)
  matrix[K_W_,Q] U_W_;  // W effort covariates (test)
  real h;  // proportion of cell sampled for Y
}

parameters {
  matrix[R,S] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b;  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega;  // genus-level correlation matrix 
  vector[Q] eta;  // W sampling effort slopes 
  vector<lower=0>[S] D;  // species bias in W
}

transformed parameters {
  matrix[R,S] b;  // species slopes
  matrix[R,G] B;  // genus slopes
  matrix<lower=0>[K_W,S] LAMBDA_W;  // latent lambda W
  matrix<lower=0>[K_Y,S] LAMBDA_Y;  // latent lambda Y
  vector<lower=0, upper=1>[K_W] E = inv_logit(U_W * eta);  // sampling effort
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega;  // B ~ mvNorm(beta, L_Omega)
  }
  b = B[,tax_i[,2]] + b_std * sigma_b; // b ~ Norm(B, sigma_b)
  LAMBDA_W = exp(X_W * b);
  LAMBDA_Y = exp(X_Y * b);
}

model {
  for(s in 1:S) {
    D[s] ~ normal(tax_i[s,3]-0.5, 1);
  }
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
  for(s in 1:S) {
    W[,s] ~ poisson(LAMBDA_W[,s].*E*D[s]);
    Y[,s] ~ poisson(LAMBDA_Y[,s]*h);
  }
}

generated quantities {
  // int<lower=0> Wpp[K_W,S];  // W counts (train)
  // int<lower=0> Ypp[K_Y,S];  // Y counts (train)
  matrix<lower=0>[K_W,S] What;  // latent lambda W
  matrix<lower=0>[K_Y,S] Yhat;  // latent lambda Y
  // int<lower=0> W_[K_W_,S];  // W counts (train)
  // int<lower=0> Y_[K_Y_,S];  // Y counts (train)
  matrix<lower=0>[K_W_,S] LAMBDA_W_ = exp(X_W_ * b);
  matrix<lower=0>[K_Y_,S] LAMBDA_Y_ = exp(X_Y_ * b);
  vector<lower=0, upper=1>[K_W_] E_ = inv_logit(U_W_ * eta);
  matrix<lower=0>[K_W_,S] W_hat;  // latent lambda W
  matrix<lower=0>[K_Y_,S] Y_hat;  // latent lambda Y
  matrix<lower=0>[K_W,S] prPres_W = 1 - exp(-LAMBDA_W);
  matrix<lower=0>[K_W_,S] prPres_W_ = 1 - exp(-LAMBDA_W_);
  matrix<lower=0>[K_Y,S] prPres_Y = 1 - exp(-LAMBDA_Y);
  matrix<lower=0>[K_Y_,S] prPres_Y_ = 1 - exp(-LAMBDA_Y_);

  for(s in 1:S) {
    What[,s] = LAMBDA_W[,s].*E*D[s];
    Yhat[,s] = LAMBDA_Y[,s]*h;
    // Wpp[,s] = poisson_rng(What[,s]);
    // Ypp[,s] = poisson_rng(Yhat[,s]);
    W_hat[,s] = LAMBDA_W_[,s].*E_*D[s];
    Y_hat[,s] = LAMBDA_Y_[,s]*h;
    // W_[,s] = poisson_rng(W_hat[,s]);
    // Y_[,s] = poisson_rng(Y_hat[,s]);
  }
}



