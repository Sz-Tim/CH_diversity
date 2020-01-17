data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> K_; // grid cells for W (test)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> J_; // grid cells for Y (test)
  int<lower=0> I; // plots (train)
  int<lower=0> I_; // plots (test)
  int<lower=0> IJ[I]; // plot to grid cell lookup (train)
  int<lower=0> IJ_[I_]; // plot to grid cell lookup (test)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> L;  // number of plot-scale covariates (no intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  
  // observed data
  int<lower=0> Y[I,S];  // Y counts (train)
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[K_+J_,R] X_;  // W grid cell covariates (test)
  matrix[I,L] V;  // Y plot covariates (train)
  matrix[I_,L] V_;  // Y plot covariates (test)
  real h;  // Y (plot area)/(cell area)
  
}



parameters {
  
  // slopes: cell level
  matrix[R,S] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b;  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega_B;  // genus-level correlation matrix 
  
  // slopes: plot level
  matrix[L,S] a_std;  // species slopes (stdNorm)
  real<lower=0> sigma_a;  // conspecific sd 
  matrix[L,G] A_std;  // genus slopes (stdNorm)
  vector[L] alpha;  // overall slopes
  cholesky_factor_corr[G] L_Omega_A;  // genus-level correlation matrix 

}



transformed parameters {
  
  // slopes
  matrix[R,S] b;  // species cell level
  matrix[R,G] B;  // genus cell level
  matrix[L,S] a;  // species plot level
  matrix[L,G] A;  // genus plot level
  
  // Lambda/lambda
  matrix<lower=0>[K+J,S] LAMBDA;  // cell level
  matrix<lower=0>[I,S] lambda;  // plot level lambda
  
  // cell level
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega_B;  // B ~ mvNorm(beta, L_Omega_B)
  }
  b = B[,tax_i[,2]] + b_std * sigma_b; // b ~ Norm(B, sigma_b)
  LAMBDA = exp(X * b);
  
  
  // plot level
  for(l in 1:L) {
    A[l,] = alpha[l] + A_std[l,] * L_Omega_A;  // A ~ mvNorm(alpha, L_Omega_A)
  }
  a = A[,tax_i[,2]] + a_std * sigma_a; // a ~ Norm(A, sigma_a)
  {
    matrix[J,S] LAMBDA_Y = block(LAMBDA, K+1, 1, J, S);  // cell level Y
    lambda = exp(V * a + log(h * LAMBDA_Y[IJ,]));
  }
  
}

model {
  
  // cell level priors
  beta[1] ~ normal(10, 1);
  beta[2:R] ~ normal(0, 1);
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    B_std[r,] ~ normal(0, 1);
  }
  sigma_b ~ normal(0, 1);
  L_Omega_B ~ lkj_corr_cholesky(2);
  
  // plot level priors
  alpha ~ normal(0, 1);
  for(l in 1:L) {
    a_std[l,] ~ normal(0, 1);
    A_std[l,] ~ normal(0, 1);
  }
  sigma_a ~ normal(0, 1);
  L_Omega_A ~ lkj_corr_cholesky(2);
  
  // likelihood
  for(s in 1:S) {
    Y[,s] ~ poisson(lambda[,s]);
  }
  
}

generated quantities {
  // matrix<lower=0>[K_W,S] What;  // latent lambda W
  // matrix<lower=0>[I_Y,S] Yhat;  // latent lambda Y
  matrix<lower=0>[K_+J_,S] LAMBDA_ = exp(X_ * b);
  matrix<lower=0>[I_,S] lambda_;
  // matrix<lower=0>[K_W_,S] W_hat;  // latent lambda W
  // matrix<lower=0>[I_Y_,S] Y_hat;  // latent lambda Y

  {
    matrix[J_,S] LAMBDA_Y_ = block(LAMBDA_, K_+1, 1, J_, S);
    lambda_ = exp(V_ * a + log(h * LAMBDA_Y_[IJ_,]));
  }
  // for(s in 1:S) {
  //   What[,s] = LAMBDA_W[,s].*E*D[s];
  //   Yhat[,s] = lambda_Y[,s];
  //   W_hat[,s] = LAMBDA_W_[,s].*E_*D[s];
  //   Y_hat[,s] = lambda_Y_[,s];
  // }
}



