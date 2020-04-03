data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> K_; // grid cells for W (test)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> J_; // grid cells for Y (test)
  int<lower=0> I_; // plots (test)
  int<lower=0> IJ_[I_]; // plot to grid cell lookup (test)
  int<lower=0> S;  // number of species
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  
  // observed data
  int<lower=0> W[K,S];  // W counts (train)
  int<lower=0> Y_[I_,S];  // Y counts (test)
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[K_+J_,R] X_;  // W grid cell covariates (test)
  matrix[K,Q] U;  // W effort covariates (train)
  matrix[K_,Q] U_;  // W effort covariates (test)
  real h;  // Y (plot area)/(cell area)
  
}



parameters {
  
  // slopes: cell level
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[S] L_Omega_b[R];  // genus-level correlation matrix 
  matrix[R,S] b_std;  // species slopes (stdNorm)
  
  // slopes: W effort
  vector[Q] eta;  // W sampling effort slopes 

}



transformed parameters {
  
  matrix[R,S] b;  // species cell level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  vector<lower=0, upper=1>[K] E = inv_logit(U * eta);  // sampling effort
  
  for(r in 1:R) {
    b[r,] = beta[r] + b_std[r,] * L_Omega_b[r];  // b ~ mvNorm(beta, L_Omega)
  }
  // cell level
  lLAMBDA = X * b;
  
}

model {
  
  // effort and species bias priors
  eta[1] ~ normal(-12, 1);
  eta[2:Q] ~ normal(0, 1);
  
  // cell level priors
  beta[1] ~ normal(6.5, 1);
  beta[2:R] ~ normal(0, 1);
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    L_Omega_b[r] ~ lkj_corr_cholesky(2);
  }
  
  // likelihood
  for(s in 1:S) {
    W[,s] ~ poisson_log(lLAMBDA[1:K,s] + log(E));
  }
  
}

generated quantities {
  matrix[K_+J_,S] lLAMBDA_ = X_ * b;
  matrix[I_,S] llambda_;
  matrix[I_,S] log_lik_lambda_;
  matrix[K+J,S] LAMBDA = exp(lLAMBDA);
  matrix[K_+J_,S] LAMBDA_ = exp(lLAMBDA_);
  matrix<lower=0, upper=1>[K+J,S] p;
  matrix<lower=0, upper=1>[K_+J_,S] p_;
  vector[K+J] ShannonH;
  vector[K_+J_] ShannonH_;

  {
    matrix[J_,S] lLAMBDA_Y_ = block(lLAMBDA_, K_+1, 1, J_, S);
    llambda_ = log(h) + lLAMBDA_Y_[IJ_,];
  }
  for(s in 1:S) {
    for(i in 1:I_) {
      log_lik_lambda_[i,s] = poisson_log_lpmf(Y_[i,s] | llambda_[i,s]);
    }
  }
  for(i in 1:(K+J)) {
    p[i,] = LAMBDA[i,] / sum(LAMBDA[i,]);
  }
  for(i in 1:(K_+J_)) {
    p_[i,] = LAMBDA_[i,] / sum(LAMBDA_[i,]);
  }
  ShannonH = - rows_dot_product(p, log(p));
  ShannonH_ = - rows_dot_product(p_, log(p_));
}



