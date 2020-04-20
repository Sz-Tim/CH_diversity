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
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  
  // observed data
  int<lower=0> W[K,S];  // W counts (train)
  int<lower=0> Y[I,S];  // Y counts (train)
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
  cholesky_factor_corr[S] L_Omega_b[R];  // species-level correlation matrix 
  vector<lower=0>[S] sigma_b[R];  // species-level correlation matrix 
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
    // b ~ mvNorm(beta, sigma_b * L_Omega_b * sigma_b);  
    b[r,] = beta[r] + b_std[r,] * diag_pre_multiply(sigma_b[r], L_Omega_b[r]);
  }
  // cell level
  lLAMBDA = X * b;
  
}

model {
  
  // effort and species bias priors
  eta[1] ~ normal(-12, 1);
  eta[2:Q] ~ normal(0, 1);
  
  // cell level priors
  beta[1] ~ normal(6, 2);
  beta[2:R] ~ normal(0, 1);
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    L_Omega_b[r] ~ lkj_corr_cholesky(2);
    sigma_b[r] ~ cauchy(0, 2);
  }
  
  // likelihood
  for(s in 1:S) {
    W[,s] ~ poisson_log(lLAMBDA[1:K,s] + log(E));
  }
  
}

generated quantities {
  
  matrix[K_+J_,S] lLAMBDA_ = X_ * b;
  matrix[I,S] llambda;
  matrix[I_,S] llambda_;
  matrix[I_,S] log_lik_lambda_;
  matrix<lower=0>[K+J,S] LAMBDA = exp(lLAMBDA);
  matrix<lower=0>[K_+J_,S] LAMBDA_ = exp(lLAMBDA_);
  matrix<lower=0, upper=1>[K+J,S] p;
  matrix<lower=0, upper=1>[K_+J_,S] p_;
  vector[K+J] ShannonH;
  vector[K_+J_] ShannonH_;
  matrix[K+J,S] prPres;
  matrix[K_+J_,S] prPres_;
  cov_matrix[S] Sigma_b[R];

  // calculated predicted LAMBDA and lambda
  {
    matrix[J,S] lLAMBDA_Y = block(lLAMBDA, K+1, 1, J, S);
    matrix[J_,S] lLAMBDA_Y_ = block(lLAMBDA_, K_+1, 1, J_, S);
    llambda = log(h) + lLAMBDA_Y[IJ,];
    llambda_ = log(h) + lLAMBDA_Y_[IJ_,];
  }
  for(s in 1:S) {
    for(i in 1:I_) {
      log_lik_lambda_[i,s] = poisson_log_lpmf(Y_[i,s] | llambda_[i,s]);
    }
  }
  
  // Shannon H: calculate p, then H
  for(i in 1:(K+J)) {
    p[i,] = LAMBDA[i,] / sum(LAMBDA[i,]);
    prPres[i,] = 1-exp(-LAMBDA[i,]);
  }
  for(i in 1:(K_+J_)) {
    p_[i,] = LAMBDA_[i,] / sum(LAMBDA_[i,]);
    prPres_[i,] = 1-exp(-LAMBDA_[i,]);
  }
  ShannonH = - rows_dot_product(p, log(p));
  ShannonH_ = - rows_dot_product(p_, log(p_));
  
  // compose correlation matrices = sigma * LL' * sigma
  for(r in 1:R) {
    Sigma_b[r] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_b[r]), 
                                sigma_b[r]);
  }
  
}



