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
  int<lower=0> Y_[I_,S];  // Y counts (test)
  
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
  real<lower=0> sigma_b[R];  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega_B;  // genus-level correlation matrix 
  
  // slopes: plot level
  matrix[L,S] a_std;  // species slopes (stdNorm)
  real<lower=0> sigma_a[L];  // conspecific sd 
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
  matrix[K+J,S] lLAMBDA;  // cell level
  matrix[I,S] llambda;  // plot level lambda
  
  // cell level
  for(r in 1:R) {
    B[r,] = beta[r] + B_std[r,] * L_Omega_B;  // B ~ mvNorm(beta, L_Omega_B)
    b[r,] = B[r,tax_i[,2]] + b_std[r,] * sigma_b[r]; // b ~ Norm(B, sigma_b)
  }
  lLAMBDA = X * b;
  
  
  // plot level
  for(l in 1:L) {
    A[l,] = alpha[l] + A_std[l,] * L_Omega_A;  // A ~ mvNorm(alpha, L_Omega_A)
    a[l,] = A[l,tax_i[,2]] + a_std[l,] * sigma_a[l]; // a ~ Norm(A, sigma_a)
  }
  {
    matrix[J,S] lLAMBDA_Y = block(lLAMBDA, K+1, 1, J, S);  // cell level Y
    llambda = V * a + log(h) + lLAMBDA_Y[IJ,];
  }
  
}

model {
  
  // cell level priors
  beta[1] ~ normal(6.5, 1);
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
    Y[,s] ~ poisson_log(llambda[,s]);
  }
  
}

generated quantities {
  matrix[K_+J_,S] lLAMBDA_ = X_ * b;
  matrix[I_,S] llambda_;
  matrix[I_,S] log_lik_lambda_;
  matrix<lower=0>[K+J,S] LAMBDA = exp(lLAMBDA);
  matrix<lower=0>[K_+J_,S] LAMBDA_ = exp(lLAMBDA_);
  matrix<lower=0, upper=1>[K+J,S] p;
  matrix<lower=0, upper=1>[K_+J_,S] p_;
  vector[K+J] ShannonH;
  vector[K_+J_] ShannonH_;

  {
    matrix[J_,S] lLAMBDA_Y_ = block(lLAMBDA_, K_+1, 1, J_, S);
    llambda_ = V_ * a + log(h) + lLAMBDA_Y_[IJ_,];
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



