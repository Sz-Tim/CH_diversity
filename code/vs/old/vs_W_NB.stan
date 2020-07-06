data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> J_; // grid cells for Y (train)
  int<lower=0> I_; // plots (test)
  int<lower=0> IJ_[I_]; // plot to grid cell lookup (test)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  
  // observed data
  int<lower=0> W[K,S];  // W counts (train)
  int<lower=0> Y_[I_,S];  // Y counts (train)
  
  // covariates
  matrix[J+K,R] X;  // W grid cell covariates (train)
  matrix[J_,R] X_;  // W grid cell covariates (test)
  real h;  // Y (plot area)/(cell area)
  
}



parameters {
  
  // slopes: cell level
  matrix[R,S] b_std;  // species slopes (stdNorm)
  real<lower=0> sigma_b[R];  // conspecific sd 
  matrix[R,G] B_std;  // genus slopes (stdNorm)
  vector[R] beta;  // overall slopes
  cholesky_factor_corr[G] L_Omega_B[R];  // genus-level correlation matrix 
  vector<lower=0>[G] sigma_B[R];  // genus-level correlation matrix 
  
  real<lower=0> disp_lam;

}



transformed parameters {
  
  // slopes
  matrix[R,S] b;  // species cell level
  matrix[R,G] B;  // genus cell level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  
  // cell level
  for(r in 1:R) {
    // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
    // b ~ Norm(B, sigma_b)
    B[r,] = beta[r] + B_std[r,] * diag_pre_multiply(sigma_B[r], L_Omega_B[r]);  
    b[r,] = B[r,tax_i[,2]] + b_std[r,] * sigma_b[r]; 
  }
  lLAMBDA = X * b;
  
}



model {

  disp_lam ~ cauchy(0, 2);
  
  // cell level priors
  beta[1] ~ normal(6, 2);
  beta[2:R] ~ normal(0, 1);
  for(r in 1:R) {
    b_std[r,] ~ normal(0, 1);
    B_std[r,] ~ normal(0, 1);
    L_Omega_B[r] ~ lkj_corr_cholesky(2);
    sigma_B[r] ~ cauchy(0, 2);
  }
  sigma_b ~ cauchy(0, 2);
  
  // likelihood
  for(k in 1:K) {
    W[k,] ~ multinomial(softmax( lLAMBDA[k,]' ));
  }
  
}



generated quantities {
  
  matrix[J_,S] lLAMBDA_ = X_ * b;
  matrix[I_,S] llambda_;
  matrix[I_,S] log_lik_;
  matrix[K+J,S] prPres;
  matrix[J_,S] prPres_;
  matrix[I_,S] prPresL_;
  vector[K+J] ShannonH;
  vector[J_] ShannonH_;
  matrix[G,G] Sigma_B[R];

    llambda_ = log(h) + lLAMBDA_[IJ_,];
    for(s in 1:S) {
      for(i in 1:I_) {
        log_lik_[i,s] = neg_binomial_2_log_lpmf(Y_[i,s] | llambda_[i,s], disp_lam);
        prPresL_[i,s] = 1 - exp(neg_binomial_2_log_lpmf(0 | llambda_[i,s], disp_lam));
      }
    }
  
  {
    matrix[K+J,S] p;
    matrix[J_,S] p_;
    matrix[K+J,S] LAMBDA = exp(lLAMBDA);
    matrix[J_,S] LAMBDA_=exp(lLAMBDA_); 
    prPres = 1 - exp(-LAMBDA);
    prPres_ = 1 - exp(-LAMBDA_);
    for(i in 1:(K+J)) {
      p[i,] = LAMBDA[i,] / sum(LAMBDA[i,]);
    }
    for(i in 1:J_) {
      p_[i,] = LAMBDA_[i,] / sum(LAMBDA_[i,]);
    }
    ShannonH = - rows_dot_product(p, log(p));
    ShannonH_ = - rows_dot_product(p_, log(p_));
  }
  
  for(r in 1:R) {
    Sigma_B[r] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_B[r]),
                                sigma_B[r]);
  }
  
}



