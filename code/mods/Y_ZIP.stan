functions {
/* zero-inflated poisson log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   lambda: mean parameter of the poisson distribution (log scale)
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_poisson_lpmf(int y, real llambda, real zi) { 
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                         poisson_log_lpmf(0 | llambda)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             poisson_log_lpmf(y | llambda); 
    } 
  }
}

data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> K_; // grid cells for W (test)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> I; // plots (train)
  int<lower=0> IJ[I]; // plot to grid cell lookup (train)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> L;  // number of plot-scale covariates (no intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  
  // observed data
  int<lower=0> Y[I,S];  // Y counts (train)
  // matrix<lower=0>[I,S] Y;  // Y counts (train)
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[K_,R] X_;  // W grid cell covariates (test)
  matrix[I,L] V;  // Y plot covariates (train)
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
  
  // slopes: plot level
  matrix[L,S] a_std;  // species slopes (stdNorm)
  real<lower=0> sigma_a[L];  // conspecific sd 
  matrix[L,G] A_std;  // genus slopes (stdNorm)
  vector[L] alpha;  // overall slopes
  cholesky_factor_corr[G] L_Omega_A[L];  // genus-level correlation matrix 
  vector<lower=0>[G] sigma_A[L];  // genus-level correlation matrix 
  
  real<lower=0, upper=1> zi;

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
    // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
    // b ~ Norm(B, sigma_b)
    B[r,] = beta[r] + B_std[r,] * diag_pre_multiply(sigma_B[r], L_Omega_B[r]);  
    b[r,] = B[r,tax_i[,2]] + b_std[r,] * sigma_b[r]; 
  }
  lLAMBDA = X * b;
  
  
  // plot level
  for(l in 1:L) {
    // A ~ mvNorm(alpha, sigma_A * L_Omega_A * sigma_A);  
    // a ~ Norm(A, sigma_a)
    A[l,] = alpha[l] + A_std[l,] * diag_pre_multiply(sigma_A[l], L_Omega_A[l]);  
    a[l,] = A[l,tax_i[,2]] + a_std[l,] * sigma_a[l]; 
  }
  {
    matrix[J,S] lLAMBDA_Y = block(lLAMBDA, K+1, 1, J, S);  // cell level Y
    llambda = V * a + log(h) + lLAMBDA_Y[IJ,];
  }
  
}

model {
  
  // disp_lam ~ normal(0, 1);
  
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
  
  
  // plot level priors
  alpha ~ normal(0, 1);
  for(l in 1:L) {
    a_std[l,] ~ normal(0, 1);
    A_std[l,] ~ normal(0, 1);
    L_Omega_A[l] ~ lkj_corr_cholesky(2);
    sigma_A[l] ~ cauchy(0, 2);
  }
  sigma_a ~ cauchy(0, 2);
  
  
  // likelihood
  for(i in 1:I) {
    for(s in 1:S) {
      Y[i,s] ~ zero_inflated_poisson(llambda[i,s], zi);
    }
  }
  
}

generated quantities {
  
  matrix[K_,S] lLAMBDA_ = X_ * b;
  // matrix[I,S] log_lik_lambda;
  matrix[K+J,S] prPres;
  matrix[K_,S] prPres_;
  vector[K+J] ShannonH;
  vector[K_] ShannonH_;
  matrix[G,G] Sigma_B[R];
  matrix[G,G] Sigma_A[L];

  // calculated predicted LAMBDA and lambda
  // log_lik_lambda = CMP_lpdf(Y | disp_lam, llambda);
  // for(s in 1:S) {
  //   for(i in 1:I) {
  //    log_lik_lambda[i,s] = poisson_log_lpmf(Y[i,s] | llambda[i,s]);  
  //   }
  // }
  
  // Shannon H: calculate p, then H
  {
    matrix[K+J,S] p;
    matrix[K_,S] p_;
    matrix[K+J,S] LAMBDA = exp(lLAMBDA);
    matrix[K_,S] LAMBDA_=exp(lLAMBDA_);
    for(i in 1:(K+J)) {
      p[i,] = LAMBDA[i,] / sum(LAMBDA[i,]);
      prPres[i,] = 1-exp(-LAMBDA[i,]);
    }
    for(i in 1:K_) {
      p_[i,] = LAMBDA_[i,] / sum(LAMBDA_[i,]);
      prPres_[i,] = 1-exp(-LAMBDA_[i,]);
    }
    ShannonH = - rows_dot_product(p, log(p));
    ShannonH_ = - rows_dot_product(p_, log(p_)); 
  }
  
  // compose correlation matrices = sigma * LL' * sigma
  for(r in 1:R) {
    Sigma_B[r] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_B[r]), 
                                sigma_B[r]);
  }
  for(l in 1:L) {
    Sigma_A[l] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_A[l]), 
                                sigma_A[l]);
  }
  
}



