data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
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
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[I,L] V;  // Y plot covariates (train)
  real h;  // Y (plot area)/(cell area)
  
}



transformed data {
  
  matrix[I,R+L] Z = append_col(X[(K+1):(K+J),][IJ,], V);

}



parameters {
  
  // slopes
  matrix[R+L,S] b_std;
  real<lower=0> sigma_b[R+L];
  matrix[R+L,G] B_std;
  vector[R+L] beta;
  cholesky_factor_corr[G] L_Omega_B[R+L];
  vector<lower=0>[G] sigma_B[R+L]; 
  vector<lower=0>[R+L-1] beta_lam;  // horseshoe prior
  real<lower=0> beta_tau;  // horseshoe prior

}



transformed parameters {
  
  // slopes
  matrix[R+L,S] b;  // species level
  matrix[R+L,G] B;  // genus level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  matrix[I,S] llambda;  // plot level
  
  // cell level
  for(m in 1:(R+L)) {
    // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
    // b ~ Norm(B, sigma_b)
    B[m,] = beta[m] + B_std[m,] * diag_pre_multiply(sigma_B[m], L_Omega_B[m]);  
    b[m,] = B[m,tax_i[,2]] + b_std[m,] * sigma_b[m]; 
  }
  lLAMBDA = X * b[1:R,] - log(h);
  llambda = Z * b;
  
}



model {
  
   // cell level priors
  beta[1] ~ normal(-4, 2);
  beta[2:(R+L)] ~ normal(0, beta_tau * beta_lam);
  beta_tau ~ cauchy(0, 2);
  beta_lam ~ student_t(4, 0, 1);
  for(m in 1:(R+L)) {
    b_std[m,] ~ normal(0, 1);
    B_std[m,] ~ normal(0, 1);
    sigma_B[m] ~ normal(0, 2);
    L_Omega_B[m] ~ lkj_corr_cholesky(2);
  }
  sigma_b ~ normal(0, 2);
  
  // likelihood
  for(s in 1:S) {
    Y[,s] ~ poisson_log(llambda[,s]);
  }
  
}



generated quantities {

  matrix[G,G] Sigma_B[R+L];

  // compose correlation matrices = sigma * LL' * sigma
  for(m in 1:(R+L)) {
    Sigma_B[m] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_B[m]),
                                sigma_B[m]);
  }
  
}



