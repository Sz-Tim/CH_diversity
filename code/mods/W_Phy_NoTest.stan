data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W (train)
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> I; // plots (train)
  int<lower=0> IJ[I]; // plot to grid cell lookup (train)
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> Q;  // number of W effort covariates
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  real<lower=0.5, upper=1.5> D_prior[S];
  
  // observed data
  int<lower=0> W[K,S];  // W counts (train)
  int<lower=0> Y[I,S];  // Y counts (train)
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates (train)
  matrix[K,Q] U;  // W effort covariates (train)
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
  
  // slopes: W effort
  vector[Q] eta;  // W sampling effort slopes 

}



transformed parameters {
  
  // slopes
  matrix[R,S] b;  // species cell level
  matrix[R,G] B;  // genus cell level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  vector<lower=0, upper=1>[K] E = inv_logit(U * eta);  // sampling effort
  
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
  
  // effort and species bias priors
  // D ~ normal(D_prior, 1);
  eta[1] ~ normal(-12, 1);
  eta[2:Q] ~ normal(0, 1);
  
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
  for(s in 1:S) {
    W[,s] ~ poisson_log(lLAMBDA[1:K,s] + log(E));
  }
  
}

generated quantities {
  
  matrix[I,S] llambda;
  matrix[I,S] log_lik_lambda;
  matrix<lower=0>[K+J,S] LAMBDA = exp(lLAMBDA);
  matrix<lower=0, upper=1>[K+J,S] p;
  vector[K+J] ShannonH;
  matrix[K+J,S] prPres;
  matrix[G,G] Sigma_B[R];

  {
    matrix[J,S] lLAMBDA_Y = block(lLAMBDA, K+1, 1, J, S);
    llambda = log(h) + lLAMBDA_Y[IJ,];
  }
  for(s in 1:S) {
    for(i in 1:I) {
      log_lik_lambda[i,s] = poisson_log_lpmf(Y[i,s] | llambda[i,s]);
    }
  }
  for(i in 1:(K+J)) {
    p[i,] = LAMBDA[i,] / sum(LAMBDA[i,]);
    prPres[i,] = 1-exp(-LAMBDA[i,]);
  }
  ShannonH = - rows_dot_product(p, log(p));
  for(r in 1:R) {
    Sigma_B[r] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_B[r]), 
                                sigma_B[r]);
  }
  
}



