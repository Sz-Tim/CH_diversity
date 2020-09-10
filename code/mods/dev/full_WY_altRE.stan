// full: All regional and local covariates  
// WY: Structured samples AND citizen science samples  
// output: Structured sample fitted values

functions {
  
  // Generalized poisson based on LaplacesDemon::dgpois
  real GP_log_lpdf(matrix y, real lambda, matrix log_theta) {
    matrix[dims(y)[1],dims(y)[2]] theta = exp(log_theta);
    matrix[dims(y)[1],dims(y)[2]] y_lambda = y*lambda;
    matrix[dims(y)[1],dims(y)[2]] theta_1mlam = theta * (1-lambda);
    return sum(log(theta_1mlam))
            + sum((y-1) .* log(theta_1mlam + y_lambda)) 
            - sum(lgamma(y+1))
            + sum(theta * (lambda-1))
            - sum(y_lambda);
  }
  
  real GP_point_log_lpdf(real y, real lambda, real log_theta) {
    real theta = exp(log_theta);
    real y_lambda = y*lambda;
    real theta_1mlam = theta * (1-lambda);
    return log(theta_1mlam) 
            + lmultiply(y-1, theta_1mlam + y_lambda)
            - lgamma(y+1) 
            + fma(theta, lambda-1, -y_lambda);
  }
  
}



data {
  
  // dimensions and indices
  int<lower=0> K; // grid cells for W 
  int<lower=0> J; // grid cells for Y 
  int<lower=0> I; // plots 
  int<lower=0> IJ[I]; // plot to grid cell lookup 
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> L;  // number of plot-scale covariates (no intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  real<lower=0.5, upper=1.5> D_prior[S];
  
  // observed data
  int<lower=0> W[K,S];  // W counts 
  matrix<lower=0>[I,S] Y;  // Y counts
  
  // covariates
  matrix[K+J,R] X;  // grid cell covariates 
  matrix[I,L] V;  // plot covariates 
  real h;  // Y (plot area)/(cell area)
  
}



transformed data {
  
  matrix[I,R+L] Z = append_col(X[(K+1):(K+J),][IJ,], V);

}



parameters {
  
  // slopes
  vector[R+L] beta;  // aggregate slopes
  vector<lower=0>[R+L-1] beta_lam;  // horseshoe prior
  real<lower=0> beta_tau;  // horseshoe prior
  
  matrix[R+L,S] b_std;  // species-specific slopes
  vector<lower=0>[R+L] sigma_b;  // sd for species-specific slopes
  cholesky_factor_corr[R+L] L_Omega_b;  // cholesky corr matrix for b's

  
  real<lower=0, upper=1> disp_lam;  // implied prior: Uniform(0,1)
  
  // bias: W species random effects
  vector<lower=0>[S] D;  // species bias in W

}



transformed parameters {
  
  // slopes
  matrix[R+L,S] b_eff;  // non-centered species-specific slope effects
  matrix[R+L,S] b; // non-centered species-specific slopes 
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  matrix[I,S] llambda;  // plot level
  
  b_eff = diag_pre_multiply(sigma_b, L_Omega_b) * b_std;
  
  for(m in 1:(R+L)) {
    // b ~ mvNorm(beta, sigma_b * L_Omega_b * sigma_b) 
    b[m,] = beta[m] + b_eff[m,]; 
  }
  
  lLAMBDA = X * b[1:R,] - log(h);
  llambda = Z * b;
  
}



model {
  
  // effort and species bias priors
  D ~ normal(D_prior, 1);
  
  // cell level priors
  beta[1] ~ normal(-4, 2);
  beta[2:(R+L)] ~ normal(0, beta_tau * beta_lam);
  beta_tau ~ cauchy(0, 2);
  beta_lam ~ student_t(4, 0, 1);
  
  to_vector(b_std) ~ normal(0, 1);
  sigma_b ~ normal(0, 2);
  L_Omega_b ~ lkj_corr_cholesky(2);
  
  // likelihood
  for(k in 1:K) {
    W[k,] ~ multinomial(softmax( lLAMBDA[k,]' .* D ));
  }
  Y ~ GP_log(disp_lam, llambda);
  
}



generated quantities {
  
  matrix[I,S] loglik;
  matrix[I,S] prPresL;
  matrix[R+L,R+L] Sigma_b;

  // calculated predicted plot-scale prPres
  for(s in 1:S) {
    for(i in 1:I) {
      loglik[i,s] = GP_point_log_lpdf(Y[i,s] | disp_lam, llambda[i,s]);
      prPresL[i,s] = 1 - exp(GP_point_log_lpdf(0 | disp_lam, llambda[i,s]));
    }
  }
  
  // compose correlation matrices = sigma * LL' * sigma
  Sigma_b = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_b), sigma_b);
  
}



