// full: All regional and local covariates  
// Y: Structured samples only  
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
  int<lower=0> J_; // grid cells for Y_
  int<lower=0> I; // plots 
  int<lower=0> I_; // plots 
  int<lower=0> IJ[I]; // plot to grid cell lookup 
  int<lower=0> IJ_[I_]; // plot to grid cell lookup 
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> L;  // number of plot-scale covariates (no intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  
  // observed data
  matrix<lower=0>[I,S] Y;  // Y counts 
  matrix<lower=0>[I_,S] Y_;  // Y_ counts 
  
  // covariates
  matrix[K+J,R] X;  // W grid cell covariates 
  matrix[J_,R] X_;  // W grid cell covariates 
  matrix[I,L] V;  // Y plot covariates 
  matrix[I_,L] V_;  // Y plot covariates 
  real h;  // Y (plot area)/(cell area)
  
}



transformed data {
  
  matrix[I,R+L] Z = append_col(X[(K+1):(K+J),][IJ,], V);
  matrix[I_,R+L] Z_ = append_col(X_[1:J_,][IJ_,], V_);

}



parameters {
  
  // aggregate slopes
  vector[R+L] beta;
  vector<lower=0>[R+L-1] beta_lam;  // horseshoe prior
  real<lower=0> beta_tau;  // horseshoe prior
  
  // genus effects
  matrix[R+L,G] B_std;
  cholesky_factor_corr[R+L] L_Omega_B;
  vector<lower=0>[R+L] sigma_B; 
  
  // species effects
  matrix[R+L,S] b_std;
  real<lower=0> sigma_b[R+L];
  
  real<lower=0, upper=1> disp_lam;  // implied prior: Uniform(0,1)

}



transformed parameters {
  
  // slopes
  matrix[R+L,S] b;  // species level
  matrix[R+L,S] b_eff;  // species level
  matrix[R+L,G] B_eff;  // genus level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  matrix[I,S] llambda;  // plot level
  
  // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
  B_eff = diag_pre_multiply(sigma_B, L_Omega_B) * B_std;
  
  // b ~ Norm(B, sigma_b)
  for(m in 1:(R+L)) {
    b_eff[m,] = B_eff[m,tax_i[,2]] + b_std[m,] * sigma_b[m]; 
    b[m,] = beta[m] + b_eff[m,]; 
  }
  
  lLAMBDA = X * b[1:R,] - log(h);
  llambda = Z * b;
  
}



model {
  
  // slope priors
  beta[1] ~ normal(-4, 2);
  beta[2:(R+L)] ~ normal(0, beta_tau * beta_lam);
  beta_tau ~ cauchy(0, 2);
  beta_lam ~ student_t(4, 0, 1);
  
  to_vector(B_std) ~ normal(0, 1);
  sigma_B ~ normal(0, 2);
  L_Omega_B ~ lkj_corr_cholesky(2);
  
  to_vector(b_std) ~ normal(0, 1);
  sigma_b ~ normal(0, 2);
  
  // likelihood
  Y ~ GP_log(disp_lam, llambda);
  
}



generated quantities {
  
  matrix[I_,S] log_lik;
  matrix[I_,S] llambda_ = Z_ * b;  // plot level

  // calculated predicted plot-scale prPres
  for(s in 1:S) {
    for(i in 1:I_) {
      log_lik[i,s] = GP_point_log_lpdf(Y_[i,s] | disp_lam, llambda_[i,s]);
    }
  }
  
}



