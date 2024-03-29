functions {  
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
  int<lower=0> K; // grid cells for W (fit)
  int<lower=0> K_; // grid cells for W (predict)
  int<lower=0> J; // grid cells for Y 
  int<lower=0> I; // plots 
  int<lower=0> IJ[I]; // plot to grid cell lookup 
  int<lower=0> S;  // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> L;  // number of plot-scale covariates (no intercept)
  
  // taxonomy
  int<lower=0> tax_i[S,2];  // species-genus lookup
  real<lower=0.5, upper=1.5> D_prior[S]; // all 1's
  
  // observed data
  int<lower=0> W[K,S];  // W counts 
  matrix<lower=0>[I,S] Y;  // Y counts 
  
  // covariates
  matrix[K+J,R] X;  // grid cell covariates 
  matrix[K_,R] X_;  // grid cell covariates (predict)
  matrix[I,L] V;  // plot covariates 
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
  
  // latent variable
  row_vector<lower=-1,upper=1>[S] gamma;
  vector[I] zeta;
  
  // Generalized Poisson dispersion parameter
  real<lower=0, upper=1> disp_lam;  
  
  // bias: W species random effects
  vector<lower=0>[S] D;  // species bias in W
}


transformed parameters {  
  // slopes
  matrix[R+L,S] b;  // species level
  matrix[R+L,G] B;  // genus level
  
  // Lambda/lambda
  matrix[K+J,S] lLAMBDA;  // cell level
  matrix[I,S] llambda;  // plot level
  
  for(m in 1:(R+L)) {
    // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
    // b ~ Norm(B, sigma_b)
    B[m,] = beta[m] + B_std[m,] * diag_pre_multiply(sigma_B[m], L_Omega_B[m]);  
    b[m,] = B[m,tax_i[,2]] + b_std[m,] * sigma_b[m]; 
  }
  lLAMBDA = X*b[1:R,] - log(h);
  llambda = Z*b + zeta*gamma;
}


model {
  // effort and species bias priors
  D ~ normal(D_prior, 1);
  
  // slope priors
  beta[1] ~ normal(-4, 2);
  beta[2:(R+L)] ~ normal(0, 1);
  for(m in 1:(R+L)) {
    b_std[m,] ~ normal(0, 1);
    B_std[m,] ~ normal(0, 1);
    sigma_B[m] ~ normal(0, 2);
    L_Omega_B[m] ~ lkj_corr_cholesky(2);
  }
  sigma_b ~ normal(0, 2);
  
  // latent variables
  gamma ~ normal(0, 1);
  zeta ~ normal(0, 1);
  
  // likelihood
  for(k in 1:K) {
    W[k,] ~ multinomial(softmax( lLAMBDA[k,]' .* D ));
  }
  Y ~ GP_log(disp_lam, llambda);
}
