// best: Optimal regional and local covariates  
// WY: Structured samples AND citizen science samples  
// output: Predicts whole canton at 1 km^2

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
  real<lower=0.5, upper=1.5> D_prior[S];
  
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
  
  matrix<lower=-1,upper=1>[2,S] LV_b;
  vector[I] LV_L;
  vector[K+J] LV_R;
  
  real<lower=0, upper=1> disp_lam;  // implied prior: Uniform(0,1)
  
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
  
  // cell level
  for(m in 1:(R+L)) {
    // B ~ mvNorm(beta, sigma_B * L_Omega_B * sigma_B);  
    // b ~ Norm(B, sigma_b)
    B[m,] = beta[m] + B_std[m,] * diag_pre_multiply(sigma_B[m], L_Omega_B[m]);  
    b[m,] = B[m,tax_i[,2]] + b_std[m,] * sigma_b[m]; 
  }
  lLAMBDA = X*b[1:R,] - log(h) + LV_R*LV_b[2,];
  llambda = Z*b + LV_L*LV_b[1,];
  
}



model {
  
  // effort and species bias priors
  D ~ normal(D_prior, 1);
  
  // cell level priors
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
  to_vector(LV_b) ~ normal(0, 1);
  LV_L ~ normal(0, 1);
  LV_R ~ normal(0, 1);
  
  // likelihood
  for(k in 1:K) {
    W[k,] ~ multinomial(softmax( lLAMBDA[k,]' .* D ));
  }
  Y ~ GP_log(disp_lam, llambda);
  
}



generated quantities {
  
  matrix[K_,S] lLAMBDA_= X_ * b[1:R,] - log(h);
  int pred_YL[I,S];
  int pred_Y[K+J,S];
  int pred_Y_[K_,S];
  matrix[I,S] prPresL;
  matrix[K+J,S] prPres;
  matrix[K_,S] prPres_;
  vector[I] ShannonH_L;
  vector[K+J] ShannonH;
  vector[K_] ShannonH_;
  int RichL[I];
  int Rich[K+J];
  int Rich_[K_];
  vector[I] tot_llam;
  vector[K+J] tot_lLAM;
  vector[K_] tot_lLAM_;
  vector[I] tot_lam;
  vector[K+J] tot_LAM;
  vector[K_] tot_LAM_;
  matrix[G,G] Sigma_B[R+L];
  real log_lik;
  matrix[2,S] LV_sig2;
  matrix[S,S] LV_Sigma[2];
  
  // calculate plot-scale quantities:
  // tot_llam, tot_lam, prPresL, RichL, ShannonH_L
  {
    // temporary variables
    matrix[I,S] p_L;
    matrix[I,S] lambda = exp(llambda);
    matrix[I,S] log_lik_is;
    
    for(i in 1:I) {
      tot_llam[i] = sum(llambda[i,]);
      tot_lam[i] = sum(exp(llambda[i,]));
      p_L[i,] = lambda[i,] / tot_lam[i];
      for(s in 1:S) {
        log_lik_is[i,s] = GP_point_log_lpdf(Y[i,s] | disp_lam, llambda[i,s]);
        prPresL[i,s] = 1 - exp(GP_point_log_lpdf(0 | disp_lam, llambda[i,s]));
        pred_YL[i,s] = bernoulli_rng(prPresL[i,s]);
      }
      RichL[i] = sum(pred_YL[i,]);
    }
    ShannonH_L = - rows_dot_product(p_L, log(p_L));
    log_lik = sum(log_lik_is);
  }
  
  // calculate cell-scale quantities:
  // prPres, tot_lLAM, tot_LAM, pred_Y, Rich, ShannonH
  {
    // temporary variables
    matrix[K+J,S] p;
    matrix[K_,S] p_;
    matrix[K+J,S] LAMBDA = exp(lLAMBDA);
    matrix[K_,S] LAMBDA_=exp(lLAMBDA_);
    prPres = 1 - exp(-LAMBDA);
    prPres_ = 1 - exp(-LAMBDA_);
    
    for(i in 1:(K+J)) {
      tot_lLAM[i] = sum(lLAMBDA[i]);
      tot_LAM[i] = sum(LAMBDA[i]);
      p[i,] = LAMBDA[i,] / tot_LAM[i];
      for(s in 1:S) {
        pred_Y[i,s] = prPres[i,s] > 0.95;
      }
      Rich[i] = sum(pred_Y[i,]);
    }
    for(i in 1:K_) {
      tot_lLAM_[i] = sum(lLAMBDA_[i]);
      tot_LAM_[i] = sum(LAMBDA_[i]);
      p_[i,] = LAMBDA_[i,] / tot_LAM_[i];
      for(s in 1:S) {
        pred_Y_[i,s] = prPres_[i,s] > 0.95;
      }
      Rich_[i] = sum(pred_Y_[i,]);
    }
    ShannonH = - rows_dot_product(p, log(p));
    ShannonH_ = - rows_dot_product(p_, log(p_)); 
  }
  
  // compose correlation matrices = sigma * LL' * sigma
  for(m in 1:(R+L)) {
    Sigma_B[m] = quad_form_diag(multiply_lower_tri_self_transpose(L_Omega_B[m]),
                                sigma_B[m]);
  }
  LV_sig2 = 1 - (LV_b .* LV_b);
  for(l in 1:2) {
    LV_Sigma[l] = LV_b[l,]' * LV_b[l,] + diag_matrix(LV_sig2[l,]');
  }
  
}



