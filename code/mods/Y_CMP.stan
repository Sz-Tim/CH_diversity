functions {
  // Conway-Maxwell Poisson distribution from brms
  // https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_com_poisson.stan
  
  // log approximate normalizing constant of the COM poisson distribuion
  // approximation based on doi:10.1007/s10463-017-0629-6
  // Args: see log_Z_com_poisson()
  real log_Z_com_poisson_approx(real log_mu, real nu) {
    real nu_mu = nu * exp(log_mu); 
    real nu2 = nu^2;
    // first 4 terms of the residual series
    real log_sum_resid = log(
      1 + nu_mu^(-1) * (nu2 - 1) / 24 + 
      nu_mu^(-2) * (nu2 - 1) / 1152 * (nu2 + 23) +
      nu_mu^(-3) * (nu2 - 1) / 414720 * (5 * nu2^2 - 298 * nu2 + 11237)
    );
    return nu_mu + log_sum_resid  - 
      ((log(2 * pi()) + log_mu) * (nu - 1) / 2 + log(nu) / 2);
  }
  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // Args:
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real log_Z_com_poisson(real log_mu, real nu) {
    real log_Z;
    real lfac;
    real term;
    real k;
    int M;
    real log_thres;
    if (nu == 1) {
      return exp(log_mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite")
    }
    if (log_mu * nu >= log(1.5) && log_mu >= log(1.5)) {
      return log_Z_com_poisson_approx(log_mu, nu);
    }
    // direct computation of the truncated series
    M = 10000;
    log_thres = log(1e-16);
    // check if the Mth term of the series is small enough
    if (nu * (M * log_mu - lgamma(M + 1)) > log_thres) {
      reject("nu is too close to zero.")
    }
    log_Z = log1p_exp(nu * log_mu);  // first 2 terms of the series
    lfac = 0;
    term = 0;
    k = 2;
    while (term > log_thres) { 
      lfac += log(k);
      term = nu * (k * log_mu - lfac);
      log_Z = log_sum_exp(log_Z, term);
      k += 1;
    }
    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args: 
  //   y: the response value 
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_mu, real nu) {
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return nu * (y * log_mu - lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu);
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
  
  real<lower=0> disp_lam;

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
  
  disp_lam ~ normal(0, 1);
  
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
      Y[i,s] ~ com_poisson_log(llambda[i,s], disp_lam);
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



