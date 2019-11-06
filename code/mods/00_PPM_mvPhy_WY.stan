data {
  int<lower=0> nCell_W; // number of grid cells
  int<lower=0> nCell_Y; // number of grid cells with data for Y
  int<lower=0> nSpp;  // number of species
  int<lower=0> nGen;  // number of species
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> tax_i[nSpp,2];  // species-genus lookup
  int<lower=0> W[nCell_W,nSpp];  // cell-scale citizen science counts
  int<lower=0> Y[nCell_Y,nSpp];  // cell-scale citizen science counts
  matrix[nCell_W,R] X_W;  // cell-scale citizen science covariates
  matrix[nCell_Y,R] X_Y;  // cell-scale citizen science covariates
  real effort_prSoil;
}

parameters {
  matrix[R,nSpp] b;  // environmental covariate species slopes
  real<lower=0> sigma_b;
  matrix[R,nGen] B;  // environmental covariate genus slopes
  vector[R] beta;
  cholesky_factor_corr[nGen] L_Omega;
  vector[R] eta;
}

transformed parameters {
  matrix<lower=0>[nCell_W,nSpp] LAMBDA_W = exp(X_W * b);
  matrix<lower=0>[nCell_Y,nSpp] LAMBDA_Y = exp(X_Y * b);
  vector<lower=0, upper=1>[nCell_W] E = inv_logit(X_W * eta);
}

model {
  L_Omega ~ lkj_corr_cholesky(2);
  sigma_b ~ normal(0, 1);
  beta[1] ~ normal(12, 1);
  eta[1] ~ normal(-10, 1);
  for(r in 2:R) {
    beta[r] ~ normal(0, 1);
    eta[r] ~ normal(0, 1);
  }
  for(r in 1:R) {
    b[r,] ~ normal(B[r,tax_i[,2]], sigma_b);
    B[r,] ~ multi_normal_cholesky(rep_vector(beta[r], nGen), L_Omega);
  }
  for(s in 1:nSpp) {
    W[,s] ~ poisson(LAMBDA_W[,s].*E);
    Y[,s] ~ poisson(LAMBDA_Y[,s]*effort_prSoil);
  }
}

// generated quantities {
//   int<lower=0,upper=1> PA_W[nCell_W,nSpp];
//   int<lower=0,upper=1> PA_Y[nCell_Y,nSpp];
//   int<lower=0,upper=nSpp> Richness_W[nCell_W];
//   int<lower=0,upper=nSpp> Richness_Y[nCell_Y];
//   
//   for(k in 1:nCell_W) {
//     for(s in 1:nSpp) {
//       PA_W[k,s] = LAMBDA_W[k,s]>0.05;
//     }
//     Richness_W[k] = sum(PA_W[k,]);
//   }
//   for(k in 1:nCell_Y) {
//     for(s in 1:nSpp) {
//       PA_Y[k,s] = LAMBDA_Y[k,s]>0.05;
//     }
//     Richness_Y[k] = sum(PA_Y[k,]);
//   }
// }



