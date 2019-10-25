data {
  int<lower=0> nCell; // number of grid cells
  int<lower=0> nSpp;  // number of species
  int<lower=0> nGen;  // number of genera
  int<lower=0> K;  // number of cell-scale covariates
  int<lower=0> tax_i[nSpp,2];  // species-genus lookup
  int<lower=0> W[nCell,nSpp];  // cell-scale citizen science counts
  int<lower=0> E[nCell];  // cell-scale citizen science total tube counts
  matrix[nCell,K+1] X_W;  // cell-scale citizen science covariates
}

parameters {
  matrix[K+1,nSpp] b;  // environmental covariate species slopes
  matrix[K+1,nGen] B;  // environmental covariate genus slopes
  vector[K+1] beta;
  real<lower=0> sigma_B;
  real<lower=0> sigma_b;
}
transformed parameters {
  matrix<lower=0>[nCell,nSpp] LAMBDA = exp(X_W * b);
}

model {
  sigma_B ~ normal(0, 1);
  sigma_b ~ normal(0, 1);
  for(k in 1:(K+1)) {
    b[k,] ~ normal(B[k,tax_i[,2]], sigma_b);
    B[k,] ~ normal(beta[k], sigma_B);
    beta[k] ~ normal(0, 1);
  }
  for(s in 1:nSpp) {
    for(j in 1:nCell) {
      W[j,s] ~ poisson(LAMBDA[j,s]*E[j]);
    }
  }
}

