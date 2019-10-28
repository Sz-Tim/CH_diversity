data {
  int<lower=0> nCell; // number of grid cells
  int<lower=0> nSpp;  // number of species
  int<lower=0> nGen;  // number of genera
  int<lower=0> R;  // number of cell-scale covariates (incl. intercept)
  int<lower=0> tax_i[nSpp,2];  // species-genus lookup
  int<lower=0> W[nCell,nSpp];  // cell-scale citizen science counts
  int<lower=0> E[nCell];  // cell-scale citizen science total tube counts
  matrix[nCell,R] X_W;  // cell-scale citizen science covariates
}

transformed data {
  int<lower=0> W_Z[nCell,nSpp];  // cell-scale citizen science counts
  for(s in 1:nSpp) {
    for(k in 1:nCell) {
      W_Z[k,s] = W[k,s]>0;
    }
  }
}

parameters {
  matrix[R,nSpp] b;  // environmental covariate species slopes
  matrix[R,nGen] B;  // environmental covariate genus slopes
  vector[R] beta;
  real<lower=0> sigma_B;
  real<lower=0> sigma_b;
  vector<lower=0, upper=1>[nSpp] p;
}
transformed parameters {
  matrix<lower=0>[nCell,nSpp] PSI = inv_logit(X_W * b);
}

model {
  sigma_B ~ normal(0, 1);
  sigma_b ~ normal(0, 1);
  for(r in 1:R) {
    b[r,] ~ normal(B[r,tax_i[,2]], sigma_b);
    B[r,] ~ normal(beta[r], sigma_B);
    beta[r] ~ normal(0, 1);
  }
  for(s in 1:nSpp) {
    for(k in 1:nCell) {
      W_Z[k,s] ~ bernoulli(PSI[k,s]*(1-(1-p[s])^E[k]));
    }
  }
}

