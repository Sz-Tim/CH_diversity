
data {
  int<lower=0> J; // grid cells for Y (train)
  int<lower=0> I; // plots (train)
  int<lower=0> IJ[I]; // plot to grid cell lookup (train)
  int<lower=0> S;  // number of species
  int P; // number of parameters
  int<lower=0> Y[I*S];  // Y counts (train)
  matrix[I*S,P] Z;
}

parameters {
  vector[P] beta;
  matrix[P,S] b;
  real<lower=0> sig_B;
}

transformed parameters {
  vector[I*S] llambda;
  
  {
    int s_i;
    for(s in 1:S) {
      s_i = (s-1)*I + 1;
      llambda[s_i:(s_i+(I-1))] = Z[s_i:(s_i+(I-1)),] * b[,s];
    }
  }
}

model {
  beta ~ normal(0, 1);
  for(p in 1:P) {
    b[p,] ~ normal(beta[p], sig_B);
  }
  Y ~ poisson_log(llambda);
}

