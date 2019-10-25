data {
  int<lower=0> nPlot; // number of structured sampling soil plots
  int<lower=0> nSite; // number of structured sampling sites 
  int<lower=0> nCell; // number of grid cells
  int<lower=0> nSpp;  // number of species
  int<lower=0> K;  // number of cell-scale covariates
  int<lower=0> L;  // number of plot-scale covariates
  matrix[nPlot,nSpp] y;  // plot-scale structured sampling counts
  matrix[nSite,nSpp] Y;  // plot-scale structured sampling counts
  matrix[nCell,nSpp] W;  // cell-scale citizen science counts
  matrix[nCell,3] E;  // cell-scale citizen science effort covariates
  matrix[nSite,K] X_Y;  // cell-scale citizen science covariates
  matrix[nCell,K] X_W;  // cell-scale citizen science covariates
  matrix[nPlot,K] V;  // plot-scale covariates
}

parameters {
  matrix<lower=0>[nCell,nSpp] LAMBDA_W;  //
  matrix<lower=0>[nSite,nSpp] LAMBDA_Y;  //
}

model {
  
}

