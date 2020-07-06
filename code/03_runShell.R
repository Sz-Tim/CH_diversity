# This script creates an rstanarm model to be used as a shell for projpred


library(rstanarm); library(tidyverse)
d.ls <- readRDS("data/stan_data/vs_no_pred_ls.rds")

# sampling parameters
nWarm=2000
nSamp=20
nChain=12


# make dataframe
Z <- cbind(d.ls$X[d.ls$IJ,], d.ls$V)
colnames(Z) <- paste0(rep(c("R_", "L_"), c(d.ls$R, d.ls$L)), colnames(Z))
Z_r <- do.call('rbind', map(1:d.ls$S, ~Z))
Y <- c(d.ls$Y)
stan_df <- as.data.frame(cbind(Y, Z_r, sp=rep(1:d.ls$S, each=d.ls$I)))


# run glmm with species random effects
fit <- stan_glmer(
  as.formula(
    paste("Y ~ ", paste(colnames(Z_r)[-1], collapse="+"), "+",
          paste("(1+", paste(colnames(Z_r)[-1], collapse="+"), "|sp)"))),
  data=stan_df,
  family=poisson("log"),
  prior=hs(),
  chains=nChain, cores=nChain, 
  warmup=nWarm, iter=nWarm+nSamp, 
  refresh=1)

# store output
saveRDS(fit, "out/shell.rds")
