# Script for munging datasets into appropriate forms
# Tim Szewczyk



##--- set up
library(rstan); library(tidyverse); library(foreach); library(doSNOW)
source("code/00_fn.R")



##--- settings
X_vars <- c("MAT", "lcH", "npp", "MAT_sq")
V_vars <- c("SoilTSt", "CnpyOpn", "CnpyMxd", "VegTot")
U_vars <- c("pop", "rdLen")



# Testing
mods <- c("W", "Y", "WY")

pars.keep <- c("sigma_b", "beta", "Sigma_B", 
               "sigma_a", "alpha", "Sigma_A",
               "b", "B", "a", "A",
               "log_lik_lambda_",
               "eta", "D",
               "lLAMBDA", "lLAMBDA_",
               "ShannonH", "ShannonH_")
pars.excl <- c("b_std", "B_std", 
               "a_std", "A_std",
               "L_Omega_B", "L_Omega_A",
               "sigma_B", "sigma_A",
               "L_Omega_a", "L_Omega_a",
               "LAMBDA", "LAMBDA_", 
               "lambda", "lambda_", 
               "llambda", "llambda_", 
               "E", "E_",
               "p", "p_")


p.c <- makeCluster(3); registerDoSNOW(p.c)
out.ls <- foreach(i=seq_along(mods), .packages="rstan") %dopar% {
  options(mc.cores=1)
  rstan_options(auto_write=TRUE)
  out <- stan(file=paste0("code/mods/PPM_mvPhy_", mods[i], "_IJK.stan"),
       sample_file=paste0("out/tests/", mods[i]),
       pars=pars.excl, include=F, chains=1,
       data=read_rdump("data/stan_data/test_realData.Rdump"), 
       warmup=50, iter=100)
  saveRDS(out, paste0("out/tests/test_", mods[i], ".rds"))
}
stopCluster(p.c)

# out.ls <- setNames(vector("list", length(mods)), mods)
# for(i in seq_along(mods)) {
#   out.ls[i] <- stan(file=paste0("code/mods/PPM_mvPhy_", mods[i], "_IJK.stan"),
#                     sample_file=paste0("out/tests/", mods[i]),
#                     pars=pars.excl, include=F, chains=2,
#                     data=read_rdump("data/stan_data/test_realData.Rdump"), 
#                     thin=5, warmup=100, iter=200)
# }

# d.ls <- read_rdump("data/stan_data/test_realData.Rdump")
# matrix(rstan::summary(out.ls$W, pars="L_Omega_A")$summary[,1], 
#        nrow=d.ls$G, byrow=T)
# matrix(rstan::summary(out.ls$W, pars="L_Omega_B")$summary[,1], 
#        nrow=d.ls$G, byrow=T)
# cbind(rstan::summary(out.ls$WY, pars=c("sigma_a", "sigma_b"))$summary[,1])
# cbind(rstan::summary(out.ls$WY, pars=c("alpha"))$summary[,1])
# cbind(rstan::summary(out.ls$WY, pars=c("beta"))$summary[,1])
# cbind(rstan::summary(out.ls$WY, pars=c("eta"))$summary[,1])
# cbind(rstan::summary(out.ls$WY, pars=c("D"))$summary[,1])





