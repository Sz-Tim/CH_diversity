# Script for munging datasets into appropriate forms
# Tim Szewczyk



##--- set up
library(rstan); library(tidyverse); library(foreach); library(doSNOW)
source("code/00_fn.R")
dataset <- c("full", "split")[1]
mod.base <- c("code/mods/%s_Phy_NoTest.stan")

data.f <- paste0("data/stan_data/opfo_", dataset)
out.f <- paste0("out/tests/", dataset, "_")


##--- settings
d.i <- readRDS(paste0(data.f, "_i.rds"))
X_vars <- colnames(d.i$X)[-(1:2)]
V_vars <- colnames(d.i$V)[-(1:2)]
U_vars <- colnames(d.i$U)[-(1:2)]



# Testing
mods <- c("W", "Y", "WY")

pars.keep <- c("sigma_b", "beta",
               "sigma_a", "alpha",
               "b", "B", "a", "A",
               "log_lik_lambda_",
               "eta", "D",
               "Y_rep",
               "lLAMBDA", "lLAMBDA_",
               "llambda", "llambda_", 
               "ShannonH", "ShannonH_")
pars.excl <- c("b_std", "B_std", 
               "a_std", "A_std",
               "Sigma_B", "Sigma_A",
               "L_Omega_B", "L_Omega_A",
               "sigma_B", "sigma_A",
               "L_Omega_a", "L_Omega_a",
               "LAMBDA", "LAMBDA_", 
               "lambda", "lambda_", 
               "E", "E_",
               "p", "p_")


cat("parallelizing\n")
p.c <- makeCluster(3); registerDoSNOW(p.c)
foreach(i=seq_along(mods), .packages="rstan") %dopar% {
  options(mc.cores=8)
  rstan_options(auto_write=TRUE)
  out <- stan(file=sprintf(mod.base, mods[i]),
              pars=pars.excl, include=F, chains=8, save_warmup=F,
              data=read_rdump(paste0(data.f, ".Rdump")), 
              warmup=1000, iter=1100)
  saveRDS(out, paste0(out.f, mods[i], ".rds"))
  saveRDS(As.mcmc.list(out), paste0(out.f, mods[i], "_MCMC.rds"))
  saveRDS(summary(out)$summary, paste0(out.f, mods[i], "_summary.rds"))
}
stopCluster(p.c)
