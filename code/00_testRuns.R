# Script for munging datasets into appropriate forms
# Tim Szewczyk



##--- set up
library(rstan); library(tidyverse); library(sf); library(googlesheets)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
gis.dir <- "../2_gis/data/VD_21781/"
ant.dir <- "../1_opfo/data/"
tax_i <- read_csv("data/tax_i.csv") #%>% filter(sNum<50)
plot_i <- read_csv(paste0(ant.dir, "opfo_envDataProcessed.csv")) %>% 
  arrange(BDM, Plot_id) %>% filter(Plot_id != "020201")
source("code/00_fn.R"); source(paste0(ant.dir, "../code/00_fn.R"))



##--- settings
X_vars <- c("MAT", "Iso", "AP", "npp", "MAT_sq", "lcH")
V_vars <- c("el", "slope", "pop")
U_vars <- c("slope", "pop", "rdLen")



# Testing
mods <- "W" #c("W", "Y", "WY")

pars.keep <- c("sigma_b", "beta", "L_Omega_B", 
               "sigma_a", "alpha", "L_Omega_A",
               "eta", "D",
               "ShannonH", "ShannonH_")
pars.excl <- c("b_std", "B_std", 
               "a_std", "A_std",
               "b", "B", "a", "A",
               "log_lik_lambda_",
               "LAMBDA", "lambda", 
               "E", "E_",
               "LAMBDA_", "lambda_", 
               "p", "p_")


out.ls <- setNames(vector("list", length(mods)), mods)
for(i in seq_along(mods)) {
  out.ls[i] <- stan(file=paste0("code/mods/PPM_mvPhy_", mods[i], "_IJK.stan"),
                    # sample_file=paste0("out/tests/", i), 
                    pars=pars.excl, include=F, chains=3,
                    data=read_rdump("data/stan_data/test_realData.Rdump"), 
                    thin=5, warmup=500, iter=1000)
}

d.ls <- read_rdump("data/stan_data/test_realData.Rdump")
matrix(rstan::summary(out.ls$W, pars="L_Omega_A")$summary[,1], 
       nrow=d.ls$G, byrow=T)
matrix(rstan::summary(out.ls$W, pars="L_Omega_B")$summary[,1], 
       nrow=d.ls$G, byrow=T)
cbind(rstan::summary(out.ls$WY, pars=c("sigma_a", "sigma_b"))$summary[,1])
cbind(rstan::summary(out.ls$WY, pars=c("alpha"))$summary[,1])
cbind(rstan::summary(out.ls$WY, pars=c("beta"))$summary[,1])
cbind(rstan::summary(out.ls$WY, pars=c("eta"))$summary[,1])
cbind(rstan::summary(out.ls$WY, pars=c("D"))$summary[,1])





