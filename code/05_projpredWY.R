# This script runs projection predictive inference using the joint model



# At the time this was run, mixed effects models are only possible on the 
# 'develop' branch of the projpred R package. Install if necessary:
# 
# devtools::install_github('stan-dev/projpred', "develop", build_vignettes=F, 
#                          force=TRUE)



library(rstan); library(rstanarm); library(projpred); library(tidyverse)
source("code/00_fn.R")
d.ls <- readRDS("data/stan_data/vs_no_pred_ls.rds")



# load model outputs
fit_shell <- readRDS("out/fit_shell.rds")
fit_WY <- read_stan_csv(dir("out/", paste0("vs_WY"), full.names=T))



# update shell with Y posteriors
fit_update <- subset_vs_data(fit_shell, fit_WY, d.ls)



# select variables
vs <- varsel(fit_update, verbose=T)
saveRDS(vs, "out/varsel_WY.rds")
