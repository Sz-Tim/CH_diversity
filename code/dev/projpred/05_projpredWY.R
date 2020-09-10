# 05_projpredWY.R
# Tim Szewczyk
#
# This script runs projection predictive inference using the structure samples



# At the time this was run, mixed effects models are only possible on the 
# 'develop' branch of the projpred R package. Install if necessary:
#
# devtools::install_github('stan-dev/projpred', "develop", build_vignettes=F, 
#                          force=TRUE)



library(rstan); library(rstanarm); library(projpred); library(tidyverse)
source("code/00_fn.R")
d.ls <- readRDS("data/stan_data/vs_no_pred_ls.rds")



# load model outputs
fit_shell <- readRDS("out/shell.rds")
fit_WY <- read_stan_csv(dir("out", "vs_WY_altRE_G", full.names=T))



# update shell with Y posteriors
fit_update <- update_rstanarm_shell(fit_shell, fit_WY, d.ls)
saveRDS(fit_update, "out/fit_update_WY.rds")



# select variables
vs <- varsel(fit_update, verbose=T, nterms_max=10, ndraws=1)
saveRDS(vs, "out/varsel_WY.rds")
