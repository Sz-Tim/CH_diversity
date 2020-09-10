# 04_projpredY.R
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
d.ls <- readRDS("data/stan_data/vs_80_ls.rds")



# load model outputs
fit_shell <- readRDS("out/vs_80/shell_80.rds")
fit_Y <- read_stan_csv(dir("out/vs_80", "vs_Y_altRE_G", full.names=T))



# update shell with Y posteriors
fit_update <- update_rstanarm_shell(fit_shell, fit_Y, d.ls)
saveRDS(fit_update, "out/vs_80/fit_update_Y.rds")



# select variables
vs <- varsel(fit_update, verbose=T, nterms_max=10, ndraws=1)
saveRDS(vs, "out/vs_80/varsel_Y.rds")
