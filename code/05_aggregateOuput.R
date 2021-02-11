# 07_aggregateOutput.R
# Tim Szewczyk
#
# This script aggregates and summarizes the output from the hierarchical models


library(tidyverse); library(rstan); library(sf); theme_set(theme_bw())
library(viridis)
source("code/00_fn.R")


### Full models ----------------------------------------------------------------
Y_opt <- "Y_4__k-2_L_SoilTSt"
# WY_opt <- "WY_5__k-2_L_Pasture"
WY_opt <- "WY_7__k-2_R_AP"
agg_vs <- aggregate_output(d.f=paste0("data/fwdSearch/", c(Y_opt, WY_opt)),
                           mods=c(Y_opt, WY_opt), 
                           out.dir="out/fwdSearch",
                           pars_save=c("lLAMBDA", "llambda",
                                       "beta", "b", "sigma_b", 
                                       "disp_lam", "D", "log_lik"))

saveRDS(agg_vs$summaries, "out/agg_vs_k-2.rds")
saveRDS(agg_vs$full, "out/agg_vs_k-2_full.rds")





### Null models ----------------------------------------------------------------

agg_null <- aggregate_output(d.f="data/opt/Y_null__opt_var_set",
                             mods="Y_null__",
                             out.dir="out/null/",
                             pars_save=c("lLAMBDA", "lLAMBDA_", "llambda",
                                         "tot_lLAM", "tot_lLAM_", "tot_llam", 
                                         "tot_LAM", "tot_LAM_", "tot_lam", 
                                         "prPres", "prPres_", "prPresL", 
                                         "ShannonH", "ShannonH_", "ShannonH_L",
                                         "Rich", "Rich_", "RichL",
                                         "pred_Y", "pred_Y_", "pred_YL",
                                         "beta", "B", "b", "sigma_b", "Sigma_B",
                                         "disp_lam", "D", 
                                         "log_lik", "log_lik_S", "log_lik_I"))

saveRDS(agg_null$summaries, "out/agg_null_Y.rds")
saveRDS(agg_null$full, "out/agg_null_full_Y.rds")





### Best models ----------------------------------------------------------------

agg_opt <- aggregate_output(d.f="data_orig/opt/WY__opt_var_set",
                            mods="WY",
                            out.dir="out/LV",
                            pars_save=c("lLAMBDA", "lLAMBDA_", "llambda",
                                        "tot_lLAM", "tot_lLAM_", "tot_llam", 
                                        "tot_LAM", "tot_LAM_", "tot_lam", 
                                        "prPres", "prPres_", "prPresL", 
                                        "ShannonH", "ShannonH_", "ShannonH_L",
                                        "Rich", "Rich_", "RichL",
                                        "pred_Y", "pred_Y_", "pred_YL",
                                        "beta", "B", "b", "sigma_b", "Sigma_B",
                                        "disp_lam", "D", 
                                        "log_lik", "log_lik_S", "log_lik_I",
                                        "zeta", "gamma", 
                                        "gamma_sig2", "gamma_Sigma"))

saveRDS(agg_opt$summaries, "out/agg_opt_Y.rds")
saveRDS(agg_opt$full, "out/agg_opt_full_Y.rds")




