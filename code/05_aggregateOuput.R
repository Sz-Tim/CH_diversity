# 07_aggregateOutput.R
# Tim Szewczyk
#
# This script aggregates and summarizes the output from the hierarchical models


library(tidyverse); library(rstan); library(sf); theme_set(theme_bw())
library(viridis)
source("code/00_fn.R")



### Null models ----------------------------------------------------------------

LV <- c("cov", "LV")[2]
mod <- c("Y", "WY")[2]

agg_null <- aggregate_output(d.f=paste0("data/opt/", LV, "_", mod, "_null__opt_var_set"),
                             mods=paste0(LV, "_", mod, "_null__"),
                             out.dir="out/null/",
                             pars_save=c("lLAMBDA", "lLAMBDA_", "llambda",
                                         paste0("tot_", c("lLAM", "lLAM_", "llam")),
                                         paste0("tot_", c("LAM", "LAM_", "lam")),
                                         paste0("prPres", c("", "_", "L")),
                                         paste0("ShannonH", c("", "_", "_L")),
                                         paste0("Rich", c("", "_", "L")),
                                         paste0("pred_Y", c("", "_", "L")),
                                         paste0("log_lik", c("", "_S", "_I")),
                                         "beta", "B", "b", "sigma_b", "Sigma_B",
                                         "disp_lam", "D", "zeta",
                                         paste0("gamma", c("", "_sig2", "_Sigma"))))

saveRDS(agg_null$summaries, paste0("out/agg_null_", LV, "_", mod, ".rds"))
saveRDS(agg_null$full, paste0("out/full_null_", LV, "_", mod, ".rds"))






### Best models ----------------------------------------------------------------

LV <- c("cov", "LV")[2]
mod <- c("Y", "WY")[2]

agg_opt <- aggregate_output(d.f=paste0("data/opt/", LV, "_", mod, "__opt_var_set"),
                            mods=paste0(LV, "_", mod, "__"),
                            out.dir="out/opt/",
                            pars_save=c("lLAMBDA", "lLAMBDA_", "llambda",
                                        paste0("tot_", c("lLAM", "lLAM_", "llam")),
                                        paste0("tot_", c("LAM", "LAM_", "lam")),
                                        paste0("prPres", c("", "_", "L")),
                                        paste0("ShannonH", c("", "_", "_L")),
                                        paste0("Rich", c("", "_", "L")),
                                        paste0("pred_Y", c("", "_", "L")),
                                        paste0("log_lik", c("", "_S", "_I")),
                                        "beta", "B", "b", "sigma_b", "Sigma_B",
                                        "disp_lam", "D", "zeta",
                                        paste0("gamma", c("", "_sig2", "_Sigma"))))

saveRDS(agg_opt$summaries, paste0("out/opt/agg_", LV, "_", mod, ".rds"))
saveRDS(agg_opt$full, paste0("out/opt/full_", LV, "_", mod, ".rds"))




