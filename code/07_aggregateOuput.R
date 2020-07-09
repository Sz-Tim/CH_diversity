# 07_aggregateOutput.R
# Tim Szewczyk
#
# This script aggregates and summarizes the output from the hierarchical models


library(tidyverse); library(rstan); library(sf); theme_set(theme_bw())
library(viridis)
source("code/00_fn.R")


### Full models ----------------------------------------------------------------

agg_vs <- aggregate_output(d.i=readRDS("data/stan_data/vs_no_pred_i.rds"),
                           mods=c("vs_Y", "vs_WY"), 
                           pars_save=c(lLAM="lLAMBDA", llam="llambda",
                                       pP_L="prPresL", 
                                       beta="beta", B="B", b="b", 
                                       sig_b="sigma_b", Sig_B="Sigma_B",
                                       disp="disp_lam", D="D"))

saveRDS(agg_vs, "out/agg_vs.rds")





### Best models ----------------------------------------------------------------

agg_best <- aggregate_output(d.i=readRDS("data/stan_data/pred_i.rds"),
                            mods=c("best_Y", "best_WY"), 
                            pars_save=c(lLAM="lLAMBDA", llam="llambda",
                                        pP_L="prPresL", pP_R="prPres",
                                        H="ShannonH",
                                        beta="beta", B="B", b="b", 
                                        sig_b="sigma_b", Sig_B="Sigma_B",
                                        disp="disp_lam", D="D"))

saveRDS(agg_best, "out/agg_best.rds")




