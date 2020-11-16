# 07_aggregateOutput.R
# Tim Szewczyk
#
# This script aggregates and summarizes the output from the hierarchical models


library(tidyverse); library(rstan); library(sf); theme_set(theme_bw())
library(viridis)
source("code/00_fn.R")


### Full models ----------------------------------------------------------------
Y_opt <- "vs_Y_altRE_G"
WY_opt <- "vs_WY_altRE_G"
# Y_opt <- "Y_7__L_aspctN"
# WY_opt <- "WY_6__R_lcH"
agg_vs <- aggregate_output(d.f=paste0("data/stan_data/", c("vs_no_pred", "vs_no_pred")),#paste0("data/fwdSearch/", c(Y_opt, WY_opt)),
                           mods=c(Y_opt, WY_opt), 
                           out.dir="out/vs_all_species",#"out/fwdSearch",
                           pars_save=c(lLAM="lLAMBDA", llam="llambda",
                                       pP_L="prPresL", 
                                       beta="beta", B="B", b="b", 
                                       sig_b="sigma_b", Sig_B="Sigma_B",
                                       disp="disp_lam", D="D", 
                                       log_lik="log_lik"))

saveRDS(agg_vs$summaries, "out/agg_vs_OLD.rds")
saveRDS(agg_vs$full, "out/agg_vs_full_OLD.rds")





### Best models ----------------------------------------------------------------

# WAIT!! RUN THIS MANUALLY WITHIN THE FUNCTION IN CASE PARNAMES FAILS!
agg_opt <- aggregate_output(d.f=paste0("data/opt/", c("WY"), "__opt_var_set"),
                             mods=c("WY"),
                             out.dir="out/opt",
                             pars_save=c(llam="llambda",
                                         lLAM="lLAMBDA", lLAM_="lLAMBDA_", 
                                         pP_L="prPresL", 
                                         pP_R="prPres", pP_R_="prPres_",
                                         H="ShannonH", H_="ShannonH_",
                                         S_R="Rich", S_R_="Rich_", S_L="RichL",
                                         pred_R="pred_Y", pred_R_="pred_Y_", 
                                         pred_L="pred_YL",
                                         beta="beta", B="B", b="b", 
                                         sig_b="sigma_b", Sig_B="Sigma_B",
                                         disp="disp_lam", D="D"))

saveRDS(agg_opt$summaries, "out/agg_opt_WY.rds")
saveRDS(agg_opt$full, "out/agg_opt_full_WY.rds")




