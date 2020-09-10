# Forward search

library(tidyverse)
library(rstan)
source("code/00_fn.R")
d.f <- "data/stan_data/vs_20"
fS_out <- "out/fwdSearch/"
fS_dat <- "data/fwdSearch/"

#--- Make full dataset with test/train split
# Run 01_makeDatasets.R with vs=T, test_prop_Y=0.2
# Creates data/stan_data/vs_20* (_ls.rds, _i.rds, .Rdump)




##### Size = 1     -------------------------------------------------------------
mod_size <- 1
#--- Load full dataset, then create a new dataset for each individual variable
make_next_datasets(full_data_base=d.f, out_base=fS_dat, opt_names="R_")


#--- Run 16 models with 1 variable each -- intercept + 1 
# Run 02-1_run_Y_0-7.sh
# Run 02-1_run_Y_8-15.sh
# Run 02-1_run_WY_0-7.sh
# Run 02-1_run_WY_8-15.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
saveRDS("R_", paste0(fS_out, "opt_Y.rds"))
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
saveRDS("R_", paste0(fS_out, "opt_WY.rds"))
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 2     -------------------------------------------------------------
mod_size <- 2
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 15 models with 2 variables each -- optimal + 1 
# Run 02-2_run_Y_0-7.sh
# Run 02-2_run_Y_8-14.sh
# Run 02-2_run_WY_0-7.sh
# Run 02-2_run_WY_8-14.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)





##### Size = 3     -------------------------------------------------------------
mod_size <- 3
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 14 models with 3 variables each -- optimal + 1 
# Run 02-3_run_Y_0-6.sh
# Run 02-3_run_Y_7-13.sh
# Run 02-3_run_WY_0-6.sh
# Run 02-3_run_WY_7-13.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 4     -------------------------------------------------------------
mod_size <- 4
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 13 models with 4 variables each -- optimal + 1 
# Run 02-4_run_Y_0-6.sh
# Run 02-4_run_Y_7-12.sh
# Run 02-4_run_WY_0-6.sh
# Run 02-4_run_WY_7-12.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 5     -------------------------------------------------------------
mod_size <- 5
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 12 models with 5 variables each -- optimal + 1 
# Run 02-5_run_Y_0-5.sh
# Run 02-5_run_Y_6-11.sh
# Run 02-5_run_WY_0-5.sh
# Run 02-5_run_WY_6-11.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 6     -------------------------------------------------------------
mod_size <- 6
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 11 models with 6 variables each -- optimal + 1 
# Run 02-6_run_Y_0-5.sh
# Run 02-6_run_Y_6-10.sh
# Run 02-6_run_WY_0-5.sh
# Run 02-6_run_WY_6-10.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 7     -------------------------------------------------------------
mod_size <- 7
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 10 models with 7 variables each -- optimal + 1 
# Run 02-7_run_Y_0-4.sh
# Run 02-7_run_Y_5-9.sh
# Run 02-7_run_WY_0-4.sh
# Run 02-7_run_WY_5-9.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 8     -------------------------------------------------------------
mod_size <- 8
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 9 models with 8 variables each -- optimal + 1 
# Run 02-8_run_Y_0-4.sh
# Run 02-8_run_Y_5-8.sh
# Run 02-8_run_WY_0-4.sh
# Run 02-8_run_WY_5-8.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 9     -------------------------------------------------------------
mod_size <- 9
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 8 models with 9 variables each -- optimal + 1 
# Run 02-9_run_Y_0-7.sh
# Run 02-9_run_WY_0-7.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 10    -------------------------------------------------------------
mod_size <- 10
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "Y_"),
                   opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size])
make_next_datasets(full_data_base=d.f, out_base=paste0(fS_dat, "WY_"), 
                   opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size])


#--- Run 7 models with 10 variables each -- optimal + 1 
# Run 02-10_run_Y_0-6.sh
# Run 02-10_run_WY_0-6.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp[1:10,]
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)
