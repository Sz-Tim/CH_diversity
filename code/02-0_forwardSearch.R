# Forward search

library(tidyverse)
library(rstan)
source("code/00_fn.R")
d.f <- "data/stan_data/vs_20"

#--- Make full dataset with test/train split
# Run 01_makeDatasets.R with vs=T, test_prop_Y=0.2
# Creates data/stan_data/vs_20* (_ls.rds, _i.rds, .Rdump)




##### Size = 1     -------------------------------------------------------------
#--- Load full dataset, then create a new dataset for each individual variable
make_next_datasets(full_data_base=d.f, out_base="data/fwdSearch/", opt_names="R_")


#--- Run 16 models with 1 variable each -- intercept + 1 
# Run 02-1_run_Y_0-7.sh
# Run 02-1_run_Y_8-15.sh
# Run 02-1_run_WY_0-7.sh
# Run 02-1_run_WY_8-15.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=1)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=1)
WY_loo$comp[1:10,]




##### Size = 2     -------------------------------------------------------------
#--- Make new datasets with selected variable + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", ""))


#--- Run 15 models with 2 variables each -- optimal + 1 
# Run 02-2_run_Y_0-7.sh
# Run 02-2_run_Y_8-14.sh
# Run 02-2_run_WY_0-7.sh
# Run 02-2_run_WY_8-14.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=2)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=2)
WY_loo$comp[1:10,]




##### Size = 3     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", ""))


#--- Run 14 models with 3 variables each -- optimal + 1 
# Run 02-3_run_Y_0-6.sh
# Run 02-3_run_Y_7-13.sh
# Run 02-3_run_WY_0-6.sh
# Run 02-3_run_WY_7-13.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=3)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=3)
WY_loo$comp[1:10,]




##### Size = 4     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", ""))


#--- Run 13 models with 4 variables each -- optimal + 1 
# Run 02-4_run_Y_0-6.sh
# Run 02-4_run_Y_7-12.sh
# Run 02-4_run_WY_0-6.sh
# Run 02-4_run_WY_7-12.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=4)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=4)
WY_loo$comp[1:10,]




##### Size = 5     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", ""))


#--- Run 12 models with 5 variables each -- optimal + 1 
# Run 02-5_run_Y_0-5.sh
# Run 02-5_run_Y_6-11.sh
# Run 02-5_run_WY_0-5.sh
# Run 02-5_run_WY_6-11.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=5)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=5)
WY_loo$comp[1:10,]




##### Size = 6     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", "", ""))


#--- Run 11 models with 6 variables each -- optimal + 1 
# Run 02-6_run_Y_0-5.sh
# Run 02-6_run_Y_6-10.sh
# Run 02-6_run_WY_0-5.sh
# Run 02-6_run_WY_6-10.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=6)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=6)
WY_loo$comp[1:10,]




##### Size = 7     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", "", "", ""))


#--- Run 10 models with 7 variables each -- optimal + 1 
# Run 02-7_run_Y_0-4.sh
# Run 02-7_run_Y_5-9.sh
# Run 02-7_run_WY_0-4.sh
# Run 02-7_run_WY_5-9.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=7)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=7)
WY_loo$comp[1:10,]




##### Size = 8     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", "", "", "", ""))


#--- Run 9 models with 8 variables each -- optimal + 1 
# Run 02-8_run_Y_0-4.sh
# Run 02-8_run_Y_5-8.sh
# Run 02-8_run_WY_0-4.sh
# Run 02-8_run_WY_5-8.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=8)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=8)
WY_loo$comp[1:10,]




##### Size = 9     -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", "", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", "", "", "", "", ""))


#--- Run 8 models with 9 variables each -- optimal + 1 
# Run 02-9_run_Y_0-7.sh
# Run 02-9_run_WY_0-7.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=9)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=9)
WY_loo$comp[1:10,]




##### Size = 10    -------------------------------------------------------------
#--- Make new datasets with selected variables + each other individually
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/Y_", 
                   opt_names=c("R_", "", "", "", "", "", "", "", "", ""))
make_next_datasets(full_data_base=d.f, 
                   out_base="data/fwdSearch/WY_", 
                   opt_names=c("R_", "", "", "", "", "", "", "", "", ""))


#--- Run 7 models with 10 variables each -- optimal + 1 
# Run 02-10_run_Y_0-6.sh
# Run 02-10_run_WY_0-6.sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir="out/fwdSearch", mod="Y", mod_size=10)
Y_loo$comp[1:10,]

WY_loo <- compare_models(fit_dir="out/fwdSearch", mod="WY", mod_size=10)
WY_loo$comp[1:10,]
