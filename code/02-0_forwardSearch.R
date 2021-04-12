# Forward search

library(tidyverse)
library(rstan)
source("code/00_fn.R")
d.dir <- "data/stan_data/"
d.f <- "cv_k_"
fS_out <- "/VOLUMES/Rocinante/opfo_div_cluster/fwd_ll/"
fS_dat <- "data/fwdSearch/"
n_folds <- length(dir(d.dir, paste0(d.f, ".*Rdump")))

#--- Make full dataset with test/train split
# Run code/01_makeDatasets.R with set_type="cv"
# Creates data/stan_data/cv_k_[1-4]* (_ls.rds, _i.rds, .Rdump)




##### Size = 0     -------------------------------------------------------------
mod_size <- 0
#--- Load full dataset, then create a new dataset for each individual variable
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x),
                        out_base=fS_dat, opt_names=NULL, type="cv", 
                        mod_size=mod_size))



#--- Run models with only intercept
# code/02-slurm/fwd_0_cov.sh

cov_Y_loo <- compare_models(fit_dir=fS_out, mod="cov_Y", mod_size=mod_size, type="cv")
saveRDS("R_", paste0(fS_out, "opt_cov_Y.rds"))

LV_Y_loo <- compare_models(fit_dir=fS_out, mod="LV_Y", mod_size=mod_size, type="cv")
saveRDS("R_", paste0(fS_out, "opt_LV_Y.rds"))

cov_WY_loo <- compare_models(fit_dir=fS_out, mod="cov_WY", mod_size=mod_size, type="cv")
saveRDS("R_", paste0(fS_out, "opt_cov_WY.rds"))

LV_WY_loo <- compare_models(fit_dir=fS_out, mod="LV_WY", mod_size=mod_size, type="cv")
saveRDS("R_", paste0(fS_out, "opt_LV_WY.rds"))




##### Size = 1     -------------------------------------------------------------
mod_size <- 1
#--- Load full dataset, then create a new dataset for each individual variable
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x),
                        out_base=fS_dat, opt_names="R_", type="cv", 
                        mod_size=mod_size))


#--- Run 14 models with 1 variable each -- intercept + 1 
# code/02-slurm/fwd_1_cov.sh
# code/02-slurm/fwd_1_LV.sh

#--- Use loo to compare models based on ability to predict test subset
cov_Y_loo <- compare_models(fS_out, "cov_Y", mod_size, type="cv")
update_opt_vars(fS_out, "cov_Y", mod_size, cov_Y_loo$comp)
cov_Y_loo$comp

LV_Y_loo <- compare_models(fS_out, "LV_Y", mod_size, type="cv")
update_opt_vars(fS_out, "LV_Y", mod_size, LV_Y_loo$comp)
LV_Y_loo$comp

cov_WY_loo <- compare_models(fS_out, "cov_WY", mod_size, type="cv")
update_opt_vars(fS_out, "cov_WY", mod_size, cov_WY_loo$comp)
cov_WY_loo$comp

LV_WY_loo <- compare_models(fS_out, "LV_WY", mod_size, type="cv")
update_opt_vars(fS_out, "LV_WY", mod_size, LV_WY_loo$comp)
LV_WY_loo$comp




##### Size = 2     -------------------------------------------------------------
mod_size <- 2
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "cov_Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_cov_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "LV_Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_LV_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "cov_WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_cov_WY.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "LV_WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_LV_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 13 models with 2 variables each -- optimal + 1 
# code/02-slurm/fwd_2_cov.sh
# code/02-slurm/fwd_2_LV.sh

#--- Use loo to compare models based on ability to predict test subset
cov_Y_loo <- compare_models(fS_out, "cov_Y", mod_size, type="cv")
update_opt_vars(fS_out, "cov_Y", mod_size, cov_Y_loo$comp)
cov_Y_loo$comp

LV_Y_loo <- compare_models(fS_out, "LV_Y", mod_size, type="cv")
update_opt_vars(fS_out, "LV_Y", mod_size, LV_Y_loo$comp)
LV_Y_loo$comp

cov_WY_loo <- compare_models(fS_out, "cov_WY", mod_size, type="cv")
update_opt_vars(fS_out, "cov_WY", mod_size, cov_WY_loo$comp)
cov_WY_loo$comp

LV_WY_loo <- compare_models(fS_out, "LV_WY", mod_size, type="cv")
update_opt_vars(fS_out, "LV_WY", mod_size, LV_WY_loo$comp)
LV_WY_loo$comp





##### Size = 3     -------------------------------------------------------------
mod_size <- 3
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 12 models with 3 variables each -- optimal + 1 
# code/02-3/02-3-[1-6]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 4     -------------------------------------------------------------
mod_size <- 4
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 11 models with 4 variables each -- optimal + 1 
# code/02-4/02-4-[1-6]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 5     -------------------------------------------------------------
mod_size <- 5
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 10 models with 5 variables each -- optimal + 1 
# code/02-5/02-5-[1-5]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 6     -------------------------------------------------------------
mod_size <- 6
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 9 models with 6 variables each -- optimal + 1 
# code/02-6/02-6-[1-5]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 7     -------------------------------------------------------------
mod_size <- 7
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 8 models with 7 variables each -- optimal + 1 
# code/02-7/02-7-[1-4]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 8     -------------------------------------------------------------
mod_size <- 8
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 7 models with 8 variables each -- optimal + 1 
# code/02-8/02-8-[1-4]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 9     -------------------------------------------------------------
mod_size <- 9
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 6 models with 9 variables each -- optimal + 1 
# code/02-9/02-9-[1-3]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 10     ------------------------------------------------------------
mod_size <- 10
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 5 models with 10 variables each -- optimal + 1 
# code/02-10/02-10-[1-3]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 11    ------------------------------------------------------------
mod_size <- 11
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 4 models with 11 variables each -- optimal + 1 
# code/02-11/02-11-[1-2]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 12     ------------------------------------------------------------
mod_size <- 12
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 3 models with 12 variables each -- optimal + 1 
# code/02-12/02-12-[1-2]_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 13     ------------------------------------------------------------
mod_size <- 13
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 2 models with 13 variables each -- optimal + 1 
# code/02-13/02-13-1_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)




##### Size = 14     ------------------------------------------------------------
mod_size <- 14
#--- Make new datasets with selected variable + each other individually
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "Y_"),
                        opt_names=readRDS(paste0(fS_out, "opt_Y.rds"))[1:mod_size],
                        mod_size=mod_size))
map(1:n_folds, 
    ~make_next_datasets(full_data_base=paste0(d.dir, d.f, .x), type="cv", 
                        out_base=paste0(fS_dat, "WY_"),
                        opt_names=readRDS(paste0(fS_out, "opt_WY.rds"))[1:mod_size],
                        mod_size=mod_size))


#--- Run 1 models with 14 variables -- optimal + 1 
# code/02-14/02-14-1_[WY-Y].sh


#--- Use loo to compare models based on ability to predict test subset
Y_loo <- compare_models(fit_dir=fS_out, mod="Y", mod_size=mod_size)
Y_loo$comp
update_opt_vars(out_dir=fS_out, mod="Y", mod_size=mod_size, comp=Y_loo$comp)

WY_loo <- compare_models(fit_dir=fS_out, mod="WY", mod_size=mod_size)
WY_loo$comp
update_opt_vars(out_dir=fS_out, mod="WY", mod_size=mod_size, comp=WY_loo$comp)















##### PERFORMANCE PLOT     -----------------------------------------------------

#--- Plot elpd improvement
talk_fonts <- theme(panel.grid=element_blank(),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    strip.text=element_text(size=16),
                    title=element_text(size=18))
mod_col <- c("Joint"="#7b3294", "Structured"="#008837")

loo.f <- dir(fS_out, "loo.*csv")
loo.df <- tibble(mod=str_split_fixed(loo.f, "_", 3)[,2], 
                 nCov=as.numeric(str_sub(str_split_fixed(loo.f, "_", 3)[,3], 1, -7)), 
                 elpd=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$elpd_loo[1]),
                 elpd_se=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$se_elpd_loo[1]),
                 looic=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$looic[1]),
                 looic_se=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$se_looic[1]),
                 v_full=map_chr(loo.f, ~read_csv(paste0(fS_out, .x))$X1[1])) %>%
  mutate(v_scale=str_sub(str_split_fixed(v_full, "__", 2)[,2], 1, 1),
         v_name=str_sub(str_split_fixed(v_full, "__", 2)[,2], 3, -1),
         model=case_when(mod=="WY" ~ "Joint",
                         mod=="Y" ~ "Structured")) %>%
  group_by(model) %>% arrange(model, looic) %>%
  mutate(elpd_diff=elpd-first(elpd),
         looic_diff=looic-first(looic))
  
ggplot(loo.df, aes(x=nCov, y=elpd, colour=model, group=model)) + 
  geom_point(size=3) + geom_line() + 
  scale_colour_manual("", values=mod_col) + 
  labs(x="Number of covariates") + 
  geom_errorbar(aes(ymin=elpd-elpd_se, ymax=elpd+elpd_se), width=0.1) +
  geom_text(aes(label=v_name), colour=1, hjust=0, vjust=1, size=3, nudge_x=0.05) +
  scale_shape_manual(values=c(1, 16)) + theme_bw()
# ggsave("eda/TEMP_varsel.jpg", width=8, height=5)

ggplot(loo.df, aes(x=nCov, y=looic, colour=model, group=model)) + 
  geom_point(size=3) + geom_line() + 
  scale_colour_manual("", values=mod_col) + 
  labs(x="Number of covariates") + 
  # geom_errorbar(aes(ymin=looic-looic_se, ymax=looic+looic_se), width=0.1) +
  geom_text(aes(label=v_name), colour=1, hjust=0, vjust=1, size=3, nudge_x=0.05) +
  scale_shape_manual(values=c(1, 16)) + theme_bw()







#--- relocate, but for curiosity...
best <- list(Y=map(1:n_folds, 
                   ~rstan::read_stan_csv(
                     dir(fS_out,
                         paste0("^Y_4__k-", .x, "_L_CnpyOpn"), full.names=T))),
             WY=map(1:n_folds, 
                    ~rstan::read_stan_csv(
                      dir(fS_out,
                          paste0("WY_5__k-", .x, "_L_VegTot"), full.names=T))))

best.beta <- map(best, ~do.call("rbind", 
                                map(.x, ~do.call("rbind", 
                                                 rstan::As.mcmc.list(.x, pars="beta")))))
best.b <- map(best, ~do.call("rbind", 
                             map(.x, ~do.call("rbind", 
                                              rstan::As.mcmc.list(.x, pars="b")))))

hdi.beta <- map(best.beta, ~HDInterval::hdi(.x) %>%
                  as_tibble %>% mutate(Limit=c("lower", "upper")) %>%
                  pivot_longer(!Limit, names_to="par", values_to="value") %>%
                  mutate(cov=str_split_fixed(par, "\\.", 2)[,2]))
hdi.b <- map(best.b, ~HDInterval::hdi(.x) %>%
               as_tibble %>% mutate(Limit=c("lower", "upper")) %>%
               pivot_longer(!Limit, names_to="par", values_to="value") %>%
               mutate(cov=str_split_fixed(par, "\\.", 3)[,2],
                      spp=str_split_fixed(par, "\\.", 3)[,3]))


p <- ggplot(hdi.b$WY, aes(value, spp)) + 
  geom_vline(xintercept=0, colour="gray80") +
  geom_line() + facet_wrap(~cov, scales="free_x") + 
  theme(panel.grid=element_blank())
ggsave("~/Desktop/WY.pdf", p)


best.i <- list(Y=readRDS("data/fwdSearch/Y_4__k-1_L_CnpyOpn_ls.rds"),
               WY=readRDS("data/fwdSearch/WY_5__k-1_L_VegTot_ls.rds"))

plot(best$Y, pars="beta") + xlim(-2.5,2.5)
c(colnames(best.i$Y$X), colnames(best.i$Y$V))

plot(best$WY, pars="beta", outer_level=0.9) + xlim(-4,4)
c(colnames(best.i$WY$X), colnames(best.i$WY$V))

stan_diag(best$WY[[1]])
stan_rhat(best$WY[[1]])
stan_ess(best$WY[[1]])
stan_mcse(best$Y)

best.df <- sjstats::tidy_stan(best$Y, prob=c(0.5, 0.9, 0.95))




