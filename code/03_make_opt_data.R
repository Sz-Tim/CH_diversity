# Forward search

library(tidyverse)
library(rstan)
source("code/00_fn.R")
d.f <- "data/stan_data/pred"
fS_out <- "out/fwd_ll/"
fS_dat <- "data/fwdSearch/"
opt_dat <- "data/opt/"

#--- Make full dataset with W prediction cells
# Run 01_makeDatasets.R with set_type="pred"
# Creates data/stan_data/pred* (_ls.rds, _i.rds, .Rdump)




##### Identify optimal variable sets     ---------------------------------------
loo.f <- dir(fS_out, "loo.*csv")
loo.df <- suppressMessages({
  tibble(mod=str_split_fixed(loo.f, "_", 4)[,3], 
         LV=str_split_fixed(loo.f, "_", 4)[,2],
         nCov=as.numeric(str_sub(str_split_fixed(loo.f, "_", 4)[,4], 1, -7)), 
         elpd=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$elpd_loo[1]),
         elpd_se=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$se_elpd_loo[1]),
         looic=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$looic[1]),
         looic_se=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$se_looic[1]),
         v_full=map_chr(loo.f, ~read_csv(paste0(fS_out, .x))$X1[1])) %>%
    group_by(mod, LV) %>% filter(abs(looic-min(looic))<4) %>% 
    filter(looic == min(looic)) %>% ungroup %>%
    mutate(file=paste(LV, mod, nCov, "_k-1", str_split_fixed(v_full, "__", 2)[,2], sep="_"))
})
opt_var <- map(loo.df$file, ~readRDS(paste0(fS_dat, .x, "_ls.rds"))) %>%
  setNames(paste(loo.df$LV, loo.df$mod, sep="_")) %>%
  map(~c(colnames(.$X), colnames(.$V)))
null_var <- map(names(opt_var), ~"R_") %>% setNames(names(opt_var))




##### Make datasets     --------------------------------------------------------
# optimal variable sets
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "cov_Y"),
                   opt_names=opt_var$cov_Y, mod_size=length(opt_var$cov_Y))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "LV_Y"),
                   opt_names=opt_var$LV_Y, mod_size=length(opt_var$LV_Y))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "cov_WY"),
                   opt_names=opt_var$cov_WY, mod_size=length(opt_var$cov_WY))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "LV_WY"),
                   opt_names=opt_var$LV_WY, mod_size=length(opt_var$LV_WY))


# intercept only null models
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "cov_Y_null"),
                   opt_names=null_var$cov_Y, mod_size=length(null_var$cov_Y))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "LV_Y_null"),
                   opt_names=null_var$LV_Y, mod_size=length(null_var$LV_Y))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "cov_WY_null"),
                   opt_names=null_var$cov_WY, mod_size=length(null_var$cov_WY))
make_next_datasets(full_data_base=d.f, type="pred", 
                   out_base=paste0(opt_dat, "LV_WY_null"),
                   opt_names=null_var$LV_WY, mod_size=length(null_var$LV_WY))
