# Forward search

library(tidyverse)
library(rstan)
source("code/00_fn.R")
d.f <- "data/stan_data/pred"
fS_out <- "out/fwdSearch/"
fS_dat <- "data/fwdSearch/"
opt_dat <- "data/opt/"

#--- Make full dataset with W prediction cells
# Run 01_makeDatasets.R with set_type="pred"
# Creates data/stan_data/pred* (_ls.rds, _i.rds, .Rdump)




##### Identify optimal variable sets     ---------------------------------------
loo.f <- dir(fS_out, "loo.*csv")
loo.df <- tibble(mod=str_split_fixed(loo.f, "_", 3)[,2], 
                 nCov=as.numeric(str_split_fixed(loo.f, "_", 4)[,3]), 
                 elpd=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$elpd_loo[1]),
                 elpd_se=map_dbl(loo.f, ~read_csv(paste0(fS_out, .x))$se_elpd_loo[1]),
                 v_full=map_chr(loo.f, ~read_csv(paste0(fS_out, .x))$X1[1])) %>%
  group_by(mod) %>% filter(abs(elpd-max(elpd))<4) %>% 
  # filter(nCov == max(nCov)) %>% ungroup %>%
  filter(elpd == max(elpd)) %>% ungroup %>%
  mutate(file=paste(mod, nCov, "_k-1", str_split_fixed(v_full, "__", 2)[,2], sep="_"))
opt_var <- map(loo.df$file, ~readRDS(paste0(fS_dat, .x, "_ls.rds"))) %>%
  setNames(loo.df$mod) %>%
  map(~c(colnames(.$X), colnames(.$V)))





##### Make datasets     --------------------------------------------------------
make_next_datasets(full_data_base=d.f, out_base=paste0(opt_dat, "Y"), 
                   opt_names=opt_var$Y, type="pred", mod_size=length(opt_var$Y))
make_next_datasets(full_data_base=d.f, out_base=paste0(opt_dat, "WY"), 
                   opt_names=opt_var$WY, type="pred", mod_size=length(opt_var$WY))

