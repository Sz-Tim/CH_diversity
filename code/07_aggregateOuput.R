
library(tidyverse); library(rstan); library(sf); theme_set(theme_bw())
library(viridis)
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

out.dir <- "cmdstan-2.23.0/opfo/out/"


d.ls <- readRDS(paste0("data/stan_data/full_ls.rds"))
d.i <- readRDS(paste0("data/stan_data/full_i.rds"))

mods <- sort(unique(str_split_fixed(dir(out.dir, "^[WY]"), "_[0-9]+", 2)[,1]))

pars.any <- c("H", "lLAMBDA", 
              "beta", "alpha", "eta", 
              "b", "a", "B", "A", "D", 
              "llambda", 
              "PRP", "prp", 
              "disp_lam")
pars.any.full <- c("ShannonH", "ShannonH_", "lLAMBDA", "lLAMBDA_",
                   "beta", "alpha", "eta",
                   "b", "a", "B", "A", "D", 
                   "llambda", 
                   "prPres", "prPres_", "prPresL",
                   "disp_lam")

out.pars <- imap(setNames(pars.any, pars.any), ~vector("list", length(mods)))
out.stan <- vector("list", length(mods)) %>% setNames(mods)



for(i in seq_along(mods)) {
  cat("\n", format(Sys.time(), "%X"), "-- Beginning model", mods[i])
  cat("\n  Reading", paste0(out.dir, mods[i]))
  out.stan[[i]] <- read_stan_csv(dir(out.dir, paste0("^", mods[i], "_[0-9]"), 
                                     full.names=T))
  pars.full <- unique(str_split_fixed(names(out.stan[[i]]), "\\[", 2)[,1])
  pars <- pars.full[pars.full %in% pars.any.full]
  out <- summary(out.stan[[i]], pars=pars)$summary
  
  out.ls <- data.frame(Parameter=rownames(out),
                       mn=out[,1],
                       se=out[,2],
                       q025=out[,4], 
                       q25=out[,5],
                       q50=out[,6],
                       q75=out[,7],
                       q975=out[,8],
                       model=mods[i],
                       Rhat=out[,10],
                       row.names=NULL) %>%
    mutate(par=str_split_fixed(Parameter, "\\[", n=2)[,1]) %>%
    split(., .$par, drop=T)
  
  # prPres (site)
  if("prPres" %in% pars) {
    out.pars$PRP[[i]] <- rbind(
      out.ls$prPres %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="fit") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
               source=ifelse(id<100000, "W", "Y")),
      out.ls$prPres_ %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="pred") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
               source=ifelse(id<100000, "W", "Y"))
    ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  
  if("disp_lam" %in% pars) {
    out.pars$disp_lam[[i]] <- out.ls$disp_lam %>%
      mutate(Parameter=as.character(Parameter),
             model=as.character(model))
  }
  
  out.pars$beta[[i]] <- out.ls$beta %>%
    mutate(ParName=c("intercept", c(colnames(d.i$X)[-(1:2)],
                                    colnames(d.i$V)[-(1:2)]))) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  
}


all.out <- map(out.pars, ~do.call('rbind', .))


left_join(filter(d.i$grd_W.sf, inbd), 
          all.out$PRP %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(mn > 0.99)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S)) + geom_sf(colour=NA) + #facet_wrap(~model) + 
  scale_fill_viridis("S mean", na.value="white", 
                     option="B", limits=c(0, d.ls$S)) +
  theme(panel.grid=element_blank(), 
        legend.position="none")
ggsave("~/Desktop/S_pred.pdf", height=6, width=6)


left_join(filter(d.i$grd_W.sf, inbd), 
          all.out$PRP %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(mn > 0.99)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(y=S)) + geom_sf(colour=NA) + #facet_wrap(~model) + 
  theme(panel.grid=element_blank(), 
        legend.position="none")


library(bayesplot)
mcmc_areas(out.stan[[1]], pars=paste0("beta[", 2:9, "]"))
