


# Testing
library(rstan); library(tidyverse); library(ggmcmc)
options(mc.cores=parallel::detectCores())
rstan::rstan_options(auto_write=TRUE)
mods <- "Y_GP"
# mods <- as.vector(outer(c("W_", "Y_", "WY_"), "GP", paste0))
# c("Pois", "GP", "NB", "ZIP"), paste0))
# pars.exc <- c("a_std", "A_std", "b_std", "B_std", 
#               "LAMBDA", "LAMBDA_",
#               #"prPres", "prPres_",# "prPresL_",
#               "lambda", "lambda_",
#               "lLAMBDA", "lLAMBDA_", "llambda", #"llambda_",
#               "L_Omega_A", "L_Omega_B", "sigma_A", "sigma_B",
#               "p", "p_")
# out.ls <- map(mods, ~stan(file="code/vs/full_Y_GP.stan",
#                           #file=paste0("code/vs/vs_", ., ".stan"),
#                           # pars=pars.exc, include=F, 
#                           init=0,
#                           data=read_rdump("data/stan_data/full.Rdump"), 
#                           # data=read_rdump("data/vs/vs_test.Rdump"), 
#                           chains=1, thin=1, warmup=50, iter=55))
out.ls <- stan(file="code/vs/full_Y_GP.stan",
               init=0,
               data=read_rdump("data/stan_data/full.Rdump"), 
               chains=1, thin=1, warmup=5, iter=10)
d.ls <- readRDS(paste0("data/stan_data/full_ls.rds"))
d.i <- readRDS(paste0("data/stan_data/full_i.rds"))


gg.llam <- map(out.ls, ~ggs(., "llambda_"))
gg.prp <- map(out.ls, ~ggs(., "prPresL_"))
out.sum <- map(out.ls, ~summary(.)$summary)

out.prp <- summary(out.ls, pars="prPres")$summary
out.prp_ <- summary(out.ls, pars="prPres_")$summary

pred.df <- rbind(
  tibble(mn=out.prp[,1],
         par=rownames(out.prp),
         set="fit") %>%
    mutate(site=str_split_fixed(str_split_fixed(par, "\\[", n=2)[,2],
                                ",", n=2)[,1],
           spp=str_sub(str_split_fixed(par, ",", n=2)[,2], 1, -2)) %>%
    arrange(as.numeric(site), as.numeric(spp)) %>%
    mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"]),
  tibble(mn=out.prp_[,1],
         par=rownames(out.prp_),
         set="fit") %>%
    mutate(site=str_split_fixed(str_split_fixed(par, "\\[", n=2)[,2],
                                ",", n=2)[,1],
           spp=str_sub(str_split_fixed(par, ",", n=2)[,2], 1, -2)) %>%
    arrange(as.numeric(site), as.numeric(spp)) %>%
    mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"])
)




for(i in seq_along(out.ls)) {
  if(grepl("Pois", names(out.ls)[i])) {
    gg.llam[[i]] <- gg.llam[[i]] %>% 
      mutate(prP=gg.prp[[i]]$value,
             sim=rpois(n(), exp(value)),
             obs=rep(c((d.i$Y_)), each=max(.$Iteration)),
             model=mods[i])
  } else if(grepl("GP", names(out.ls)[i])) {
    gg.d <- ggs(out.ls[[i]], "disp_lam")
    gg.llam[[i]] <- gg.llam[[i]] %>% 
      mutate(prP=gg.prp[[i]]$value,
             phi=gg.d$value[match(Iteration, gg.d$Iteration)],
             # sim=LaplacesDemon::rgpois(n(), exp(value), phi), # only dgpois...
             sim=rpois(n(), exp(value)),
             obs=rep(c((d.i$Y_)), each=max(.$Iteration)),
             model=mods[i]) %>% 
      select(-phi)
  } else if(grepl("NB", names(out.ls)[i])) {
    gg.d <- ggs(out.ls[[i]], "disp_lam")
    gg.llam[[i]] <- gg.llam[[i]] %>% 
      mutate(prP=gg.prp[[i]]$value,
             phi=gg.d$value[match(Iteration, gg.d$Iteration)],
             sim=rnbinom(n(), size=phi, mu=exp(value)),
             obs=rep(c((d.i$Y_)), each=max(.$Iteration)),
             model=mods[i]) %>% 
      select(-phi)
  } else if(grepl("ZIP", names(out.ls)[i])) {
    gg.z <- ggs(out.ls[[i]], "zi")
    gg.llam[[i]] <- gg.llam[[i]] %>% 
      mutate(prP=gg.prp[[i]]$value,
             z=gg.z$value[match(Iteration, gg.z$Iteration)],
             sim=(1-rbinom(n(), 1, z))*rpois(n(), exp(value)),
             obs=rep(c((d.i$Y_)), each=max(.$Iteration)),
             model=mods[i]) %>% 
      select(-z)
  }
}

out.all <- do.call('rbind', gg.llam)

out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim)) %>%
  ggplot(aes(x=prP, y=as.numeric(obs>0))) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, method.args=list(family="binomial")) +
  facet_wrap(~model) + xlim(0, 1)

out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim)) %>%
  ggplot(aes(x=prP, y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth() + 
  facet_wrap(~model) + xlim(0, 1)

out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim)) %>%
  ggplot(aes(x=sim, y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  facet_wrap(~model) + xlim(0, 10)

out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim)) %>%
  ggplot(aes(x=as.factor(obs), y=prP, fill=model)) + 
  geom_boxplot(outlier.colour=NA)
out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim)) %>%
  ggplot(aes(x=obs>0, y=prP, fill=model)) + 
  geom_boxplot(outlier.colour=NA)
out.all %>% group_by(model, Parameter) %>%
  summarise(obs=mean(obs), prP=mean(prP), sim=mean(sim), lam=mean(value)) %>%
  ggplot(aes(x=as.factor(obs), y=exp(lam), fill=model)) + 
  geom_boxplot(outlier.colour=NA) + ylim(0, 5)



# lpd (Y_ | lambda_)
library(tidyverse); library(ggmcmc); theme_set(theme_bw())
gg.ll <- map_dfr(out.ls, ~ggs(., "log_lik_"), .id="model")
gg.ll.lpd <- gg.ll %>%   # raw values: more negative = worse
  filter(!is.infinite(value)) %>%
  group_by(model, Parameter) %>% 
  summarise(mn=log(mean(exp(value)))) %>%  # more negative = worse
  group_by(model) %>% 
  summarise(lpd=-2*sum(mn)/n())  # larger = worse
gg.ll.lpd %>% arrange(lpd)  # smallest = best, largest = worse

# waic, loo
library(loo)
map(out.ls, ~loo::waic(loo::extract_log_lik(., "log_lik_"))) %>%
  loo::loo_compare()
map(out.ls[-c(2,4)], ~loo::loo(loo::extract_log_lik(., "log_lik_"))) %>%
  loo::loo_compare()

# parameter estimates
mod_cols <- c(W="#e41a1c", Y="#377eb8", WY="#984ea3")

beta.out <- aggregate_aggSlopes(out.ls, list(agg=NA), "beta")
ggplot(beta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_rug(data=beta.out$post_mns, aes(x=mean), sides="b") +
  geom_density(aes(group=paste(model, Chain))) + 
  scale_colour_manual(values=mod_cols) 

alpha.out <- aggregate_aggSlopes(out.ls, list(agg=NA), "alpha")
ggplot(alpha.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_rug(data=alpha.out$post_mns, aes(x=mean), sides="b") +
  geom_density(aes(group=paste(model, Chain))) + 
  scale_colour_manual(values=mod_cols) 

eta.out <- aggregate_eta(out.ls, NA)
ggplot(eta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_density(aes(group=paste(model, Chain))) + 
  scale_colour_manual(values=mod_cols) 

b.out <- aggregate_spSlopes(out.ls, 
                            list(sp=matrix(nrow=out.ls$W@par_dims$b[1], 
                                           ncol=out.ls$W@par_dims$b[2])), 
                            tax_i[,3:5], "b")
b.out$sum.gg %>% select(model, Parameter, mn, R) %>% 
  spread(model, mn) %>%
  ggplot() +
  geom_point(aes(x=W, y=WY), colour="black") + 
  geom_point(aes(x=Y, y=WY), colour="red") + 
  facet_wrap(~R, scales="free") +
  geom_abline(linetype=2) 
b.out$sum.gg %>% select(model, Parameter, mn, R) %>% 
  spread(model, mn) %>%
  ggplot() + geom_point(aes(x=WY, y=W-WY), colour="black") +
  geom_point(aes(x=WY, y=Y-WY), colour="red") +
  facet_wrap(~R, scales="free") +
  geom_hline(yintercept=0)

d.ls <- read_rdump("data/stan_data/vs_30.Rdump")
lam.out <- map_dfr(out.ls, ~ggs(., "^lambda_"), .id="model") %>%
  group_by(model, Parameter) %>%
  summarise(mn=exp(mean(log(value)))) %>%
  mutate(I=str_split_fixed(str_remove(Parameter, "lambda_\\["), ",", 2)[,1],
         S=str_split_fixed(str_remove(Parameter, "\\]"), ",", 2)[,2]) %>%
  arrange(model, as.numeric(I), as.numeric(S))
lam.out$obs <- rep(c(t(d.ls$Y_)), times=3)
lam.out$p <- dpois(lam.out$obs, lam.out$mn)
ggplot(lam.out, aes(mn, obs, colour=model)) + geom_point(alpha=0.2) +
  facet_wrap(~model) + 
  scale_colour_manual(values=mod_cols) 
ggplot(lam.out, aes(1-exp(-mn), obs, colour=model)) + 
  geom_point(alpha=0.2) +facet_wrap(~model) + 
  scale_colour_manual(values=mod_cols) 

I.sum <- lam.out %>% group_by(model, I) %>%
  summarise(nObs=sum(obs), nLam0_5=sum(mn>0.1), npP0_5=sum((1-exp(-mn))>0.05))

ggplot(I.sum, aes(x=nLam0_5, y=nObs, colour=model)) + 
  geom_point(alpha=0.3) + facet_wrap(~model) + 
  scale_colour_manual(values=mod_cols) 
ggplot(I.sum, aes(x=npP0_5, y=nObs, colour=model)) + 
  geom_point(alpha=0.3) + facet_wrap(~model) + 
  scale_colour_manual(values=mod_cols) 

tax_i %>% mutate(D=rstan::summary(out.ls$WY, pars="D")$summary[,1]) %>% 
  arrange(D) %>% print.AsIs

matrix(rstan::summary(out.ls$WY, pars="L_Omega_A")$summary[,1], nrow=8, byrow=T)
matrix(rstan::summary(out.ls$WY, pars="L_Omega_B")$summary[,1], nrow=8, byrow=T)
rstan::summary(out.ls$WY, pars=c("sigma_a", "sigma_b"))$summary[,1]
rstan::summary(out.ls$WY, pars=c("alpha"))$summary[,1]
rstan::summary(out.ls$WY, pars=c("beta"))$summary[,1]
rstan::summary(out.ls$WY, pars=c("eta"))$summary[,1]
