
library(tidyverse); library(rstan); theme_set(theme_bw())
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
out.dir <- "out/vs/"

d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))

mods <- sort(unique(str_split_fixed(dir(out.dir, "^[WY]"), "_[0-9]+", 2)[,1]))
pars.any <- c("LL",
              "H", "lLAMBDA", 
              "beta", "alpha", "eta", 
              "b", "a", "B", "A", "D", 
              "llambda", 
              "PRP", "prp", 
              "disp_lam", "zi")
pars.any.full <- c("log_lik_", 
                   "ShannonH", "ShannonH_", "lLAMBDA", "lLAMBDA_",
                   "beta", "alpha", "eta",
                   "b", "a", "B", "A", "D", 
                   "llambda", "llambda_", 
                   "prPres", "prPres_", "prPresL_",
                   "disp_lam", "zi")
out.pars <- imap(setNames(pars.any, pars.any), ~vector("list", length(mods)))
out.loo <- out.stan <- vector("list", length(mods)) %>% setNames(mods)

Y.J <- matrix(0, ncol=d.ls$S, nrow=d.ls$J)
Y.J_ <- matrix(0, ncol=d.ls$S, nrow=d.ls$J_)
for(i in 1:d.ls$J) Y.J[i,] <- colSums(d.ls$Y[d.ls$IJ==i,])
for(i in 1:d.ls$J_) Y.J_[i,] <- colSums(d.ls$Y_[d.ls$IJ_==i,])
obs <- data.frame(el=c(d.i$X[,2], d.i$X_[,2]), 
                  H.obs=vegan::diversity(rbind(d.ls$W, Y.J, Y.J_)),
                  N.obs=rowSums(rbind(d.ls$W, Y.J, Y.J_)),
                  S.obs=rowSums(rbind(d.ls$W, Y.J, Y.J_)>0),
                  site=c(1:nrow(d.i$X), 1:nrow(d.i$X_)),
                  set=rep(c("train", "test"), 
                          times=c(nrow(d.i$X), nrow(d.i$X_))),
                  source=rep(c("W", "Y", "Y"), 
                             times=c(d.ls$K, d.ls$J, d.ls$J_))) %>%
  mutate(set=as.character(set), source=as.character(source))


for(i in seq_along(mods)) {
  cat("\n", format(Sys.time(), "%X"), "-- Beginning model", mods[i])
  cat("\n  Reading", paste0(out.dir, mods[i]))
  out.stan[[i]] <- read_stan_csv(dir(out.dir, paste0("^", mods[i], "_[0-9]"), 
                                     full.names=T))
  out.loo[[i]] <- loo::extract_log_lik(out.stan[[i]], parameter_name="log_lik_")
  # out.loo[[i]] <- loo::loo(out.stan[[i]], pars="log_lik_")
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
  
  # log likelihood
  if("log_lik_" %in% pars) {
    out.pars$LL[[i]] <- out.ls$log_lik_ %>%
      mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             set="test") %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      arrange(plot, spp) %>%
      mutate(Plot_id=d.i$V_[plot,"Plot_id"], el=d.i$V_[plot,"el"],
             Parameter=as.character(Parameter), 
             model=as.character(model),
             obs=c(t(d.i$Y_)))
  }
  
  # prPres (plot)
  if("prPresL_" %in% pars) {
    out.pars$prp[[i]] <- out.ls$prPresL_ %>%
      mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             set="test") %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      arrange(plot, spp) %>%
      mutate(Plot_id=d.i$V_[plot,"Plot_id"], el=d.i$V_[plot,"el"],
             Parameter=as.character(Parameter), 
             model=as.character(model),
             obs=c(t(d.i$Y_)))
  }
  
  # prPres (site)
  if("prPres" %in% pars) {
    out.pars$PRP[[i]] <- rbind(
      out.ls$prPres %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="train") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
               source=ifelse(id<100000, "W", "Y"),
               obs=c(t(rbind(d.i$W, Y.J)))),
      out.ls$prPres_ %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="test") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
               source=ifelse(id<100000, "W", "Y"),
               obs=c(t(Y.J_)))
    ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  
  # llambda
  if("llambda_" %in% pars) {
    out.pars$llambda[[i]] <- out.ls$llambda_ %>%
      mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             set="test") %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      arrange(plot, spp) %>%
      mutate(Plot_id=d.i$V_[plot,"Plot_id"], el=d.i$V_[plot,"el"],
             Parameter=as.character(Parameter), 
             model=as.character(model),
             obs=c(t(d.i$Y_)))
  }
  
  if("zi" %in% pars) {
    out.pars$zi[[i]] <- out.ls$zi %>%
      mutate(Parameter=as.character(Parameter),
             model=as.character(model))
  }
  
  if("disp_lam" %in% pars) {
    out.pars$disp_lam[[i]] <- out.ls$disp_lam %>%
      mutate(Parameter=as.character(Parameter),
             model=as.character(model))
  }
  out.pars$beta[[i]] <- out.ls$beta %>%
    # mutate(ParName=c("intercept", colnames(d.i$X)[-(1:2)])) %>% 
    mutate(ParName=c("intercept", c(colnames(d.i$X)[-(1:2)],
                                    colnames(d.i$V)[-(1:2)]))) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  if("alpha" %in% pars) {
    out.pars$alpha[[i]] <- out.ls$alpha %>%
      mutate(ParName=colnames(d.i$V)[-(1:2)]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model)) 
  }
}


out.loo.clean <- map(out.loo, ~replace(., is.na(.), -1e20) %>%
                       replace(., is.infinite(.), -1e20)) %>%
  map(loo::loo) 

all.out <- map(out.pars, ~do.call('rbind', .))
# saveRDS(all.out, "out/vs/allPars_temp_HS.rds")

mod_cols <- c(W_GP="#a50f15", Y_GP="#08519c", WY_GP="#54278f",
              W_NB="#de2d26", Y_NB="#3182bd", WY_NB="#756bb1", 
              W_Pois="#fb6a4a", Y_Pois="#6baed6", WY_Pois="#9e9ac8", 
              W_ZIP="#fc9272", Y_ZIP="#9ecae1", WY_ZIP="#bcbddc",
              W_HS="#fcbba1", Y_HS="#c6dbef", WY_HS="#dadaeb")



loo::loo_compare(out.loo.clean)

all.out$LL %>% group_by(model) %>% 
  # filter(obs>0) %>%
  summarise(lpd=-2*sum(mn)/n()) %>% arrange(lpd)



ggplot(all.out$prp, aes(x=q50, y=as.numeric(obs>0), colour=model)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, 
              method.args=list(family="binomial"), fullrange=T) + 
  stat_smooth(aes(x=q025), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q25), method="glm", size=0.5, linetype=2, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q75), method="glm", size=0.5, linetype=2, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q975), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  scale_colour_manual(values=mod_cols, guide=F) +
  facet_wrap(~model) + xlim(0,1)

ggplot(all.out$prp, aes(x=mn, y=as.numeric(obs>0), colour=model)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", method.args=list(family="binomial"), 
              fullrange=T, size=0.5, se=F) + 
  scale_colour_manual(values=mod_cols) + xlim(0,1) + facet_wrap(~spp)

ggplot(all.out$prp, aes(x=factor(obs), y=mn, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.25, outlier.alpha=0.25) + 
  scale_fill_manual(values=mod_cols) + 
  labs(x="Observed abundance", y="Predicted pr(presence)",
       title="Plot-level out-of-sample predictions")

ggplot(all.out$prp, aes(x=factor(as.numeric(obs>0)), y=mn, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.25, outlier.alpha=0.25) + 
  scale_fill_manual(values=mod_cols) + 
  labs(x="Observed absence or presence", y="Predicted pr(presence)",
       title="Plot-level out-of-sample predictions")


ggplot(filter(all.out$PRP, source=="Y" & set=="test"), 
       aes(x=factor(obs), y=mn, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.25, outlier.alpha=0.25) + 
  scale_fill_manual(values=mod_cols) + 
  labs(x="Observed abundance", y="Predicted pr(presence)",
       title="Site-level out-of-sample predictions")

ggplot(filter(all.out$PRP, source=="Y" & set=="test"), 
       aes(x=factor(as.numeric(obs>0)), y=mn, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.25, outlier.alpha=0.25) + 
  scale_fill_manual(values=mod_cols) + 
  labs(x="Observed absence or presence", y="Predicted pr(presence)",
       title="Site-level out-of-sample predictions")


ggplot(filter(all.out$PRP, source=="Y" & set=="test"), 
       aes(x=model, alpha=factor(mn>0.95), fill=model)) +
  geom_bar(position="fill") + facet_wrap(~obs>0) +
  scale_alpha_manual("prP > 0.95", values=c(0.5, 1)) + 
  scale_fill_manual(values=mod_cols) + 
  theme(panel.grid=element_blank()) 

ggplot(filter(all.out$PRP, source=="Y" & set=="test"), 
       aes(x=model, alpha=factor(mn>0.95), fill=model)) +
  geom_bar(position="fill") + facet_wrap(~obs) +
  scale_alpha_manual("prP > 0.95", values=c(0.5, 1)) + 
  scale_fill_manual(values=mod_cols) + 
  theme(panel.grid=element_blank()) 

ggplot(all.out$prp, aes(x=model, alpha=factor(mn>0.05), fill=model)) +
  geom_bar(position="fill") + facet_wrap(~obs>0) +
  scale_alpha_manual("prP > 0.05", values=c(0.5, 1)) + 
  scale_fill_manual(values=mod_cols) + 
  theme(panel.grid=element_blank()) 

ggplot(all.out$prp, aes(x=model, alpha=factor(mn>0.05), fill=model)) +
  geom_bar(position="fill") + facet_wrap(~obs) +
  scale_alpha_manual("prP > 0.05", values=c(0.5, 1)) + 
  scale_fill_manual(values=mod_cols) + 
  theme(panel.grid=element_blank()) 



ggplot(filter(all.out$PRP, source=="Y" & set=="test"), 
       aes(x=mn, y=as.numeric(obs>0))) + 
  geom_point(alpha=0.5, shape=1, size=0.75) + 
  stat_smooth(method="glm", size=0.5, method.args=list(family="binomial")) +
  facet_wrap(~model)



ggplot(all.out$llambda, aes(x=exp(q50), y=as.numeric(obs>0))) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, method.args=list(family="binomial")) + 
  stat_smooth(aes(x=exp(q025)), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial")) +
  stat_smooth(aes(x=exp(q25)), method="glm", size=0.5, linetype=2, se=F,
              method.args=list(family="binomial")) +
  stat_smooth(aes(x=exp(q75)), method="glm", size=0.5, linetype=2, se=F,
              method.args=list(family="binomial")) +
  stat_smooth(aes(x=exp(q975)), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial")) +
  facet_wrap(~model) + xlim(0,10)#xlim(0,max(exp(all.out$llambda$q50)))

ggplot(all.out$llambda, aes(x=mn, y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, method.args=list(family="poisson")) +
  facet_wrap(~model) + xlim(0,10)#xlim(0,max(exp(all.out$llambda$mn)))

ggplot(all.out$llambda, aes(x=exp(mn), y=obs)) + 
  geom_abline(linetype=3, colour="gray30") +
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="lm", size=0.5) +
  # stat_smooth(aes(x=exp(q025)), method="lm", size=0.5, linetype=3, se=F) +
  # stat_smooth(aes(x=exp(q25)), method="lm", size=0.5, linetype=2, se=F) +
  # stat_smooth(aes(x=exp(q75)), method="lm", size=0.5, linetype=2, se=F) +
  # stat_smooth(aes(x=exp(q975)), method="lm", size=0.5, linetype=3, se=F) +
  facet_wrap(~model, scales="free_x")# + xlim(0,10)

ggplot(all.out$llambda, aes(x=sim, y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="lm", size=0.5) +
  facet_wrap(~model) + xlim(0,10)#xlim(0,max(sim))

ggplot(filter(all.out$llambda, obs>0), aes(x=mn, y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, method.args=list(family="poisson")) +
  facet_wrap(~model)

ggplot(filter(all.out$llambda, obs>0), aes(x=exp(mn), y=obs)) + 
  geom_point(alpha=0.5, shape=1, size=0.75) +
  stat_smooth(method="glm", size=0.5, method.args=list(family="poisson"), se=F) +
  stat_smooth(aes(x=exp(mn+2*se)), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="poisson")) +
  stat_smooth(aes(x=exp(mn-2*se)), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="poisson")) +
  facet_wrap(~model, scales="free_x") #+ xlim(0,max(exp(all.out$llambda$mn)))



all.out$LL %>% filter(obs>0) %>%
  ggplot(aes(factor(obs), exp(mn), fill=model)) + 
  geom_boxplot(outlier.size=0.5, outlier.alpha=0.3, outlier.shape=1) +
  scale_fill_manual(values=mod_cols) +
  labs(x=expression("Observed colony abundance"~(Y[is])), 
       y=expression(Pr(Y[is]~'|'~lambda,phi,z))) 



all.out$LL %>% filter(!grepl("W_", model)) %>%
  # filter(obs>0) %>%
  mutate(obs.f=as.factor(obs)) %>%
  aov(exp(mn) ~ model, data=.) %>% TukeyHSD %>% plot(cex.axis=0.5, las=2)



ggplot(all.out$beta, aes(x=ParName, y=q50, ymin=q025, ymax=q975, 
                         colour=model, shape=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=3) +
  geom_point(position=position_dodge(width=0.8)) + 
  geom_linerange(position=position_dodge(width=0.8)) + 
  scale_colour_manual(values=mod_cols) + 
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  coord_flip()

ggplot(filter(all.out$beta, ParName != "intercept" & !grepl("W_", model)), 
       aes(x=ParName, y=q50, ymin=q025, ymax=q975, 
                         colour=model, shape=sign(q25)==sign(q75))) + 
  geom_hline(yintercept=0, linetype=3) +
  geom_point(position=position_dodge(width=0.8), size=3) + 
  geom_linerange(aes(ymin=q25, ymax=q75), size=1,
                 position=position_dodge(width=0.8)) + 
  geom_linerange(position=position_dodge(width=0.8), size=0.25) + 
  scale_colour_manual(values=mod_cols) + 
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  coord_flip()

ggplot(all.out$b, aes(x=ParName, y=q50, ymin=q025, ymax=q975, 
                      colour=model, shape=sign(q025)==sign(q975))) + 
  geom_point(position=position_dodge(width=0.8)) + 
  geom_linerange(position=position_dodge(width=0.8)) + 
  scale_colour_manual(values=mod_cols) + 
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  coord_flip()



all.out$zi
all.out$disp_lam

Y_sim <- rbind(filter(all.out$llambda, grepl("ZIP", model)) %>%
                 mutate(sim=(1-rbinom(n(), 1, all.out$zi$mn[match(model, all.out$zi$model)]))*rpois(n(), exp(mn))),
               filter(all.out$llambda, !grepl("ZIP", model)) %>%
                 mutate(sim=rpois(n(), exp(mn))))

ggplot(Y_sim, aes(x=sim, y=obs)) + geom_abline() + 
  geom_point(alpha=0.3) + stat_smooth(method="lm") +
  facet_wrap(~model) + xlim(0, 10) + ylim(0, 10)

