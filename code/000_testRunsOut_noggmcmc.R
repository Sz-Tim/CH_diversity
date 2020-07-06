

library(tidyverse); library(rstan); theme_set(theme_bw())
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

# vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
vs_d <- str_sub(dir("data/stan_data/", ".Rdump"), 1, -7)[1]
out.dir <- "out/samp_files/"

# d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
# d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))
d.ls <- readRDS(paste0("data/stan_data/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/stan_data/", vs_d, "_i.rds"))

mods <- sort(unique(str_split_fixed(dir(out.dir), "_[0-9]+", 2)[,1]))
pars.any <- c("H", "lLAMBDA", "beta", "alpha", "eta", "b", "a", "B", "A", "D", 
              "llambda", "prPres", "disp_lam", "zi")
pars.any.full <- c("ShannonH", "ShannonH_", "lLAMBDA", "lLAMBDA_", 
                   "beta", "alpha", "eta", 
                   "b", "a", "B", "A", "D", 
                   "llambda", "llambda_", "prPres", "prPres_", "disp_lam", "zi")
out.pars <- imap(setNames(pars.any, pars.any), ~vector("list", length(mods)))
out.raw <- vector("list", length(mods))

Y.J <- matrix(0, ncol=d.ls$S, nrow=d.ls$J)
Y.J_ <- matrix(0, ncol=d.ls$S, nrow=d.ls$J_)
for(i in 1:d.ls$J) {
  Y.J[i,] <- colSums(d.ls$Y[d.ls$IJ==i,])
}
for(i in 1:d.ls$J_) {
  Y.J_[i,] <- colSums(d.ls$Y_[d.ls$IJ_==i,])
}
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
  out <- read_stan_csv(dir(out.dir, paste0("^", mods[i], "_[0-9]"), full.names=T))
  out.raw[[i]] <- out
  pars.full <- unique(str_split_fixed(names(out), "\\[", 2)[,1])
  pars <- pars.full[pars.full %in% pars.any.full]
  out <- summary(out, pars=pars)$summary
  
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
  
  out.pars$H[[i]] <- rbind(
    out.ls$ShannonH %>% 
      mutate(site=str_remove(str_split_fixed(Parameter,  "\\[", 2)[,2], "]"),
             site=as.numeric(site),
             set="train") %>%
      arrange(site) %>%
      mutate(id=d.i$X[,"id"], el=d.i$X[,"el"], 
             source=ifelse(id<100000, "W", "Y")),
    out.ls$ShannonH_ %>% 
      mutate(site=str_remove(str_split_fixed(Parameter, "\\[", 2)[,2], "]"),
             site=as.numeric(site),
             set="test") %>%
      arrange(site) %>%
      mutate(id=d.i$X_[,"id"], el=d.i$X_[,"el"], 
             source=ifelse(id<100000, "W", "Y"))
  ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$lLAMBDA[[i]] <- rbind(
    out.ls$lLAMBDA %>%
      mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             set="train") %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
             source=ifelse(id<100000, "W", "Y")),
    out.ls$lLAMBDA_ %>%
      mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             set="test") %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
             source=ifelse(id<100000, "W", "Y"))
  ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$beta[[i]] <- out.ls$beta %>%
    # mutate(ParName=c("intercept", colnames(d.i$X)[-(1:2)])) %>% 
    mutate(ParName=c("intercept", colnames(d.i$X)[-(1:2)], 
                     colnames(d.i$V)[-(1:2)])) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$B[[i]] <- out.ls$B %>%
    mutate(ParNum=str_split_fixed(str_split_fixed(Parameter, ",", 2)[,1], 
                                  "\\[", 2)[,2],
           GenNum=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]")) %>%
    mutate(ParNum=as.numeric(ParNum),
           GenNum=as.numeric(GenNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)], 
                     colnames(d.i$V)[-(1:2)])[ParNum],
           GenName=unique(d.i$tax_i$genus)[GenNum]) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$b[[i]] <- out.ls$b %>%
    mutate(ParNum=str_split_fixed(str_split_fixed(Parameter, ",", 2)[,1], 
                                  "\\[", 2)[,2],
           SpNum=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]")) %>%
    mutate(ParNum=as.numeric(ParNum),
           SpNum=as.numeric(SpNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)], 
                     colnames(d.i$V)[-(1:2)])[ParNum],
           SpName=d.i$tax_i$species[SpNum],
           GenName=str_split_fixed(SpName, "_", 2)[,1]) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  if("alpha" %in% pars) {
    out.pars$alpha[[i]] <- out.ls$alpha %>%
      mutate(ParName=colnames(d.i$V)[-(1:2)]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    out.pars$A[[i]] <- out.ls$A %>%
      mutate(ParNum=str_split_fixed(str_split_fixed(Parameter, ",", 2)[,1], 
                                    "\\[", 2)[,2],
             GenNum=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]")) %>%
      mutate(ParNum=as.numeric(ParNum),
             GenNum=as.numeric(GenNum),
             ParName=c(colnames(d.i$V)[-(1:2)])[ParNum],
             GenName=unique(d.i$tax_i$genus)[GenNum]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    out.pars$a[[i]] <- out.ls$a %>%
      mutate(ParNum=str_split_fixed(str_split_fixed(Parameter, ",", 2)[,1], 
                                    "\\[", 2)[,2],
             SpNum=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]")) %>%
      mutate(ParNum=as.numeric(ParNum),
             SpNum=as.numeric(SpNum),
             ParName=c(colnames(d.i$V)[-(1:2)])[ParNum],
             SpName=d.i$tax_i$species[SpNum],
             GenName=str_split_fixed(SpName, "_", 2)[,1]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  if("eta" %in% pars) {
    out.pars$eta[[i]] <- out.ls$eta %>%
      mutate(ParName=c("intercept", colnames(d.i$U)[-1])) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  if("D" %in% pars) {
    out.pars$D[[i]] <- out.ls$D %>%
      mutate(SpNum=str_sub(str_split_fixed(as.character(Parameter), "\\[", 2)[,2],
                           1, -2) %>% as.numeric) %>% 
      mutate(SpName=d.i$tax_i$species[SpNum],
             GenName=str_split_fixed(SpName, "_", 2)[,1],
             Parameter=as.character(Parameter), 
             model=as.character(model))
  }
  if("llambda" %in% pars) {
    out.pars$llambda[[i]] <- rbind(
      out.ls$llambda %>%
        mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="train") %>%
        mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
        arrange(plot, spp) %>%
        mutate(Plot_id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               source="Y")
    ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  if("prPres" %in% pars) {
    out.pars$prPres[[i]] <- rbind(
      out.ls$prPres %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="train") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
               source=ifelse(id<100000, "W", "Y")),
      out.ls$prPres_ %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
               set="test") %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
               source=ifelse(id<100000, "W", "Y"))
    ) %>% mutate(Parameter=as.character(Parameter), model=as.character(model))
  }
  if("zi" %in% pars) {
    out.pars$zi[[i]] <- out.ls$zi
  }
  if("disp_lam" %in% pars) {
    out.pars$disp_lam[[i]] <- out.ls$disp_lam 
  }
}

all.pars <- map(out.pars, ~do.call('rbind', .))
saveRDS(all.pars, "out/tests/allPars_temp.rds")


mod_cols <- c(W_GP="#a50f15", Y_GP="#08519c", WY_GP="#54278f",
              W_NB="#de2d26", Y_NB="#3182bd", WY_NB="#756bb1", 
              W_Pois="#fb6a4a", Y_Pois="#6baed6", WY_Pois="#9e9ac8", 
              W_ZIP="#fc9272", Y_ZIP="#9ecae1", WY_ZIP="#bcbddc",
              W_HS="#fcbba1", Y_HS="#c6dbef", WY_HS="#dadaeb")

ggplot(all.pars$H, aes(x=el, y=mn, ymin=mn-2*se, ymax=mn+2*se, colour=model)) + 
  geom_point(alpha=0.5, aes(shape=set)) + 
  geom_linerange(size=0.2, alpha=0.75) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + ylim(0, NA) + 
  labs(x="Elevation (m)", y=expression(Predicted~Shannon~H~(1~km^2)))
ggsave("eda/all_H.pdf", width=8.5, height=4)

all.pars$lLAMBDA %>% filter(model!="W_Pois") %>%
  group_by(model, site, set, id, el, source) %>%
  summarise(tot=sum(exp(mn)), tot_q025=sum(q025), tot_q975=sum(q975)) %>%
  ggplot(aes(el, tot, ymin=tot_q025, ymax=tot_q975, colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  # geom_linerange(size=0.2, alpha=0.75) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_wrap(~source, scales="free") + #coord_trans(y="exp") + 
  labs(x="Elevation (m)", y=expression(Predicted~total~log(Lambda)~(1~km^2)))

all.pars$lLAMBDA %>% filter(model!="W_Pois") %>%
  group_by(model, site, set, id, el, source) %>%
  summarise(med_mn=median(mn), min_mn=min(mn), max_mn=max(mn)) %>%
  ggplot(aes(el, exp(med_mn), colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_wrap(~source, scales="free") +  
  labs(x="Elevation (m)", y=expression(Predicted~median~Lambda~(1~km^2)))
ggsave("eda/all_LAMBDA.pdf", width=8.5, height=4)

ggplot(filter(all.pars$beta, ParName!="intercept" & model != "W_Pois"), 
       aes(x=ParName, y=mn, colour=model, group=model,
           # shape=sign(mn+2*se)==sign(mn-2*se))) +
           shape=sign(q025)==sign(q975))) +
  geom_hline(yintercept=0, linetype=2) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=q025, ymax=q975),
  # geom_linerange(aes(ymin=mn-2*se, ymax=mn+2*se), 
                 position=position_dodge(width=0.5)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(beta)))
ggsave("eda/all_beta.pdf", width=6, height=5)
ggplot(filter(all.pars$B, ParName!="intercept"), 
       aes(x=ParName, y=mn, colour=model, group=model, 
           shape=sign(q025)==sign(q975),
           size=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=2) +
  geom_point(alpha=0.7, position=position_dodge(width=0.5)) +
  scale_size_manual("Non-zero", values=c(0.5, 2)) +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(B)))
ggsave("eda/all_B.pdf", width=6, height=5)
ggplot(filter(all.pars$b, ParName!="intercept"), 
       aes(x=ParName, y=mn, colour=model, group=model, 
           shape=sign(q025)==sign(q975),
           size=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=2) +
  geom_point(alpha=0.7, position=position_dodge(width=0.5)) +
  scale_size_manual("Non-zero", values=c(0.5, 2)) +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(b)))
ggsave("eda/all_b.pdf", width=6, height=5)
ggplot(filter(all.pars$b, ParName!="intercept"),
       aes(x=ParName, y=mn, colour=model, group=paste0(model, ParName))) +
  geom_boxplot(outlier.size=0.25) +
  geom_hline(yintercept=0, colour="gray", size=0.25) +
  scale_colour_manual(values=mod_cols) + coord_flip()
all.pars$b %>% filter(ParName!="intercept") %>% 
  group_by(ParName, model) %>% 
  summarise(nSig=sum(sign(q025)==sign(q975))) %>% 
  arrange(model, nSig) %>% ungroup %>%
  mutate(ParNameOrd=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(x=ParNameOrd, y=nSig, fill=model)) + 
  geom_bar(stat='identity', position="dodge") + 
  scale_fill_manual(values=mod_cols) + coord_flip() + 
  facet_grid(.~model) + theme(legend.position='none') +
  labs(x="", y="Number of species with significant relationship")
all.pars$b %>% filter(ParName!="intercept" & model!="W_Pois") %>% 
  group_by(ParName, model) %>% 
  summarise(mnEff=mean(abs(mn)), sd=sd(mn)) %>%
  ggplot(aes(x=mnEff, y=sd, label=ParName)) + 
  geom_point() + geom_text(size=3, hjust=0, nudge_x=0.05) +
  stat_smooth(method="lm", se=F, size=0.2) + facet_wrap(~model)


ggplot(filter(all.pars$eta, ParName!="intercept"), 
       aes(x=ParName, y=mn, colour=model, shape=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=2) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=q025, ymax=q975), 
                 position=position_dodge(width=0.5)) +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(eta)))
ggsave("eda/all_eta.pdf", width=6, height=5)


ggplot(all.pars$alpha, 
       aes(x=ParName, y=mn, colour=model, group=model,
           shape=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=2) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=q025, ymax=q975), 
                 position=position_dodge(width=0.5)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(alpha)))
ggsave("eda/all_alpha.pdf", width=6, height=5)
ggplot(all.pars$A, 
       aes(x=ParName, y=mn, colour=model, group=model, 
           shape=sign(q025)==sign(q975))) + 
  geom_hline(yintercept=0, linetype=2) +
  geom_point(size=1, alpha=0.7, position=position_dodge(width=0.5)) +
  scale_shape_manual("Non-zero", values=c(1, 19)) +
  scale_colour_manual(values=mod_cols) + coord_flip() +
  labs(x="Variable", y=expression(Slope~posterior~distribution~(A)))
ggsave("eda/all_A.pdf", width=6, height=5)





D.cols <- c("red3", "gray30")[(sign(all.pars$D$q025-1) != sign(all.pars$D$q975-1))+1]
ggplot(all.pars$D, aes(x=SpName, y=mn, ymin=q025, ymax=q975,
                       colour=sign(q025-1) != sign(q975-1))) + 
  geom_hline(yintercept=1, linetype=2) +
  geom_linerange() + geom_point() + 
  scale_colour_manual(values=c("red3", "gray50"), guide=F) +
  labs(x="", y="Proportional taxonomic bias (D)") +
  theme(axis.text.x=element_text(angle=270, vjust=0.5, hjust=0, size=7,
                                 colour=D.cols)) +
  facet_wrap(~model, nrow=2)
ggsave("eda/all_D.pdf", width=7, height=9)



all.pars$lLAMBDA %>% group_by(model, set, source, site, el) %>% 
  summarise(S=sum(1-exp(-exp(mn)) > 0.95)) %>% 
  ggplot(aes(el, S, colour=model)) + ylim(0, d.ls$S) +
  geom_point(aes(shape=set), alpha=0.2) + facet_grid(.~source) + 
  labs(x="Elevation (m)", y="Predicted richness (1km2)") + 
  scale_colour_manual(values=mod_cols) + stat_smooth(method="loess", se=F) + 
  scale_shape_manual(values=c(1, 19))

all.pars$prPres %>% #filter(model!="W_Pois") %>%
  mutate(se=replace(se, is.nan(se), 0)) %>%
  group_by(model, set, source, site, el) %>% 
  summarise(S.lo=sum(q25 > 0.95),
            S.mn=sum(q50 > 0.95),
            S.hi=sum(q75 > 0.95)) %>%
  # summarise(S.lo=sum((mn-2*se) > 0.95),
  #           S.mn=sum(mn > 0.95),
  #           S.hi=sum((mn+2*se) > 0.95)) %>%
  ggplot(aes(el, S.mn, colour=model, fill=model)) + 
  geom_point(size=0.75, shape=1) + ylim(0, d.ls$S) +
  facet_wrap(.~model, ncol=2) + 
  stat_smooth(method="loess", se=F, size=0.75) + 
  stat_smooth(aes(y=S.lo), method="loess", se=F, linetype=2, size=0.5) + 
  stat_smooth(aes(y=S.hi), method="loess", se=F, linetype=2, size=0.5) + 
  # geom_linerange(aes(ymin=S.lo, ymax=S.hi), alpha=0.3) +
  # geom_point(data=filter(obs, source=="Y"), aes(el, S.obs), 
  #            colour="black", alpha=0.5, fill="black") +
  labs(x="Elevation (m)", y="Mean ± 2SE predicted richness (1km2)") + 
  scale_colour_manual(values=mod_cols) 
ggsave("eda/all_S_CI.pdf", width=8, height=10)





library(viridis)
left_join(filter(d.i$grd_W.sf, inbd), all.pars$H, by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=mn)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("H", na.value="white", option="B", limits=c(0, NA))

left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            mutate(se=replace(se, is.nan(se), 0)) %>%
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(mn-2*se > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("S 2.5%", na.value="white", option="B", limits=c(0, d.ls$S))
left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            mutate(se=replace(se, is.nan(se), 0)) %>%
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(mn+2*se > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("S 97.5%", na.value="white", option="B", limits=c(0, d.ls$S))
left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            mutate(se=replace(se, is.nan(se), 0)) %>%
            group_by(model, set, source, site, el, id) %>% 
            summarise(S_lo=sum(mn-2*se > 0.95),
                      S_hi=sum(mn+2*se > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S_hi-S_lo)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("S CI Width", na.value="white", option="D", limits=c(0, NA))
left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(mn > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("S mean", na.value="white", option="B", limits=c(0, d.ls$S))

left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(S=sum(q50 > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=S)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("S (pPr)", na.value="white", option="B", limits=c(0, d.ls$S))

left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$lLAMBDA %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(N=median(mn)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=exp(N))) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("med(N)", na.value="white", option="B")

left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$lLAMBDA %>% 
            group_by(model, set, source, site, el, id) %>% 
            summarise(N=sum(mn)), 
          by="id") %>% 
  filter(!is.na(model)) %>%
  ggplot(aes(fill=N)) + geom_sf(colour=NA) + facet_wrap(~model) + 
  scale_fill_viridis("sum(N)", na.value="white", option="B")

all.pars$lLAMBDA %>% filter(model!="W_Pois") %>%
  group_by(model, set, source, site, el, id) %>% 
  summarise(N=median(mn)) %>%
  ggplot(aes(el, N, colour=model)) + geom_point(alpha=0.6, shape=1) +
  scale_colour_manual(values=mod_cols) 

left_join(filter(d.i$grd_W.sf, inbd), 
          all.pars$prPres %>% 
            mutate(se=replace(se, is.nan(se), 0)) %>%
            group_by(model, set, source, site, el, id) %>% 
            summarise(S_lo=sum(mn-2*se > 0.95),
                      S_hi=sum(mn+2*se > 0.95)), 
          by="id") %>% 
  filter(!is.na(model)) %>% filter(model!="W_Pois") %>%
  ggplot(aes(el, S_hi-S_lo, colour=model)) + geom_point(alpha=0.6, shape=1) +
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~model)





all.pars$prPres %>% filter(el>2000 & source=="Y" & mn>0.975) %>% 
  mutate(SpName=d.i$tax_i$species[as.numeric(spp)]) %>%
  group_by(model) %>% select(model, SpName) %>% print.AsIs




l.p <- ggplot(filter(all.pars$lLAMBDA, source=="Y"),
              aes(x=el, y=exp(mn), colour=model)) + 
  geom_point(alpha=0.5) + facet_wrap(~spp, scales="free_y") +
  stat_smooth(method="loess", se=F, size=1) +
  scale_colour_manual(values=mod_cols)
ggsave("eda/lLAMBDA_spp.jpg", l.p, width=17, height=15)





ggplot(filter(obs, S.obs>0), aes(el, S.obs)) + 
  geom_point(alpha=0.5) + 
  stat_smooth(method="loess", se=F, colour="gray30") + 
  facet_grid(.~source) + labs(x="Elevation (m)", y="Richness (observed)")
ggsave("eda/obs_S.pdf", width=8.5, height=4)

ggplot(filter(obs, S.obs>0), aes(el, H.obs)) + 
  geom_point(alpha=0.5) + 
  stat_smooth(method="loess", se=F, colour="gray30") + 
  facet_grid(.~source) + labs(x="Elevation (m)", y="Shannon H (observed)")
ggsave("eda/obs_H.pdf", width=8.5, height=4)

ggplot(filter(obs, S.obs>0), aes(el, N.obs)) + 
  geom_point(alpha=0.5) + 
  stat_smooth(method="loess", se=F, colour="gray30") + 
  facet_grid(.~source) + labs(x="Elevation (m)", y="Total detections (observed)")
ggsave("eda/obs_N.pdf", width=8.5, height=4)



all.pars$llambda %>% group_by(model, set, source, plot, el) %>% 
  summarise(S=sum(exp(mn) > 0.5)) %>% 
  ggplot(aes(el, S, colour=model)) + 
  geom_point(aes(shape=set)) + facet_wrap(~model, scales="free") + 
  labs(x="Elevation (m)", y="Predicted richness (0.75m2)") + 
  scale_colour_manual(values=mod_cols) + stat_smooth(method="loess", se=F) + 
  scale_shape_manual(values=c(1, 19))

all.pars$llambda %>% group_by(model, set, source, plot, el) %>% 
  summarise(S.mn=sum(1-exp(-exp(mn)) > 0.5),
            S.025=sum(1-exp(-exp(q025)) > 0.5),
            S.25=sum(1-exp(-exp(q25)) > 0.5),
            S.50=sum(1-exp(-exp(q50)) > 0.5),
            S.75=sum(1-exp(-exp(q75)) > 0.5),
            S.975=sum(1-exp(-exp(q975)) > 0.5)) %>% 
  ggplot(aes(el, S.50, colour=model, fill=model)) + 
  # ylim(0, 7) +
  geom_point(aes(shape=set)) + 
  stat_smooth(method="loess", se=F, size=0.75) + 
  stat_smooth(aes(y=S.025), method="loess", se=F, linetype=2, size=0.5) + 
  stat_smooth(aes(y=S.975), method="loess", se=F, linetype=2, size=0.5) + 
  facet_wrap(~model, scales="free") + 
  labs(x="Elevation (m)", y="Predicted richness (.75m2)") + 
  scale_colour_manual(values=mod_cols) + stat_smooth(method="loess", se=F) + 
  scale_shape_manual(values=c(1, 19))

all.pars$llambda %>% filter(!grepl("W_", model)) %>%
  group_by(model, plot, set, el, source) %>%
  summarise(tot=median(mn), tot_q025=sum(q025), tot_q975=sum(q975)) %>%
  ggplot(aes(el, exp(tot), colour=model)) + 
  geom_point(alpha=0.3, aes(shape=set)) + 
  # geom_linerange(aes(ymin=tot_q025, ymax=tot_q975), size=0.2, alpha=0.75) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + #coord_trans(y="exp") + 
  labs(x="Elevation (m)", y=expression(Predicted~total~log(lambda)~(.75~m^2)))


plotS.obs <- tibble(model="obs", el=d.i$V[,2], obs=rowSums(d.i$Y>0))
all.pars$llambda %>% group_by(model, plot, el) %>% 
  summarise(S.mn=sum(1-exp(-exp(mn)) > 0.5)) %>% ungroup %>%
  select(model, el, S.mn) %>% 
  mutate(obs=rep(plotS.obs$obs, times=length(mods))) %>%
  ggplot(aes(obs, S.mn)) + geom_abline() + 
  geom_point(alpha=0.5, shape=1) +
  facet_wrap(~model, scales="free") + 
  labs(x="Observed richness (.75m2)", y="Predicted richness (.75m2)") + 
  stat_smooth(method="lm", se=F) 


all.pars$llambda %>% 
  mutate(obs=rep(c(t(d.i$Y)), times=length(mods))) %>%
  ggplot(aes(exp(mn), obs)) + geom_abline() + 
  geom_point(alpha=0.3, size=0.75, shape=1) + stat_smooth(method="lm") +
  facet_wrap(~model, scales="free_x")
all.pars$llambda %>% 
  mutate(obs=rep(c(t(d.i$Y)), times=length(mods))) %>%
  ggplot(aes(factor(obs), exp(mn))) + 
  geom_boxplot() +
  facet_wrap(~model, scales="free_y")


all.pars$llambda %>% 
  mutate(obs=rep(c(t(d.i$Y)), times=length(mods)),
         prPois=dpois(obs, exp(mn), log=F)) %>%
  group_by(model) %>% filter(obs>0) %>%
  # summarise(prMod=mean(prPois)) %>% arrange(prMod)
  ggplot(aes(prPois, colour=model)) + geom_density() + 
  scale_colour_manual(values=mod_cols) 

all.pars$llambda %>% 
  mutate(obs=rep(c(t(d.i$Y)), times=length(mods)),
         prPois=dpois(obs, exp(mn), log=F)) %>%
  filter(obs > 0) %>%
  group_by(model) %>% 
  ggplot(aes(factor(obs), prPois, colour=model)) + 
  geom_boxplot(outlier.size=0.5,
               outlier.alpha=0.3, outlier.shape=1) +
  scale_colour_manual(values=mod_cols) +
  labs(x=expression("Observed colony abundance"~(Y[is])), 
       y=expression(Pr(Y[is]~'|'~lambda)))

plot(d.i$V[,2], rowSums(d.i$Y))
plot(d.i$V[,2], rowSums(d.i$Y>0))

all.pars$llambda %>% 
  mutate(obs=rep(c(t(d.i$Y)), times=length(mods)),
         prPres=1-exp(-exp(mn)),
         prP.lo=1-exp(-exp(q025)),
         prP.hi=1-exp(-exp(q975))) %>%
  filter(model != "Y_CMP") %>% filter(model != "W_Pois") %>%
  # filter(model=="Y_ZIP") %>%
  ggplot(aes(x=prPres, y=as.numeric(obs>0))) + 
  # geom_point(alpha=0.25, shape=1, size=0.75) + 
  stat_smooth(method="glm", size=0.5, 
              method.args=list(family="binomial"), se=F) + 
  stat_smooth(aes(x=prP.lo), method="glm", size=0.5, 
              method.args=list(family="binomial"), linetype=2, se=F) + 
  stat_smooth(aes(x=prP.hi), method="glm", size=0.5, 
              method.args=list(family="binomial"), linetype=2, se=F) + 
  facet_wrap(~model) + ylim(0,1) + xlim(0,1)

all.pars$lLAMBDA %>% 
  mutate(obs=rep(c(t(rbind(d.i$W, Y.J, d.i$W_))), times=length(mods)),
         prPres=1-exp(-exp(mn))) %>%
  filter(model != "Y_CMP") %>%
  filter(source=="Y") %>%
  ggplot(aes(x=prPres, y=as.numeric(obs>0))) + 
  geom_point(alpha=0.25, shape=1, size=0.75) + 
  stat_smooth(method="glm", size=0.5, method.args=list(family="binomial")) + 
  facet_wrap(~model) + ylim(0,1)






iter.S <- map(out.raw, 
              ~t(as.matrix(.)[,grepl("prPres\\[", colnames(as.matrix(.)))]>0.95))
iter.S.df <- map_dfr(seq_along(iter.S),
                 ~tibble(Pres=c(iter.S[[.x]]),
                         iter=rep(1:ncol(iter.S[[.x]]), 
                                  each=nrow(iter.S[[.x]])),
                         Param=rep(rownames(iter.S[[.x]]), 
                                   times=ncol(iter.S[[.x]])), 
                         par=str_split_fixed(Param, "\\[", 2)[,1],
                         site=str_split_fixed(
                           str_split_fixed(Param, "\\[", 2)[,2], 
                           ",", 2)[,1],
                         spp=str_sub(
                           str_split_fixed(Param, ",", 2)[,2], 
                           1L, -2L),
                         mod=mods[.x]) %>%
                   group_by(mod, iter, par, site) %>%
                   summarise(S=sum(Pres)) %>%
                   mutate(el=d.i$X[as.numeric(site),"el"])) 


iter.S_ <- map(out.raw, 
              ~t(as.matrix(.)[,grepl("prPres_", colnames(as.matrix(.)))]>0.95))
iter.S_.df <- map_dfr(seq_along(iter.S_),
                     ~tibble(Pres=c(iter.S_[[.x]]),
                             iter=rep(1:ncol(iter.S_[[.x]]), 
                                      each=nrow(iter.S_[[.x]])),
                             Param=rep(rownames(iter.S_[[.x]]), 
                                       times=ncol(iter.S_[[.x]])), 
                             par=str_split_fixed(Param, "\\[", 2)[,1],
                             site=str_split_fixed(
                               str_split_fixed(Param, "\\[", 2)[,2], 
                               ",", 2)[,1],
                             spp=str_sub(
                               str_split_fixed(Param, ",", 2)[,2], 
                               1L, -2L),
                             mod=mods[.x]) %>%
                       group_by(mod, iter, par, site) %>%
                       summarise(S=sum(Pres)) %>%
                       mutate(el=d.i$X_[as.numeric(site),"el"])) 



ggplot(iter.S.df, aes(x=el, y=S, group=iter, colour=mod)) + ylim(0, 80) +
  # geom_point(alpha=0.01) + 
  facet_wrap(~mod) +
  scale_colour_manual(values=mod_cols) +
  geom_line(stat="smooth", alpha=0.05, method="lm", formula=y~x+I(x^2)+I(x^3)) +
  labs(x="Elevation (m)", y="Predicted richness (1 km2)")

ggplot(iter.S_.df, aes(x=el, y=S, group=iter, colour=mod)) + ylim(0, 80) +
  # geom_point(alpha=0.01) + 
  facet_wrap(~mod) +
  scale_colour_manual(values=mod_cols) +
  geom_line(stat="smooth", alpha=0.05, method="lm", formula=y~x+I(x^2)+I(x^3)) +
  labs(x="Elevation (m)", y="Predicted richness (1 km2)")









l_glm <- map(mods, 
    ~glm(obsPres ~ prPres, 
         data=filter(all.pars$llambda %>% 
                       mutate(obsPres=rep(c(t(d.i$Y)), times=length(mods))>0,
                              prPres=1-exp(-exp(mn))), 
                     model==.x), 
         family="binomial"))
tibble(model=mods, 
       aic=map_dbl(l_glm, ~.$aic),
       deviance=map_dbl(l_glm, ~.$deviance)) %>% arrange(aic)

r_glm <- map(mods, 
             ~glm(obsPres ~ prPres, 
                  data=filter(all.pars$lLAMBDA %>% 
                                mutate(obsPres=rep(c(t(rbind(d.i$W, Y.J, Y.J_))), times=length(mods))>0,
                                       prPres=1-exp(-exp(mn))), 
                              model==.x) %>% filter(source=="Y"), 
                  family="binomial"))
tibble(model=mods, 
       aic=map_dbl(r_glm, ~.$aic),
       deviance=map_dbl(r_glm, ~.$deviance)) %>% arrange(aic)





















# generate some predicted W & Y...
site.sf <- agg_str_site_data()
library(viridis)
all.pars$lLAMBDA %>% filter(source=="Y") %>% 
  rename(BDM=id) %>%
  group_by(BDM, model) %>%
  summarise(S.pr=sum(1-exp(-exp(mn)) > 0.95),
            S.ct=sum(mn>4)) %>%
  full_join(site.sf, ., by="BDM") %>%
  ggplot(aes(fill=S.ct, colour=S.ct)) + geom_sf() + 
  scale_fill_viridis(option="B") + scale_colour_viridis(option="B") + 
  facet_wrap(~model)




S = cA^z
log(S) = log(c) + z*log(A)