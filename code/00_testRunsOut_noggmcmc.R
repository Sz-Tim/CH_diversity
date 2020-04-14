

library(tidyverse); library(rstan)
theme_set(theme_bw())

d.ls <- readRDS("data/stan_data/test_realData_ls.rds")
d.i <- readRDS("data/stan_data/test_realData_i.rds")

mods <- c("W", "WY")
pars.any <- c("H", "lLAMBDA", "beta", "alpha", "eta", "b", "a", "B", "A", "D")
out.pars <- imap(setNames(pars.any, pars.any), ~vector("list", length(mods)))

Y.J <- matrix(0, ncol=d.ls$S, nrow=d.ls$J)
Y.J_ <- matrix(0, ncol=d.ls$S, nrow=d.ls$J_)
for(i in 1:d.ls$J) {
  Y.J[i,] <- colSums(d.ls$Y[d.ls$IJ==i,])
}
for(i in 1:d.ls$J_) {
  Y.J_[i,] <- colSums(d.ls$Y_[d.ls$IJ_==i,])
}
obs <- data.frame(el=c(d.i$X[,2], d.i$X_[,2]), 
                  H.obs=vegan::diversity(rbind(d.ls$W, Y.J, d.ls$W_, Y.J_)),
                  tot.obs=rowSums(rbind(d.ls$W, Y.J, d.ls$W_, Y.J_)),
                  site=c(1:nrow(d.i$X), 1:nrow(d.i$X_)),
                  set=rep(c("train", "test"), 
                          times=c(nrow(d.i$X), nrow(d.i$X_))),
                  source=rep(c("W", "Y", "W", "Y"), 
                             times=c(d.ls$K, d.ls$J, d.ls$K_, d.ls$J_))) %>%
  mutate(set=as.character(set), source=as.character(source))

for(i in seq_along(mods)) {
  out <- readRDS(paste0("out/tests/testSum_2_", mods[i], ".rds"))
  pars <- unique(str_split_fixed(rownames(out), "\\[", 2)[,1])
  pars <- pars[grep("lp__", pars, value=F, invert=T)]
  
  out.ls <- data.frame(Parameter=rownames(out),
                       mn=out[,1],
                       se=out[,2],
                       q025=out[,4], 
                       q25=out[,5],
                       q75=out[,7],
                       q975=out[,8],
                       model=mods[i],
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
    mutate(ParName=c("intercept", colnames(d.i$X)[-(1:2)])) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$B[[i]] <- out.ls$B %>%
    mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
           GenNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
    mutate(ParNum=as.numeric(ParNum),
           GenNum=as.numeric(GenNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)])[ParNum],
           GenName=unique(d.i$tax_i$genus)[GenNum]) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  out.pars$b[[i]] <- out.ls$b %>%
    mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
           SpNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
    mutate(ParNum=as.numeric(ParNum),
           SpNum=as.numeric(SpNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)])[ParNum],
           SpName=d.i$tax_i$species[SpNum],
           GenName=str_split_fixed(SpName, "_", 2)[,1]) %>% 
    mutate(Parameter=as.character(Parameter), model=as.character(model))
  if("alpha" %in% pars) {
    out.pars$alpha[[i]] <- out.ls$alpha %>%
      mutate(ParName=colnames(d.i$V)[-(1:2)]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    out.pars$A[[i]] <- out.ls$A %>%
      mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
             GenNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
      mutate(ParNum=as.numeric(ParNum),
             GenNum=as.numeric(GenNum),
             ParName=c(colnames(d.i$V)[-(1:2)])[ParNum],
             GenName=unique(d.i$tax_i$genus)[GenNum]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    out.pars$a[[i]] <- out.ls$a %>%
      mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
             SpNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
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
}

all.pars <- map(out.pars, ~do.call('rbind', .))
saveRDS(all.pars, "out/tests/allPars.rds")


mod_cols <- c(W="#33a02c", Y="#1f78b4", WY="#e31a1c")

ggplot(all.pars$H, aes(x=el, y=mn, ymin=q025, ymax=q975, colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  geom_linerange(size=0.2, alpha=0.75) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + 
  labs(x="Elevation (m)", y=expression(Predicted~Shannon~H~(1~km^2)))
ggsave("eda/all_H.pdf", width=8.5, height=4)

all.pars$lLAMBDA %>% group_by(model, site, set, id, el, source) %>%
  summarise(tot=sum(mn), tot_q025=sum(q025), tot_q975=sum(q975)) %>%
  ggplot(aes(el, tot, ymin=tot_q025, ymax=tot_q975, colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  geom_linerange(size=0.2, alpha=0.75) + 
  stat_smooth(method="loess", se=F) +
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + #coord_trans(y="exp") + 
  labs(x="Elevation (m)", y=expression(Predicted~total~log(Lambda)~(1~km^2)))
ggsave("eda/all_LAMBDA.pdf", width=8.5, height=4)

ggplot(filter(all.pars$beta, ParName!="intercept"), 
       aes(x=model, y=mn, colour=model)) + 
  # geom_point(data=filter(all.pars$b, ParName !="intercept"), 
  #            size=0.5, alpha=0.3) +
  # geom_point(data=filter(all.pars$B, ParName !="intercept"), 
  #            size=1.5, shape=1, alpha=0.5) +
  geom_point(size=3.5, shape=1) + 
  geom_linerange(aes(ymin=q025, ymax=q975)) + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName) +
  labs(x="Model", y=expression(Slope~posterior~distribution~(beta~B~b)))
ggsave("eda/all_beta.pdf", width=5, height=4)

ggplot(filter(all.pars$eta, ParName!="intercept"), 
       aes(x=model, y=mn, ymin=q025, ymax=q975, colour=model)) + 
  geom_point(size=3.5, shape=1) + 
  geom_linerange(colour="gray30") + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName) +
  labs(x="Model", y=expression(Slope~posterior~distribution~(eta)))
ggsave("eda/all_eta.pdf", width=4, height=4)

ggplot(filter(all.pars$alpha, ParName!="intercept"), 
       aes(x=model, y=mn, ymin=q025, ymax=q975, colour=model)) + 
  # geom_point(data=filter(all.pars$a, ParName !="intercept"), 
             # size=0.5, alpha=0.3) +
  # geom_point(data=filter(all.pars$A, ParName !="intercept"), 
             # size=1.5, shape=1, alpha=0.5) +
  geom_point(size=3.5, shape=1) + 
  geom_linerange(colour="gray30") + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName, scales="free_y") +
  labs(x="Model", y=expression(Slope~posterior~distribution~(alpha~A~a)))
ggsave("eda/all_alpha.pdf", width=4, height=4)


D.cols <- c("red3", "gray30")[(sign(all.pars$D$q025-1) != sign(all.pars$D$q975-1))+1]
ggplot(all.pars$D, aes(x=SpName, y=mn, ymin=q025, ymax=q975,
                       colour=sign(q025-1) != sign(q975-1))) + 
  geom_hline(yintercept=1, linetype=2) +
  geom_linerange() + geom_point() + 
  scale_colour_manual(values=c("red3", "gray50"), guide=F) +
  labs(x="", y="Proportional taxonomic bias (D)") +
  theme(axis.text.x=element_text(angle=270, vjust=0.5, hjust=0, size=7,
                                 colour=D.cols))
ggsave("eda/all_D.pdf", width=7, height=5)


all.pars$H %>% left_join(., obs, by=c("site", "el", "set", "source")) %>%
  ggplot(aes(x=H.obs, y=mn, ymin=q025, ymax=q975, colour=model)) + 
  stat_smooth(method="lm") +
  geom_linerange() + geom_point() + facet_grid(.~source)


ggplot(all.pars$lLAMBDA, aes(x=mn, colour=el, group=el)) + 
  geom_density() + facet_wrap(~model, scales="free")

all.pars$lLAMBDA %>% group_by(model, set, source, site, el) %>% 
  summarise(S=sum(mn > 0)) %>% 
  ggplot(aes(el, S, colour=model)) + 
  geom_point(aes(shape=set)) + facet_grid(.~source) + 
  labs(x="Elevation (m)", y="Predicted richness (1km2)") + 
  scale_colour_manual(values=mod_cols) + stat_smooth(method="loess", se=F) + 
  scale_shape_manual(values=c(1, 19))


all.pars$lLAMBDA %>% group_by(model, set, source, site, el) %>% 
  summarise(S=sum(1-exp(-exp(mn)) > 0.95)) %>% 
  ggplot(aes(el, S, colour=model)) + ylim(0, d.ls$S) + 
  geom_point(aes(shape=set)) + facet_grid(.~source) + 
  labs(x="Elevation (m)", y="Predicted richness (1km2)") + 
  scale_colour_manual(values=mod_cols) + stat_smooth(method="loess", se=F) + 
  scale_shape_manual(values=c(1, 19))
