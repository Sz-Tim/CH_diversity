
library(tidyverse); library(rstan); library(ggmcmc)
theme_set(theme_bw())

d.ls <- readRDS("data/stan_data/test_realData_ls.rds")
d.i <- readRDS("data/stan_data/test_realData_i.rds")

mods <- c("W", "Y", "WY")
pars.any <- c("H", "lLAMBDA", "beta", "alpha", "eta", "b", "a", "B", "A")
out.pars <- imap(setNames(pars.any, pars.any), ~vector("list", length(mods)))

for(i in seq_along(mods)) {
  out <- readRDS(paste0("out/tests/test_", mods[i], ".rds"))
  cat("loaded output for", mods[i], "\n")
  # out.mx <- as.matrix(out)
  pars <- unique(str_split_fixed(names(out), "\\[", 2)[,1])
  pars <- pars[grep("lp__", pars, value=F, invert=T)]
  
  out.gg.sum <- map(pars, ~ggs(out, paste0("^", .x, "\\[")) %>%
                      group_by(Parameter) %>% 
                      summarise(mn=mean(value), 
                                q025=quantile(value, 0.025),
                                q25=quantile(value, 0.25),
                                q75=quantile(value, 0.75),
                                q975=quantile(value, 0.975)) %>%
                      mutate(model=mods[i])) %>% setNames(pars)
  rm(out)
  saveRDS(out.gg.sum, paste0("out/tests/ggsum_", mods[i], ".rds"))
  # out.gg <- map(pars, ~ggs(out, paste0("^", .x, "\\["))) %>% setNames(pars)
  # rm(out)
  # out.gg.sum <- map(out.gg, ~.x %>% group_by(Parameter) %>% 
  #                     summarise(mn=mean(value), 
  #                               q025=quantile(value, 0.025),
  #                               q25=quantile(value, 0.25),
  #                               q75=quantile(value, 0.75),
  #                               q975=quantile(value, 0.975)) %>%
  #                     mutate(model=mods[i]))
  
  
  # add site info, need to align 'elevation' with shannon h, lambdas
  out.pars$H[[i]] <- rbind(
    out.gg.sum$ShannonH %>% 
      mutate(site=str_remove(str_split_fixed(Parameter,  "\\[", 2)[,2], "]"),
             site=as.numeric(site),
             par=str_split_fixed(Parameter, "\\[", n=2)[,1],
             set="train") %>%
      arrange(site) %>%
      mutate(id=d.i$X[,"id"], el=d.i$X[,"el"], 
             source=ifelse(id<100000, "W", "Y")),
    out.gg.sum$ShannonH_ %>% 
      mutate(site=str_remove(str_split_fixed(Parameter, "\\[", 2)[,2], "]"),
             site=as.numeric(site),
             par=str_split_fixed(Parameter, "\\[", n=2)[,1],
             set="test") %>%
      arrange(site) %>%
      mutate(id=d.i$X_[,"id"], el=d.i$X_[,"el"], 
             source=ifelse(id<100000, "W", "Y"))
  )
  out.pars$lLAMBDA[[i]] <- rbind(
    out.gg.sum$lLAMBDA %>%
      mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             par=str_split_fixed(Parameter, "\\[", n=2)[,1],
             set="train") %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
             source=ifelse(id<100000, "W", "Y")),
    out.gg.sum$lLAMBDA_ %>%
      mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2),
             par=str_split_fixed(Parameter, "\\[", n=2)[,1],
             set="test") %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
             source=ifelse(id<100000, "W", "Y"))
  )
  out.pars$beta[[i]] <- out.gg.sum$beta %>%
    mutate(ParName=c("intercept", colnames(d.i$X)[-(1:2)]))
  out.pars$B[[i]] <- out.gg.sum$B %>%
    mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
           GenNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
    mutate(ParNum=as.numeric(ParNum),
           GenNum=as.numeric(GenNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)])[ParNum],
           GenName=unique(d.i$tax_i$genus)[GenNum])
  out.pars$b[[i]] <- out.gg.sum$b %>%
    mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
           SpNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
    mutate(ParNum=as.numeric(ParNum),
           SpNum=as.numeric(SpNum),
           ParName=c("intercept", colnames(d.i$X)[-(1:2)])[ParNum],
           SpName=d.i$tax_i$species[SpNum],
           GenName=str_split_fixed(SpName, "_", 2)[,1])
  if("alpha" %in% pars) {
    out.pars$alpha[[i]] <- out.gg.sum$alpha %>%
      mutate(ParName=colnames(d.i$V)[-(1:2)])
    out.pars$A[[i]] <- out.gg.sum$A %>%
      mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
             GenNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
      mutate(ParNum=as.numeric(ParNum),
             GenNum=as.numeric(GenNum),
             ParName=c(colnames(d.i$V)[-(1:2)])[ParNum],
             GenName=unique(d.i$tax_i$genus)[GenNum])
    out.pars$a[[i]] <- out.gg.sum$a %>%
      mutate(ParNum=str_sub(str_split_fixed(Parameter, ",", 2)[,1], -1, -1),
             SpNum=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
      mutate(ParNum=as.numeric(ParNum),
             SpNum=as.numeric(SpNum),
             ParName=c(colnames(d.i$V)[-(1:2)])[ParNum],
             SpName=d.i$tax_i$species[SpNum],
             GenName=str_split_fixed(SpName, "_", 2)[,1])
  }
  if("eta" %in% pars) {
    out.pars$eta[[i]] <- out.gg.sum$eta %>%
      mutate(ParName=c("intercept", colnames(d.i$U)[-1]))
  }
  cat("Finished", mods[i], "\n")
}

all.pars <- map(out.pars, ~do.call('rbind', .))
saveRDS(all.pars, "out/tests/allPars2.rds")


mod_cols <- c(W="#33a02c", Y="#1f78b4", WY="#e31a1c")

ggplot(all.pars$H, aes(x=el, y=mn, ymin=q25, ymax=q75, colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  geom_linerange(size=0.2, alpha=0.75) + 
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + 
  labs(x="Elevation (m)", y=expression(Predicted~Shannon~H~(1~km^2)))
ggsave("eda/all_H.pdf", width=8.5, height=4)

all.pars$lLAMBDA %>% group_by(model, site, set, id, el, source) %>%
  summarise(tot=sum(mn), tot_q25=sum(q25), tot_q75=sum(q75)) %>%
  ggplot(aes(el, tot, ymin=tot_q25, ymax=tot_q75, colour=model)) + 
  geom_point(alpha=0.75, aes(shape=set)) + 
  geom_linerange(size=0.2, alpha=0.75) + 
  scale_colour_manual(values=mod_cols) + scale_shape_manual(values=c(1,19)) + 
  facet_grid(.~source) + #coord_trans(y="exp") + 
  labs(x="Elevation (m)", y=expression(Predicted~total~log(Lambda)~(1~km^2)))
ggsave("eda/all_LAMBDA.pdf", width=8.5, height=4)

ggplot(filter(all.pars$beta, ParName!="intercept"), 
       aes(x=model, y=mn, ymin=q25, ymax=q75, colour=model)) + 
  geom_point(data=filter(all.pars$b, ParName !="intercept"), 
             size=0.5, alpha=0.5) +
  geom_point(data=filter(all.pars$B, ParName !="intercept"), 
             size=1.5, shape=1) +
  geom_point(size=3.5, shape=1) + 
  geom_linerange() + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName) +
  labs(x="Model", y=expression(Slope~posterior~distribution~(beta~B~b)))
ggsave("eda/all_beta.pdf", width=5, height=4)

ggplot(filter(all.pars$eta, ParName!="intercept"), 
       aes(x=model, y=mn, ymin=q25, ymax=q75, colour=model)) + 
  geom_point(size=3.5, shape=1) + 
  geom_linerange() + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName) +
  labs(x="Model", y=expression(Slope~posterior~distribution~(eta)))
ggsave("eda/all_eta.pdf", width=4, height=4)

ggplot(filter(all.pars$alpha, ParName!="intercept"), 
       aes(x=model, y=mn, ymin=q25, ymax=q75, colour=model)) + 
  geom_point(data=filter(all.pars$a, ParName !="intercept"), 
             size=0.5, alpha=0.5) +
  geom_point(data=filter(all.pars$A, ParName !="intercept"), 
             size=1.5, shape=1) +
  geom_point(size=3.5, shape=1) + 
  geom_linerange() + 
  geom_hline(yintercept=0, linetype=2) +
  scale_colour_manual(values=mod_cols) + facet_wrap(~ParName) +
  labs(x="Model", y=expression(Slope~posterior~distribution~(alpha~A~a)))
ggsave("eda/all_alpha.pdf", width=4, height=4)



Y.J <- matrix(0, ncol=d.ls$S, nrow=d.ls$J)
Y.J_ <- matrix(0, ncol=d.ls$S, nrow=d.ls$J_)
for(i in 1:d.ls$J) {
  Y.J[i,] <- colSums(d.ls$Y[d.ls$IJ==i,])
}
for(i in 1:d.ls$J_) {
  Y.J_[i,] <- colSums(d.ls$Y_[d.ls$IJ_==i,])
}


obs <- data.frame(el=c(d.i$X[,2], d.i$X_[,2]), 
                  H=vegan::diversity(rbind(d.ls$W, Y.J, d.ls$W_, Y.J_)),
                  tot=rowSums(rbind(d.ls$W, Y.J, d.ls$W_, Y.J_)),
                  set=rep(c("train", "test"), times=c(nrow(d.i$X), nrow(d.i$X_))),
                  source=rep(c("W", "Y", "W", "Y"), 
                             times=c(d.ls$K, d.ls$J, d.ls$K_, d.ls$J_))) %>%
  filter(H>0)
ggplot(obs, aes(el, H)) + geom_point(alpha=0.5) + facet_grid(.~source) +
  stat_smooth(method="loess", se=F, colour="gray30", span=1) +
  labs(x="Elevation (m)", y="Shannon H (observed)")
ggsave("eda/all_obsH.pdf", width=8.5, height=4)
ggplot(obs, aes(el, tot)) + geom_point(alpha=0.5) + facet_grid(.~source) +
  labs(x="Elevation (m)", y="Total abundance (observed)")
ggsave("eda/all_obsLAMBDA.pdf", width=8.5, height=4)
