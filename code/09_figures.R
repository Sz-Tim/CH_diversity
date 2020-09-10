library(tidyverse); theme_set(theme_bw() + theme(panel.grid=element_blank()))

d.i <- readRDS("data/stan_data/vs_30_i.rds")
agg_vs <- readRDS("out/agg_vs.rds")


mod_col <- c("WY"="#7b3294", "Y"="#008837")

filter(agg_vs$beta, ParName != "intercept") %>%
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  mutate(Scale=str_sub(ParName, -1L, -1L),
         Par=str_sub(ParName, 1L, -3L),
         model=str_split_fixed(model, "_", 4)[,2]) %>%
  mutate(Scale=factor(Scale, levels=c("R", "L"), 
                      labels=c("Regional", "Local"))) %>%
  ggplot(aes(x=ParName, y=mn, ymin=q05, ymax=q95, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(position=position_dodge(width=0.5)) + 
  facet_grid(Scale~., scales="free_y") +
  coord_flip() + theme(panel.grid=element_blank()) + 
  labs(x="", y="beta") + 
  scale_colour_manual(values=mod_col)
ggsave("eda/beta_fullModel.png", width=6, height=6)

ggplot(filter(agg_vs$B, ParName != "intercept"), 
       aes(x=ParName, y=q50, shape=sign(q05)==sign(q95), colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_point(alpha=0.8) + 
  scale_shape_manual(values=c(1, 19)) +
  coord_flip() + facet_wrap(~genName)

ggplot(filter(agg_vs$b, ParName != "intercept") %>%
         filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")), 
       aes(x=ParName, y=mn, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_boxplot() + 
  # geom_point(aes(shape=sign(q05)==sign(q95)), alpha=0.5) + 
  scale_shape_manual(values=c(1, 19)) +
  coord_flip() + facet_wrap(~model)

agg_vs$B %>% filter(ParName != "intercept") %>%
  mutate(effect=case_when(q05<0 & q95<0 ~ "negative",
                          q05>0 & q95>0 ~ "positive",
                          q05<0 & q95>0 ~ "no effect",
                          q05>0 & q95<0 ~ "no effect"),
         sig=sign(q05)==sign(q95)) %>%
  group_by(model, ParName, effect) %>% summarise(nGen=n(), nSig=sum(sig)) %>%
  group_by(model, ParName) %>%
  mutate(effect=factor(effect, levels=c("no effect", "negative", "positive")),
         nSig=sum(nSig)) %>%
  arrange(model, nSig) %>% ungroup %>%
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nGen/21, fill=effect)) + 
  geom_bar(stat="identity", colour="gray30", size=0.25) +
  facet_wrap(~model) + coord_flip() + 
  scale_fill_manual(values=c("gray90", "blue", "red")) + 
  labs(x="", y="Proportion of genera")

filter(agg_vs$B, ParName != "intercept") %>%
  filter(sign(q05)==sign(q95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/21, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of genera", x="")

agg_vs$b %>% filter(ParName != "intercept") %>%
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  mutate(effect=case_when(q05<0 & q95<0 ~ "negative",
                          q05>0 & q95>0 ~ "positive",
                          q05<0 & q95>0 ~ "ns",
                          q05>0 & q95<0 ~ "ns"),
         sig=sign(q05)==sign(q95)) %>%
  mutate(Scale=str_sub(ParName, -1L, -1L),
         Par=str_sub(ParName, 1L, -3L),
         model=str_split_fixed(model, "_", 4)[,2]) %>%
  group_by(model, Par, Scale, ParName, effect) %>% 
  summarise(nSpp=n(), nSig=sum(sig)) %>%
  group_by(model, Par, Scale, ParName) %>%
  mutate(effect=factor(effect, levels=c("ns", "negative", "positive")),
         nSig=sum(nSig)) %>%
  arrange(model, Scale, nSig) %>% ungroup %>%
  mutate(ParName=factor(ParName, levels=unique(ParName)),
         Scale=factor(Scale, levels=c("R", "L"), 
                      labels=c("Regional", "Local"))) %>%
  ggplot(aes(ParName, nSpp/80, fill=effect)) + 
  geom_bar(stat="identity", colour="gray30", size=0.25) +
  facet_grid(Scale~model, scales="free") + coord_flip() + 
  scale_fill_manual("", values=c("gray90", "blue", "red")) + 
  labs(x="", y="Proportion of species")
ggsave("eda/b_bar_fullModel.png", width=8, height=8)

filter(agg_vs$b, ParName != "intercept") %>%
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  filter(sign(q05)==sign(q95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/80, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of species", x="")

agg_vs$b %>% filter(ParName != "intercept") %>%
  mutate(effect=case_when(q05<0 & q95<0 ~ "negative",
                          q05>0 & q95>0 ~ "positive",
                          q05<0 & q95>0 ~ "no effect",
                          q05>0 & q95<0 ~ "no effect"),
         sig=sign(q05)==sign(q95)) %>%
  group_by(model, sppName, sig) %>% summarise(nVar=n()) %>% 
  filter(sig) %>% summary # 1-12 significant variables
  ggplot(aes(sig, nVar, colour=model)) + geom_boxplot()
  
agg_vs$b %>% select(model, sppName, ParName, mn) %>%
  filter(ParName != "intercept") %>%
  pivot_wider(names_from="model", values_from="mn") %>%
  ggplot(aes(vs_Y_altRE_G, vs_WY_altRE_G)) + 
  geom_hline(yintercept=0, colour="gray90") + 
  geom_vline(xintercept=0, colour="gray90") + 
  geom_point(alpha=0.5) + geom_abline() + 
  facet_wrap(~ParName) + ylim(-4,4) + xlim(-4,4)
  
  




agg_vs$lLAM %>% filter(id>1e5) %>% 
    filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
    mutate(model=str_split_fixed(model, "_", 4)[,2]) %>%
    group_by(site, model, el) %>% 
  summarise(n=sum(exp(mn))) %>%
  ggplot(aes(el, n, group=model, colour=model)) + 
  geom_point(alpha=0.5) + stat_smooth(span=2, se=F) +
  labs(x="Elevation (m)", y="Predicted colonies per km2") + 
  scale_colour_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, 1.7e6)) + 
  theme(panel.grid=element_blank())
ggsave("eda/lLAMBDA_fullModel.png", width=6, height=5)

agg_vs$lLAM %>% filter(id>1e5) %>% 
  group_by(site, model, el) %>% 
  summarise(S=sum(1-exp(-exp(q025))>0.95)) %>%
  ggplot(aes(el, S, group=model, colour=model)) + 
  geom_point(alpha=0.5) + stat_smooth(span=2) + ylim(0,80) +
  labs("Elevation (m)", "Predicted richness per km2")

agg_vs$lLAM %>% filter(id>1e5) %>% 
  ggplot(aes(el, q025, colour=model)) + geom_hline(yintercept=0, linetype=3) +
  stat_smooth(se=F) +
  facet_wrap(~sppName, scales="free_y")

agg_vs$lLAM %>% mutate(pP_R=1-exp(-exp(q025))) %>%
  ggplot(aes(x=el, y=pP_R, colour=model)) + ylim(0,1) +
  stat_smooth(method="loess", se=F) + facet_wrap(~sppName)


agg_vs$llam %>% 
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  mutate(model=str_split_fixed(model, "_", 4)[,2]) %>%
  group_by(plot, model, el) %>% 
  summarise(N_hi=sum(exp(q95)),
            N_mn=sum(exp(q50)),
            N_lo=sum(exp(q05))) %>%
  ggplot(aes(el)) + ylim(0, NA) +
  stat_smooth(aes(y=N_mn), size=1, se=F) + 
  stat_smooth(aes(y=N_lo), size=0.5, linetype=2, se=F) + 
  stat_smooth(aes(y=N_hi), size=0.5, linetype=2, se=F) + 
  facet_wrap(~model) + 
  theme(panel.grid=element_blank()) + 
  labs(x="Elevation (m)", y="Mean abundance per 0.75m2")
ggsave("eda/llambda_fullModel.png", width=9, height=4)

agg_vs$pP_L %>% group_by(plot, model, el) %>% 
  summarise(S_hi=sum(q95 > 0.5),
            S_mn=sum(q50 > 0.5),
            S_lo=sum(q05 > 0.5)) %>%
  ggplot(aes(el)) + ylim(0, NA) +
  stat_smooth(aes(y=S_mn), size=0.5, se=F) + 
  stat_smooth(aes(y=S_lo), size=0.25, se=F) + 
  stat_smooth(aes(y=S_hi), size=0.25, se=F) + 
  facet_wrap(~model) + labs("Elevation (m)", "Predicted richness per 0.75m2")

ggplot(agg_vs$pP_L, aes(x=el, y=mn, colour=model)) + 
  stat_smooth(method="loess") + facet_wrap(~sppName, scales="free_y")

agg_vs$D %>%
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  mutate(model=str_split_fixed(model, "_", 4)[,2]) %>%
  ggplot(aes(sppName, q50, ymin=q05, ymax=q95, 
                     colour=sign(q05-1)==sign(q95-1))) + 
  geom_hline(yintercept=1, colour="gray") +
  geom_point() + geom_linerange() + labs(x="", y="Proportional Bias") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        panel.grid=element_blank()) + 
  scale_colour_manual(values=c("gray", "darkred"), guide=F)
ggsave("eda/D_fullModel.png", width=9, height=4)




ggplot(agg_vs$sig_b, aes(ParName, q50, ymin=q05, ymax=q95, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(position=position_dodge(width=0.5)) + 
  coord_flip() + labs(x="Standard deviation among congeners", y="")


agg_vs$Sig_B %>% filter(gen1 != gen2) %>% 
  filter(ParName != "intercept") %>%
  group_by(model, ParName) %>%
  summarise(mnAbsCor=mean(abs(mn)), 
            minAbsCor=min(abs(mn)),
            maxAbsCor=max(abs(mn))) %>%
  ggplot(aes(ParName, mnAbsCor, ymin=minAbsCor, ymax=maxAbsCor, colour=model)) +
  geom_point(position=position_dodge(width=0.5)) + 
  geom_linerange(position=position_dodge(width=0.5))


agg_vs$pP_L$obs <- rep(c(t(d.i$Y)), n_distinct(agg_vs$pP_L$model))
agg_vs$pP_L$disp <- agg_vs$disp$mn[match(agg_vs$pP_L$model, agg_vs$disp$model)]


ggplot(agg_vs$pP_L, aes(mn, obs)) + geom_point(alpha=0.1) + 
  stat_smooth() + facet_wrap(~model)
ggplot(agg_vs$pP_L, aes(mn, as.numeric(obs>0), group=model, colour=model)) + 
  #geom_point(alpha=0.1) + 
  stat_smooth(aes(x=q05), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q50), method="glm", size=0.5, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q95), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  # facet_wrap(~model) +
  ylim(0,1)
ggplot(agg_vs$pP_L, aes(factor(obs), mn, fill=model)) + geom_boxplot()


agg_vs$pP_L <- agg_vs$pP_L %>%
  mutate(ll=LaplacesDemon::dgpois(obs, mn, disp, log=T))
agg_vs$pP_L %>% group_by(model) %>% 
  summarise(nll=-sum(ll)*2) %>% 
  arrange(nll)





lLAMBDA.ls <- map(agg_vs$full, 
                  ~1-exp(-exp(rstan::extract(., pars="lLAMBDA")[[1]]))) %>%
  map(~array(., dim=dim(.), 
             dimnames=list("iter"=1:240, 
                           "id"=d.i$X[,"id"], 
                           "spp"=d.i$tax_i$species)))
lLAMBDA.tb <- map(lLAMBDA.ls, ~cubelyr::as.tbl_cube(.) %>% as_tibble) 
LAM.df <- map(1:length(lLAMBDA.tb), ~lLAMBDA.tb[[.x]] %>% 
                mutate(model=names(lLAMBDA.tb)[[.x]])) %>%
  do.call("rbind", .)

LAM.df %>% filter(id > 1e5) %>%
  filter(model %in% c("vs_WY_altRE_G", "vs_Y_altRE_G")) %>%
  mutate(model=str_split_fixed(model, "_", 4)[,2]) %>%
  group_by(model, id, spp) %>%
  summarise(mnP=mean(.)) %>% 
  mutate(el=d.i$X[,'el'][match(id, d.i$X[,'id'])]) %>%
  group_by(model, el) %>%
  summarise(S=sum(mnP>0.99)) %>%
  ggplot(aes(el, S, colour=model)) + geom_point() + ylim(0, nrow(d.i$tax_i)) +
  stat_smooth(method="loess", se=F, formula=y~x, size=0.5) +
  theme(panel.grid=element_blank()) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Predicted richness (1km2)")
ggsave("eda/S_fullModel.png", width=6, height=5)

loo.all <- map(agg_vs$full, ~loo::loo(., par="loglik"))
loo::loo_compare(loo.all)
