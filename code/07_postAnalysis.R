
library(tidyverse); theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

# model inputs and posteriors
d.ls <- map(paste0("data/opt_noCnpy_rdExc/LV_", c("WY", "Y"), 
                   "__opt_var_set_ls.rds"), readRDS)
d.i <- map(paste0("data/opt_noCnpy_rdExc/LV_", c("WY", "Y"), 
                  "__opt_var_set_i.rds"), readRDS)
# aggregate summarised output for optimal models (WY/Y, cov/LV)
agg_opt.ls <- c(dir("out/", "agg_[cov,LV].*_Y", full.names=T),
                dir("out/opt_noCnpy_rdExc", 
                    "agg_[cov,LV].*_WY", full.names=T)) %>% 
  map(readRDS) %>% map(., ~discard(.x, is.null))
agg_opt <- vector("list", n_distinct(unlist(map(agg_opt.ls, names)))) %>%
  setNames(sort(unique(unlist(map(agg_opt.ls, names)))))
for(i in names(agg_opt)) {
  agg_opt[[i]] <- map(agg_opt.ls, ~.[[i]]) %>% 
    do.call('rbind', .) %>%
    mutate(dataset=if_else(grepl("WY", model), "Joint", "Structured"),
           LV=if_else(grepl("cov", model), "", " + LV"),
           model=paste0(dataset, LV))
}
rm(agg_opt.ls)

# aggregate summarised output for null models (WY/Y, cov/LV)
agg_null.ls <- dir("out/", "agg_null_", full.names=T) %>% 
  map(readRDS) %>% map(., ~discard(.x, is.null))
agg_null <- vector("list", n_distinct(unlist(map(agg_null.ls, names)))) %>%
  setNames(sort(unique(unlist(map(agg_null.ls, names)))))
for(i in names(agg_null)) {
  agg_null[[i]] <- map(agg_null.ls, ~.[[i]]) %>% 
    do.call('rbind', .) %>%
    mutate(dataset=if_else(grepl("WY", model), "Joint", "Structured"),
           LV=if_else(grepl("cov", model), "", " + LV"),
           model=paste0(dataset, LV))
}
rm(agg_null.ls)



site_i <- read_csv("../1_opfo/data/opfo_siteSummaryProcessed.csv")

ants <- load_ant_data(str_type="soil", clean_spp=T, 
                      DNA_dir="../1_opfo/data/DNA_ID_clean")
tax_i <- read_csv("data/tax_i.csv") %>% 
  mutate(across(contains("Full"), as.factor))
  

det_Y <- which(colSums(d.ls$Y)>0)

site.mns <- ants$str %>% filter(TypeOfSample=="soil") %>%
  st_set_geometry(NULL) %>% group_by(BDM, Plot_id) %>%
  summarise(nTubes=n(), S=n_distinct(SPECIESID), mnt25=mean(mnt25)) %>%
  group_by(BDM) %>% 
  summarise(mnTubes=sum(nTubes)/25, seTubes=sd(nTubes)/5,
            mnS=sum(S)/25, seS=sd(S)/5,
            el=mean(mnt25))

det.base <- ants$all %>% sf::st_set_geometry(NULL) %>%
  rename(sppName=SPECIESID) %>%
  group_by(sppName, source) %>% summarise(nDet=n()) %>%
  group_by(sppName) %>% mutate(pDet=nDet/sum(nDet))




########------------------------------------------------------------------------
## DATA SUMMARIES
########------------------------------------------------------------------------







########------------------------------------------------------------------------
## SLOPES
########------------------------------------------------------------------------
agg_opt$beta %>% select(ParName, model, L025, mean, median, L975)







########------------------------------------------------------------------------
## LAMBDAS
########------------------------------------------------------------------------

spp_obs <- which(colSums(d.i$Y)>0)

agg_opt$lam %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(resid=mean - obs) %>%
  filter(obs > 0) %>%
  group_by(dataset, LV, el) %>%
  summarise(RMSE=sqrt(mean(resid^2))) %>%
  ggplot(aes(el, RMSE)) + stat_smooth() + geom_point(shape=1, alpha=0.5) + 
  facet_grid(dataset~LV)

ggplot(agg_opt$LL_I, aes(exp(median), fill=model)) + geom_boxplot()
ggplot(agg_opt$LL_I, aes(el, exp(median))) + 
  stat_smooth() + geom_point(shape=1, alpha=0.5) + 
  facet_grid(dataset~LV)

ggplot(agg_opt$LL_S, aes(exp(median), fill=model)) + geom_boxplot()
ggplot(agg_opt$LL_S, aes(exp(median), xmin=exp(L025), xmax=exp(L975),
                         y=sppName, colour=dataset)) + 
  geom_point(position=position_dodge(width=1)) + 
  geom_linerange(position=position_dodge(width=1)) + 
  facet_grid(.~LV)

agg_opt$lam %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(resid=mean - obs) %>%
  group_by(model) %>%
  summarise(RMSE=sqrt(mean(resid^2)),
            SS_tot=sum((obs - mean(obs))^2),
            SS_resid=sum(resid^2),
            R2=1 - SS_resid/SS_tot)

agg_opt$pP_L %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(obs=as.numeric(obs>0),
         resid=mean - obs) %>%
  group_by(model) %>%
  summarise(RMSE=sqrt(mean(resid^2)),
            SS_tot=sum((obs - mean(obs))^2),
            SS_resid=sum(resid^2),
            R2=1 - SS_resid/SS_tot)

agg_opt$pP_L %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(obs=as.numeric(obs>0), 
         resid=mean - obs) %>%
  filter(spp %in% spp_obs) %>%
  group_by(model, sppName) %>%
  summarise(RMSE=sqrt(mean(resid^2)),
            SS_tot=sum((obs - mean(obs))^2),
            SS_resid=sum(resid^2),
            R2=1 - SS_resid/SS_tot) %>% arrange(desc(R2)) %>%
  ggplot(aes(model, R2)) + geom_boxplot()

agg_opt$lam %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(resid=mean - obs) %>%
  filter(spp %in% spp_obs) %>%
  group_by(model, sppName) %>%
  summarise(RMSE=sqrt(mean(resid^2)),
            SS_tot=sum((obs - mean(obs))^2),
            SS_resid=sum(resid^2),
            R2=1 - SS_resid/SS_tot) %>% 
  group_by(sppName) %>% arrange(model) %>% 
  summarise(pctDiff=(last(RMSE)-first(RMSE))/last(RMSE)) %>% ungroup %>% summary()

agg_opt$lam %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  filter(spp %in% spp_obs) %>%
  # filter(obs > 0) %>%
  ggplot(aes(mean, obs, colour=model)) + geom_point(shape=1, alpha=0.5) +
  stat_smooth(method="lm") 
  # facet_wrap(~sppName) 

agg_opt$pP_L %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  filter(spp %in% spp_obs) %>%
  ggplot(aes(mean, as.numeric(obs>0), colour=model)) + 
  geom_point(shape=1, alpha=0.5) +
  stat_smooth(method="glm", size=0.5, se=F,
              method.args=list(family="binomial"), fullrange=T) + 
  facet_wrap(~sppName) +
  xlim(0, 1) + ylim(0, 1)

agg_opt$pP_L %>% mutate(obs=rep(c(t(d.i$Y)), times=4)) %>%
  mutate(resid=median - (obs>0)) %>%
  ggplot(aes(x=obs>0, y=L975, fill=model)) + geom_boxplot() + ylim(0,1)





########------------------------------------------------------------------------
## PSEUDO-R2
########------------------------------------------------------------------------

mod_col <- c("Joint"="#7b3294", "Structured"="#008837")

# REDO WITH RE-FIT FINAL MODELS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
obs_lam <- agg_opt$lam %>% 
  mutate(obs=rep(c(t(d.ls$Y)), times=4),
         null=rep(filter(agg_null$lam, model=="Structured: None")$mean, times=4),
         LL=LaplacesDemon::dgpois(obs, mean, agg_opt$disp$mean, log=T),
         LL_null=LaplacesDemon::dgpois(obs, null, filter(agg_null$disp, model=="Structured: None")$mean, log=T))
obs_lam %>% group_by(model) %>% 
  summarise(LL=sum(LL), LL_null=sum(LL_null)) %>%
  mutate(R2_mcf=1-(LL/LL_null),
         D2=((-2*LL_null) - (-2*LL))/(-2*LL_null))

agg_opt$LL %>% 
  mutate(R2=1-mean/filter(agg_null$LL, model=="Structured: None")$mean) %>%
  select(model, R2)

R2_spp <- obs_lam %>% group_by(model, sppName) %>% 
  summarise(LL=sum(LL), LL_null=sum(LL_null), obs=any(obs>0), N_Y=sum(obs)) %>%
  mutate(Pr_Y=N_Y/sum(N_Y),
         N_W=colSums(d.ls$W), 
         Pr_W=N_W/sum(N_W), 
         R2_mcf=1-(LL/LL_null),
         D2=((-2*LL_null) - (-2*LL))/(-2*LL_null))

ggplot(R2_spp, aes(R2_mcf, fill=model)) + 
  geom_density(alpha=0.75) + 
  labs(x=expression('Pseudo-R'^2~'by species'), y="Density")

ggplot(R2_spp, aes(R2_mcf, obs, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75) + 
  labs(x=expression('Pseudo-R'^2), y="Species observed in Y")

ggplot(R2_spp, aes(Pr_Y - Pr_W, R2_mcf, colour=model)) + geom_point() 

R2_spp %>% filter(R2_mcf < 0)
R2_spp %>% group_by(model) %>% 
  summarise(PrImprove=sum(R2_mcf>0)/n(),
            mn=mean(R2_mcf),
            se=sd(R2_mcf)/sqrt(n()))
# Mean pseudo-R2 among species:
#  - Structured: 0.16 ± 0.01
#  - Joint: 0.24 ± 0.02


ggplot(R2_spp, aes(R2_mcf, sppName, colour=model)) + 
  geom_point() +
  labs(x="", y=expression('Pseudo-R'^2))

R2_spp %>% arrange(sppName, model) %>% group_by(sppName) %>%
  summarise(diff=first(R2_mcf)-last(R2_mcf), 
            obs=first(obs)) %>%
  ggplot(aes(diff, fill=obs)) + 
  geom_density(alpha=0.7) +
  geom_vline(xintercept=0)

(R2_spp %>% arrange(sppName, model) %>% group_by(sppName) %>%
    summarise(diff=first(R2_mcf)-last(R2_mcf), 
              obs=as.factor(first(obs))) %>%
    group_by(obs) %>%
    summarise(t=list(t.test(diff))))$t









agg_opt$gamma_Sigma %>% filter(spp1 != spp2) %>% 
  filter(model=="Joint: Local LV") %>%
  # filter(sppName1 %in% names(det_Y) & sppName2 %in% names(det_Y)) %>%
  mutate(sppName1=factor(sppName1, levels=sort(unique(sppName1), decreasing=F)),
         sppName2=factor(sppName2, levels=sort(unique(sppName2), decreasing=T))) %>%
  ggplot(aes(sppName1, sppName2, fill=mean)) + 
  geom_point(size=3, shape=22, colour="gray") + 
  scale_fill_gradient2("Mean\nresidual\ncorrelation", 
                       low=scales::muted("blue"), 
                       high=scales::muted("red"), limits=c(-1,1)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) + 
  labs(x="", y="", title="Joint model") +
  facet_grid(genName2~genName1, scales="free", space="free", switch="y")

agg$gamma_Sigma %>% filter(spp1 != spp2) %>%
  ggplot(aes(genName1, mean, colour=genName2)) + 
  geom_point() + facet_wrap(~genName2, scales="free")
agg$gamma_Sigma %>% filter(spp1 != spp2) %>%
  ggplot(aes(mean, sppName1, fill=genName1==genName2)) + 
  geom_vline(xintercept=0) +
  ggridges::geom_density_ridges(rel_min_height=0.01, alpha=0.5, size=0.1)
agg$gamma_Sigma %>% filter(genName1 != genName2) %>%
  ggplot(aes(mean, genName1)) + 
  geom_vline(xintercept=0) +
  ggridges::geom_density_ridges()


ggplot(agg$zeta, aes(el, mean, colour=model)) + geom_point(shape=1, size=0.5) +
  scale_colour_manual(values=mod_col)
  





########------------------------------------------------------------------------
## RICHNESS
########------------------------------------------------------------------------

library(iNEXT)
J_obs <- cbind(index=d.ls[[1]]$IJ, d.ls[[1]]$Y) %>% as_tibble %>%
  group_by(index) %>% summarise_all(sum) %>% ungroup %>% select(-index)
chao <- do.call('rbind', apply(J_obs, 1, ChaoRichness)) %>%
  mutate(BDM=d.i[[1]]$X[(d.ls[[1]]$K+1):(d.ls[[1]]$K+d.ls[[1]]$J),1][unique(d.ls[[1]]$IJ)],
         el=d.i[[1]]$X[(d.ls[[1]]$K+1):(d.ls[[1]]$K+d.ls[[1]]$J),2][unique(d.ls[[1]]$IJ)])

ggplot(chao, aes(el, Estimator, ymin=`95% Lower`, ymax=`95% Upper`)) + 
  geom_point() + geom_linerange() + 
  geom_point(aes(y=Observed), shape=1, size=2)

ggplot(chao, aes(Observed, Estimator)) + geom_point()

ggplot(chao, aes(el, `95% Upper`)) + geom_point()

inext.out <- iNEXT(as.data.frame(t(as.matrix(J_obs))))

Y_0 <- which(rowSums(d.ls[[1]]$Y)==0)
chao_Y <- do.call('rbind', apply(d.ls[[1]]$Y[-Y_0,], 1, ChaoRichness)) %>%
  mutate(el=d.i[[1]]$V[-Y_0,2])

ggplot(chao_Y, aes(el, Estimator)) + geom_point()





########------------------------------------------------------------------------
## DIVERSITY
########------------------------------------------------------------------------









########------------------------------------------------------------------------
## BETA DIVERSITY
########------------------------------------------------------------------------
library(betapart)
lam.site.ls <- agg_opt$lam %>% 
  mutate(site=str_pad(str_sub(as.character(id), 1, -5), 2, "left", "0"),
         BDM=arrange(site_i, BDM_id)$BDM[as.numeric(site)]) %>%
  select(model, sppName, site, BDM, id, median) %>%
  pivot_wider(names_from="sppName", values_from="median")
beta.lam.df <- bind_rows(site.mns %>% mutate(model="Joint: None"),
                         site.mns %>% mutate(model="Structured: None"),
                         site.mns %>% mutate(model="Joint: Local LV"),
                         site.mns %>% mutate(model="Structured: Local LV")) %>% 
  filter(BDM %in% lam.site.ls$BDM) %>%
  mutate(beta.BRAY.BAL=NA, 
         beta.BRAY.GRA=NA,
         beta.BRAY=NA)
for(i in 1:nrow(beta.lam.df)) {
  BDM_i <- beta.lam.df$BDM[i]
  mod_i <- beta.lam.df$model[i]
  beta_i <- beta.multi.abund(filter(lam.site.ls, BDM==BDM_i & model==mod_i)[,5:23])
  beta.lam.df$beta.BRAY.BAL[i] <- beta_i$beta.BRAY.BAL
  beta.lam.df$beta.BRAY.GRA[i] <- beta_i$beta.BRAY.GRA
  beta.lam.df$beta.BRAY[i] <- beta_i$beta.BRAY
}

# total
AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model=="Joint: None")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: None"))),
                   paste0("jN.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el, data=filter(beta.lam.df, model=="Joint: None")))

AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model=="Joint: Local LV")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: Local LV"))),
                   paste0("jL.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el, data=filter(beta.lam.df, model!="Joint: Local LV")))

AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model=="Structured: None")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: None"))),
                   paste0("sN.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el + I(el^2), data=filter(beta.lam.df, model=="Structured: None")))

AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model=="Structured: Local LV")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: Local LV"))),
                   paste0("sL.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el, data=filter(beta.lam.df, model!="Structured: Local LV")))


# balanced variation
AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model=="Joint: None")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: None"))),
                   paste0("jN.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el + I(el^2), data=filter(beta.lam.df, model=="Joint: None")))

AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model=="Joint: Local LV")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: Local LV"))),
                   paste0("jL.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el + I(el^2), data=filter(beta.lam.df, model!="Joint: Local LV")))

AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model=="Structured: None")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: None"))),
                   paste0("sN.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el, data=filter(beta.lam.df, model=="Structured: None")))

AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model=="Structured: Local LV")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: Local LV"))),
                   paste0("sL.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el, data=filter(beta.lam.df, model!="Structured: Local LV")))


# abundance gradient
AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model=="Joint: None")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: None"))),
                   paste0("jN.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el + I(el^2), data=filter(beta.lam.df, model=="Joint: None")))

AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model=="Joint: Local LV")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint: Local LV"))),
                   paste0("jL.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el + I(el^2), data=filter(beta.lam.df, model!="Joint: Local LV")))

AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model=="Structured: None")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: None"))),
                   paste0("sN.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el, data=filter(beta.lam.df, model=="Structured: None")))

AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model=="Structured: Local LV")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Structured: Local LV"))),
                   paste0("sL.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el, data=filter(beta.lam.df, model!="Structured: Local LV")))







########------------------------------------------------------------------------
## DPCoA
########------------------------------------------------------------------------

tax_dist <- tax_i %>% select(contains("Full")) %>% as.data.frame
rownames(tax_dist) <- paste(tax_dist$FullGen, tax_dist$FullSpp, sep="_")

lam.plot.mx <- as.matrix(lam.site.ls[,4:83])
rownames(lam.plot.mx) <- lam.site.ls$id

env.plot <- data.frame(el=as.factor((d.i$V[,2] %/% 100) * 100),
                       el_cont=d.i$V[,2],
                       region=site_i$region[match(lam.site.ls$BDM, site_i$BDM)])

LAM.site.mx <- as.matrix(LAM.site.ls[,3:82])
rownames(LAM.site.mx) <- LAM.site.ls$id

env.site <- data.frame(el=as.factor((LAM.site.ls$el %/% 100) * 100),
                       el_cont=LAM.site.ls$el,
                       region=site_i$region[match(LAM.site.ls$id, site_i$BDM)])

library(ade4)
library(adegraphics) 
adegpar(pbackground.col = "grey", 
        pgrid.col = "white", 
        psub.cex = 1.5, 
        ppoints = list(cex = 0.4, alpha = 0.8), 
        pellipse = list(alpha = 0.4, axes = list(draw = FALSE)))

lam.dpcoa <- dpcoa(as.data.frame(lam.plot.mx), 
                   vegan::taxa2dist(tax_dist, varstep=T), scannf=F, nf=4)
lam.bc <- bca(lam.dpcoa, env.plot$el, scannf=F, nf=5)

lam.bc$ratio




p <- list(
  a=s.class(lam.dpcoa$li, env.plot$el, 
            col=terrain.colors(n_distinct(env.plot$el))),
  b=s.class(lam.dpcoa$dls, tax_dist$FullSpp, 
            col=scales::viridis_pal()(n_distinct(tax_dist$FullSpp))),
  c=s.class(lam.dpcoa$dls, as.factor(tax_dist$FullGen), 
            col=scales::viridis_pal(option="E")(n_distinct(tax_dist$FullGen))),
  d=s.class(lam.dpcoa$dls, tax_dist$FullSF, 
            col=scales::viridis_pal(option="E")(n_distinct(tax_dist$FullSF))),
  e=s.class(lam.bc$ls, env.plot$el, 
            col=terrain.colors(n_distinct(env.plot$el))),
  f=s.class(lam.bc$dls, tax_dist$FullSpp, 
            col=scales::viridis_pal()(n_distinct(tax_dist$FullSpp))),
  g=s.class(lam.bc$dls, tax_dist$FullGen,
            col=scales::viridis_pal(option="E")(n_distinct(tax_dist$FullGen))),
  h=s.class(lam.bc$dls, tax_dist$FullSF, 
            col=scales::viridis_pal(option="E")(n_distinct(tax_dist$FullSF)))
)

Fig2 <- ADEgS(p, layout = c(2,4))

