
library(tidyverse); theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

agg.ls <- map(paste0("out/agg_opt_", c("WY", "Y"), ".rds"), readRDS) %>%
  map(., ~discard(.x, is.null))
agg <- vector("list", n_distinct(c(names(agg.ls[[1]]), names(agg.ls[[2]])))) %>%
  setNames(unique(names(agg.ls[[1]]), names(agg.ls[[2]])))
for(i in names(agg)) {
  agg[[i]] <- rbind(agg.ls[[1]][[i]], agg.ls[[2]][[i]]) %>%
    mutate(., model=case_when(model=="Y" ~ "Structured", 
                              model=="WY" ~ "Joint"))
}
rm(agg.ls)

d.ls <- readRDS("data/opt/Y__opt_var_set_ls.rds")
d.i <- readRDS("data/opt/Y__opt_var_set_i.rds")
site_i <- read_csv("../1_opfo/data/opfo_siteSummaryProcessed.csv")

ants <- load_ant_data(str_type="soil", clean_spp=T)
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
agg$beta







########------------------------------------------------------------------------
## LAMBDAS
########------------------------------------------------------------------------







########------------------------------------------------------------------------
## RICHNESS
########------------------------------------------------------------------------









########------------------------------------------------------------------------
## DIVERSITY
########------------------------------------------------------------------------









########------------------------------------------------------------------------
## BETA DIVERSITY
########------------------------------------------------------------------------
library(betapart)
lam.site.ls <- agg$lam %>% 
  mutate(site=str_pad(str_sub(as.character(id), 1, -5), 2, "left", "0"),
         BDM=arrange(site_i, BDM_id)$BDM[as.numeric(site)]) %>%
  filter(model=="Joint") %>%
  select(sppName, site, BDM, id, L025) %>%
  pivot_wider(names_from="sppName", values_from="L025")
beta.lam.df <- site.mns %>% 
  filter(BDM %in% lam.site.ls$BDM) %>%
  mutate(beta.BRAY.BAL=NA, 
         beta.BRAY.GRA=NA,
         beta.BRAY=NA)
for(i in 1:n_distinct(beta.lam.df$BDM)) {
  BDM_i <- unique(beta.lam.df$BDM)[i]
  beta_i <- beta.multi.abund(filter(lam.site.ls, BDM==BDM_i)[,4:83])
  beta.lam.df$beta.BRAY.BAL[i] <- beta_i$beta.BRAY.BAL
  beta.lam.df$beta.BRAY.GRA[i] <- beta_i$beta.BRAY.GRA
  beta.lam.df$beta.BRAY[i] <- beta_i$beta.BRAY
}

# total
AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model=="Joint")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint"))),
                   paste0("j.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el, data=filter(beta.lam.df, model=="Joint")))

AICcmodavg::aictab(list(lm(beta.BRAY ~ el, 
                           data=filter(beta.lam.df, model!="Joint")),
                        lm(beta.BRAY ~ el + I(el^2), 
                           data=filter(beta.lam.df, model!="Joint"))),
                   paste0("s.tot.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY ~ el + I(el^2), data=filter(beta.lam.df, model!="Joint")))


# balanced variation
AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model=="Joint")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint"))),
                   paste0("j.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el + I(el^2), data=filter(beta.lam.df, model=="Joint")))

AICcmodavg::aictab(list(lm(beta.BRAY.BAL ~ el, 
                           data=filter(beta.lam.df, model!="Joint")),
                        lm(beta.BRAY.BAL ~ el + I(el^2), 
                           data=filter(beta.lam.df, model!="Joint"))),
                   paste0("s.bal.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.BAL ~ el + I(el^2), data=filter(beta.lam.df, model!="Joint")))


# abundance gradient
AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model=="Joint")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model=="Joint"))),
                   paste0("j.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el, data=filter(beta.lam.df, model=="Joint")))

AICcmodavg::aictab(list(lm(beta.BRAY.GRA ~ el, 
                           data=filter(beta.lam.df, model!="Joint")),
                        lm(beta.BRAY.GRA ~ el + I(el^2), 
                           data=filter(beta.lam.df, model!="Joint"))),
                   paste0("s.gra.", c("linear", "quadr")), sort=T)
summary(lm(beta.BRAY.GRA ~ el + I(el^2), data=filter(beta.lam.df, model!="Joint")))







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

