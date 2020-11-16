
library(tidyverse); theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

# d.ls <- map(paste0("data/opt/", c("WY", "Y"), "__opt_var_set_ls.rds"), readRDS) 
# d.i <- map(paste0("data/opt/", c("WY", "Y"), "__opt_var_set_i.rds"), readRDS) 
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

ants <- load_ant_data(str_type="soil", clean_spp=T)
tax_i <- read_csv("data/tax_i.csv") %>%
  mutate(subf=as.factor(case_when(genus=="Apha" ~ "Myrmicinae",
                                  genus=="Camp" ~ "Formicinae",
                                  genus=="Colo" ~ "Formicinae",
                                  genus=="Crem" ~ "Myrmicinae",
                                  genus=="Doli" ~ "Dolichoderinae",
                                  genus=="Form" ~ "Formicinae",
                                  genus=="Formx" ~ "Myrmicinae",
                                  genus=="Harp" ~ "Myrmicinae",
                                  genus=="Lasi" ~ "Formicinae",
                                  genus=="Lept" ~ "Myrmicinae",
                                  genus=="Mani" ~ "Myrmicinae",
                                  genus=="Myrm" ~ "Myrmicinae",
                                  genus=="Myrme" ~ "Myrmicinae",
                                  genus=="Plag" ~ "Formicinae",
                                  genus=="Poly" ~ "Formicinae",
                                  genus=="Pone" ~ "Ponerinae",
                                  genus=="Sole" ~ "Myrmicinae",
                                  genus=="Sten" ~ "Myrmicinae",
                                  genus=="Tapi" ~ "Dolichoderinae",
                                  genus=="Temn" ~ "Myrmicinae",
                                  genus=="Tetr" ~ "Myrmicinae")),
         fam=as.factor("Formicidae"))

mod_col <- c("Joint"="#7b3294", "Structured"="#008837")

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

talk_fonts <- theme(panel.grid=element_blank(),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    strip.text=element_text(size=16),
                    title=element_text(size=18))






########------------------------------------------------------------------------
## DATA SUMMARIES
########------------------------------------------------------------------------

spp.det_byEl <- bind_rows(as_tibble(d.ls$Y) %>% 
                            mutate(el=d.i$V[,"el"], 
                                   Dataset="Structured"), 
                          as_tibble(d.ls$W) %>%
                            mutate(el=d.i$X[1:d.ls$K, "el"],
                                   Dataset="Cit. Sci.")) %>%
  pivot_longer(1:80, names_to="Species", values_to="nDet") %>%
  uncount(nDet) %>%
  mutate(Genus=tax_i$genus[match(Species, tax_i$species)],
         Subf=tax_i$subf[match(Species, tax_i$species)])

spp.det_byEl %>%
  ggplot(aes(x=el, fill=Dataset, colour=Dataset)) +
  geom_hline(yintercept=0, size=0.25, colour="gray") +
  geom_density(alpha=0.5, size=0.25) + 
  geom_text(aes(y=0), size=1.5, label="-", alpha=0.4, hjust=1.2, angle=90,
            position=position_dodge(width=0.15)) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) +
  facet_wrap(~Species, scales="free_y", ncol=10) + 
  labs(x="Elevation (m)", y="", 
       title="Detections by dataset by species") 
ggsave("eda/obs_elRange_Spp.pdf", height=12, width=17, units="in")

spp.det_byEl %>%
  ggplot(aes(x=el, y=Species, fill=Dataset, colour=Dataset)) +
  ggridges::geom_density_ridges(alpha=0.5, size=0.25, 
                                rel_min_height=0.01, scale=0.9) + 
  geom_text(size=1.5, label="-", alpha=0.4, hjust=1.6, angle=90,
            position=position_dodge(width=0.15)) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) +
  facet_wrap(~Genus, scales="free_y", ncol=7) + 
  labs(x="Elevation (m)", y="", 
       title="Detections by dataset by species") +
  theme(panel.grid.major.y=element_line(size=0.25, colour="gray"))
ggsave("eda/obs_elRange_Spp_byGen.pdf", height=12, width=17, units="in")

spp.det_byEl %>%
  ggplot(aes(x=el, y=Genus, fill=Dataset, colour=Dataset)) +
  ggridges::geom_density_ridges(alpha=0.5, size=0.25, 
                                rel_min_height=0.01, scale=0.8) + 
  geom_text(size=1.5, label="-", alpha=0.4, hjust=1.6, angle=90,
            position=position_dodge(width=0.15)) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) +
  facet_grid(Subf~., scales="free_y", space="free_y") + 
  labs(x="Elevation (m)", y="", 
       title="Detections by dataset by genus") +
  theme(panel.grid.major.y=element_line(size=0.25, colour="gray"),
        panel.spacing=unit(0.1, "cm"))
ggsave("eda/obs_elRange_Gen_tall.pdf", height=10, width=5, units="in")

spp.det_byEl %>%
  ggplot(aes(x=el, y=Genus, fill=Dataset, colour=Dataset)) +
  ggridges::geom_density_ridges(alpha=0.5, size=0.25, 
                                rel_min_height=0.01, scale=0.8) + 
  geom_text(size=1.5, label="-", alpha=0.4, hjust=1.6, angle=90,
            position=position_dodge(width=0.15)) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) +
  facet_wrap(~Subf, scales="free_y") + 
  labs(x="Elevation (m)", y="", 
       title="Detections by dataset by genus") +
  theme(panel.grid.major.y=element_line(size=0.25, colour="gray"))
ggsave("eda/obs_elRange_Gen_square.pdf", height=6, width=8, units="in")

spp.det_byEl %>%
  ggplot(aes(x=el, y=Subf, fill=Dataset, colour=Dataset)) +
  ggridges::geom_density_ridges(alpha=0.5, size=0.25, 
                                rel_min_height=0.01, scale=1) + 
  geom_text(size=2.5, label="-", alpha=0.3, hjust=1.75, angle=90,
             position=position_dodge(width=0.08)) + 
  scale_colour_brewer(type="qual", palette=2) +
  scale_fill_brewer(type="qual", palette=2) +
  labs(x="Elevation (m)", y="", 
       title="Detections by dataset by subfamily") +
  theme(panel.grid.major.y=element_line(size=0.25, colour="gray"),
        panel.spacing=unit(0.1, "cm"))
ggsave("eda/obs_elRange_SF.png", height=5, width=5, units="in", dpi=300)










########------------------------------------------------------------------------
## SLOPES
########------------------------------------------------------------------------

## Aggregate ---------------------------

agg$beta %>% filter(Parameter != "beta[1]") %>%
  mutate(cov=as.numeric(str_sub(Parameter, 6L, -2L))) %>%
  mutate(Scale=factor(str_sub(ParName, 1L, 1L), levels=c("R", "L"), 
                      labels=c("Regional", "Local")),
         Par=str_sub(ParName, 3L, -1L)) %>%
  ggplot(aes(x=Par, y=mn, colour=model)) + 
  geom_hline(yintercept=0, linetype=1, colour="gray", size=0.35) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_linerange(aes(ymin=q25, ymax=q75), 
                 position=position_dodge(width=0.5), size=1.5) + 
  geom_linerange(aes(ymin=q05, ymax=q95),
                 position=position_dodge(width=0.5), size=0.5) + 
  facet_grid(Scale~., scales="free_y", space="free_y") +
  coord_flip() + 
  labs(x="", y=expression(beta), title="Aggregate slopes") + 
  scale_colour_manual("", values=mod_col) + talk_fonts 
ggsave("eda/beta_opt.png", width=9, height=8)



## Genus -------------------------------

ggplot(filter(agg$B, ParName != "intercept"), 
       aes(x=ParName, y=q50, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(q05)==sign(q95)), 
             position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=q25, ymax=q75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=q05, ymax=q95),
                 position=position_dodge(width=0.5), size=0.2) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  labs(y="Genus-level response (B)", x="") +
  facet_wrap(~genName, ncol=7) 
ggsave("eda/Bg_opt_byGenus.png", width=15, height=6)

ggplot(filter(agg$B, ParName != "intercept"), 
       aes(x=genName, y=q50, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(q05)==sign(q95)), 
             position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=q25, ymax=q75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=q05, ymax=q95),
                 position=position_dodge(width=0.5), size=0.2) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  labs(y="Genus-level response (B)", x="") +
  facet_wrap(~ParName, scales="free_x")
ggsave("eda/Bg_opt_byParam.png", width=9, height=8)

agg$B %>% filter(ParName != "intercept") %>%
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
ggsave("eda/Bg_opt_bar.png", width=7, height=8)

filter(agg$B, ParName != "intercept") %>%
  filter(sign(q05)==sign(q95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/21, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  scale_fill_manual(values=mod_col) +
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of genera", x="")




## Species -----------------------------

ggplot(agg$b, aes(x=abs(q95-q05), y=ParName, fill=spp %in% det_Y)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) + 
  scale_fill_brewer("Species detected in\nstructured samples", type="qual") +
  labs(x="90% CI width", y="", 
       title="Distribution of species-level response CI widths") +
  facet_wrap(~model)
ggsave("eda/b_opt_obsY_CI_ridges.png", width=7, height=6)

ggplot(agg$b, aes(x=mn, y=ParName, fill=spp %in% det_Y)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_fill_brewer("Species detected in\nstructured samples", type="qual") +
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~model)
ggsave("eda/b_opt_obsY_mn_ridges.png", width=7, height=6)

agg$b %>% mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=mn, y=genName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, 
                                rel_min_height=0.01, size=0.25) +
  geom_point(aes(colour=model), size=1.5, shape="|") +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("", values=mod_col) + 
  scale_fill_manual("", values=mod_col) + 
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~ParName, scales="free_x") +
  theme(panel.grid.major.y=element_line(colour="gray30", size=0.2))
ggsave("eda/b_opt_aggGenus_byParam_ridges.png", width=7, height=8)

agg$b %>% mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=mn, y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, 
                                rel_min_height=0.01, size=0.25) +
  geom_point(aes(colour=model), size=1.5, shape="|") +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("", values=mod_col) + 
  scale_fill_manual("", values=mod_col) + 
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~genName, scales="free_x", ncol=7) + 
  theme(panel.grid.major.y=element_line(colour="gray30", size=0.2))
ggsave("eda/b_opt_aggGenus_byGenus_ridges.png", width=13, height=6)

agg$b %>% mutate(sfName=tax_i$subf[match(sppName, tax_i$species)]) %>%
  ggplot(aes(x=mn, y=sfName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, 
                                rel_min_height=0.01, size=0.25) +
  geom_point(aes(colour=model), size=1.5, shape="|") +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("", values=mod_col) + 
  scale_fill_manual("", values=mod_col) + 
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~ParName, scales="free_x") + 
  theme(panel.grid.major.y=element_line(colour="gray30", size=0.2))
ggsave("eda/b_opt_aggSF_byParam_ridges.png", width=7, height=8)

agg$b %>% mutate(sfName=tax_i$subf[match(sppName, tax_i$species)]) %>%
  ggplot(aes(x=mn, y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, 
                                rel_min_height=0.01, size=0.25) +
  geom_point(aes(colour=model), size=1.5, shape="|") +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("", values=mod_col) + 
  scale_fill_manual("", values=mod_col) + 
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~sfName, scales="free_x") + 
  theme(panel.grid.major.y=element_line(colour="gray30", size=0.2))
ggsave("eda/b_opt_aggSF_bySF_ridges.png", width=6, height=5)

agg$b %>% mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=ParName, y=q50, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(q05)==sign(q95)),
             position=position_dodge(width=0.5)) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  geom_linerange(aes(ymin=q25, ymax=q75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=q05, ymax=q95),
                 position=position_dodge(width=0.5), size=0.2) + 
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  facet_wrap(~sppName, ncol=10) + 
  labs(y="Species-level responses (b)", x="")
ggsave("eda/b_opt_bySpecies.png", width=14, height=10)

agg$b %>% mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=sppName, y=q50, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(q05)==sign(q95)),
             position=position_dodge(width=0.5)) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  geom_linerange(aes(ymin=q25, ymax=q75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=q05, ymax=q95),
                 position=position_dodge(width=0.5), size=0.2) + 
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  facet_grid(genName~ParName, scales="free", space="free_y") +
  labs(y="Species-level responses (b)", x="") + 
  theme(strip.text.y=element_blank(), strip.background.y=element_blank(),
        panel.spacing.y=unit(0.05, 'cm'))
ggsave("eda/b_opt_byParam.png", width=9, height=12)

agg$b %>% filter(cov != "1") %>% 
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Par=str_sub(ParName, 3L, -1L)) %>%
  mutate(effect=case_when(q05<0 & q95<0 ~ "negative",
                          q05>0 & q95>0 ~ "positive",
                          q05<0 & q95>0 ~ "ns",
                          q05>0 & q95<0 ~ "ns"),
         sig=sign(q05)==sign(q95)) %>%
  group_by(model, Par, Scale, ParName, effect) %>% 
  summarise(nSpp=n(), nSig=sum(sig)) %>%
  group_by(model, Par, Scale, ParName) %>%
  mutate(effect=factor(effect, levels=c("ns", "negative", "positive")),
         nSig=sum(nSig)) %>%
  arrange(model, Scale, nSig) %>% ungroup %>%
  mutate(Par=factor(Par, levels=unique(Par)),
         Scale=factor(Scale, levels=c("R", "L"), 
                      labels=c("Regional", "Local"))) %>%
  ggplot(aes(Par, nSpp/80, fill=effect)) + 
  geom_bar(stat="identity", colour="gray30", size=0.25) +
  facet_grid(Scale~model, scales="free", space="free_y") + coord_flip() + 
  scale_fill_manual("", values=c("gray90", "blue", "red")) + 
  labs(x="", y="Proportion of species", title="Species Responses") + talk_fonts
ggsave("eda/b_opt_bar.png", width=10, height=8)

agg$b %>% filter(cov != "1") %>%
  filter(sign(q05)==sign(q95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/80, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  scale_fill_manual(values=mod_col) +
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of species", x="")

agg$b %>% filter(cov != "1") %>%
  mutate(effect=case_when(q05<0 & q95<0 ~ "negative",
                          q05>0 & q95>0 ~ "positive",
                          q05<0 & q95>0 ~ "no effect",
                          q05>0 & q95<0 ~ "no effect"),
         sig=sign(q05)==sign(q95)) %>%
  group_by(model, sppName, sig) %>% summarise(nVar=n()) %>% 
  # filter(sig) %>% summary # 1-12 significant variables
  ggplot(aes(sig, nVar, fill=sppName %in% names(det_Y))) + 
  geom_boxplot() + ylim(0, NA) + 
  scale_fill_brewer("Species detected in\nstructured samples", type="qual")

agg$b %>% filter(cov != "1") %>%
  select(model, sppName, ParName, mn) %>%
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  pivot_wider(names_from="model", values_from="mn") %>%
  filter(!is.na(Joint) & !is.na(Structured)) %>%
  ggplot(aes(Structured, Joint, colour=ParName)) + 
  geom_hline(yintercept=0, colour="gray90") + 
  geom_vline(xintercept=0, colour="gray90") + 
  geom_abline(colour="gray30", size=0.25) + geom_point(shape=1) + 
  scale_colour_brewer(type="qual", palette=2) + 
  facet_wrap(~genName, ncol=7) + ylim(-4,4) + xlim(-4,4)
ggsave("eda/b_opt_Y_WY_scatter.png", width=13, height=6)









########------------------------------------------------------------------------
## LAMBDAS
########------------------------------------------------------------------------

## Regional ----------------------------

## total ----------
lLAM.loess <- agg$lLAM %>%
  group_by(model, id, el) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100),
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)")
ggsave("eda/R_lLAMBDA_opt.png", width=6, height=5)

ggplot(lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100),
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare")
ggsave("eda/R_LAMBDA_opt.png", width=6, height=5)


## by species -----
spp.lLAM.loess <- agg$lLAM %>%
  group_by(model, id, el, sppName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

det.spp <- det.base %>% select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", str_pad(round(p*100), 3, "left", " "), "%", 
                         "\ns: ", str_pad(round(s*100), 3, "left", " "), "%"),
         model="Joint")

ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/R_lLAMBDA_opt_el_Species.png", width=14, height=10)

ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species.png", width=14, height=10)

ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/R_lLAMBDA_opt_el_Species_free.png", width=16, height=10)

ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species_free.png", width=16, height=10)


## by genus -----
det.gen <- det.base %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(genName, source) %>% summarise(nDet=sum(nDet)) %>%
  group_by(genName) %>% mutate(pDet=nDet/sum(nDet)) %>%
  select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", round(p*100), "%", 
                         "\ns: ", round(s*100), "%"),
         model="Joint")

gen.lLAM.loess <- agg$lLAM %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, id, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=8.5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/R_lLAMBDA_opt_el_Genus.png", width=15, height=6)

ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus.png", width=15, height=6)

ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=8.5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)") + 
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/R_lLAMBDA_opt_el_Genus_free.png", width=15, height=6)

ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus_free.png", width=15, height=6)




## Local -------------------------------

llam.loess <- agg$llam %>%
  group_by(model, plot, el) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

ggplot(llam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi), fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md), colour=model)) + 
  geom_point(data=site.mns, aes(y=mnTubes)) +
  geom_linerange(data=site.mns,
                 aes(ymin=mnTubes-seTubes, ymax=mnTubes+seTubes), size=0.25) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Mean abundance per plot")
ggsave("eda/L_llambda_opt.png", width=9, height=6)


## by species -----
spp.llam.loess <- agg$llam %>%
  group_by(model, plot, el, sppName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

det.spp <- det.base %>% select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", str_pad(round(p*100), 3, "left", " "), "%", 
                         "\ns: ", str_pad(round(s*100), 3, "left", " "), "%"),
         model="Joint")

ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_llambda_opt_el_Species.png", width=14, height=10)

ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_llambda_opt_el_Species_pts.png", width=14, height=10)

ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_lambda_opt_el_Species.png", width=14, height=10)

ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_lambda_opt_el_Species_pts.png", width=14, height=10)

ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_llambda_opt_el_Species_free.png", width=16, height=10)

ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_llambda_opt_el_Species_free_pts.png", width=16, height=10)

ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_lambda_opt_el_Species_free.png", width=16, height=10)

ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_lambda_opt_el_Species_free_pts.png", width=16, height=10)


## by genus -----
det.gen <- det.base %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(genName, source) %>% summarise(nDet=sum(nDet)) %>%
  group_by(genName) %>% mutate(pDet=nDet/sum(nDet)) %>%
  select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", round(p*100), "%", 
                         "\ns: ", round(s*100), "%"),
         model="Joint")

gen.llam.loess <- agg$llam %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q10 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(q50 ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q90 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=8.5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_llambda_opt_el_Genus.png", width=15, height=6)

agg$llam %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  ggplot(aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_llambda_opt_el_Genus_pts.png", width=15, height=6)

ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_lambda_opt_el_Genus.png", width=15, height=6)

agg$llam %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  ggplot(aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_pts.png", width=15, height=6)

ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=1, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=8.5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") + 
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_llambda_opt_el_Genus_free.png", width=15, height=6)

agg$llam %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  ggplot(aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot  (log-axis)") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_llambda_opt_el_Genus_free_pts.png", width=15, height=6)

ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_free.png", width=15, height=6)

agg$llam %>% 
  mutate(genName=tax_i$genus[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("q10", "q50", "q90"), ~log(sum(exp(.))))) %>%
  ggplot(aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(q50), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_free_pts.png", width=15, height=6)










########------------------------------------------------------------------------
## RICHNESS
########------------------------------------------------------------------------

## Regional ----------------------------
S.loess <- agg$pP_R %>%
  group_by(model, id, el) %>%
  summarise(across(one_of("q05", "mn", "q95"), ~sum(.>0.95))) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(q05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(mn ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(q95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(S.loess, aes(x=el, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), alpha=0.25, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  ylim(0, NA) +
  labs(x="Elevation (m)", y="Richness")
ggsave("eda/S_opt.png", width=6, height=5)


agg$lLAM %>% #filter(id>1e5) %>% 
  group_by(site, model, el) %>% 
  summarise(S_lo=sum(1-exp(-exp(q05))>0.995),
            S_md=sum(1-exp(-exp(q50))>0.995),
            S_hi=sum(1-exp(-exp(q95))>0.995)) %>%
  ggplot(aes(el, group=model, colour=model)) + 
  stat_smooth(aes(y=S_lo), se=F, size=0.5, span=2) + 
  stat_smooth(aes(y=S_md), se=F, span=2) + 
  stat_smooth(aes(y=S_hi), se=F, size=0.5, span=2) + 
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2")

agg$pP_R %>% #filter(id>1e5) %>% 
  group_by(id, model, el) %>% 
  summarise(across(one_of("q025", "q05", "q10", "q25", "q50", "mn",
                          "q75", "q90", "q95", "q975"), ~sum(.>0.95))) %>%
  ggplot(aes(el, group=model, colour=model)) + 
  stat_smooth(aes(y=q025), size=0.1, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q05), size=0.2, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q10), size=0.5, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q25), size=0.75, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q50), size=1, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=mn), size=1, se=F, linetype=2, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q75), size=0.75, se=F, method="loess", formula=y~x) +  
  stat_smooth(aes(y=q90), size=0.5, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q95), size=0.2, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=q975), size=0.1, se=F, method="loess", formula=y~x) + 
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2") + facet_wrap(~model)

agg$S_R %>%
  ggplot(aes(el, mn, colour=model)) + 
  geom_point() +
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2")

agg$pred_R %>% group_by(model, id, el) %>% 
  summarise(S=sum(mn>0.95)) %>% 
  ggplot(aes(el, S, colour=model)) + 
  geom_point(alpha=0.5) + 
  scale_colour_manual(values=mod_col) + ylim(0, 80)


## Local -------------------------------

agg$S_L %>% 
  ggplot(aes(el, colour=model)) + ylim(0, NA) +
  geom_point(data=site.mns, aes(y=mnS), colour="black") +
  geom_linerange(data=site.mns,
                 aes(ymin=mnS-seS, ymax=mnS+seS), size=0.25, colour="black") +
  stat_smooth(aes(y=mn), size=0.5, se=F) + 
  stat_smooth(aes(y=q05), size=0.25, se=F) + 
  stat_smooth(aes(y=q95), size=0.25, se=F) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~model) + labs("Elevation (m)", "Predicted richness per 0.75m2")








########------------------------------------------------------------------------
## DIVERSITY
########------------------------------------------------------------------------

agg$H %>%
  ggplot(aes(el, mn, colour=model)) + 
  geom_point() +
  scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted diversity per km2")







########------------------------------------------------------------------------
## PROBABILITY OF PRESENCE
########------------------------------------------------------------------------

## Regional ----------------------------

agg$lLAM %>% mutate(pP_R=1-exp(-exp(q50))) %>%
  filter(spp %in% det_Y) %>%
  ggplot(aes(x=el, y=pP_R, colour=model)) + ylim(0,1) +
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="pr(Presence): site") +
  stat_smooth(method="loess", se=F) + facet_wrap(~sppName)

agg$pP_R %>% 
  # filter(spp %in% det_Y) %>%
  ggplot(aes(x=el, y=mn, colour=model)) + ylim(0,1) +
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="pr(Presence): site") +
  facet_wrap(~sppName)

agg$pP_R %>% filter(id>1e5) %>%
  ggplot(aes(x=el, colour=model)) + 
  # geom_point(alpha=0.5, size=0.5, shape=1) + 
  stat_smooth(aes(y=q95), method="loess", se=F, size=0.1) + 
  stat_smooth(aes(y=q75), method="loess", se=F, size=0.25) + 
  stat_smooth(aes(y=q50), method="loess", se=F, size=0.5) + 
  stat_smooth(aes(y=q25), method="loess", se=F, size=0.25) + 
  stat_smooth(aes(y=q05), method="loess", se=F, size=0.1) + 
  ylim(0, 1) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~sppName)



## Local -------------------------------

agg$pP_L %>% 
  ggplot(aes(x=el, y=mn, colour=model)) + 
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~sppName, scales="free_y") + 
  labs(x="Elevation (m)", y="pr(Presence): plot")














########------------------------------------------------------------------------
## TAXONOMIC BIAS
########------------------------------------------------------------------------

agg$D %>%
  ggplot(aes(sppName, q50, ymin=q05, ymax=q95, 
             colour=sign(q05-1)==sign(q95-1))) + 
  geom_hline(yintercept=1, colour="gray") +
  geom_point() + geom_linerange() + labs(x="", y="Proportional Bias") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        panel.grid=element_blank()) + 
  scale_colour_manual(values=c("gray", "darkred"), guide=F)
ggsave("eda/D_opt.png", width=9, height=4)










########------------------------------------------------------------------------
## SIGMAS
########------------------------------------------------------------------------

## Among genera ------------------------

agg$Sig_B %>% filter(gen1 != gen2) %>%
  mutate(sf1=tax_i$subf[match(gen1Name, tax_i$genus)],
         sf2=tax_i$subf[match(gen2Name, tax_i$genus)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(x=gen1Name, y=gen2Name, colour=mn)) + 
  geom_point(size=3) + facet_grid(model~ParName) + 
  scale_colour_gradient2(midpoint=0, low=scales::muted("blue"), high=scales::muted("red")) + 
  theme(panel.grid=element_line(colour="gray", size=0.25),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

agg$Sig_B %>% filter(gen1 != gen2) %>%
  mutate(sf1=tax_i$subf[match(gen1Name, tax_i$genus)],
         sf2=tax_i$subf[match(gen2Name, tax_i$genus)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(ParName, mn, fill=sf2)) + 
  geom_hline(yintercept=0) + 
  geom_boxplot() + facet_grid(model~sf1)

agg$Sig_B %>% filter(gen1 != gen2) %>%
  mutate(sf1=tax_i$subf[match(gen1Name, tax_i$genus)],
         sf2=tax_i$subf[match(gen2Name, tax_i$genus)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(mn, ParName, fill=sf2)) + 
  geom_vline(xintercept=0) +
  ggridges::geom_density_ridges(alpha=0.5) +
  facet_grid(model~sf1)

agg$Sig_B %>% filter(gen1 != gen2) %>%
  ggplot(aes(x=mn, y=ParName, fill=model)) + 
  geom_vline(xintercept=0) + 
  ggridges::geom_density_ridges(alpha=0.5) +
  scale_fill_manual("", values=mod_col)


## Within genera -----------------------

agg$sig_b %>%
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Par=str_sub(ParName, 3L, -1L)) %>%
  mutate(Scale=factor(Scale, levels=c("R", "L"), 
                      labels=c("Regional", "Local"))) %>%
  ggplot(aes(ParName, q50, ymin=q05, ymax=q95, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(position=position_dodge(width=0.5)) + 
  scale_colour_manual("", values=mod_col) +
  coord_flip() + labs(x="Standard deviation among congeners", y="")










agg$pP_L$obs <- rep(c(t(d.ls$Y)), n_distinct(agg$pP_L$model))
agg$pP_L$disp <- agg$disp$mn[match(agg$pP_L$model, agg$disp$model)]


ggplot(agg$pP_L, aes(mn, obs)) + geom_point(alpha=0.1) + 
  stat_smooth() + facet_wrap(~model)
ggplot(agg$pP_L, aes(mn, as.numeric(obs>0), group=model, colour=model)) + 
  #geom_point(alpha=0.1) + 
  stat_smooth(aes(x=q05), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q50), method="glm", size=0.5, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=q95), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  # facet_wrap(~model) +
  ylim(0,1)
ggplot(agg$pP_L, aes(factor(obs), mn, fill=model)) + geom_boxplot() +
  labs(x="Number of detections at plot", y=expression(lambda))
ggplot(agg$pP_L, aes(factor(obs>0), mn, fill=model)) + geom_boxplot() +
  labs(x="Number of detections at plot", y=expression(lambda))


agg$pP_L <- agg$pP_L %>%
  mutate(ll=LaplacesDemon::dgpois(obs, mn, disp, log=T))
agg$pP_L %>% group_by(model) %>% 
  summarise(nll=-sum(ll)*2) %>% 
  arrange(nll)




agg$pP_R  %>% group_by(site, model, el) %>% 
  filter(id>1e5) %>%
  summarise(S_mn=sum(mn > 0.95)) %>%
  ggplot(aes(el)) + ylim(0, 80) + geom_point(aes(y=S_mn)) + 
  facet_wrap(~model) + labs("Elevation (m)", "Predicted richness per km2")

agg$S_R %>% filter(id>1e5) %>%
  ggplot(aes(el, q025)) + geom_point() + ylim(0,80)


agg_full <- readRDS("out/agg_opt_full_Y.rds")

S_R.ls <- map(agg_full, 
              ~rstan::extract(., pars="Rich")[[1]]) %>% 
  map(~matrix(., nrow=3000, dimnames=list(1:3000, d.i$X[,"id"])) %>%
        as_tibble %>% 
        pivot_longer(cols=everything(), names_to="id", values_to="S") %>%
        mutate(el=d.i$X[match(id, d.i$X[,"id"]),"el"],
               iter=rep(1:3000, each=nrow(d.i$X)))) 

S_L.ls <- map(agg_full, 
              ~rstan::extract(., pars="RichL")[[1]]) %>% 
  map(~matrix(., nrow=3000, dimnames=list(1:3000, d.i$V[,"Plot_id"])) %>%
        as_tibble %>% 
        pivot_longer(cols=everything(), names_to="Plot_id", values_to="S") %>%
        mutate(el=d.i$V[match(Plot_id, d.i$V[,"Plot_id"]),"el"],
               iter=rep(1:3000, each=nrow(d.i$V)))) 

lLAMBDA.ls <- map(agg_full, 
                  ~rstan::extract(., pars="lLAMBDA")[[1]]) %>%
  map(~array(., dim=dim(.), 
             dimnames=list("iter"=1:3000, 
                           "id"=d.i$X[,"id"], 
                           "spp"=tax_i$species)))
lLAMBDA.ls <- map(lLAMBDA.ls, ~.[,which(d.i$X[,"id"] > 1e5),])
lLAMBDA.tb <- map(lLAMBDA.ls, ~cubelyr::as.tbl_cube(.) %>% as_tibble) 
LAM.df <- map(1:length(lLAMBDA.tb), ~lLAMBDA.tb[[.x]] %>% 
                mutate(model=names(lLAMBDA.tb)[[.x]])) %>%
  do.call("rbind", .)  %>% 
  rename(lLAM=".")



LAM.df %>% group_by(id, spp) %>% 
  summarise(prop_pos=sum((lLAM-log(100))>0)/n()) %>%
  mutate(el=d.i$X[match(as.numeric(id), d.i$X[,"id"]),"el"]) %>% 
  ggplot(aes(el, prop_pos)) + 
  geom_point(alpha=0.5, shape=1, size=0.5) + 
  facet_wrap(~spp)

LAM.df %>% group_by(id, spp) %>% 
  mutate(el=d.i$X[match(as.numeric(id), d.i$X[,"id"]),"el"],
         elBin=(el%/%200)*200) %>% 
  # filter(spp %in% c("Apha_subt", "Lasi_emar")) %>%
  ggplot(aes(x=lLAM, y=as.factor(elBin))) + 
  ggridges::geom_density_ridges(alpha=0.5) + 
  facet_wrap(~spp) + xlim(0,1)





lLAMBDA.ls <- map(agg_full, 
                  ~1-exp(-exp(rstan::extract(., pars="lLAMBDA")[[1]]))) %>%
  map(~array(., dim=dim(.), 
             dimnames=list("iter"=1:3000, #1:210, 
                           "id"=d.i$X[,"id"], 
                           "spp"=d.ls$tax_i$species)))
lLAMBDA.tb <- map(lLAMBDA.ls, ~cubelyr::as.tbl_cube(.) %>% as_tibble) 
LAM.df <- map(1:length(lLAMBDA.tb), ~lLAMBDA.tb[[.x]] %>% 
                mutate(model=names(lLAMBDA.tb)[[.x]])) %>%
  do.call("rbind", .) %>%
  filter(id > 1e5)

LAM.df %>%
  mutate(sim=rbinom(n(), 1, .)) %>%
  # group_by(model, id, iter) %>% summarise(S=sum(sim)) %>%
  group_by(model, id, spp) %>% summarise(mnP=mean(.)) %>%
  mutate(el=d.ls$X[,'el'][match(id, d.ls$X[,'id'])]) %>%
  group_by(model, el) %>%
  summarise(S=sum(mnP>0.995)) %>%
  ggplot(aes(el, S, colour=model)) + geom_point() + ylim(0, nrow(d.ls$tax_i)) +
  stat_smooth(method="loess", se=F, formula=y~x, size=0.5) +
  theme(panel.grid=element_blank()) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Predicted richness (1km2)")
ggsave("eda/S_opt.png", width=6, height=5)

loo.all <- map(agg$full, ~loo::loo(., par="loglik"))
loo::loo_compare(loo.all)



library(sf); library(viridis)
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp") %>%
  filter(!grepl("Lac ", NAME)) 
site.sf <- agg_str_site_data() %>% arrange(BDM)
X_all <- rbind(d.i$X, d.i$X_) %>% as.data.frame
lLAM.site <- agg$lLAM %>% filter(id<1e5) %>% 
  group_by(id, model) %>% 
  summarise(N_lo=log(sum(exp(q05))),
            N_md=log(sum(exp(q50))),
            N_mn=log(sum(exp(mn))),
            N_hi=log(sum(exp(q95)))) %>%
  mutate(model=case_when(model=="Joint" ~ "Joint",
                         model=="Structured" ~ "Str")) %>%
  pivot_wider(names_from="model", values_from=3:6)
out.sf <- d.i$grd_W.sf %>% 
  left_join(., select(X_all, -el), by="id") %>%
  left_join(., lLAM.site, by="id") %>%
  left_join(., agg$pP_R %>% 
              group_by(id, model) %>% 
              summarise(across(one_of("q05", "q50", "mn","q95"), 
                               ~sum(.>0.95))) %>%
              rename(S_lo=q05, S_md=q50, S_mn=mn, S_hi=q95) %>%
              mutate(model=case_when(model=="Joint" ~ "Joint",
                                     model=="Structured" ~ "Str")) %>%
              pivot_wider(names_from="model", values_from=3:6), by="id")
  
ggplot(out.sf, aes(fill=(S_hi_Str-S_lo_Str)/S_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[S[Str]]), option="B", limits=c(0, 1.7))

ggplot(out.sf, aes(fill=(S_hi_Joint-S_lo_Joint)/S_md_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[S[Joint]]), option="B", limits=c(0, 1.7))

ggplot(out.sf, aes(fill=(exp(N_hi_Str)-exp(N_lo_Str))/exp(N_md_Str))) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[N[Str]]), option="B", limits=c(0, 5.3))

ggplot(out.sf, aes(fill=(exp(N_hi_Joint)-exp(N_lo_Joint))/exp(N_md_Joint))) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[N[Str]]), option="B", limits=c(0, 5.3))

ggplot(out.sf, aes(fill=S_mn_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(S[Joint]), option="B", limits=c(0, 80))

ggplot(out.sf, aes(fill=S_mn_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(S[Str]), option="B", limits=c(0, 80))


map_theme <- theme(legend.position=c(0.35, 0.075), 
                   legend.direction="horizontal",
                   legend.key.height=unit(0.3, "cm"), 
                   legend.title=element_text(size=8, hjust=1), 
                   legend.text=element_text(size=5), 
                   axis.text=element_blank(), 
                   axis.ticks=element_blank(), 
                   axis.title=element_blank(), 
                   panel.border=element_rect(size=0.2, colour="grey30", fill=NA)) 

for(m in 1:2) {
  m.full <- unique(agg$lLAM$model)[m]
  m.abb <- ifelse(m.full=="Joint", "WY", "Y")
  
  sp.N <- sp.pP <- sp.PA <- sp.LamCV <- sp.PrCI <- vector("list", n_distinct(tax_i$species))
  gen.N <- gen.S <- vector("list", n_distinct(tax_i$genus))
  sf.N <- sf.S <- vector("list", n_distinct(tax_i$subf))
  for(i in seq_along(tax_i$species)) {
    sp.i <- tax_i$species[i]
    sp.f <- str_replace(sp.i, "/", "_")
    lLAM.i <- filter(agg$lLAM, sppName==sp.i & model==m.full)
    pP.i <- filter(agg$pP_R, sppName==sp.i & model==m.full)
    obs.i <- filter(ants$all, SPECIESID==sp.i) 
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    
    sp.N[[i]] <- d.i$grd_W.sf %>%
      mutate(N_sp=lLAM.i$mn[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_sp)) %>%
      ggplot() + geom_sf(aes(fill=exp(N_sp)), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(Lambda), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_LAM_", m.abb, ".jpg"), 
           sp.N[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.pP[[i]] <- d.i$grd_W.sf %>%
      mutate(P_sp=pP.i$mn[match(id, pP.i$id)]) %>%
      filter(!is.na(P_sp)) %>%
      ggplot() + geom_sf(aes(fill=P_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) +
      scale_fill_viridis("pr(Pres)", option="B", limits=c(0, 1)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_pP_", m.abb, ".jpg"), 
           sp.pP[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.PA[[i]] <- d.i$grd_W.sf %>%
      mutate(P_sp=pP.i$mn[match(id, pP.i$id)]>0.99) %>%
      filter(!is.na(P_sp)) %>%
      ggplot() + geom_sf(aes(fill=P_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_manual("pr(Pres) > 99%", 
                        values=c("FALSE"="gray70", "TRUE"="#2171b5")) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_PA_", m.abb, ".jpg"), 
           sp.PA[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.PrCI[[i]] <- d.i$grd_W.sf %>%
      mutate(P_CI=pP.i$q95[match(id, pP.i$id)]-pP.i$q05[match(id, pP.i$id)]) %>%
      filter(!is.na(P_CI)) %>%
      ggplot() + geom_sf(aes(fill=P_CI), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) +
      scale_fill_viridis(expression(CI[90](PrP)), option="E", limits=c(0,1)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_prCI_", m.abb, ".jpg"), 
           sp.PrCI[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.LamCV[[i]] <- d.i$grd_W.sf %>%
      mutate(CV_sp=log(exp(lLAM.i$sd[match(id, lLAM.i$id)])/
                         exp(lLAM.i$mn[match(id, lLAM.i$id)]))) %>%
      filter(!is.na(CV_sp)) %>%
      ggplot() + geom_sf(aes(fill=CV_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(log(CV[Lambda])), option="E") + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_LamCV_", m.abb, ".jpg"), 
           sp.LamCV[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  for(i in 1:n_distinct(tax_i$genus)) {
    gen.i <- unique(tax_i$genus)[i]
    gen.sppNum <- which(tax_i$genus == gen.i)
    obs.i <- filter(ants$all, SPECIESID %in% tax_i$species[gen.sppNum])
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    lLAM.i <- filter(agg$lLAM, spp %in% gen.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mn=sum(exp(mn)))
    S.i <- filter(agg$pP_R, spp %in% gen.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mn=sum(mn>0.99))
    gen.S[[i]] <- d.i$grd_W.sf %>%
      mutate(S_gen=S.i$mn[match(id, S.i$id)]) %>%
      filter(!is.na(S_gen)) %>%
      ggplot() + geom_sf(aes(fill=S_gen), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", gen.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis("Pred.\nRich.", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/gen/", gen.i, "_S_", m.abb, ".jpg"), 
           gen.S[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    gen.N[[i]] <- d.i$grd_W.sf %>%
      mutate(N_gen=lLAM.i$mn[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_gen)) %>%
      ggplot() + geom_sf(aes(fill=N_gen), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", gen.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(Lambda), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/gen/", gen.i, "_LAM_", m.abb, ".jpg"),
           gen.N[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  for(i in 1:n_distinct(tax_i$subf)) {
    sf.i <- unique(tax_i$subf)[i]
    sf.sppNum <- which(tax_i$subf == sf.i)
    obs.i <- filter(ants$all, SPECIESID %in% tax_i$species[sf.sppNum])
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    lLAM.i <- filter(agg$lLAM, spp %in% sf.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mn=sum(exp(mn)))
    S.i <- filter(agg$pP_R, spp %in% sf.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mn=sum(mn>0.99))
    sf.S[[i]] <- d.i$grd_W.sf %>%
      mutate(S_sf=S.i$mn[match(id, S.i$id)]) %>%
      filter(!is.na(S_sf)) %>%
      ggplot() + geom_sf(aes(fill=S_sf), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sf.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis("Pred.\nRich.", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/sf/", sf.i, "_S_", m.abb, ".jpg"), 
           sf.S[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sf.N[[i]] <- d.i$grd_W.sf %>%
      mutate(N_sf=lLAM.i$mn[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_sf)) %>%
      ggplot() + geom_sf(aes(fill=N_sf), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sf.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(Lambda), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/sf/", sf.i, "_LAM_", m.abb, ".jpg"), 
           sf.N[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  
  ggpubr::ggarrange(plotlist=sp.N, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_N_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.pP, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_pP_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.PA, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_PA_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.PrCI, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_PrCI_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.LamCV, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_LamCV_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  
  ggpubr::ggarrange(plotlist=gen.N, ncol=7, nrow=3) %>%
    ggsave(paste0("eda/maps/gen_N_", m.abb, ".jpg"), ., 
           width=21, height=9, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=gen.S, ncol=7, nrow=3) %>%
    ggsave(paste0("eda/maps/gen_S_", m.abb, ".jpg"), ., 
           width=21, height=9, units="in", dpi=300)
  
  ggpubr::ggarrange(plotlist=sf.N, ncol=2, nrow=2) %>%
    ggsave(paste0("eda/maps/sf_N_", m.abb, ".jpg"), ., 
           width=6, height=6, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sf.S, ncol=2, nrow=2) %>%
    ggsave(paste0("eda/maps/sf_S_", m.abb, ".jpg"), ., 
           width=6, height=6, units="in", dpi=300)
  
}



agg$lLAM %>% select(sppName, id, el, model, mn, q05, q95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mn", "q05", "q95")) %>%
  ggplot(aes(el, mn_WY-mn_Y)) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_hline(yintercept=0) + facet_wrap(~sppName)


agg$lLAM %>% select(sppName, id, el, model, mn, q05, q95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mn", "q05", "q95")) %>%
  ggplot(aes(mn_WY, mn_Y)) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_abline() + facet_wrap(~sppName)

agg$llam %>% select(sppName, id, el, model, mn, q05, q95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mn", "q05", "q95")) %>%
  ggplot(aes(mn_WY, mn_Y)) + 
  geom_point(alpha=0.5, size=0.5, shape=1) + 
  geom_abline() + facet_wrap(~sppName, scales="free")

agg$lLAM %>% select(sppName, id, el, model, mn, q05, q95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mn", "q05", "q95")) %>%
  ggplot(aes(el, (q95_WY-q05_WY)-(q95_Y-q05_Y))) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_hline(yintercept=0) + facet_wrap(~sppName)





d.i$grd_W.sf %>% 
  left_join(., X_all, by="id") %>%
  mutate(N_lo=lLAM.i$q05[match(id, lLAM.i$id)],
         N_md=lLAM.i$q50[match(id, lLAM.i$id)],
         N_hi=lLAM.i$q95[match(id, lLAM.i$id)]) %>%
  ggplot(aes(R_grwnDD0, exp(N_md))) + geom_point() 




bdm_27 <- map_dfr(1:80, ~tibble(llam=filter(agg$llam, spp==.x & 
                                              plot %in% which(d.i$IJ$Y==1) &
                                              model=="Joint")$mn,
                                lLAM=filter(agg$lLAM, spp==.x & 
                                              id==527158 & 
                                              model=="Joint")$mn,
                                spp=.x))

ggplot(bdm_27, aes(exp(lLAM), exp(llam))) + geom_point()

exp(filter(agg$lLAM, spp==1 & id==527158 & model=="Joint")$mn)/
  exp(filter(agg$llam, spp==1 & plot %in% which(d.i$IJ$Y==1) & model=="Joint")$mn)





out <- stan("code/mods/Y_fwdSearch_PSI_test.stan", data=read_rdump("data/stan_data/test_PSI.Rdump"), chains=1, iter=500)