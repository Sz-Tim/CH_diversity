# 06_figures
# Tim Szewczyk
#
# Script for ms figures and eda 



########------------------------------------------------------------------------
## SET UP
########------------------------------------------------------------------------

# libraries etc
library(tidyverse); library(sf)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
source("code/00_fn.R"); source("../1_opfo/code/00_fn.R")

ms_dir <- "ms/1_Ecography/1/"
mod_col <- c("Joint"="#7b3294", "Structured"="#008837")


# auxilliary datasets
site_i <- read_csv("../1_opfo/data/opfo_siteSummaryProcessed.csv")
lc_rast <- raster::raster("../2_gis/data/VD_21781/lc_21781.tif")
#   table(raster::values(.)) %>%
#   setNames(., readxl::read_xlsx("../1_opfo/data/landcover_id.xlsx", 1)$LC)
lc_sum <- readxl::read_xlsx("../1_opfo/data/landcover_id.xlsx", 1) %>%
  mutate(VD=c(52381568, 24099, 2361482, 14810092, 8926571, 7856914, 95030, 
              101983, 6954192, 571985, 1114617, 960859, 5151839, 603243, 8666373), 
         VD_prop=VD/sum(VD))
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp") 
VD <- st_union(VD_raw %>% filter(!grepl("Lac ", NAME)))
dem <- raster::raster("../2_gis/data/VD_21781/dem_VD_21781.tif") %>%
  raster::mask(., st_zm(VD_raw %>% filter(!grepl("Lac ", NAME))))
world <- st_read("../2_gis/data/world/World_Countries__Generalized_.shp") %>%
  st_transform(st_crs(VD_raw)) %>%
  select(COUNTRY)

site.sf <- st_read("../2_gis/data/VD_21781/site_env_sf.shp")
names(site.sf) <- c(read_csv("../2_gis/data/VD_21781/site_env_sf_names.csv")$full,
                    "geometry")
ants_raw <- load_ant_data(str_type="all", clean_spp=T)
ants <- list(pub=bind_rows(ants_raw$pub, 
                      ants_raw$str %>% filter(TypeOfSample != "soil") %>%
                        select(TubeNo, SPECIESID, SampleDate)), 
             str=ants_raw$str %>% filter(TypeOfSample == "soil") %>%
               select(TubeNo, SPECIESID, SampleDate))
ants$all <- bind_rows(ants$pub %>% mutate(source="p"),
                      ants$str %>% mutate(source="s"))
tax_i <- read_csv("data/tax_i.csv")


# model inputs and posteriors
d.ls <- map(paste0("data/opt/", c("WY", "Y"), "__opt_var_set_ls.rds"), readRDS)
d.i <- map(paste0("data/opt/", c("WY", "Y"), "__opt_var_set_i.rds"), readRDS)
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


# summaries
det_Y <- which(colSums(d.ls[[1]]$Y)>0)

site.mns <- ants_raw$str %>% filter(TypeOfSample=="soil") %>%
  st_set_geometry(NULL) %>% 
  mutate(spp_fac=factor(SPECIESID)) %>%
  group_by(BDM, Plot_id) %>%
  summarise(nTubes=n(), S=n_distinct(SPECIESID), 
            H=vegan::diversity(tabulate(spp_fac)), mnt25=mean(mnt25)) %>%
  group_by(BDM) %>% 
  summarise(mnTubes=sum(nTubes)/25, seTubes=sd(nTubes)/5,
            mnS=sum(S)/25, seS=sd(S)/5,
            mnH=sum(H)/25, seH=sd(H)/5,
            el=mean(mnt25))

det.base <- ants$all %>% sf::st_set_geometry(NULL) %>%
  rename(sppName=SPECIESID) %>%
  group_by(sppName, source) %>% summarise(nDet=n()) %>%
  group_by(sppName) %>% mutate(pDet=nDet/sum(nDet))


# ggplot theme details
talk_fonts <- theme(panel.grid=element_blank(),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    strip.text=element_text(size=16),
                    title=element_text(size=18))
ms_fonts <- theme(panel.grid=element_blank(),
                  axis.text=element_text(size=11),
                  axis.title=element_text(size=12),
                  legend.text=element_text(size=10),
                  legend.title=element_text(size=12),
                  strip.text=element_text(size=12),
                  title=element_text(size=13))





################################################################################
############-------- MANUSCRIPT FIGURES
################################################################################

########------------------------------------------------------------------------
## VD sampling map
########------------------------------------------------------------------------

png(paste0(ms_dir, "figs/map_VD.png"), 
    height=7, width=7, res=400, units="in")
par(mar=c(0.5, 0.5, 0, 0), fig=c(0, 1, 0, 1))
raster::plot(dem, legend=F, axes=F, box=F, 
             col=colorRampPalette(c("gray45", "white"))(255))
raster::scalebar(d=10000, xy=c(570000, 162000), below="km", 
                 label=c(0, 5, 10), type="bar")
plot(VD, add=TRUE, lwd=0.5)
plot(select(d.i[[1]]$grd_W.sf, inbd), add=TRUE, 
     col=NA, fill=NA, border="gray10", lwd=0.2)
plot(select(ants$pub, TubeNo), add=TRUE, 
     col=rgb(5/256,113/256,176/256,0.75), cex=0.4)
plot(select(site.sf, el), add=TRUE, col=NA, fill=NA, border="#b2182b", lwd=1.5)
legend(565000, 182000, pch=c(0,1), pt.cex=c(1.2, 1), y.intersp=1.25,
       col=c("#b2182b", rgb(5/256,113/256,176/256,0.75)),
       legend=c("Structured\nsampling site", "Public sample"), 
       bty="n")
dev.off()

p <- ggplot(world) + geom_sf() + xlim(-3e5, 1.5e6) + ylim(-4e5, 6e5) + 
  geom_sf(data=st_union(VD_raw), fill="#fdd49e", colour="black", size=0.5) +
  theme(axis.line=element_blank(), axis.ticks=element_blank(), 
        axis.text=element_blank(), panel.border=element_rect(size=2))
ggsave(paste0(ms_dir, "figs/map_inset.png"), p, height=3, width=5, units="in")





########------------------------------------------------------------------------
## Posterior richness, abundance, diversity
########------------------------------------------------------------------------

el.pred <- data.frame(el=seq(300, 3000, by=10))
S.loess <- bind_rows(
  agg$S_R %>%
    group_by(model, id, el) %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    select(-data) %>%
    mutate(el=el.pred$el,
           Scale="Regional"), 
  agg$S_L %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    mutate(el=el.pred$el,
           Scale="Local",
           loess_lo=if_else(loess_lo < 0, 0, loess_lo)) 
) %>% mutate(Scale=factor(Scale, levels=c("Regional", "Local")))
H.loess <- bind_rows(
  agg$H %>%
    group_by(model, id, el) %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    select(-data) %>%
    mutate(el=el.pred$el,
           Scale="Regional"), 
  agg$H_L %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    mutate(el=el.pred$el,
           Scale="Local") 
) %>% mutate(Scale=factor(Scale, levels=c("Regional", "Local")))
lam.loess <- bind_rows(
  agg$tot_LAM %>%
    group_by(model, id, el) %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025/100 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median/100 ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975/100 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    select(-data) %>%
    mutate(el=el.pred$el,
           Scale="Regional"), 
  agg$tot_lam %>% group_by(model) %>% nest() %>%
    mutate(loess_lo=map(data, ~predict(loess(L025 ~ el, data=.x), el.pred)),
           loess_md=map(data, ~predict(loess(median ~ el, data=.x), el.pred)),
           loess_hi=map(data, ~predict(loess(L975 ~ el, data=.x), el.pred))) %>%
    unnest(cols=c(loess_lo, loess_md, loess_hi)) %>%
    mutate(el=el.pred$el,
           Scale="Local") 
) %>% mutate(Scale=factor(Scale, levels=c("Regional", "Local")))

S.p <- ggplot(S.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi, fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md, colour=model)) + 
  scale_colour_manual("Model", values=mod_col) + 
  scale_fill_manual("Model", values=mod_col) + 
  xlim(350, 3000) + 
  facet_grid(Scale~., scales="free_y") + 
  labs(x="Elevation (m)", y="Richness") + 
  ms_fonts + theme(panel.spacing=unit(3, "mm"),
                   axis.title.y=element_text(margin=margin(r=6)))
H.p <- ggplot(H.loess, aes(el)) + 
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi, fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md, colour=model)) + 
  scale_colour_manual("Model", values=mod_col) + 
  scale_fill_manual("Model", values=mod_col) + 
  xlim(350, 3000) + 
  facet_grid(Scale~., scales="free_y") + 
  labs(x="Elevation (m)", y="Diversity") + 
  ms_fonts + theme(panel.spacing=unit(3, "mm"),
                   axis.title.y=element_text(margin=margin(l=12, r=6)))
lam.p <- ggplot(lam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi, fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md, colour=model)) + 
  scale_colour_manual("Model", values=mod_col) + 
  scale_fill_manual("Model", values=mod_col) + 
  xlim(350, 3000) + scale_y_continuous(labels=scales::comma, limits=c(0, NA)) + 
  facet_grid(Scale~., scales="free_y") + 
  labs(x="Elevation (m)", y="Colony intensity") + 
  ms_fonts + theme(panel.spacing=unit(3, "mm"), 
                   axis.title.y=element_text(margin=margin(l=12, r=6)))

p <- ggpubr::ggarrange(S.p, H.p, lam.p, nrow=1, widths=c(.9, .9, 1),
                       labels=c("a.", "b.", "c."), label.x=0.02,
                       common.legend=T, legend="bottom") 
ggsave(paste0(ms_dir, "figs/el_patterns.png"), p, width=10, height=5, units="in")





########------------------------------------------------------------------------
## Posterior species responses
########------------------------------------------------------------------------

parName.df <- tibble(par=c("intercept", 
                           "R_grwnDD0", "R_grwnDD0_sq", 
                           "R_AP", 
                           "R_Forest",
                           "R_rdLen", 
                           "L_SoilTSt", 
                           "L_Pasture", "L_Crop", 
                           "L_CnpyOpn", "L_CnpyMxd"),
                     full=c("Intercept", 
                            "Growing deg. days",
                            "Growing deg. days (sq)", 
                            "Annual precipation",
                            "Forest proportion", 
                            "Road length", 
                            "Soil temperature", 
                            "Pasture", "Crop", 
                            "Open canopy", "Mixed canopy"))
b_post <- agg$b %>% filter(ParName != "intercept") %>%
  mutate(Scale=factor(str_sub(ParName, 1L, 1L), levels=c("R", "L"), 
                      labels=c("Regional", "Local")),
         Par=parName.df$full[match(ParName, parName.df$par)],
         Par=factor(Par, levels=rev(unique(parName.df$full))))
beta_post <- agg$beta %>% filter(ParName != "intercept") %>%
  mutate(Scale=factor(str_sub(ParName, 1L, 1L), levels=c("R", "L"), 
                      labels=c("Regional", "Local")),
         Par=parName.df$full[match(ParName, parName.df$par)],
         Par=factor(Par, levels=rev(unique(parName.df$full))))

p.A <- ggplot(b_post, aes(x=mean, y=Par, fill=model, colour=model)) +
  ggridges::geom_density_ridges(colour="gray30", alpha=0.75, scale=0.7, 
                                size=0.25, rel_min_height=0.001) +
  geom_point(data=filter(beta_post, model=="Joint"), shape=1, 
             position=position_nudge(y=-0.075)) + 
  geom_linerange(data=filter(beta_post, model=="Joint"), aes(xmin=L10, xmax=L90), 
                 size=0.75, position=position_nudge(y=-0.075)) + 
  geom_linerange(data=filter(beta_post, model=="Joint"), aes(xmin=L025, xmax=L975), 
                 size=0.3, position=position_nudge(y=-0.075)) + 
  geom_point(data=filter(beta_post, model!="Joint"), shape=1, 
             position=position_nudge(y=-0.2)) + 
  geom_linerange(data=filter(beta_post, model!="Joint"), aes(xmin=L10, xmax=L90), 
                 size=0.75, position=position_nudge(y=-0.2)) + 
  geom_linerange(data=filter(beta_post, model!="Joint"), aes(xmin=L025, xmax=L975), 
                 size=0.3, position=position_nudge(y=-0.2)) + 
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("Model", values=mod_col) + 
  scale_fill_manual("Model", values=mod_col) + 
  facet_grid(Scale~., scales="free_y", space="free_y") +
  labs(x="Posterior slope", y="") +
  ms_fonts + 
  theme(panel.grid.major.y=element_line(size=0.1, colour="gray30"),
        legend.position="bottom")

# p.B <- agg$b %>% 
#   mutate(Scale=str_sub(ParName, 1L, 1L),
#          Scale=factor(if_else(Scale=="i", "R", Scale),
#                       levels=c("R", "L"), labels=c("Regional", "Local")), 
#          Par=parName.df$full[match(ParName, parName.df$par)],
#          Par=factor(Par, levels=rev(unique(parName.df$full))),
#          SppInY=c("Species not in Y", 
#                   "Species in Y")[(spp %in% det_Y)+1]) %>%
#   ggplot(aes(x=abs(L975-L025), y=Par, fill=model)) + 
#   geom_vline(xintercept=0, colour="gray30", size=0.1) +
#   ggridges::geom_density_ridges(colour="gray30", alpha=0.75, scale=1, 
#                                 size=0.25, rel_min_height=0.001) + 
#   scale_fill_manual("Model", values=mod_col) +
#   labs(x="95% HPDI width", y="") +
#   facet_grid(Scale~SppInY, scales="free_y", space="free_y") +
#   ms_fonts + 
#   theme(panel.grid.major.y=element_line(size=0.1, colour="gray30"),
#         legend.position="bottom")
# 
# p <- ggpubr::ggarrange(p.A, p.B, nrow=1, labels=c("a.", "b."), 
#                        widths=c(0.7, 1), label.x=c(0.05, 0.07),
#                        common.legend=T, legend="bottom") 
# ggsave(paste0(ms_dir, "figs/slope_means+HDI.png"), p, width=10, height=5, units="in")

hdi_width.df <- agg$b %>% 
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Scale=factor(if_else(Scale=="i", "R", Scale),
                      levels=c("R", "L"), labels=c("Regional", "Local")), 
         Par=parName.df$full[match(ParName, parName.df$par)],
         Par=factor(Par, levels=rev(unique(parName.df$full))),
         SppInY=c("Species not in Y", 
                  "Species in Y")[(spp %in% det_Y)+1],
         HDI_w=L975-L025) %>%
  select(Scale, Par, sppName, SppInY, model, HDI_w) %>%
  pivot_wider(names_from="model", values_from="HDI_w") %>%
  mutate(HDI_diff=Joint-Structured) %>%
  filter(!is.na(HDI_diff))
hdi_width.sum <- hdi_width.df %>%
  group_by(Scale, Par, SppInY) %>%
  summarise(mn=mean(HDI_diff),
            se=sd(HDI_diff)/sqrt(n())) 

p.B <- ggplot(hdi_width.df, aes(y=Par, fill=SppInY, colour=SppInY)) + 
  ggridges::geom_density_ridges(aes(x=HDI_diff),
                                colour="gray30", alpha=0.75, scale=0.7, 
                                size=0.25, rel_min_height=0.001) + 
  geom_point(data=filter(hdi_width.sum, SppInY=="Species in Y"), 
             aes(x=mn), shape=1, position=position_nudge(y=-0.075)) + 
  geom_linerange(data=filter(hdi_width.sum, SppInY=="Species in Y"), 
                 aes(xmin=mn-2*se, xmax=mn+2*se), 
                 size=0.5, position=position_nudge(y=-0.075)) + 
  geom_point(data=filter(hdi_width.sum, SppInY!="Species in Y"), 
             aes(x=mn), shape=1, position=position_nudge(y=-0.2)) + 
  geom_linerange(data=filter(hdi_width.sum, SppInY!="Species in Y"), 
                 aes(xmin=mn-2*se, xmax=mn+2*se), 
                 size=0.5, position=position_nudge(y=-0.2)) + 
  geom_vline(xintercept=0, colour="gray30", size=0.1) +
  labs(x="Change in 95% HPDI width\n(Joint - Structured)", y="") +
  facet_grid(Scale~., scales="free_y", space="free_y") +
  ms_fonts + 
  scale_fill_brewer("", type="qual", palette=3, direction=-1) + 
  scale_colour_brewer("", type="qual", palette=3, direction=-1) + 
  theme(panel.grid.major.y=element_line(size=0.1, colour="gray30"),
        legend.position="bottom")

p <- ggpubr::ggarrange(p.A, p.B, nrow=1, labels=c("a.", "b."), 
                       label.x=c(0.05, 0.07),
                       common.legend=F, legend="bottom") 
ggsave(paste0(ms_dir, "figs/slope_means+HDI.png"), p, width=10, height=5, units="in")





########------------------------------------------------------------------------
## Beta diversity
########------------------------------------------------------------------------

library(betapart)
lam.site.ls <- agg$lam %>% 
  mutate(site=str_pad(str_sub(as.character(id), 1, -5), 2, "left", "0"),
         BDM=arrange(site_i, BDM_id)$BDM[as.numeric(site)]) %>%
  select(model, sppName, site, BDM, id, median) %>%
  pivot_wider(names_from="sppName", values_from="median")
beta.lam.df <- bind_rows(site.mns %>% mutate(model="Joint"),
                         site.mns %>% mutate(model="Structured")) %>% 
  filter(BDM %in% lam.site.ls$BDM) %>%
  mutate(beta.BRAY.BAL=NA, 
         beta.BRAY.GRA=NA,
         beta.BRAY=NA)
for(i in 1:nrow(beta.lam.df)) {
  BDM_i <- beta.lam.df$BDM[i]
  mod_i <- beta.lam.df$model[i]
  beta_i <- beta.multi.abund(filter(lam.site.ls, BDM==BDM_i & model==mod_i)[,5:84])
  beta.lam.df$beta.BRAY.BAL[i] <- beta_i$beta.BRAY.BAL
  beta.lam.df$beta.BRAY.GRA[i] <- beta_i$beta.BRAY.GRA
  beta.lam.df$beta.BRAY[i] <- beta_i$beta.BRAY
}

beta.df <- beta.lam.df %>% 
  pivot_longer(10:12, names_to="BetaPart", values_to="Beta") %>%
  mutate(BetaPart=factor(BetaPart, 
                         levels=paste0("beta.BRAY", c("", ".BAL", ".GRA")), 
                         labels=c("Overall", "Balanced\nvariation",
                                  "Abundance\ngradient"))) 
p <- ggplot(beta.df, aes(el, Beta, colour=model, fill=model,
                         linetype=BetaPart, shape=BetaPart)) + 
  geom_point(size=0.9) +
  stat_smooth(data=filter(beta.df, BetaPart=="Overall" &
                            model=="Joint"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x) + 
  stat_smooth(data=filter(beta.df, BetaPart == "Balanced\nvariation" &
                            model=="Joint"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x+I(x^2)) + 
  stat_smooth(data=filter(beta.df, BetaPart == "Abundance\ngradient" &
                            model=="Joint"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x) + 
  stat_smooth(data=filter(beta.df, BetaPart=="Overall" &
                            model=="Structured"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x+I(x^2)) + 
  stat_smooth(data=filter(beta.df, BetaPart == "Balanced\nvariation" &
                            model=="Structured"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x+I(x^2)) + 
  stat_smooth(data=filter(beta.df, BetaPart == "Abundance\ngradient" &
                            model=="Structured"), 
              alpha=0.25, size=0.5, se=F, method="lm", formula=y~x+I(x^2)) + 
  scale_colour_manual("Model", values=mod_col) + 
  scale_fill_manual("Model", values=mod_col) + 
  scale_linetype_manual(expression(beta~partition), values=c(1,5,3)) +
  scale_shape_manual(expression(beta~partition), values=c(1,4,0)) +
  labs(x="Elevation (m)", y=expression('Within-site'~beta~diversity)) + 
  ylim(0,1) + ms_fonts + 
  guides(linetype=guide_legend(override.aes=list(colour="gray30"))) +
  theme(legend.key.height=unit(9, "mm"))
ggsave(paste0(ms_dir, "figs/beta_diversity.png"), p, width=5, height=3.5)




########------------------------------------------------------------------------
## Elevational assemblages
########------------------------------------------------------------------------

# DPCoA on plot-level posteriors
library(ade4); library(adegraphics)
adegpar(psub.cex=1.5, 
        ppoints=list(cex=0.4, alpha=0.8), 
        plines=list(lwd=0.25),
        pellipse=list(alpha=0.4, axes=list(draw=FALSE)))

col_region <- c(Alps="#56B4E9", Jura="#E69F00", Plateau="#666666")

tax_dist <- tax_i %>% select(contains("Full")) %>% as.data.frame
rownames(tax_dist) <- paste(tax_dist$FullGen, tax_dist$FullSpp, sep="_")


lam.site.ls <- agg$lam %>% 
  mutate(site=str_pad(str_sub(as.character(id), 1, -5), 2, "left", "0"),
         BDM=arrange(site_i, BDM_id)$BDM[as.numeric(site)]) %>%
  filter(model=="Joint") %>%
  select(sppName, site, BDM, id, median) %>%
  pivot_wider(names_from="sppName", values_from="median")
lam.plot.mx <- as.matrix(lam.site.ls[,4:83])

env.plot <- tibble(el=as.factor((d.i[[1]]$V[,2] %/% 100) * 100),
                   el_cont=d.i[[1]]$V[,2],
                   region=site_i$region[match(lam.site.ls$BDM, site_i$BDM)]) %>%
  mutate(region=factor(ifelse(region=="Low", "Plateau", region)))

lam.dpcoa <- dpcoa(as.data.frame(lam.plot.mx), 
                   vegan::taxa2dist(tax_dist, varstep=T), scannf=F, nf=4)
p.A1 <- s.class(lam.dpcoa$li, env.plot$el, 
               psub=list(text="a.", position="topleft"),
               col=terrain.colors(n_distinct(env.plot$el)))
p.A2 <- plotEig(lam.dpcoa$eig[1:30], lam.dpcoa$nf, pbackground=list(box=TRUE), 
                psub=list(text="Eigenvalues", cex=1), plot=FALSE)
p.A <- insert(p.A2, p.A1, posi=c(0.02, 0.02), ratio=0.25)
p.B <- s.class(lam.dpcoa$li, env.plot$region, 
               psub=list(text="b.", position="topleft"), col=col_region)

png(paste0(ms_dir, "figs/DPCoA.png"), height=4, width=8, res=400, units="in")
ADEgS(list(p.A, p.B), layout=c(1,2))
dev.off()


# Raw samples: genus proportions
genProp.df <- ants$all %>%
  mutate(el=raster::extract(dem, .),
         elCat=c("Montane", "Plateau")[(el < 1000) + 1],
         elCat=factor(elCat, levels=c("Plateau", "Montane")),
         genFull=tax_i$FullGen[match(SPECIESID, tax_i$species)],
         source=factor(source, labels=c("Presence\nonly", "Structured\nabundance"))) %>%
  st_set_geometry(NULL) %>% filter(!is.na(elCat))
gen_1pct <- unique((genProp.df %>% 
                      group_by(genFull, elCat, source) %>% summarise(N=n()) %>%
                      group_by(elCat, source) %>% mutate(prop=N/sum(N)) %>%
                      ungroup %>% filter(prop > 0.01))$genFull)
genProp.df <- genProp.df %>%
  mutate(genus_1pct=if_else(genFull %in% gen_1pct, genFull, "Other"), 
         genus_1pct=factor(genus_1pct, levels=c(gen_1pct, "Other")))

p <- ggplot(genProp.df, aes(elCat, fill=genus_1pct)) + 
  geom_hline(yintercept=0, colour="gray30", size=0.25) + 
  geom_bar(position="fill", colour="black", size=0.25) +
  scale_fill_viridis_d("Genus", option="D") + 
  facet_grid(.~source) + 
  ms_fonts + 
  labs(x="", y="Proportion of samples")

ggsave(paste0(ms_dir, "figs/genus_assemblages.png"), p, width=5, height=5)





########------------------------------------------------------------------------
## Taxonomic bias
########------------------------------------------------------------------------

p <- agg$D %>%
  mutate(genFull=tax_i$FullGen[match(sppName, tax_i$species)],
         sppName=str_replace(str_remove(sppName, "-GR"), "_", ".")) %>%
  ggplot(aes(x=median, xmin=L05, xmax=L95, y=sppName,
             colour=sign(L05-1)==sign(L95-1))) + 
  geom_vline(xintercept=1, colour="gray") +
  geom_point() + geom_linerange() + 
  facet_grid(genFull~., scales="free_y", space="free_y") +
  scale_colour_manual(values=c("gray", "darkred"), guide=F) +
  ms_fonts  + theme(axis.text.y=element_text(hjust=1, vjust=0.5),
                    strip.background=element_blank(),
                    strip.text=element_blank(), 
                    panel.spacing=unit(0.5, "mm")) + 
  labs(y="", x="Proportional Bias")
ggsave(paste0(ms_dir, "figs/D_opt.png"), p, width=4, height=13)











################################################################################
############-------- DRAFTS AND EDA
################################################################################

########------------------------------------------------------------------------
## DATA SUMMARIES
########------------------------------------------------------------------------

spp.det_byEl <- bind_rows(as_tibble(d.ls[[1]]$Y) %>% 
                            mutate(el=d.i[[1]]$V[,"el"], 
                                   Dataset="Structured"), 
                          as_tibble(d.ls[[1]]$W) %>%
                            mutate(el=d.i[[1]]$X[1:d.ls[[1]]$K, "el"],
                                   Dataset="Cit. Sci.")) %>%
  pivot_longer(1:80, names_to="Species", values_to="nDet") %>%
  uncount(nDet) %>%
  mutate(Genus=tax_i$FullGen[match(Species, tax_i$species)],
         Subf=tax_i$FullSF[match(Species, tax_i$species)])

p <- spp.det_byEl %>%
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
ggsave("eda/obs_elRange_Spp.pdf", p, height=12, width=17, units="in")

p <- spp.det_byEl %>%
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
ggsave("eda/obs_elRange_Spp_byGen.pdf", p, height=12, width=17, units="in")

p <- spp.det_byEl %>%
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
ggsave("eda/obs_elRange_Gen_tall.pdf", p, height=10, width=5, units="in")

p <- spp.det_byEl %>%
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
ggsave("eda/obs_elRange_Gen_square.pdf", p, height=6, width=8, units="in")

p <- spp.det_byEl %>%
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
ggsave("eda/obs_elRange_SF.png", p, height=5, width=5, units="in", dpi=300)

ants$all$LC <- raster::extract(lc_rast, ants$all)
lc.prop.df <- lc_sum %>%
  mutate(W=c(table(filter(ants$all, source=="p")$LC)),
         W_prop=W/sum(W),
         Y=c(table(filter(ants$all, source=="s")$LC)),
         Y_prop=Y/sum(Y),
         LC=factor(LC, levels=unique(LC), 
                   labels=c("Open Canopy", "Quarry", "Orchard/Vineyard",
                            "Evergreen Forest", "Deciduous Forest", 
                            "Mixed Forest", "Gravel", "Wetland", "Edge", 
                            "Marsh", "Scree", "Dry Prairie", "Roadside", 
                            "Alluvial Zone", "Built Zone"))) 
lc.prop.long <- lc.prop.df %>%
  select(1,5,7,9) %>% pivot_longer(2:4, names_to="Source", values_to="Prop") %>%
  mutate(Source=factor(Source, levels=c("VD_prop", "Y_prop", "W_prop"),
                       labels=c("Region", "Structured samples", "Public samples")))
p <- ggplot(lc.prop.long, aes(x=Prop, y=LC, fill=Source)) + 
  geom_vline(xintercept=0, colour="grey30", size=0.25) +
  geom_bar(stat="identity", position="dodge", colour="grey30", size=0.25) + 
  scale_fill_brewer(name="", type="qual", palette=3, 
                    direction=-1, guide=guide_legend(reverse=T)) +
  labs(title="Distribution of habitat types", x="Proportion", y="") + 
  theme(legend.position=c(0.8, 0.5)) + talk_fonts
ggsave("eda/obs_LC_distr.png", p, height=8, width=7, units="in", dpi=300)

lc.prop.df %>% select(1, 4, 6) %>%
  column_to_rownames("LC") %>% as.matrix %>%
  chisq.posthoc.test::chisq.posthoc.test() %>%
  filter(Value=="p values") %>% 
  summarise(nDiff=sum(W<0.05), 
            pDiff=nDiff/n())
lc.prop.df %>% select(1, 4, 8) %>%
  column_to_rownames("LC") %>% as.matrix %>%
  chisq.posthoc.test::chisq.posthoc.test() %>%
  filter(Value=="p values") %>% 
  summarise(nDiff=sum(Y<0.05), 
            pDiff=nDiff/n())











########------------------------------------------------------------------------
## SLOPES
########------------------------------------------------------------------------

## Aggregate ---------------------------

p <- agg$beta %>% filter(ParName != "intercept") %>%
  mutate(cov=as.numeric(str_sub(Parameter, 6L, -1L))) %>%
  mutate(Scale=factor(str_sub(ParName, 1L, 1L), levels=c("R", "L"), 
                      labels=c("Regional", "Local")),
         Par=str_sub(ParName, 3L, -1L)) %>%
  ggplot(aes(x=Par, y=mean, colour=model)) + 
  geom_hline(yintercept=0, linetype=1, colour="gray", size=0.35) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_linerange(aes(ymin=L25, ymax=L75), 
                 position=position_dodge(width=0.5), size=1.5) + 
  geom_linerange(aes(ymin=L05, ymax=L95),
                 position=position_dodge(width=0.5), size=0.5) + 
  facet_grid(Scale~., scales="free_y", space="free_y") +
  coord_flip() + 
  labs(x="", y=expression(beta), title="Aggregate slopes") + 
  scale_colour_manual("", values=mod_col) + talk_fonts 
ggsave("eda/beta_opt.png", p, width=9, height=8)



## Genus -------------------------------

p <- ggplot(agg$B, aes(x=mean, y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_fill_manual(values=mod_col) +
  labs(x="Genus-level responses (B)", y="", 
       title="Distribution of genus-level responses") 
ggsave("eda/Bg_opt_mn_ridges_byModel.png", p, width=7, height=6)

p <- ggplot(filter(agg$B, ParName != "intercept"), 
            aes(x=ParName, y=median, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_jitter(aes(shape=sign(L05)==sign(L95)), width=0.2, height=0) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  scale_colour_manual("", values=mod_col) + 
  labs(y="Genus-level response (B)", x="") +
  facet_wrap(~ParName, scales="free_x", nrow=1) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("eda/Bg_opt_scatter.png", p, width=12, height=4)

p <- ggplot(filter(agg$B, ParName != "intercept"), 
       aes(x=ParName, y=mean, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(L05)==sign(L95)), 
             position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=L25, ymax=L75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=L05, ymax=L95),
                 position=position_dodge(width=0.5), size=0.2) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  labs(y="Genus-level response (B)", x="") +
  facet_wrap(~genName, ncol=7) 
ggsave("eda/Bg_opt_byGenus.png", p, width=15, height=6)

p <- ggplot(filter(agg$B, ParName != "intercept"), 
       aes(x=genName, y=median, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(L05)==sign(L95)), 
             position=position_dodge(width=0.5)) + 
  geom_linerange(aes(ymin=L25, ymax=L75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=L05, ymax=L95),
                 position=position_dodge(width=0.5), size=0.2) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  labs(y="Genus-level response (B)", x="") +
  facet_wrap(~ParName, scales="free_x", nrow=2)
ggsave("eda/Bg_opt_byParam.png", p, width=14, height=8)

p <- agg$B %>% filter(ParName != "intercept") %>%
  mutate(effect=case_when(L05<0 & L95<0 ~ "negative",
                          L05>0 & L95>0 ~ "positive",
                          L05<0 & L95>0 ~ "no effect",
                          L05>0 & L95<0 ~ "no effect"),
         sig=sign(L05)==sign(L95)) %>%
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
ggsave("eda/Bg_opt_bar.png", p, width=7, height=8)

p <- filter(agg$B, ParName != "intercept") %>%
  filter(sign(L05)==sign(L95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/21, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  scale_fill_manual(values=mod_col) +
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of genera", x="")

agg$B %>% 
  filter(model=="Joint") %>%
  select(mean, ParName, genName) %>% 
  filter(ParName != "intercept") %>% 
  pivot_wider(names_from=ParName, values_from=mean) %>% 
  column_to_rownames("genName") %>% 
  data.matrix %>% 
  gplots::heatmap.2(col=colorRampPalette(c("blue4", "gray80", "red4")))
agg$B %>% 
  filter(model=="Structured") %>%
  select(mean, ParName, genName) %>% 
  filter(ParName != "intercept") %>% 
  pivot_wider(names_from=ParName, values_from=mean) %>% 
  column_to_rownames("genName") %>% 
  data.matrix %>% 
  gplots::heatmap.2(col=colorRampPalette(c("blue4", "gray80", "red4")))



## Species -----------------------------

p <- ggplot(agg$b, aes(x=mean, y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_fill_manual(values=mod_col) +
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") 
ggsave("eda/b_opt_mn_ridges_byModel.png", p, width=7, height=6)

p <- ggplot(agg$b, aes(x=abs(L95-L05), y=ParName, fill=spp %in% det_Y)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) + 
  scale_fill_brewer("Species detected in\nstructured samples", type="qual", palette=3) +
  labs(x="90% HDI width", y="", 
       title="Distribution of species-level response HDI widths") +
  facet_wrap(~model)
ggsave("eda/b_opt_HDI_ridges_byObsY.png", p, width=7, height=6)

p <- ggplot(agg$b, aes(x=abs(L95-L05), y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) + 
  scale_fill_manual(values=mod_col) +
  labs(x="90% HDI width", y="", 
       title="Distribution of species-level response HDI widths") +
  facet_wrap(~(spp %in% det_Y))
ggsave("eda/b_opt_obsY_HDI_ridges_byModel.png", p, width=7, height=6)

p <- ggplot(agg$b, aes(x=mean, y=ParName, fill=spp %in% det_Y)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_fill_brewer("Species detected in\nstructured samples", type="qual", palette=3) +
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~model)
ggsave("eda/b_opt_mn_ridges_byObsY.png", p, width=7, height=6)

p <- ggplot(agg$b, aes(x=mean, y=ParName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, size=0.25) +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_fill_manual(values=mod_col) +
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~(spp %in% det_Y))
ggsave("eda/b_opt_obsY_mn_ridges_byModel.png", p, width=7, height=6)

p <- agg$b %>% mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=mean, y=genName, fill=model)) + 
  ggridges::geom_density_ridges(alpha=0.75, scale=1, 
                                rel_min_height=0.01, size=0.25) +
  geom_point(aes(colour=model), size=1.5, shape="|") +
  geom_vline(xintercept=0, linetype=3, colour="gray30", size=0.5) +
  scale_colour_manual("", values=mod_col) + 
  scale_fill_manual("", values=mod_col) + 
  labs(x="Species-level responses (b)", y="", 
       title="Distribution of species-level responses") +
  facet_wrap(~ParName, scales="free_x", ncol=5) +
  theme(panel.grid.major.y=element_line(colour="gray30", size=0.2))
ggsave("eda/b_opt_aggGenus_byParam_ridges.png", p, width=12, height=7)

p <- agg$b %>% mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=mean, y=ParName, fill=model)) + 
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
ggsave("eda/b_opt_aggGenus_byGenus_ridges.png", p, width=13, height=7)

p <- agg$b %>% mutate(sfName=tax_i$FullSF[match(sppName, tax_i$species)]) %>%
  ggplot(aes(x=mean, y=sfName, fill=model)) + 
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
ggsave("eda/b_opt_aggSF_byParam_ridges.png", p, width=7, height=8)

p <- agg$b %>% mutate(sfName=tax_i$FullSF[match(sppName, tax_i$species)]) %>%
  ggplot(aes(x=mean, y=ParName, fill=model)) + 
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
ggsave("eda/b_opt_aggSF_bySF_ridges.png", p, width=6, height=5)

p <- agg$b %>% mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=ParName, y=mean, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(L05)==sign(L95)),
             position=position_dodge(width=0.5)) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  geom_linerange(aes(ymin=L25, ymax=L75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=L05, ymax=L95),
                 position=position_dodge(width=0.5), size=0.2) + 
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  facet_wrap(~sppName, ncol=10) + 
  labs(y="Species-level responses (b)", x="")
ggsave("eda/b_opt_bySpecies.png", p, width=14, height=12)

p <- agg$b %>% mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  filter(ParName != "intercept") %>%
  ggplot(aes(x=sppName, y=median, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray30", size=0.5) +
  geom_point(aes(shape=sign(L05)==sign(L95)),
             position=position_dodge(width=0.5)) + 
  scale_shape_manual("sig", values=c(1, 19)) +
  geom_linerange(aes(ymin=L25, ymax=L75), 
                 position=position_dodge(width=0.5), size=0.5) + 
  geom_linerange(aes(ymin=L05, ymax=L95),
                 position=position_dodge(width=0.5), size=0.2) + 
  coord_flip() + scale_colour_manual("", values=mod_col) + 
  facet_grid(genName~ParName, scales="free", space="free_y") +
  labs(y="Species-level responses (b)", x="") + 
  theme(strip.text.y=element_blank(), strip.background.y=element_blank(),
        panel.spacing.y=unit(0.05, 'cm'))
ggsave("eda/b_opt_byParam.png", p, width=12, height=12)

p <- agg$b %>% filter(cov != "1") %>% 
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Par=str_sub(ParName, 3L, -1L)) %>%
  mutate(effect=case_when(L05<0 & L95<0 ~ "negative",
                          L05>0 & L95>0 ~ "positive",
                          L05<0 & L95>0 ~ "ns",
                          L05>0 & L95<0 ~ "ns"),
         sig=sign(L05)==sign(L95)) %>%
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
ggsave("eda/b_opt_bar.png", p, width=10, height=8)

p <- agg$b %>% filter(cov != "1") %>% 
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Par=str_sub(ParName, 3L, -1L)) %>%
  mutate(effect=case_when(L25<0 & L75<0 ~ "negative",
                          L25>0 & L75>0 ~ "positive",
                          L25<0 & L75>0 ~ "ns",
                          L25>0 & L75<0 ~ "ns"),
         sig=sign(L25)==sign(L75)) %>%
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

p <- agg$b %>% filter(cov != "1") %>%
  filter(sign(L05)==sign(L95)) %>%
  group_by(model, ParName) %>% summarise(nSig=n()) %>%
  arrange(nSig) %>% arrange(model, nSig) %>% 
  mutate(ParName=factor(ParName, levels=unique(ParName))) %>%
  ggplot(aes(ParName, nSig/80, fill=model)) + 
  geom_hline(yintercept=0, colour="gray30") + 
  geom_bar(stat="identity", colour="gray30") + 
  scale_fill_manual(values=mod_col) +
  coord_flip() + facet_wrap(~model) + labs(y="Proportion of species", x="")

p <- agg$b %>% filter(cov != "1") %>%
  mutate(effect=case_when(L05<0 & L95<0 ~ "negative",
                          L05>0 & L95>0 ~ "positive",
                          L05<0 & L95>0 ~ "no effect",
                          L05>0 & L95<0 ~ "no effect"),
         sig=sign(L05)==sign(L95)) %>%
  group_by(model, sppName, sig) %>% summarise(nVar=n()) %>% 
  # filter(sig) %>% summary # 1-12 significant variables
  ggplot(aes(sig, nVar, fill=sppName %in% names(det_Y))) + 
  geom_boxplot() + ylim(0, NA) + 
  scale_fill_brewer("Species detected in\nstructured samples", type="qual")

p <- agg$b %>% filter(cov != "1") %>%
  select(model, sppName, ParName, mean) %>%
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  pivot_wider(names_from="model", values_from="mean") %>%
  filter(!is.na(Joint) & !is.na(Structured)) %>%
  ggplot(aes(Structured, Joint, colour=ParName)) + 
  geom_hline(yintercept=0, colour="gray90") + 
  geom_vline(xintercept=0, colour="gray90") + 
  geom_abline(colour="gray30", size=0.25) + geom_point(shape=1) + 
  scale_colour_brewer(type="qual", palette=2) + 
  facet_wrap(~genName, ncol=7) + ylim(-4,4) + xlim(-4,4)
ggsave("eda/b_opt_Y_WY_scatter_byGenus.png", p, width=13, height=6)

p <- agg$b %>% filter(cov != "1") %>%
  select(model, sppName, ParName, mean) %>%
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  pivot_wider(names_from="model", values_from="mean") %>%
  filter(!is.na(Joint) & !is.na(Structured)) %>%
  ggplot(aes(Structured, Joint)) + 
  geom_hline(yintercept=0, colour="gray90") + 
  geom_vline(xintercept=0, colour="gray90") + 
  geom_abline(colour="gray30", size=0.25) + geom_point(shape=1) + 
  # scale_colour_brewer(type="qual", palette=2) + 
  facet_wrap(~ParName) + ylim(-4,4) + xlim(-4,4)
ggsave("eda/b_opt_Y_WY_scatter_byParam.png", p, width=8, height=8)

agg$b %>% 
  filter(model=="Joint") %>%
  filter(ParName != "intercept") %>% 
  # filter(spp %in% det_Y) %>%
  select(mean, ParName, sppName) %>% 
  pivot_wider(names_from=ParName, values_from=mean) %>% 
  column_to_rownames("sppName") %>% 
  data.matrix %>% 
  gplots::heatmap.2(col=colorRampPalette(c("blue4", "gray80", "red4")), 
                    cexRow=0.5, cexCol=0.5)
agg$b %>% 
  filter(model=="Structured") %>%
  filter(ParName != "intercept") %>% 
  filter(spp %in% det_Y) %>%
  select(mean, ParName, sppName) %>% 
  pivot_wider(names_from=ParName, values_from=mean) %>% 
  column_to_rownames("sppName") %>% 
  data.matrix %>% 
  gplots::heatmap.2(col=colorRampPalette(c("blue4", "gray80", "red4")), 
                    cexRow=0.5, cexCol=0.5)







########------------------------------------------------------------------------
## LAMBDAS
########------------------------------------------------------------------------

## Regional ----------------------------

## total ----------
lLAM.loess <- agg$tot_lLAM %>%
  # filter(id %in% d.i[[1]]$X[,1]) %>%
  group_by(model, id, el) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L025 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L975 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
LAM.loess <- agg$tot_LAM %>%
  # filter(id %in% d.i[[1]]$X[,1]) %>%
  # filter(id > 1e5) %>%
  group_by(model, id, el) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L025 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L975 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

p <- ggplot(lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo)/100, ymax=exp(loess_hi)/100),
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md)/100)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare (log-axis)")
ggsave("eda/R_lLAMBDA_opt.png", p, width=6, height=5)

p <- ggplot(LAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo/100, ymax=loess_hi/100),
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md/100)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  # scale_y_continuous(labels=scales::comma, limits=c(0, NA)) +
  labs(x="Elevation (m)", y="Colonies per hectare")
ggsave("eda/R_LAMBDA_opt.png", p, width=6, height=5)


## by species -----
spp.lLAM.loess <- agg$lLAM %>%
  filter(id %in% d.i[[1]]$X[,1]) %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
spp.LAM.loess <- agg$LAM %>%
  filter(id %in% d.i[[1]]$X[,1]) %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

det.spp <- det.base %>% select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", str_pad(round(p*100), 3, "left", " "), "%", 
                         "\ns: ", str_pad(round(s*100), 3, "left", " "), "%"),
         model="Joint")

p <- ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/R_lLAMBDA_opt_el_Species.png", p, width=14, height=10)

p <- ggplot(spp.LAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo/100, ymax=loess_hi/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species.png", p, width=14, height=10)

p <- ggplot(spp.LAM.loess, aes(el, group=model, colour=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median/100), size=0.5, shape=1, alpha=0.25) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species_pts.png", p, width=14, height=10)

p <- ggplot(spp.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/R_lLAMBDA_opt_el_Species_free.png", p, width=16, height=10)

p <- ggplot(spp.LAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo/100, ymax=loess_hi/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md/100)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species_free.png", p, width=16, height=10)

p <- ggplot(spp.LAM.loess, aes(el, group=model, colour=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median/100), size=0.5, shape=1, alpha=0.25) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~sppName,scales="free_y", ncol=10)
ggsave("eda/R_LAMBDA_opt_el_Species_free_pts.png", p, width=16, height=10)


## by genus -----
det.gen <- det.base %>% 
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(genName, source) %>% summarise(nDet=sum(nDet)) %>%
  group_by(genName) %>% mutate(pDet=nDet/sum(nDet)) %>%
  select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", round(p*100), "%", 
                         "\ns: ", round(s*100), "%"),
         model="Joint")

gen.lLAM.loess <- agg$lLAM %>% 
  filter(id %in% d.i[[1]]$X[,1]) %>%
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(model, id, el, genName) %>%
  summarise(across(one_of("L05", "median", "L95"), ~log(sum(exp(.))))) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
gen.LAM.loess <- agg$LAM %>% 
  filter(id %in% d.i[[1]]$X[,1]) %>%
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(model, id, el, genName) %>%
  summarise(across(one_of("L05", "median", "L95"), sum)) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

p <- ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/R_lLAMBDA_opt_el_Genus.png", p, width=15, height=6)

p <- ggplot(gen.LAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo/100, ymax=loess_hi/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus.png", p, width=15, height=6)

p <- ggplot(gen.LAM.loess, aes(el, group=model, colour=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median/100), size=0.5, shape=1, alpha=0.25) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus_pts.png", p, width=15, height=6)

p <- ggplot(gen.lLAM.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/R_lLAMBDA_opt_el_Genus_free.png", p, width=15, height=6)

p <- ggplot(gen.LAM.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo/100, ymax=loess_hi/100), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md/100)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus_free.png", p, width=15, height=6)

p <- ggplot(gen.LAM.loess, aes(el, group=model, colour=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median/100), size=0.5, shape=1, alpha=0.25) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_y_continuous(labels=scales::comma) +
  labs(x="Elevation (m)", y="Colonies per hectare") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/R_LAMBDA_opt_el_Genus_free_pts.png", p, width=15, height=6)




## Local -------------------------------

agg$lam %>% 
  group_by(model, plot, el) %>%
  summarise(across(one_of("L025", "median", "L975"), sum)) %>%
  ggplot(aes(el, median, ymin=L025, ymax=L975, colour=model)) + 
  geom_point(size=0.5, shape=1) + geom_linerange(size=0.25) +
  scale_colour_manual(values=mod_col)

llam.loess <- agg$tot_llam %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
lam.loess <- agg$tot_lam %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

p <- ggplot(llam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi), fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md), colour=model)) + 
  geom_point(data=site.mns, aes(y=mnTubes)) +
  geom_linerange(data=site.mns,
                 aes(ymin=mnTubes-seTubes, ymax=mnTubes+seTubes), size=0.25) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Mean abundance per plot")
ggsave("eda/L_llambda_opt.png", p, width=6, height=5)

p <- ggplot(lam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi, fill=model), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md, colour=model)) + 
  geom_point(data=site.mns, aes(y=mnTubes)) +
  geom_linerange(data=site.mns,
                 aes(ymin=mnTubes-seTubes, ymax=mnTubes+seTubes), size=0.25) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Mean abundance per plot")
ggsave("eda/L_lambda_opt.png", p, width=6, height=5)


## by species -----
spp.llam.loess <- agg$llam %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
spp.lam.loess <- agg$lam %>%
  group_by(model, sppName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

det.spp <- det.base %>% select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", str_pad(round(p*100), 3, "left", " "), "%", 
                         "\ns: ", str_pad(round(s*100), 3, "left", " "), "%"),
         model="Joint")

p <- ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=exp(loess_lo), ymax=exp(loess_hi)), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=exp(loess_md))) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=5) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log)") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_llambda_opt_el_Species.png", p, width=14, height=10)

p <- ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(median), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_llambda_opt_el_Species_pts.png", p, width=14, height=10)

p <- ggplot(spp.lam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_lambda_opt_el_Species.png", p, width=14, height=10)

p <- ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(median), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, ncol=10)
ggsave("eda/L_lambda_opt_el_Species_pts.png", p, width=14, height=10)

p <- ggplot(spp.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/L_llambda_opt_el_Species_free.png", p, width=16, height=10)

p <- ggplot(agg$llam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=exp(median), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_llambda_opt_el_Species_free_pts.png", p, width=16, height=10)

p <- ggplot(spp.lam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  geom_text(data=det.spp, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_lambda_opt_el_Species_free.png", p, width=16, height=10)

p <- ggplot(agg$lam, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median, colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~sppName, scales="free_y", ncol=10)
ggsave("eda/L_lambda_opt_el_Species_free_pts.png", p, width=16, height=10)


## by genus -----
det.gen <- det.base %>% 
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(genName, source) %>% summarise(nDet=sum(nDet)) %>%
  group_by(genName) %>% mutate(pDet=nDet/sum(nDet)) %>%
  select(-nDet) %>%
  pivot_wider(names_from=source, values_from=pDet,
              values_fill=0) %>%
  mutate(prop_lab=paste0("p: ", round(p*100), "%", 
                         "\ns: ", round(s*100), "%"),
         model="Joint")

gen.llam.loess <- agg$llam %>% 
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("L05", "median", "L95"), ~log(sum(exp(.))))) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

gen.lam.loess <- agg$lam %>% 
  mutate(genName=tax_i$FullGen[match(sppName, tax_i$species)]) %>%
  group_by(model, plot, el, genName) %>%
  summarise(across(one_of("L05", "median", "L95"), sum)) %>%
  group_by(model, genName) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))

p <- ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/L_llambda_opt_el_Genus.png", p, width=15, height=6)

p <- ggplot(gen.llam.loess, aes(el)) + 
  geom_point(aes(y=exp(median), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot (log-axis)") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_llambda_opt_el_Genus_pts.png", p, width=15, height=6)

p <- ggplot(gen.lam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_lambda_opt_el_Genus.png", p, width=15, height=6)

p <- ggplot(gen.lam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median, colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_pts.png", p, width=15, height=6)

p <- ggplot(gen.llam.loess, aes(el, group=model, colour=model, fill=model)) + 
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
ggsave("eda/L_llambda_opt_el_Genus_free.png", p, width=15, height=6)

p <- ggplot(gen.llam.loess, aes(el)) + 
  geom_point(aes(y=exp(median), colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  scale_y_log10() +
  labs(x="Elevation (m)", y="Colonies per plot  (log-axis)") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_llambda_opt_el_Genus_free_pts.png", width=15, height=6)

p <- ggplot(gen.lam.loess, aes(el, group=model, colour=model, fill=model)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), 
              alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  geom_text(data=det.gen, aes(x=-Inf, y=Inf, label=prop_lab),
            colour=1, size=2, hjust=-0.1, vjust=1.2) +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_free.png", p, width=15, height=6)

p <- ggplot(gen.lam.loess, aes(el)) + 
  geom_hline(yintercept=0, size=0.25, colour="gray30") +
  geom_point(aes(y=median, colour=model), alpha=0.2, shape=1, size=0.3) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Colonies per plot") +
  facet_wrap(~genName, scales="free_y", ncol=7)
ggsave("eda/L_lambda_opt_el_Genus_free_pts.png", p, width=15, height=6)










########------------------------------------------------------------------------
## RICHNESS
########------------------------------------------------------------------------

## Regional ----------------------------
S.loess <- agg$pP_R %>%
  filter(id %in% d.i[[1]]$X[,1]) %>%
  group_by(model, id, el) %>%
  summarise(across(one_of("L05", "median", "L95"), ~sum(.>0.95))) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x, span=0.7)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x, span=0.7)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x, span=0.7)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(S.loess, aes(x=el, colour=model, fill=model)) + 
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), alpha=0.25, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  ylim(0, NA) +
  labs(x="Elevation (m)", y="Richness")

S_R.loess <- agg$S_R %>%
  # filter(id %in% d.i[[1]]$X[,1]) %>%
  # filter(id > 1e5) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L025 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L975 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(S_R.loess, aes(x=el, colour=model, fill=model)) + 
  # geom_point(aes(y=median), shape=1, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  ylim(0, NA) +
  labs(x="Elevation (m)", y="Richness")
ggsave("eda/R_S_opt.png", width=6, height=5)
agg$S_R %>%
  ggplot(aes(el, median, ymin=L025, ymax=L975, colour=model)) + 
  geom_errorbar(size=0.25, width=10) + 
  geom_point(shape=1, size=1) +
  ylim(0,80) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2")

agg$pP_R %>% group_by(model, id, el) %>% 
  summarise(S=sum(L025>0.95)) %>% 
  ggplot(aes(el, S, colour=model)) + 
  geom_point(alpha=0.5) + 
  scale_colour_manual(values=mod_col) + ylim(0, 80)

agg$pP_R %>% group_by(model, id, el) %>% 
  summarise(across(one_of("L05", "median", "L95"), ~sum(.>0.95))) %>%
  ggplot(aes(el, median, ymin=L05, ymax=L95, colour=model)) + 
  geom_linerange(alpha=0.5, size=0.25) +
  geom_point(alpha=0.5, shape=1) + 
  scale_colour_manual(values=mod_col) + ylim(0, 80)


agg$LAM %>% #filter(id>1e5) %>% 
  group_by(site, model, el) %>% 
  summarise(S_lo=sum(1-exp(-L025)>0.95),
            S_md=sum(1-exp(-median)>0.95),
            S_hi=sum(1-exp(-L975)>0.95)) %>%
  ggplot(aes(el, group=model, colour=model)) + 
  stat_smooth(aes(y=S_lo), se=F, size=0.5) + 
  stat_smooth(aes(y=S_md), se=F) + 
  stat_smooth(aes(y=S_hi), se=F, size=0.5) + 
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2")

agg$pP_R %>% #filter(id>1e5) %>% 
  group_by(id, model, el) %>% 
  summarise(across(one_of("L025", "L05", "L10", "L25", "median", "mean",
                          "L75", "L90", "L95", "L975"), ~sum(.>0.95))) %>%
  ggplot(aes(el, group=model, colour=model)) + 
  stat_smooth(aes(y=L025), size=0.1, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L05), size=0.2, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L10), size=0.5, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L25), size=0.75, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=median), size=1, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=mean), size=1, se=F, linetype=2, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L75), size=0.75, se=F, method="loess", formula=y~x) +  
  stat_smooth(aes(y=L90), size=0.5, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L95), size=0.2, se=F, method="loess", formula=y~x) + 
  stat_smooth(aes(y=L975), size=0.1, se=F, method="loess", formula=y~x) + 
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2") + facet_wrap(~model)

agg$S_R %>%
  ggplot(aes(el, L025, colour=model)) + 
  geom_point(shape=1, size=1) +
  ylim(0,80) + scale_colour_manual(values=mod_col) + 
  labs("Elevation (m)", "Predicted richness per km2")



## Local -------------------------------

S_L.loess <- agg$S_L %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L05 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L95 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(S_L.loess, aes(x=el, colour=model)) + 
  geom_linerange(data=site.mns, aes(ymin=mnS-seS, ymax=mnS+seS), 
                 size=0.25, colour="black") +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi, fill=model), alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  geom_point(data=site.mns, aes(y=mnS), colour="black") +
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="Predicted richness per 0.75m2")
ggsave("eda/L_S_opt.png", width=6, height=5)
agg$S_L %>% 
  ggplot(aes(el, colour=model)) + ylim(0, NA) +
  geom_point(data=site.mns, aes(y=mnS), colour="black") +
  geom_linerange(data=site.mns,
                 aes(ymin=mnS-seS, ymax=mnS+seS), size=0.25, colour="black") +
  stat_smooth(aes(y=mean), size=0.5, se=F) + 
  stat_smooth(aes(y=L05), size=0.25, se=F) + 
  stat_smooth(aes(y=L95), size=0.25, se=F) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~model) + labs(x="Elevation (m)", y="Predicted richness per 0.75m2")








########------------------------------------------------------------------------
## DIVERSITY
########------------------------------------------------------------------------

H_R.loess <- agg$H %>%
  # filter(id %in% d.i[[1]]$X[,1]) %>%
  # filter(id > 1e5) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L025 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L975 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(H_R.loess, aes(x=el, colour=model, fill=model)) + 
  # geom_point(aes(y=median), shape=1, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  # ylim(0, NA) +
  labs(x="Elevation (m)", y="Predicted diversity per km2")
ggsave("eda/R_H_opt.png", width=6, height=5)
H_L.loess <- agg$H_L %>%
  # filter(id %in% d.i[[1]]$X[,1]) %>%
  # filter(id > 1e5) %>%
  group_by(model) %>%
  nest() %>%
  mutate(loess_lo=map(data, ~loess(L025 ~ el, data=.x)[['fitted']]),
         loess_md=map(data, ~loess(median ~ el, data=.x)[['fitted']]),
         loess_hi=map(data, ~loess(L975 ~ el, data=.x)[['fitted']])) %>%
  unnest(cols=c(data, loess_lo, loess_md, loess_hi))
ggplot(H_L.loess, aes(x=el, colour=model, fill=model)) + 
  # geom_point(aes(y=median), shape=1, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin=loess_lo, ymax=loess_hi), alpha=0.5, colour=NA) +
  geom_line(aes(y=loess_md)) + 
  scale_colour_manual(values=mod_col) + 
  scale_fill_manual(values=mod_col) + 
  # ylim(0, NA) +
  labs(x="Elevation (m)", y="Predicted diversity per 0.75 m2")
ggsave("eda/L_H_opt.png", width=6, height=5)







########------------------------------------------------------------------------
## PROBABILITY OF PRESENCE
########------------------------------------------------------------------------

## Regional ----------------------------

agg$lLAM %>% mutate(pP_R=1-exp(-exp(median))) %>%
  filter(spp %in% det_Y) %>%
  ggplot(aes(x=el, y=pP_R, colour=model)) + ylim(0,1) +
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="pr(Presence): site") +
  facet_wrap(~sppName)

agg$pP_R %>% 
  # filter(spp %in% det_Y) %>%
  ggplot(aes(x=el, y=L025, colour=model)) + ylim(0,1) +
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual(values=mod_col) + 
  labs(x="Elevation (m)", y="pr(Presence): site") +
  facet_wrap(~sppName)

agg$pP_R %>% filter(id>1e5) %>%
  ggplot(aes(x=el, colour=model)) + 
  # geom_point(alpha=0.5, size=0.5, shape=1) + 
  stat_smooth(aes(y=L95), method="loess", se=F, size=0.1) + 
  stat_smooth(aes(y=L75), method="loess", se=F, size=0.25) + 
  stat_smooth(aes(y=median), method="loess", se=F, size=0.5) + 
  stat_smooth(aes(y=L25), method="loess", se=F, size=0.25) + 
  stat_smooth(aes(y=L05), method="loess", se=F, size=0.1) + 
  ylim(0, 1) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~sppName)



## Local -------------------------------

agg$pP_L %>% 
  ggplot(aes(x=el, y=mean, colour=model)) + 
  geom_point(alpha=0.3, size=0.5, shape=1) + 
  scale_colour_manual("", values=mod_col) +
  facet_wrap(~sppName, scales="free_y") + 
  labs(x="Elevation (m)", y="pr(Presence): plot")














########------------------------------------------------------------------------
## TAXONOMIC BIAS
########------------------------------------------------------------------------

agg$D %>%
  ggplot(aes(sppName, median, ymin=L05, ymax=L95, 
             colour=sign(L05-1)==sign(L95-1))) + 
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
  mutate(sf1=tax_i$FullSF[match(gen1Name, tax_i$FullGen)],
         sf2=tax_i$FullSF[match(gen2Name, tax_i$FullGen)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(x=gen1Name, y=gen2Name, colour=mean)) + 
  geom_point(size=3) + facet_grid(model~ParName) + 
  scale_colour_gradient2(midpoint=0, low=scales::muted("blue"), high=scales::muted("red")) + 
  theme(panel.grid=element_line(colour="gray", size=0.25),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

agg$Sig_B %>% filter(gen1 != gen2) %>%
  mutate(sf1=tax_i$FullSF[match(gen1Name, tax_i$FullGen)],
         sf2=tax_i$FullSF[match(gen2Name, tax_i$FullGen)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(ParName, mean, fill=sf2)) + 
  geom_hline(yintercept=0) + 
  geom_boxplot() + facet_grid(model~sf1)

agg$Sig_B %>% filter(gen1 != gen2) %>%
  mutate(sf1=tax_i$FullSF[match(gen1Name, tax_i$FullGen)],
         sf2=tax_i$FullSF[match(gen2Name, tax_i$FullGen)]) %>%
  arrange(sf1, sf2, gen1Name, gen2Name) %>%
  mutate(gen1Name=factor(gen1Name, levels=unique(gen1Name)),
         gen2Name=factor(gen2Name, levels=unique(gen2Name))) %>%
  ggplot(aes(mean, ParName, fill=sf2)) + 
  geom_vline(xintercept=0) +
  ggridges::geom_density_ridges(alpha=0.5) +
  facet_grid(model~sf1)

agg$Sig_B %>% filter(gen1 != gen2) %>%
  ggplot(aes(x=abs(mean), y=ParName, fill=model)) + 
  geom_vline(xintercept=0) + 
  ggridges::geom_density_ridges(alpha=0.5) +
  scale_fill_manual("", values=mod_col)


## Within genera -----------------------

agg$sig_b %>%
  mutate(Scale=str_sub(ParName, 1L, 1L),
         Par=str_sub(ParName, 3L, -1L)) %>%
  mutate(Scale=factor(Scale, levels=c("R", "L"), 
                      labels=c("Regional", "Local"))) %>%
  ggplot(aes(ParName, median, ymin=L05, ymax=L95, colour=model)) + 
  geom_hline(yintercept=0, linetype=3, colour="gray", size=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(position=position_dodge(width=0.5)) + 
  scale_colour_manual("", values=mod_col) +
  coord_flip() + labs(x="Standard deviation among congeners", y="")










########------------------------------------------------------------------------
## ASSEMBLAGES
########------------------------------------------------------------------------

# I'm not sure this makes any sense.... It's modelled and the species aren't independent
md_LAM <- agg$LAM %>% select(model, median, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=median)
lo_LAM <- agg$LAM %>% select(model, L025, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=L025)
hi_LAM <- agg$LAM %>% select(model, L975, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=L975)

md_lam <- agg$lam %>% select(model, median, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=median)
lo_lam <- agg$lam %>% select(model, L025, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=L025)
hi_lam <- agg$lam %>% select(model, L975, sppName, id, el) %>%
  pivot_wider(names_from=sppName, values_from=L975)

library(vegan)

pca.lam.joint <- filter(lo_LAM, model=="Joint") %>% 
  select(-model, -id, -el) %>% rda() 
pca.lam.str <- filter(lo_LAM, model=="Structured") %>% 
  select(-model, -id, -el) %>% rda() 

biplot(pca.lam.joint, type=c("text", "points"), cex=1.5)
biplot(pca.lam.str, type=c("text", "points"), cex=1.5)









agg$pP_L$obs <- rep(c(t(d.ls[[1]]$Y)), n_distinct(agg$pP_L$model))
agg$lam$obs <- agg$llam$obs <- rep(c(t(d.ls[[1]]$Y)), n_distinct(agg$pP_L$model))
agg$pP_L$disp <- agg$disp$mean[match(agg$pP_L$model, agg$disp$model)]


ggplot(agg$pP_L, aes(mean, obs)) + geom_point(alpha=0.1) + 
  stat_smooth() + facet_wrap(~model)
ggplot(agg$pP_L, aes(mean, as.numeric(obs>0), group=model, colour=model)) + 
  #geom_point(alpha=0.1) + 
  stat_smooth(aes(x=L05), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=median), method="glm", size=0.5, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  stat_smooth(aes(x=L95), method="glm", size=0.5, linetype=3, se=F,
              method.args=list(family="binomial"), fullrange=T) +
  # facet_wrap(~model) +
  ylim(0,1)
ggplot(agg$pP_L, aes(factor(obs), mean, fill=model)) + geom_boxplot() +
  labs(x="Number of detections at plot", y="prPresence")
ggplot(agg$pP_L, aes(mean, factor(obs), fill=model)) + 
  geom_vline(xintercept=0, colour="gray") + 
  ggridges::geom_density_ridges(alpha=0.5) + 
  scale_fill_manual(values=mod_col) +
  labs(y="Number of detections at plot", x="prPresence")
ggplot(agg$lam, aes(mean, factor(obs), fill=model)) + 
  geom_vline(xintercept=0, colour="gray") + 
  ggridges::geom_density_ridges(alpha=0.5) + 
  scale_fill_manual(values=mod_col) +
  labs(y="Number of detections at plot", x=expression(lambda))
ggplot(agg$llam, aes(factor(obs>0), mean, fill=model)) + geom_boxplot() +
  labs(x="Number of detections at plot", y=expression(log(lambda)))


agg$pP_L <- agg$pP_L %>%
  mutate(ll=LaplacesDemon::dgpois(obs, mean, disp, log=T))
agg$pP_L %>% group_by(model) %>% 
  summarise(nll=-sum(ll)*2) %>% 
  arrange(nll)




agg$pP_R  %>% group_by(site, model, el) %>% 
  filter(id>1e5) %>%
  summarise(S_mn=sum(mean > 0.95)) %>%
  ggplot(aes(el)) + ylim(0, 80) + geom_point(aes(y=S_mn)) + 
  facet_wrap(~model) + labs("Elevation (m)", "Predicted richness per km2")

agg$S_R %>% filter(id>1e5) %>%
  ggplot(aes(el, L025)) + geom_point() + ylim(0,80)


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
                           "id"=d.i[[1]]$X[,"id"], 
                           "spp"=tax_i$species)))
lLAMBDA.ls <- map(lLAMBDA.ls, ~.[,which(d.i[[1]]$X[,"id"] > 1e5),])
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
X_all <- rbind(d.i[[1]]$X, d.i[[1]]$X_) %>% as.data.frame
lLAM.site <- agg$tot_LAM %>% filter(id<1e5) %>% 
  group_by(id, model) %>% 
  summarise(N_lo=sum(L05),
            N_md=sum(median),
            N_mn=sum(mean),
            N_hi=sum(L95)) %>%
  mutate(model=case_when(model=="Joint" ~ "Joint",
                         model=="Structured" ~ "Str")) %>%
  pivot_wider(names_from="model", values_from=3:6)
S.site <- agg$S_R %>% filter(id<1e5) %>% 
  group_by(id, model) %>% 
  summarise(S_lo=sum(L05),
            S_md=sum(median),
            S_mn=sum(mean),
            S_hi=sum(L95)) %>%
  mutate(model=case_when(model=="Joint" ~ "Joint",
                         model=="Structured" ~ "Str")) %>%
  pivot_wider(names_from="model", values_from=3:6)
out.sf <- d.i[[1]]$grd_W.sf %>% 
  left_join(., select(X_all, -el), by="id") %>%
  # filter(R_rdLen > 5) %>%
  left_join(., lLAM.site, by="id") %>%
  left_join(., S.site) 
  
ggplot(out.sf, aes(fill=(S_hi_Str-S_lo_Str)/S_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[S[Str]]), option="E", limits=c(0, 0.4))
ggplot(out.sf, aes(fill=(S_hi_Joint-S_lo_Joint)/S_md_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[S[Joint]]), option="E", limits=c(0, 0.4))

ggplot(out.sf, aes(fill=(N_hi_Str-N_lo_Str)/N_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[N[Str]]), option="B", limits=c(0, 3))
ggplot(out.sf, aes(fill=(N_hi_Joint-N_lo_Joint)/N_md_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(CV[N[Joint]]), option="B", limits=c(0, 3))

ggplot(out.sf, aes(fill=S_md_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(S[Joint]), limits=c(0, 80))

ggplot(out.sf, aes(fill=S_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(S[Str]), limits=c(0, 80))

ggplot(out.sf, aes(fill=S_md_Joint - S_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_gradient2()

ggplot(out.sf, aes(fill=N_md_Joint)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(N[Joint]), option="B", limits=c(0, 1800000))

ggplot(out.sf, aes(fill=N_md_Str)) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_viridis(expression(N[Str]), option="B", limits=c(0, 1800000))

ggplot(out.sf, aes(fill=log(N_md_Joint) - log(N_md_Str))) + 
  geom_sf(colour="black", size=0.1) + 
  scale_fill_gradient2()


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
  
  sp.pP <- sp.PA <- sp.LamCV <- vector("list", n_distinct(tax_i$species))
  sp.N_mn <- sp.N_L95 <- sp.lN_mn <- sp.lN_L95 <- sp.pP
  gen.N <- gen.lN <- gen.S <- vector("list", n_distinct(tax_i$FullGen))
  sf.N <- sf.lN <- sf.S <- vector("list", n_distinct(tax_i$FullSF))
  for(i in seq_along(tax_i$species)) {
    sp.i <- tax_i$species[i]
    sp.f <- str_replace(sp.i, "/", "_")
    LAM.i <- filter(agg$LAM, sppName==sp.i & model==m.full)
    lLAM.i <- filter(agg$lLAM, sppName==sp.i & model==m.full)
    pP.i <- filter(agg$pP_R, sppName==sp.i & model==m.full)
    obs.i <- filter(ants$all, SPECIESID==sp.i) 
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    
    sp.N_mn[[i]] <- out.sf %>%
      mutate(N_sp=LAM.i$median[match(id, LAM.i$id)]) %>%
      filter(!is.na(N_sp)) %>%
      ggplot() + geom_sf(aes(fill=N_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(Lambda)), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_LAM_mn_", m.abb, ".jpg"), 
           sp.N_mn[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.N_L95[[i]] <- out.sf %>%
      mutate(N_sp=LAM.i$L025[match(id, LAM.i$id)]) %>%
      filter(!is.na(N_sp)) %>%
      ggplot() + geom_sf(aes(fill=N_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(Lambda[L95]), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_LAM_L95_", m.abb, ".jpg"), 
           sp.N_L95[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.lN_mn[[i]] <- out.sf %>%
      mutate(N_sp=lLAM.i$median[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_sp)) %>%
      ggplot() + geom_sf(aes(fill=N_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(log(Lambda))), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_lLAM_mn_", m.abb, ".jpg"), 
           sp.lN_mn[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.lN_L95[[i]] <- out.sf %>%
      mutate(N_sp=lLAM.i$L025[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_sp)) %>%
      ggplot() + geom_sf(aes(fill=N_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(log(Lambda)[L95]), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_lLAM_L95_", m.abb, ".jpg"), 
           sp.lN_L95[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.pP[[i]] <- out.sf %>%
      mutate(P_sp=pP.i$L025[match(id, pP.i$id)]) %>%
      filter(!is.na(P_sp)) %>%
      ggplot() + geom_sf(aes(fill=P_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) +
      scale_fill_viridis(expression(pr(Pres)[L95]), option="B", limits=c(0, 1)) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_pP_", m.abb, ".jpg"), 
           sp.pP[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.PA[[i]] <- out.sf %>%
      mutate(P_sp=pP.i$L025[match(id, pP.i$id)]>0.99) %>%
      filter(!is.na(P_sp)) %>%
      ggplot() + geom_sf(aes(fill=P_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_manual(expression(pr(Pres)[L95]~'> 99%'), 
                        values=c("FALSE"="gray70", "TRUE"="#2171b5")) + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_PA_", m.abb, ".jpg"), 
           sp.PA[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sp.LamCV[[i]] <- out.sf %>%
      mutate(CV_sp=LAM.i$sd[match(id, LAM.i$id)]/
                         LAM.i$mean[match(id, LAM.i$id)]) %>%
      filter(!is.na(CV_sp)) %>%
      ggplot() + geom_sf(aes(fill=CV_sp), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sp.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(CV[Lambda]), option="E") + 
      map_theme
    ggsave(paste0("eda/maps/spp/", sp.f, "_LamCV_", m.abb, ".jpg"), 
           sp.LamCV[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  for(i in 1:n_distinct(tax_i$FullGen)) {
    gen.i <- unique(tax_i$FullGen)[i]
    gen.sppNum <- which(tax_i$FullGen == gen.i)
    obs.i <- filter(ants$all, SPECIESID %in% tax_i$species[gen.sppNum])
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    LAM.i <- filter(agg$LAM, spp %in% gen.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=sum(median))
    lLAM.i <- filter(agg$lLAM, spp %in% gen.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=log(sum(exp(median))))
    S.i <- filter(agg$pP_R, spp %in% gen.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=sum(L025>0.95))
    
    gen.S[[i]] <- out.sf %>%
      mutate(S_gen=S.i$mean[match(id, S.i$id)]) %>%
      filter(!is.na(S_gen)) %>%
      ggplot() + geom_sf(aes(fill=S_gen), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", gen.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis("Pred.\nRich.", limits=c(0, length(gen.sppNum))) + 
      map_theme
    ggsave(paste0("eda/maps/gen/", gen.i, "_S_", m.abb, ".jpg"), 
           gen.S[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    gen.N[[i]] <- out.sf %>%
      mutate(N_gen=LAM.i$mean[match(id, LAM.i$id)]) %>%
      filter(!is.na(N_gen)) %>%
      ggplot() + geom_sf(aes(fill=N_gen), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", gen.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(Lambda)), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/gen/", gen.i, "_LAM_", m.abb, ".jpg"),
           gen.N[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    gen.lN[[i]] <- out.sf %>%
      mutate(N_gen=lLAM.i$mean[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_gen)) %>%
      ggplot() + geom_sf(aes(fill=N_gen), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", gen.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(log(Lambda))), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/gen/", gen.i, "_lLAM_", m.abb, ".jpg"),
           gen.lN[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  for(i in 1:n_distinct(tax_i$FullSF)) {
    sf.i <- unique(tax_i$FullSF)[i]
    sf.sppNum <- which(tax_i$FullSF == sf.i)
    obs.i <- filter(ants$all, SPECIESID %in% tax_i$species[sf.sppNum])
    if(m.abb=="Y") obs.i <- filter(obs.i, source=="s")
    LAM.i <- filter(agg$LAM, spp %in% sf.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=sum(median))
    lLAM.i <- filter(agg$lLAM, spp %in% sf.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=log(sum(exp(median))))
    S.i <- filter(agg$pP_R, spp %in% sf.sppNum & model==m.full) %>% 
      group_by(id) %>%
      summarise(mean=sum(L025>0.95))
    
    sf.S[[i]] <- out.sf %>%
      mutate(S_sf=S.i$mean[match(id, S.i$id)]) %>%
      filter(!is.na(S_sf)) %>%
      ggplot() + geom_sf(aes(fill=S_sf), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sf.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis("Pred.\nRich.", limits=c(0, length(sf.sppNum))) + 
      map_theme
    ggsave(paste0("eda/maps/sf/", sf.i, "_S_", m.abb, ".jpg"), 
           sf.S[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sf.N[[i]] <- out.sf %>%
      mutate(N_sf=LAM.i$mean[match(id, LAM.i$id)]) %>%
      filter(!is.na(N_sf)) %>%
      ggplot() + geom_sf(aes(fill=N_sf), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sf.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(Lambda)), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/sf/", sf.i, "_LAM_", m.abb, ".jpg"), 
           sf.N[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
    
    sf.lN[[i]] <- out.sf %>%
      mutate(N_sf=lLAM.i$mean[match(id, lLAM.i$id)]) %>%
      filter(!is.na(N_sf)) %>%
      ggplot() + geom_sf(aes(fill=N_sf), colour="black", size=0.1) + 
      annotate("text", label=paste0(m.full, "\n", sf.i), size=3,
               x=495000, y=200000, hjust=0, vjust=0.5) + 
      scale_fill_viridis(expression(hat(log(Lambda))), option="B", limits=c(0, NA)) + 
      map_theme
    ggsave(paste0("eda/maps/sf/", sf.i, "_lLAM_", m.abb, ".jpg"), 
           sf.lN[[i]] + 
             geom_sf(data=obs.i, colour="deepskyblue1", shape=1, size=0.5), 
           width=4, height=4, units="in", dpi=300)
  }
  
  
  ggpubr::ggarrange(plotlist=sp.N_mn, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_N_mn_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.N_L95, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_N_L95_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.lN_mn, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_lN_mn_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.lN_L95, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_lN_L95_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.pP, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_pP_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.PA, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_PA_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sp.LamCV, ncol=10, nrow=8) %>%
    ggsave(paste0("eda/maps/spp_LamCV_", m.abb, ".jpg"), ., 
           width=30, height=24, units="in", dpi=300)
  
  ggpubr::ggarrange(plotlist=gen.N, ncol=7, nrow=3) %>%
    ggsave(paste0("eda/maps/gen_N_", m.abb, ".jpg"), ., 
           width=21, height=9, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=gen.lN, ncol=7, nrow=3) %>%
    ggsave(paste0("eda/maps/gen_lN_", m.abb, ".jpg"), ., 
           width=21, height=9, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=gen.S, ncol=7, nrow=3) %>%
    ggsave(paste0("eda/maps/gen_S_", m.abb, ".jpg"), ., 
           width=21, height=9, units="in", dpi=300)
  
  ggpubr::ggarrange(plotlist=sf.N, ncol=2, nrow=2) %>%
    ggsave(paste0("eda/maps/sf_N_", m.abb, ".jpg"), ., 
           width=6, height=6, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sf.lN, ncol=2, nrow=2) %>%
    ggsave(paste0("eda/maps/sf_lN_", m.abb, ".jpg"), ., 
           width=6, height=6, units="in", dpi=300)
  ggpubr::ggarrange(plotlist=sf.S, ncol=2, nrow=2) %>%
    ggsave(paste0("eda/maps/sf_S_", m.abb, ".jpg"), ., 
           width=6, height=6, units="in", dpi=300)
  
}



agg$lLAM %>% select(sppName, id, el, model, median, L05, L95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("median", "L05", "L95")) %>%
  ggplot(aes(el, median_WY-median_Y)) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_hline(yintercept=0) + facet_wrap(~sppName)


agg$lLAM %>% select(sppName, id, el, model, mean, L05, L95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mean", "L05", "L95")) %>%
  ggplot(aes(mean_WY, mean_Y)) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_abline() + facet_wrap(~sppName)

agg$llam %>% select(sppName, id, el, model, mean, L05, L95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mean", "L05", "L95")) %>%
  ggplot(aes(mn_WY, mn_Y)) + 
  geom_point(alpha=0.5, size=0.5, shape=1) + 
  geom_abline() + facet_wrap(~sppName, scales="free")

agg$lLAM %>% select(sppName, id, el, model, mean, L05, L95) %>% 
  mutate(model=case_when(model=="Joint" ~ "WY",
                         model!="Joint" ~ "Y")) %>%
  pivot_wider(names_from="model", values_from=c("mean", "L05", "L95")) %>%
  ggplot(aes(el, (q95_WY-q05_WY)-(q95_Y-q05_Y))) + 
  geom_point(alpha=0.1, size=0.5, shape=1) + 
  geom_hline(yintercept=0) + facet_wrap(~sppName)





ggplot(out.sf, aes(R_rdLen, N_mn_Joint)) + geom_point(alpha=0.5)



bdm_27 <- map_dfr(1:80, ~tibble(llam=filter(agg$llam, spp==.x & 
                                              plot %in% which(d.i[[1]]$IJ$Y==1) &
                                              model=="Joint")$mean,
                                lLAM=filter(agg$lLAM, spp==.x & 
                                              id==527158 & 
                                              model=="Joint")$mean,
                                spp=.x))

ggplot(bdm_27, aes(exp(lLAM), exp(llam))) + geom_point()

exp(filter(agg$lLAM, spp==1 & id==527158 & model=="Joint")$mean)/
  exp(filter(agg$llam, spp==1 & plot %in% which(d.i[[1]]$IJ$Y==1) & model=="Joint")$mean)



