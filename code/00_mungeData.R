# Script for munging datasets into appropriate forms
# Tim Szewczyk



##--- set up
library(tidyverse); library(sf); library(googlesheets)
gis.dir <- "../2_gis/data/VD_21781/"
ant.dir <- "../1_opfo/data/"
tax_i <- read_csv("data/tax_i.csv") #%>% filter(sNum<50)
veg_i <- read_csv(paste0(ant.dir, "vegcover_id.csv")) # van der Maarel 2007
lc_i <- readxl::read_xlsx(paste0(ant.dir, "landcover_id.xlsx"), 1)
plot_i <- read_csv(paste0(ant.dir, "opfo_envDataProcessed.csv")) %>% 
  arrange(BDM, Plot_id) %>% filter(Plot_id != "020201") %>%
  group_by(BDM) %>% 
  mutate(SoilTSt=(SoilTemp-mean(SoilTemp,na.rm=T))/sd(SoilTemp,na.rm=T),
         CnpyOpn=(lc_i$Canopy[match(Categorie, lc_i$LC)]=="Open")*1,
         CnpyMxd=(lc_i$Canopy[match(Categorie, lc_i$LC)]=="Mixed")*1) %>%
  mutate_at(c("Grass", "Forb", "Shrub", "Bare", "Litter", "Moss"), 
            ~veg_i$Pct[match(., veg_i$class)]) %>%
  mutate(VegTot=Grass + Forb + Shrub)
source("code/00_fn.R"); source(paste0(ant.dir, "../code/00_fn.R"))



##--- settings
test_prop_W <- 0.1 # proportion of W data to use for testing models
test_prop_Y <- 0.1 # proportion of Y data to use for testing models
X_vars <- c("MAT", "Iso", "AP", "lcH", "npp", "MAT_sq", "AP_sq")
V_vars <- c("SoilTSt", "CnpyOpn", "CnpyMxd", "VegTot")
U_vars <- c("pop", "rdLen")



##--- grids
# W: 1km2 grid of Vaud
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp")
VD <- st_union(VD_raw)
VD_ext <- raster::extent(matrix(st_bbox(VD), ncol=2))
grd_W <- raster::raster(ext=VD_ext, crs=st_crs(VD), resolution=1000) %>%
  raster::rasterize(VD_ext, .) %>% 
  raster::rasterToPolygons(., n=4)
grd_W@data$layer <- 1:raster::ncell(grd_W)
grd_W.sf <- st_as_sf(grd_W) %>% st_set_crs(st_crs(VD)) %>% rename(id=layer) %>%
  mutate(inbd=c(st_intersects(., VD, sparse=F)))
# Y: 44 x 1km2 sampling sites
site.sf <- agg_str_site_data() %>% arrange(BDM)



##--- testing / training
grd_W.sf <- grd_W.sf %>% 
  mutate(K_orig=row_number(), K_rand=sample(1:n())) %>%
  arrange(K_rand)
site.sf <- site.sf %>%
  mutate(J_orig=row_number(), J_rand=sample(1:n())) %>%
  arrange(J_rand)
plot_i <- plot_i %>% ungroup %>%
  mutate(J_orig=site.sf$J_orig[match(plot_i$BDM, site.sf$BDM)],
         J_rand=site.sf$J_rand[match(plot_i$BDM, site.sf$BDM)]) %>%
  arrange(J_rand, Plot_id) %>%
  mutate(I_rand=row_number())



##--- W
ants <- load_ant_data(str_type="soil", clean_spp=T)
grid.W <- st_join(ants$pub, grd_W.sf) %>%
  group_by(SPECIESID, id, inbd, K_orig, K_rand) %>% st_set_geometry(NULL) %>%
  summarise(nObs=n()) %>% ungroup %>%
  filter(inbd) %>%
  arrange(K_rand)
W <- matrix(0, nrow=nrow(grd_W.sf), ncol=max(tax_i$sNum))
for(s in 1:max(tax_i$sNum)) {
  sp.occ <- filter(grid.W, SPECIESID==tax_i$species[s])
  if(nrow(sp.occ)>0) W[sp.occ$K_rand,s] <- sp.occ$nObs
}
W_id <- which(rowSums(W)>0)



##--- Y
box.Y <- ants$str %>% group_by(BDM, Plot_id, SPECIESID) %>%
  summarise(nObs=n()) %>% 
  left_join(., select(plot_i, BDM, Plot_id, J_orig, J_rand, I_rand), 
            by=c("BDM", "Plot_id")) %>%
  arrange(I_rand)
# Y: ordered by plot_i$Plot_id = arrange(BDM, Plot_id)
Y <- matrix(0, nrow=sum(site.sf$nPlot_Total), ncol=max(tax_i$sNum))
for(s in 1:max(tax_i$sNum)) {
  sp.occ <- filter(box.Y, SPECIESID==tax_i$species[s])
  if(nrow(sp.occ)>0) Y[which(plot_i$Plot_id %in% sp.occ$Plot_id),s] <- sp.occ$nObs
}



##--- indexes
K <- list(W=floor(length(W_id)*(1-test_prop_W)), 
          W_=ceiling(length(W_id)*test_prop_W))
J <- list(Y=floor(nrow(site.sf)*(1-test_prop_Y)), 
          Y_=ceiling(nrow(site.sf)*test_prop_Y))
IJ <- list(Y=rep(1:J$Y, times=site.sf$nPlot_Total[1:J$Y]), 
           Y_=rep(1:J$Y_, times=site.sf$nPlot_Total[(1:J$Y_)+J$Y]))
I <- list(Y=length(IJ$Y),
          Y_=length(IJ$Y_))


# ##--- testing / training randomization
# K_ord <- tibble(orig=1:length(W_id), rand=1:length(W_id)) %>%#rand=sample(orig)) %>%
#   arrange(rand) %>% mutate(train=row_number()<=K$W)
# J_ord <- tibble(orig=1:nrow(site.sf), rand=1:nrow(site.sf)) %>%#rand=sample(orig)) %>%
#   arrange(rand) %>% mutate(train=row_number()<=J$Y)
# IJ_ord <- list(Y_ord=rep(J_ord$orig[J_ord$train], 
#                          times=site.sf$nPlot_Total[J_ord$orig[J_ord$train]]),
#                Y_ord_=rep(J_ord$orig[!J_ord$train], 
#                           times=site.sf$nPlot_Total[J_ord$orig[!J_ord$train]]),
#                Y=rep(1:J$Y, 
#                      times=site.sf$nPlot_Total[J_ord$orig[J_ord$train]]),
#                Y_=rep(1:J$Y_, 
#                       times=site.sf$nPlot_Total[J_ord$orig[!J_ord$train]]))
# I_ord <- list(Y=length(IJ_ord$Y),
#               Y_=length(IJ_ord$Y_))




##--- tax_i: temporary
if(is.null(tax_i)) {
  tax_i <- tibble(species=sort(unique(ants$all$SPECIESID)),
                  genus=str_split_fixed(species, "_", 2)[,1],
                  sNum=as.numeric(factor(species)),
                  gNum=as.numeric(factor(genus)),
                  Dprior=rbinom(length(species), 1, 0.3)+0.5)
  write_csv(tax_i, "data/tax_i.csv")
}



##--- covariates
clim <- dir(gis.dir, "chelsa") %>%
  setNames(str_remove(., "_chelsa_VD_21781.tif")) %>%
  magrittr::extract(names(.) %in% X_vars) %>%
  map(., ~raster::raster(paste0(gis.dir, .x)))
dem <- raster::raster(paste0(gis.dir, "dem_VD_21781.tif"))
slope <- raster::raster(paste0(gis.dir, "slope_VD_21781.tif"))
npp <- raster::raster(paste0(gis.dir, "MODIS_2010-2019_VD_21781.tif"))
lc.W <- st_read(paste0(gis.dir, "lc_1kmZones_21781.shp")) %>%
  filter(id %in% W_id) %>% st_set_geometry(NULL) %>% select(4:18) 
lc.Y <- st_read(paste0(gis.dir, "site_lcZones_21781.shp")) %>% 
  arrange(BDM) %>% st_set_geometry(NULL) %>% select(3:17)
pop <- raster::raster(paste0(gis.dir, "popDensity_VD_21781.tif")) 
roads <- st_read(paste0(gis.dir, "roads_VD_21781.shp")) %>% st_set_crs(21781)
# bldgs <- st_read(paste0(gis.dir, "buildings_VD_21781.shp")) %>% st_set_crs(21781)



##--- X
# !!! Need to store all covariate df's, including _sq terms (XW, XY, V, U)
# calculate means within W cells
# X_W.df <- filter(grd_W.sf, id %in% W_id) %>%
#   mutate(el=raster::extract(dem, ., fun=mean),
#          slope=raster::extract(slope, ., fun=mean),
#          npp=raster::extract(npp, ., fun=mean, na.rm=T),
#          lcH=vegan::diversity(as.matrix(lc.W))) %>%
#   bind_cols(., map_dfc(clim, ~raster::extract(., filter(grd_W.sf, id%in%W_id), 
#                                                 fun=mean))) %>%
#   arrange(id)
# write_sf(X_W.df, "data/cov/X_W-df.shp")
X_W.df <- st_read("data/cov/X_W-df.shp") %>%
  mutate(J_rand=match(id, grd_W.sf$id)) %>% 
  arrange(J_rand)
# select covariates and add quadratic terms
X_W.mx <- as.matrix(X_W.df %>% st_set_geometry(NULL) %>% 
                      select(-inbd) %>% 
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      select("J_rand", "el", one_of(X_vars)))
X_Y.mx <- as.matrix(site.sf %>% st_set_geometry(NULL) %>% 
                      select(-region, -BDM_id, -Sample_date,
                             -contains("nPlot")) %>%
                      mutate(lcH=vegan::diversity(as.matrix(lc.Y))) %>%
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      select("BDM", "el", one_of(X_vars)))
# scale and carefully combine
X.mx <- rbind(X_W.mx[1:K$W,], 
              X_Y.mx[1:J$Y,],
              X_W.mx[K$W+(1:K$W_),], 
              X_Y.mx[J$Y+(1:J$Y_),])
X.scale <- scale(X.mx[,-(1:2)])
X.all <- cbind(1, X.scale)




##--- V
# V.df <- st_read("../2_gis/data/opfo/opfo_soil_25per.shp") %>%
#   st_transform(st_crs(site.sf)) %>%
#   mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
#   select(Plot_id) %>%
#   left_join(., select(plot_i, BDM, Plot_id, all_of(V_vars)), by="Plot_id") %>%
#   mutate(el=raster::extract(dem, .),
#          slope=raster::extract(slope, .),
#          pop=replace_na(raster::extract(pop, .), 0)) %>%
#   arrange(BDM, Plot_id)
# write_sf(V.df, "data/cov/V-df.shp")
V.df <- st_read("data/cov/V-df.shp") %>% 
  mutate(J_rand=site.sf$J_rand[match(BDM, site.sf$BDM)],
         I_rand=plot_i$I_rand[match(Plot_id, plot_i$Plot_id)]) %>% 
  arrange(I_rand)
V_Y.mx <- as.matrix(V.df %>% st_set_geometry(NULL) %>% 
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      mutate(Plot_id=as.numeric(as.character(Plot_id))) %>%
                      select("Plot_id", "el", all_of(V_vars)))
V.mx <- rbind(V_Y.mx[1:I$Y,], 
              V_Y.mx[I$Y+(1:I$Y_),])
V.scale <- scale(V.mx[,-(1:2)])
if("CnpyOpn" %in% V_vars) {
  V.scale[,"CnpyOpn"] <- V.scale[,"CnpyOpn"] * 
    attr(V.scale, "scaled:scale")["CnpyOpn"] + 
    attr(V.scale, "scaled:center")["CnpyOpn"]
}
if("CnpyMxd" %in% V_vars) {
  V.scale[,"CnpyMxd"] <- V.scale[,"CnpyMxd"] * 
    attr(V.scale, "scaled:scale")["CnpyMxd"] + 
    attr(V.scale, "scaled:center")["CnpyMxd"]
}



##--- U
# # load covariates
# rdLen.grd_W <- roads %>%
#   st_intersection(., filter(grd_W.sf, id %in% W_id)) %>%
#   mutate(length=st_length(.)) %>%
#   group_by(id) %>%
#   summarise(rdLen=as.numeric(sum(length)))
# bldgs.grd_W <- bldgs %>%
#   st_intersection(., filter(grd_W.sf, id %in% W_id)) %>%
#   mutate(area=st_area(.)) %>%
#   group_by(id) %>%
#   summarise(bldgArea=as.numeric(sum(area)))
# 
# # calculate means within W cells
# U_W.df <- grd_W.sf %>% filter(id %in% W_id) %>%
#   mutate(slope=raster::extract(slope, ., fun=mean),
#          pop=raster::extract(pop, ., fun=mean, na.rm=T)) %>%
#   arrange(id)
# U_W.df <- left_join(U_W.df, st_set_geometry(rdLen.grd_W, NULL), by="id")
# U_W.df[is.na(U_W.df)] <- 0
# write_sf(U_W.df, "data/cov/U_W-df.shp")
U_W.df <- st_read("data/cov/U_W-df.shp") %>%
  mutate(J_rand=match(id, grd_W.sf$id)) %>% 
  arrange(J_rand)
U_W.mx <- as.matrix(st_set_geometry(U_W.df, NULL) %>% 
                      select("J_rand", all_of(U_vars)))
U.scale <- scale(U_W.mx[,-1])
U.all <- cbind(1, U.scale)



##--- NA's
# identify NA cells
na.W <- which(is.na(rowSums(X_W.mx)[1:K$W]))
na.W_ <- which(is.na(rowSums(X_W.mx)[K$W+(1:K$W_)]))
na.Y <- which(is.na(rowSums(V.mx)[1:I$Y]))
na.Y_ <- which(is.na(rowSums(V.mx)[I$Y+(1:I$Y_)]))



# d.ls <- read_rdump("data/stan_data/test_realData.Rdump")
# d.ls$V <- V.scale[1:I$Y,][-na.Y,]
# d.ls$V_ <- V.scale[I$Y+(1:I$Y_),][-na.Y_,]
# d.ls$L <- dim(V.scale)[2]

##--- export stan data
d.ls <- list(K=K$W-length(na.W), K_=K$W_-length(na.W_), 
             J=J$Y, J_=J$Y_, 
             IJ=IJ$Y[-na.Y], IJ_=IJ$Y_[-na.Y_], 
             I=I$Y-length(na.Y), I_=I$Y_-length(na.Y_),
             S=max(tax_i$sNum), 
             G=max(tax_i$gNum), 
             tax_i=tax_i[,c("sNum", "gNum")], 
             D_prior=tax_i$Dprior, 
             R=dim(X.all)[2], 
             L=dim(V.scale)[2],
             Q=dim(U.all)[2], 
             W=W[(1:K$W)[-na.W],], 
             W_=W[(K$W+(1:K$W_))[-na.W_],], 
             Y=Y[1:I$Y,][-na.Y,],
             Y_=Y[I$Y+(1:I$Y_),][-na.Y_,],
             X=X.all[(1:(K$W+J$Y))[-na.W],], 
             X_=X.all[((K$W+J$Y)+(1:(K$W_+J$Y_)))[-na.W_],], 
             V=V.scale[1:I$Y,][-na.Y,], 
             V_=V.scale[I$Y+(1:I$Y_),][-na.Y_,],
             U=U.all[(1:K$W)[-na.W],], 
             U_=U.all[(K$W+(1:K$W_))[-na.W_],], 
             h=7.5e-7)
d.i <- list(X=X.mx[(1:(K$W+J$Y))[-na.W],],
            X_=X.mx[((K$W+J$Y)+(1:(K$W_+J$Y_)))[-na.W_],],
            V=V.mx[1:I$Y,][-na.Y,],
            V_=V.mx[I$Y+(1:I$Y_),][-na.Y_,],
            U=U_W.mx[(1:K$W)[-na.W],],
            U_=U_W.mx[(K$W+(1:K$W_))[-na.W_],],
            tax_i=tax_i)
# d.ls <- list(K=K$W, K_=K$W_,
#              J=J$Y, J_=J$Y_,
#              IJ=IJ$Y, IJ_=IJ$Y_,
#              I=I$Y, I_=I$Y_,
#              S=max(tax_i$sNum),
#              G=max(tax_i$gNum),
#              tax_i=tax_i[,c("sNum", "gNum")],
#              D_prior=tax_i$Dprior,
#              R=dim(X.all)[2],
#              L=dim(V.scale)[2],
#              Q=dim(U.all)[2],
#              W=W[1:K$W,],
#              W_=W[K$W+(1:K$W_),],
#              Y=Y[1:I$Y,],
#              Y_=Y[I$Y+(1:I$Y_),],
#              X=X.all[1:(K$W+J$Y),],
#              X_=X.all[(K$W+J$Y)+(1:(K$W_+J$Y_)),],
#              V=V.scale[1:I$Y,],
#              V_=V.scale[I$Y+(1:I$Y_),],
#              U=U.all[1:K$W,],
#              U_=U.all[K$W+(1:K$W_),],
#              h=7.5e-7)
d.f <- "data/stan_data/test_realData"
saveRDS(d.ls, paste0(d.f, "_ls.rds"))
saveRDS(d.i, paste0(d.f, "_i.rds"))
rstan::stan_rdump(ls(d.ls),
                  file="data/stan_data/test_realData.Rdump",
                  envir=list2env(d.ls))



# Testing
library(rstan)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
mods <- "W" #c("W", "Y", "WY")
pars.exc <- c("a_std", "A_std", "b_std", "B_std", 
              "LAMBDA", "LAMBDA_", "lambda", "lambda_",
              "lLAMBDA", "lLAMBDA_", "llambda", "llambda_",
              "L_Omega_A", "L_Omega_B", "sigma_A", "sigma_B",
              "p", "p_")
out.ls <- setNames(vector("list", length(mods)), mods)
for(i in seq_along(mods)) {
  out.ls[i] <- stan(file=paste0("code/mods/PPM_mvPhy_", mods[i], "_IJK.stan"),
                    # sample_file=paste0("out/tests/", i), 
                    pars=pars.exc, include=F, chains=1,
                    data=read_rdump("data/stan_data/test_realData.Rdump"), 
                    thin=1, warmup=500, iter=750)
}
# out.ls <- map(mods, ~read_stan_csv(dir("out/tests", ., full.names=T)))

# lpd (Y_ | lambda_)
library(tidyverse); library(ggmcmc); theme_set(theme_bw())
gg.ll <- map_dfr(out.ls, ~ggs(., "log_lik_lambda_"), .id="model")
gg.ll.lpd <- gg.ll %>% group_by(model, Parameter) %>% 
  summarise(mn=log(mean(exp(value)))) %>% 
  group_by(model) %>% 
  summarise(lpd=-2*sum(mn)/n())
gg.ll.lpd %>% arrange(lpd)

# waic, loo
library(loo)
map(out.ls, ~loo::waic(loo::extract_log_lik(., "log_lik_lambda_"))) %>%
  loo::loo_compare()
map(out.ls, ~loo::loo(loo::extract_log_lik(., "log_lik_lambda_"))) %>%
  loo::loo_compare()

# parameter estimates
mod_cols <- c(W="#e41a1c", Y="#377eb8", WY="#984ea3", "WY_con"="black")

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

d.ls <- read_rdump("data/stan_data/test_realData.Rdump")
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
