# Script for munging datasets into appropriate forms
# Tim Szewczyk



##--- set up
library(tidyverse); library(sf); library(googlesheets); library(rdrop2)
gis.dir <- "~/Documents/unil/GIS/data/VD_21781/"
ant.dir <- "~/Documents/unil/opfo_str_sampling/data/"
tax_i <- read_csv("data/tax_i.csv") %>% filter(sNum<50)
plot_i <- read_csv(paste0(ant.dir, "opfo_envDataProcessed.csv")) %>% 
  arrange(BDM, Plot_id) %>% filter(Plot_id != "020201")
source("code/00_fn.R"); source(paste0(ant.dir, "../code/opfo_fn.R"))



##--- settings
test_prop_W <- 0.3 # proportion of W data to use for testing models
test_prop_Y <- 0.2 # proportion of Y data to use for testing models
X_vars <- c("MAT", "Iso", "AP", "NPP")
V_vars <- c("el", "slope", "pop")
U_vars <- c("slope", "pop", "rdLen")



##--- grids
# W: 1km2 grid of Vaud
VD_raw <- st_read(paste0(gis.dir, "Vaud_boundaries.shp")) 
VD <- st_union(VD_raw)
VD_ext <- raster::extent(matrix(st_bbox(VD), ncol=2))
grd_W <- raster::raster(ext=VD_ext, crs=st_crs(VD), resolution=1000) %>%
  raster::rasterize(VD_ext, .) %>% 
  raster::rasterToPolygons(., n=4)
grd_W@data$layer <- 1:raster::ncell(grd_W)
grd_W.sf <- st_as_sf(grd_W) %>% st_set_crs(st_crs(VD)) %>% rename(id=layer) %>%
  mutate(inbd=c(st_intersects(., VD, sparse=F)))
# Y: 44 x 1km2 sampling sites
site.sf <- agg_str_site_data(gis.dir) %>% arrange(BDM)



##--- W
ants <- load_ant_data(pub.dir=ant.dir, str_type="soil", clean_spp=T)
grid.W <- st_join(ants$pub, grd_W.sf) %>%
  group_by(SPECIESID, id) %>% st_set_geometry(NULL) %>%
  summarise(nObs=n()) %>% 
  filter(id %in% filter(grd_W.sf, inbd)$id) %>%
  arrange(id)
W <- matrix(0, nrow=nrow(grd_W.sf), ncol=max(tax_i$sNum))
for(s in 1:max(tax_i$sNum)) {
  sp.occ <- filter(grid.W, SPECIESID==tax_i$species[s])
  if(nrow(sp.occ)>0) W[sp.occ$id,s] <- sp.occ$nObs
}
W_id <- which(rowSums(W)>0)



##--- Y
box.Y <- ants$str %>% group_by(BDM, Plot_id, SPECIESID) %>%
  summarise(nObs=n()) %>%
  arrange(BDM)
Y <- matrix(0, nrow=nrow(plot_i), ncol=max(tax_i$sNum))
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



##--- tax_i: temporary
if(is.null(tax_i)) {
  tax_i <- tibble(species=sort(unique(c(ants$pub$SPECIESID, ants$str$SPECIESID))),
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
NPP <- raster::raster(paste0(gis.dir, "MODIS_2010-2019_VD_21781.tif"))
pop <- raster::raster(paste0(gis.dir, "popDensity_VD_21781.tif")) 
roads <- st_read(paste0(gis.dir, "roads_VD_21781.shp")) %>% st_set_crs(21781)



##--- X
# calculate means within W cells
X_W.df <- filter(grd_W.sf, id %in% W_id) %>%
  mutate(el=raster::extract(dem, ., fun=mean),
         slope=raster::extract(slope, ., fun=mean),
         NPP=raster::extract(NPP, ., fun=mean)) %>%
  bind_cols(., map_dfc(clim, ~raster::extract(., filter(grd_W.sf, id%in%W_id), 
                                                fun=mean))) %>%
  arrange(id)
# select covariates and add quadratic terms
X_W.mx <- as.matrix(X_W.df %>% st_set_geometry(NULL) %>% 
                      select(-id, -inbd) %>% 
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      select(one_of(X_vars)))
X_Y.mx <- as.matrix(site.sf %>% st_set_geometry(NULL) %>% 
                      select(-BDM, -region, -BDM_id, -Sample_date,
                             -contains("nPlot")) %>% 
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      select(one_of(X_vars)))
# scale and carefully combine
X.scale <- scale(rbind(X_W.mx[1:K$W,], 
                       X_Y.mx[1:J$Y,],
                       X_W.mx[K$W+(1:K$W_),], 
                       X_Y.mx[J$Y+(1:J$Y_),]))
X.all <- cbind(1, X.scale)



##--- V
V.df <- st_read(paste0(ant.dir, "gis/opfo_soil_25per.shp")) %>% 
  st_transform(st_crs(site.sf)) %>%
  mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
  select(Plot_id) %>% 
  left_join(., select(plot_i, BDM, Plot_id), by="Plot_id") %>% 
  mutate(el=raster::extract(dem, .),
         slope=raster::extract(slope, .),
         pop=replace_na(raster::extract(pop, .), 0)) %>%
  arrange(BDM, Plot_id) 
V.mx <- as.matrix(V.df %>% st_set_geometry(NULL) %>% 
                    select(-BDM, -Plot_id) %>% 
                    mutate_all(.funs=list(sq=~.^2)) %>%
                    select(V_vars))
V.scale <- scale(rbind(V.mx[1:I$Y,], 
                       V.mx[I$Y+(1:I$Y_),]))


##--- U
# load covariates
rdLen.grd_W <- roads %>%
  st_intersection(., filter(grd_W.sf, id %in% W_id)) %>%
  mutate(length=st_length(.)) %>%
  group_by(id) %>%
  summarise(rdLen=as.numeric(sum(length)))

# calculate means within W cells
U_W.df <- grd_W.sf %>% filter(id %in% W_id) %>%
  mutate(slope=raster::extract(slope, ., fun=mean),
         pop=raster::extract(pop, ., fun=mean, na.rm=T)) %>% 
  arrange(id)
U_W.df <- left_join(U_W.df, st_set_geometry(rdLen.grd_W, NULL), by="id")
U_W.df[is.na(U_W.df)] <- 0
U_W.mx <- as.matrix(st_set_geometry(U_W.df, NULL) %>% select(U_vars))
U.scale <- scale(U_W.mx)
U.all <- cbind(1, U.scale)



##--- export stan data
d.ls <- list(K=K$W, K_=K$W_, 
             J=J$Y, J_=J$Y_, 
             IJ=IJ$Y, IJ_=IJ$Y_, 
             I=I$Y, I_=I$Y_,
             S=max(tax_i$sNum), 
             G=max(tax_i$gNum), 
             tax_i=tax_i[,c("sNum", "gNum")], 
             D_prior=tax_i$Dprior, 
             R=dim(X.all)[2], 
             L=dim(V.scale)[2],
             Q=dim(U.all)[2], 
             W=W[1:K$W,], 
             W_=W[K$W+(1:K$W_),], 
             Y=Y[1:I$Y,],
             Y_=Y[I$Y+(1:I$Y_),],
             X=X.all[1:(K$W+J$Y),], 
             X_=X.all[(K$W+J$Y)+(1:(K$W_+J$Y_)),], 
             V=V.scale[1:I$Y,], 
             V_=V.scale[I$Y+(1:I$Y_),],
             U=U.all[1:K$W,], 
             U_=U.all[K$W+(1:K$W_),], 
             h=7.5e-7)
rstan::stan_rdump(ls(d.ls),
                  file=paste0("data/stan_data/test_realData.Rdump"),
                  envir=list2env(d.ls))



# Testing
library(rstan)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
mods <- c("W")#, "Y", "WY")
pars.exc <- c("a.std", "A.std", "b.std", "B.std", "LAMBDA", "LAMBDA_", "lambda")
out.ls <- vector("list", length(mods))
for(i in seq_along(mods)) {
  out.ls[i] <- stan(file=paste0("code/mods/PPM_mvPhy_", mods[i], "_IJK.stan"),
                    # sample_file=paste0("out/tests/", i), 
                    pars=pars.exc, include=F,
                    data=read_rdump("data/stan_data/test_realData.Rdump"), 
                    thin=5, warmup=1000, iter=2000)
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
  geom_density() + scale_colour_manual(values=mod_cols) 

alpha.out <- aggregate_aggSlopes(out.ls, list(agg=NA), "alpha")
ggplot(alpha.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_rug(data=alpha.out$post_mns, aes(x=mean), sides="b") +
  geom_density() + scale_colour_manual(values=mod_cols) 

eta.out <- aggregate_eta(out.ls, NA)
ggplot(eta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_density() + scale_colour_manual(values=mod_cols) 

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
