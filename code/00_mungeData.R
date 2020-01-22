# Script for munging datasets into appropriate forms
# Tim Szewczyk

##--- set up
library(tidyverse); library(sf); library(googlesheets); library(rdrop2)
gis.dir <- "~/Documents/unil/GIS/data/"
ant.dir <- "~/Documents/unil/opfo_str_sampling/data/"
tax_i <- read_csv("data/tax_i.csv") %>% filter(sNum<40)
site_i <- read_csv(paste0(ant.dir, "opfo_siteSummaryProcessed.csv"))
plot_i <- read_csv(paste0(ant.dir, "opfo_envDataProcessed.csv")) %>% 
  arrange(BDM, Plot_id) %>% filter(Plot_id != "020201")
source("code/00_fn.R")

##--- settings
test_prop_W <- 0.3 # proportion of W data to use for testing models
test_prop_Y <- 0.2 # proportion of Y data to use for testing models
X_vars <- c("MAT", "MAT_sq", "Iso", "AP")
V_vars <- c("el", "slope")
U_vars <- c("slope")



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
box_Y.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/BDMplus.shp")) %>%
  select(Id) %>% rename(BDM=Id) %>% st_transform(21781) %>%
  left_join(., select(site_i, BDM, region, BDM_id), by="BDM") %>%
  arrange(BDM) %>% 
  mutate(n_plots=(plot_i %>% group_by(BDM) %>% 
           summarise(n_plots=n()) %>% arrange(BDM))$n_plots)



##--- W
dropbox_i <- drop_dir("Inventaire fourmis/5. BASE DE DONNEES/") %>% 
  select(name, path_lower, path_display) %>%
  filter(grepl("FourmisData", name))
drop_download(dropbox_i$path_lower[1], paste0(ant.dir), overwrite=T)
ant.W <- paste0(ant.dir, dropbox_i$name[1]) %>%
  readxl::read_xlsx(., 1) %>%
  rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID) %>%
  clean_species_names()
# align coordinate systems
ant.W <- rbind(ant.W %>% filter(is.na(LATITUDE)) %>%
                 filter(!is.na(SWISSCOORDINATE_X)) %>%
                 st_as_sf(coords=c("SWISSCOORDINATE_X", "SWISSCOORDINATE_Y")) %>%
                 st_set_crs(21781) %>%
                 select(TubeNo, SPECIESID), 
               ant.W %>% filter(!is.na(LONGITUDE)) %>%
                 st_as_sf(coords=c("LONGITUDE", "LATITUDE")) %>% 
                 st_set_crs(4326) %>% st_transform(21781) %>%
                 select(TubeNo, SPECIESID))
# add grid id
grid.W <- st_join(ant.W, grd_W.sf) %>%
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
ant.Y <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer",
                 col_types=cols(PLOT=col_character())) %>%
  filter(!is.na(lon)) %>% filter(TypeOfSample=="soil") %>%
  rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
  mutate(Plot_id=as.character(Plot_id)) %>%
  clean_species_names()
box.Y <- left_join(select(ant.Y, BDM, Plot_id, TubeNo, SPECIESID), 
                   box_Y.sf, by="BDM") %>%
  group_by(BDM, Plot_id, SPECIESID) %>%
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
J <- list(Y=floor(nrow(box_Y.sf)*(1-test_prop_Y)), 
          Y_=ceiling(nrow(box_Y.sf)*test_prop_Y))
IJ <- list(Y=rep(1:J$Y, times=box_Y.sf$n_plots[1:J$Y]), 
           Y_=rep(1:J$Y_, times=box_Y.sf$n_plots[(1:J$Y_)+J$Y]))
I <- list(Y=length(IJ$Y),
          Y_=length(IJ$Y_))



##--- tax_i: temporary
if(is.null(tax_i)) {
  tax_i <- tibble(species=sort(unique(c(ant.W$SPECIESID, ant.Y$SPECIESID))),
                  genus=str_split_fixed(species, "_", 2)[,1],
                  sNum=as.numeric(factor(species)),
                  gNum=as.numeric(factor(genus)),
                  Dprior=rbinom(length(species), 1, 0.3)+0.5)
  write_csv(tax_i, "data/tax_i.csv")
}



##--- X
# load covariates
clim_i <- read_csv(paste0(gis.dir, "bioclim_variable_names.csv"))
clim.f <- dir(paste0(gis.dir, "chelsa/")) %>%
  setNames(., clim_i$abbrev[match(str_remove(str_remove(., "CHELSA_"), ".tif"),
                               clim_i$code)])
dem <- raster::raster(paste0(gis.dir, "aster/CH_21781.tif")) %>%
  raster::crop(., VD_raw) 
slope <- raster::raster(paste0(gis.dir, "aster/CH_slope_21781.tif")) %>%
  raster::crop(., VD_raw) 
clim.r <- map(clim.f[c(1,3,12)], 
              ~raster::raster(paste0(gis.dir, "chelsa/", .x)) %>%
                raster::projectRaster(., dem) %>%
                raster::crop(., VD_raw))
# calculate means within W cells
X_W.df <- filter(grd_W.sf, id %in% W_id) %>%
  mutate(el=raster::extract(dem, ., fun=mean),
         slope=raster::extract(slope, ., fun=mean)) %>%
  bind_cols(., map_dfc(clim.r, ~raster::extract(., filter(grd_W.sf, id%in%W_id), 
                                                fun=mean))) %>%
  arrange(id)
# calculate means within Y boxes
X_Y.df <- box_Y.sf %>%
  mutate(el=raster::extract(dem, ., fun=mean),
         slope=raster::extract(slope, ., fun=mean)) %>%
  bind_cols(., map_dfc(clim.r, ~raster::extract(., box_Y.sf, fun=mean))) %>%
  arrange(BDM)
# select covariates and add quadratic terms
X_W.mx <- as.matrix(X_W.df %>% st_set_geometry(NULL) %>% 
                      select(-id, -inbd) %>% 
                      mutate_all(.funs=list(sq=~.^2)) %>%
                      select(X_vars)) 
X_Y.mx <- as.matrix(X_Y.df %>% st_set_geometry(NULL) %>% 
                      select(-BDM, -region, -BDM_id, -n_plots) %>% 
                      mutate_all(.funs=list(sq=~.^2)) %>%
                      select(X_vars))
# scale and carefully combine
X.scale <- scale(rbind(X_W.mx[1:K$W,], 
                       X_Y.mx[1:J$Y,],
                       X_W.mx[K$W+(1:K$W_),], 
                       X_Y.mx[J$Y+(1:J$Y_),]))
X.all <- cbind(1, X.scale)



##--- V
V.df <- st_read(paste0(ant.dir, "gis/opfo_soil_25per.shp")) %>% 
  st_transform(st_crs(box_Y.sf)) %>%
  mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
  select(Plot_id) %>% 
  left_join(., select(plot_i, BDM, Plot_id), by="Plot_id") %>% 
  mutate(el=raster::extract(dem, .),
         slope=raster::extract(slope, .)) %>%
  arrange(BDM, Plot_id) 
V.mx <- as.matrix(V.df %>% st_set_geometry(NULL) %>% 
                    select(-BDM, -Plot_id) %>% 
                    mutate_all(.funs=list(sq=~.^2)) %>%
                    select(V_vars))
V.scale <- scale(rbind(V.mx[1:I$Y,], 
                       V.mx[I$Y+(1:I$Y_),]))


##--- U
# load covariates
slope <- raster::raster(paste0(gis.dir, "aster/CH_slope_21781.tif")) %>%
  raster::crop(., VD_raw) 
# calculate means within W cells
U_W.df <- grd_W.sf %>% filter(id %in% W_id) %>%
  mutate(slope=raster::extract(slope, ., fun=mean)) %>% 
  st_set_geometry(NULL) %>% arrange(id)
U_W.mx <- as.matrix(U_W.df %>% select(U_vars))
U.scale <- scale(U_W.mx)
U.all <- cbind(1, U.scale)



##--- export stan data
stan_data <- list(K=K$W, K_=K$W_, 
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



# Testing
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
out.ls <- list(W=stan(file="code/mods/PPM_mvPhy_W_IJK.stan",
                 data=stan_data, thin=5, warmup=1000, iter=2000),
               Y=stan(file="code/mods/PPM_mvPhy_Y_IJK.stan",
                       data=stan_data, thin=5, warmup=1000, iter=2000),
               WY=stan(file="code/mods/PPM_mvPhy_WY_IJK.stan",
                      data=stan_data, thin=5, warmup=1000, iter=2000))

# lpd (Y_ | lambda_)
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
library(ggmcmc); theme_set(theme_bw())
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
                            list(sp=matrix(nrow=stan_data$R, ncol=stan_data$S)), 
                            tax_i[,3:5], "b")
b.out$sum.gg %>% select(model, Parameter, mn, R) %>% 
  spread(model, mn) %>%
  ggplot(aes(x=W, y=WY)) + facet_wrap(~R, scales="free") +
  geom_abline(linetype=2) +
  geom_point()
  # labs(title="b", x="Model", y="RMSE among species") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/b_RMSE.pdf", width=8, height=3, units="in")