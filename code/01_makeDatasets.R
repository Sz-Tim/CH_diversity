# 01_makeDatasets
# Tim Szewczyk
#
# Script for munging datasets into appropriate forms



##--- set up
set_type <- c("vs", "pred")[1]  # vs = variable selection; pred = all of VD
library(tidyverse); library(sf); library(googlesheets)
gis.dir <- "../2_gis/data/VD_21781/"
ant.dir <- "../1_opfo/data/"
tax_i <- read_csv("data/tax_i.csv") 
veg_i <- read_csv(paste0(ant.dir, "vegcover_id.csv")) # van der Maarel 2007
lc_i <- readxl::read_xlsx(paste0(ant.dir, "landcover_id.xlsx"), 1)
plot_i <- read_csv(paste0(ant.dir, "opfo_envDataProcessed.csv")) %>% 
  arrange(BDM, Plot_id) %>% filter(Plot_id != "020201") %>%
  group_by(BDM) %>% 
  mutate(SoilTSt=(SoilTemp-mean(SoilTemp,na.rm=T))/sd(SoilTemp,na.rm=T),
         CnpyOpn=(lc_i$Canopy[match(Categorie, lc_i$LC)]=="Open")*1,
         CnpyMxd=(lc_i$Canopy[match(Categorie, lc_i$LC)]=="Mixed")*1,
         Pasture=(TypeOfOpen=="pasture")*1,
         Crop=(TypeOfOpen=="crop")*1,
         Meadow=(TypeOfOpen=="meadow")*1,
         Built=(Categorie %in% c("ZoneConstruite", "transport"))) %>%
  mutate_at(c("Grass", "Forb", "Shrub", "Bare", "Litter", "Moss"), 
            ~veg_i$Pct[match(., veg_i$class)]) %>%
  mutate(VegTot=Grass + Forb + Shrub)
plot_i$Pasture[plot_i$Categorie=="PrairieSeche"] <- 1
plot_i$Pasture[is.na(plot_i$Pasture)] <- 0
plot_i$Crop[is.na(plot_i$Crop)] <- 0
plot_i$Built[is.na(plot_i$Built)] <- 0
source("code/00_fn.R"); source(paste0(ant.dir, "../code/00_fn.R"))



##--- settings
test_prop_Y <- ifelse(set_type=="vs", 0.2, 0)  # prop BDM to withold
X_vars <- c("grwnDD0", "grwnDD0_sq", 
            "AP",
            "npp",
            "lcH",
            "Edge",
            "bldgPer", "rdLen",
            "aspctN")
V_vars <- c("SoilTSt", 
            "VegTot",
            "CnpyOpn", "CnpyMxd",
            "Pasture", "Crop",
            "aspctN")




##--- grids
# W: 1km2 grid of Vaud
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp") %>%
  filter(!grepl("Lac ", NAME)) 
VD <- st_union(VD_raw)
VD_ext <- raster::extent(matrix(st_bbox(VD), ncol=2))
grd_W <- raster::raster(ext=VD_ext, crs=sp::CRS("+init=epsg:21781"), 
                        resolution=1000) %>%
  raster::rasterize(VD_ext, .) %>% 
  raster::rasterToPolygons(., n=4)
grd_W@data$layer <- 1:raster::ncell(grd_W)
grd_W.sf <- st_as_sf(grd_W) %>% st_transform(st_crs(VD)) %>% rename(id=layer) %>%
  mutate(inbd=c(st_intersects(., VD, sparse=F)))
# Y: 44 x 1km2 sampling sites
site.sf <- agg_str_site_data() %>% arrange(BDM)




##--- W
ants <- load_ant_data(str_type="soil", clean_spp=T)
grid.W <- st_join(ants$pub, grd_W.sf) %>%
  group_by(SPECIESID, id, inbd) %>% st_set_geometry(NULL) %>%
  summarise(nObs=n()) %>% ungroup %>%
  filter(inbd) %>% arrange(id)
W <- matrix(0, nrow=nrow(grd_W.sf), ncol=max(tax_i$sNum), 
            dimnames=list(grd_W.sf$id, tax_i$species))
for(s in 1:max(tax_i$sNum)) {
  sp.occ <- filter(grid.W, SPECIESID==tax_i$species[s])
  if(nrow(sp.occ)>0) W[sp.occ$id,s] <- sp.occ$nObs
}
W_id <- which(rowSums(W)>0)
W_inbd <- which(grd_W.sf$inbd)
grd_W.sf$W_obs <- rowSums(W)>0




##--- Y
box.Y <- ants$str %>% group_by(BDM, Plot_id, SPECIESID) %>%
  summarise(nObs=n()) %>% 
  left_join(., select(plot_i, BDM, Plot_id), 
            by=c("BDM", "Plot_id")) %>%
  arrange(BDM, Plot_id)
Y <- matrix(0, nrow=nrow(plot_i), ncol=max(tax_i$sNum),
            dimnames=list(plot_i$Plot_id, tax_i$species))
for(s in 1:max(tax_i$sNum)) {
  sp.occ <- filter(box.Y, SPECIESID==tax_i$species[s])
  if(nrow(sp.occ)>0) Y[which(plot_i$Plot_id %in% sp.occ$Plot_id),s] <- sp.occ$nObs
}





##--- indexes
if(set_type=="pred") {
  K <- list(id_W=W_id,  # row of W, row/id of grd_W.sf
            id_W_=W_inbd[-which(W_inbd %in% W_id)])
  K$W <- length(K$id_W)
  K$W_ <- length(K$id_W_) 
} else {
  K <- list(id_W=W_id,
            W=length(W_id))
}

J <- list(id_Y=sample(1:nrow(site.sf), floor(nrow(site.sf)*(1-test_prop_Y))))
J$BDM_Y <- site.sf$BDM[J$id_Y]
J$Y <- length(J$id_Y)
if(test_prop_Y>0) {
  J$id_Y_ <- (1:nrow(site.sf))[-J$id_Y]
  J$BDM_Y_ <- site.sf$BDM[J$id_Y_]
  J$Y_ <- length(J$id_Y_)
}

if(test_prop_Y>0) {  # test/train
  IJ.ls <- list(Y=map(J$BDM_Y, ~which(plot_i$BDM==.x)),
                Y_=map(J$BDM_Y_, ~which(plot_i$BDM==.x)))
  IJ <- list(id_Y=unlist(IJ.ls$Y),   # row of plot_i
             id_Y_=unlist(IJ.ls$Y_),
             Y=unlist(imap(IJ.ls$Y, ~rep(.y, length(.x)))),
             Y_=unlist(imap(IJ.ls$Y_, ~rep(.y, length(.x)))))
  I <- list(Y=length(IJ$Y),
            Y_=length(IJ$Y_))
} else {  # no test
  IJ.ls <- list(Y=map(J$BDM_Y, ~which(plot_i$BDM==.x)))
  IJ <- list(id_Y=unlist(IJ.ls$Y),   # row of plot_i
             Y=unlist(imap(IJ.ls$Y, ~rep(.y, length(.x)))))
  I <- list(Y=length(IJ$Y))
}


# 
# 
# ##--- use only species in Y
# if(set_type=="vs") {
#   in_Y <- which(colSums(Y[IJ$id_Y,])>0)
#   W <- W[,in_Y]
#   W_id <- which(rowSums(W)>0)
#   W_inbd <- which(grd_W.sf$inbd)
#   grd_W.sf$W_obs <- rowSums(W)>0
#   Y <- Y[,in_Y]
#   tax_i <- tax_i %>% filter(species %in% names(in_Y)) %>%
#     mutate(sNum=row_number(), gNum=as.numeric(factor(genus)))
# }




##--- tax_i: temporary
if(is.null(tax_i)) {
  tax_i <- tibble(species=sort(unique(ants$all$SPECIESID)),
                  genus=str_split_fixed(species, "_", 2)[,1],
                  sNum=as.numeric(factor(species)),
                  gNum=as.numeric(factor(genus)),
                  Dprior=1)
  write_csv(tax_i, "data/tax_i.csv")
}




##--- covariates
clim <- dir(gis.dir, "chelsa") %>%
  setNames(str_remove(., "_chelsa_VD_21781.tif$")) %>%
  magrittr::extract(names(.) %in% X_vars) %>%
  map(., ~raster::raster(paste0(gis.dir, .x)))
envirem <- dir(gis.dir, "envirem") %>%
  setNames(str_remove(., "_envirem_VD_21781.tif$")) %>%
  # magrittr::extract(names(.) %in% X_vars) %>%
  map(., ~raster::raster(paste0(gis.dir, .x)))
dem <- raster::raster(paste0(gis.dir, "dem_VD_21781.tif"))
aspect <- raster::raster(paste0(gis.dir, "aspect_VD_21781.tif"))
slope <- raster::raster(paste0(gis.dir, "slope_VD_21781.tif"))
npp <- raster::raster(paste0(gis.dir, "MODIS_2010-2019_VD_21781.tif"))
lc.W <- st_read(paste0(gis.dir, "lc_1kmZones_21781.shp")) %>%
  st_set_geometry(NULL) %>%
  filter(inbd==1) %>% arrange(id) %>% 
  mutate(LC_Forest=LC_4+LC_5+LC_6) %>%
  select(3:18)
lc.Y <- st_read(paste0(gis.dir, "site_lcZones_21781.shp")) %>% 
  arrange(BDM) %>% st_set_geometry(NULL) %>% select(1,3:17) %>%
  mutate(LC_Forest=LC_4+LC_5+LC_6)
pop <- raster::raster(paste0(gis.dir, "popDensity_VD_21781.tif")) 
roads <- st_read(paste0(gis.dir, "roads_VD_21781.shp")) %>% st_set_crs(21781) 
# rdLen.W <- roads %>% st_intersection(., filter(grd_W.sf, inbd)) %>%
#   mutate(length=st_length(.)) %>% st_set_geometry(NULL) %>%
#   group_by(id) %>% summarise(rdLen=as.numeric(sum(length)))
rdLen.Y <- roads %>% st_intersection(., site.sf %>% st_set_crs(21781)) %>%
  mutate(length=st_length(.)) %>% st_set_geometry(NULL) %>%
  group_by(BDM) %>% summarise(rdLen=as.numeric(sum(length)))
# bldgs.W <- st_read(paste0(gis.dir, "bldgs_1km_21781.shp")) %>% 
#   mutate(area=st_area(.), perimeter=lwgeom::st_perimeter(.)) %>% 
#   st_set_geometry(NULL) %>% group_by(id) %>%
#   summarise(area=as.numeric(sum(area)), 
#             perimeter=as.numeric(sum(perimeter)))
bldgs.Y <- st_read(paste0(gis.dir, "bldgs_site_21781.shp")) %>% 
  mutate(area=st_area(.), perimeter=lwgeom::st_perimeter(.)) %>% 
  st_set_geometry(NULL) %>% group_by(BDM) %>%
  summarise(area=as.numeric(sum(area)), 
            perimeter=as.numeric(sum(perimeter)))



##--- X
# calculate means within W cells
# X_W.df <- filter(grd_W.sf, inbd) %>% 
#   mutate(lcH=vegan::diversity(as.matrix(lc.W)[,-16]),
#          Forest=lc.W$LC_Forest,
#          Edge=lc.W$LC_9,
#          bldgArea=bldgs.W$area[match(id, bldgs.W$id)],
#          bldgArea=log(replace(bldgArea, is.na(bldgArea), 0)+1),
#          bldgPer=bldgs.W$perimeter[match(id, bldgs.W$id)],
#          bldgPer=log(replace(bldgPeri, is.na(bldgPeri), 0)+1),
#          rdLen=rdLen.W$rdLen[match(id, rdLen.W$id)],
#          rdLen=log(replace(rdLen, is.na(rdLen), 0)+1),
#          el=raster::extract(dem, ., fun=mean),
#          aspect=raster::extract(aspect, ., fun=mean, na.rm=T),
#          slope=raster::extract(slope, ., fun=mean),
#          npp=raster::extract(npp, ., fun=mean, na.rm=T)) %>%
#   bind_cols(., map_dfc(clim, ~raster::extract(., filter(grd_W.sf, inbd),
#                                                 fun=mean))) %>%
#   arrange(id)
# write_sf(X_W.df, "data/cov/X_W-df.shp")

# cos(aspect*pi/180): N=1, S=-1, E/W=0
# sin(aspect*pi/180): E=1, W=-1, N/S=0
X_W.df <- st_read("data/cov/X_W-df.shp") %>% arrange(id)
# select covariates and add quadratic terms
X_W.mx <- as.matrix(
  X_W.df %>% st_set_geometry(NULL) %>% 
    select(-inbd) %>% 
    mutate(aspctN=cos(aspect*pi/180),
           aspctE=sin(aspect*pi/180)) %>%
    mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
    select("id", "el", one_of(X_vars))
)
X_Y.mx <- as.matrix(
  site.sf %>% select(-region, -Sample_date, -contains("nPlot"), -npp) %>%
    mutate(lcH=vegan::diversity(as.matrix(lc.Y[,-c(1, 17)])),
           Forest=lc.Y$LC_Forest,
           Edge=lc.Y$LC_9,
           grwnDD0=raster::extract(envirem$growingDegDays0, ., fun=mean),
           bldgPer=bldgs.Y$perimeter[match(BDM, bldgs.Y$BDM)],
           bldgPer=replace(bldgPer, is.na(bldgPer), 0),
           bldgPer=log(bldgPer+1),
           rdLen=rdLen.Y$rdLen[match(BDM, rdLen.Y$BDM)],
           rdLen=log(rdLen+1),
           npp=raster::extract(npp, ., fun=mean, na.rm=T),
           aspect=raster::extract(aspect, ., fun=mean, na.rm=T)) %>%
    mutate(aspctN=cos(aspect*pi/180),
           aspctE=sin(aspect*pi/180)) %>%
    mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
    st_set_geometry(NULL) %>% 
    select("BDM", "el", one_of(X_vars))
  )
# scale and carefully combine
X.mx <- rbind(X_W.mx[match(K$id_W, X_W.mx[,'id']),], 
              X_Y.mx[match(J$BDM_Y, X_Y.mx[,'BDM']),])
if(set_type=="pred") {
  X.mx <- rbind(X.mx,
                X_W.mx[match(K$id_W_, X_W.mx[,'id']),])
} 
if(test_prop_Y > 0) {  # test/train
  X.mx <- rbind(X.mx,
                X_Y.mx[match(J$BDM_Y_, X_Y.mx[,'BDM']),])
}
X.scale <- scale(X.mx[,-(1:2)])
if(any(grep("_sq", X_vars))) {
  sqX <- grep("_sq", X_vars, value=T) %>% list(base=str_sub(., 1, -4), sq=.)
  X.scale[,sqX$sq] <- X.scale[,sqX$base]^2 
}
X.all <- cbind(1, X.scale)




##--- V
# V.df <- st_read("../2_gis/data/opfo/opfo_soil_25per.shp") %>%
#   st_transform(st_crs(site.sf)) %>%
#   mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
#   select(Plot_id) %>%
#   left_join(., select(plot_i, BDM, Plot_id, any_of(V_vars)), by="Plot_id") %>%
#   mutate(el=raster::extract(dem, .),
#          slope=raster::extract(slope, .),
#          pop=replace_na(raster::extract(pop, .), 0),
#          aspect=raster::extract(aspect, ., fun=mean, na.rm=T)) %>%
#   mutate(aspctN=cos(aspect*pi/180),
#          aspctE=sin(aspect*pi/180)) %>%
#   arrange(BDM, Plot_id)
# write_sf(V.df, "data/cov/V-df.shp")
V.df <- st_read("data/cov/V-df.shp") %>% 
  arrange(BDM, Plot_id)
V_Y.mx <- as.matrix(V.df %>% st_set_geometry(NULL) %>% 
                      mutate_if(is.numeric, .funs=list(sq=~.^2)) %>%
                      mutate(Plot_id=as.numeric(as.character(Plot_id))) %>%
                      select("Plot_id", "el", all_of(V_vars)))
V.mx <- V_Y.mx[IJ$id_Y,]
if(test_prop_Y > 0) {  # test/train
  V.mx <- rbind(V.mx, 
                V_Y.mx[IJ$id_Y_,]) 
}
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
if("Pasture" %in% V_vars) {
  V.scale[,"Pasture"] <- V.scale[,"Pasture"] * 
    attr(V.scale, "scaled:scale")["Pasture"] + 
    attr(V.scale, "scaled:center")["Pasture"]
}
if("Crop" %in% V_vars) {
  V.scale[,"Crop"] <- V.scale[,"Crop"] * 
    attr(V.scale, "scaled:scale")["Crop"] + 
    attr(V.scale, "scaled:center")["Crop"]
}
if("Built" %in% V_vars) {
  V.scale[,"Built"] <- V.scale[,"Built"] * 
    attr(V.scale, "scaled:scale")["Built"] + 
    attr(V.scale, "scaled:center")["Built"]
}




##--- NA's
# identify NA cells
na.W <- which(is.na(rowSums(X.mx[1:K$W,]))) # rows of X with NA
na.Y <- which(is.na(rowSums(V.mx[1:I$Y,]))) # rows of V with NA
if(set_type=="pred") {
  na.W_ <- which(is.na(rowSums(X.mx[K$W+J$Y+(1:K$W_),]))) # rows of X_ with NA
} else {
  na.W_ <- numeric(0)
}
if(test_prop_Y > 0) {  # test/train
  na.Y_ <- which(is.na(rowSums(V.mx[I$Y+(1:I$Y_),]))) # rows of V_ with NA
}



##--- export stan data
if(set_type=="vs") {
  if(test_prop_Y > 0) {  # test/train
    d.ls <- list(K=K$W-length(na.W),
                 J=J$Y, J_=J$Y_,
                 IJ=IJ$Y[-na.Y], IJ_=IJ$Y_[-na.Y_],
                 I=I$Y-length(na.Y), I_=I$Y_-length(na.Y_),
                 S=max(tax_i$sNum), G=max(tax_i$gNum), 
                 tax_i=tax_i[,c("sNum", "gNum")], 
                 D_prior=tax_i$Dprior, 
                 R=dim(X.all)[2], 
                 L=dim(V.scale)[2],
                 W=W[K$id_W,][-na.W,],
                 Y=Y[IJ$id_Y,][-na.Y,],
                 Y_=Y[IJ$id_Y_,][-na.Y_,],
                 Y_int=Y[IJ$id_Y_,][-na.Y_,],
                 X=X.all[1:(K$W+J$Y),][-na.W,], 
                 X_=X.all[K$W+J$Y+(1:J$Y_),], 
                 V=V.scale[1:I$Y,][-na.Y,], 
                 V_=V.scale[I$Y+(1:I$Y_),][-na.Y_,],
                 h=7.5e-7)
    d.i <- list(X=X.mx[1:(K$W+J$Y),][-na.W,],
                X_=X.mx[K$W+J$Y+(1:J$Y_),],
                V=V.mx[1:I$Y,][-na.Y,],
                V_=V.mx[I$Y+(1:I$Y_),][-na.Y_,],
                W=W[K$id_W,][-na.W,], 
                Y=Y[IJ$id_Y,][-na.Y,],
                Y_=Y[IJ$id_Y_,][-na.Y_,],
                I=I, J=J, IJ=IJ, K=K, 
                na.W=na.W, na.Y=na.Y, na.Y_=na.Y_,
                tax_i=tax_i,
                grd_W.sf=grd_W.sf)
  d.f <- paste0("data/stan_data/", set_type, "_", test_prop_Y*100)
  saveRDS(d.ls, paste0(d.f, "_ls.rds"))
  saveRDS(d.i, paste0(d.f, "_i.rds"))
  rstan::stan_rdump(ls(d.ls),
                    file=paste0(d.f, ".Rdump"),
                    envir=list2env(d.ls))
  } else {
    d.ls <- list(K=K$W-length(na.W), 
                 J=J$Y, 
                 IJ=IJ$Y[-na.Y],
                 I=I$Y-length(na.Y),
                 S=max(tax_i$sNum), G=max(tax_i$gNum), 
                 tax_i=tax_i[,c("sNum", "gNum")], 
                 D_prior=tax_i$Dprior, 
                 R=dim(X.all)[2], 
                 L=dim(V.scale)[2],
                 Q=dim(U.all)[2], 
                 W=W[K$id_W,],
                 Y=Y[IJ$id_Y,][-na.Y,],
                 X=X.all[1:(K$W+J$Y),], 
                 V=V.scale[1:I$Y,][-na.Y,], 
                 U=U.all[1:K$W,], 
                 h=7.5e-7)
    d.i <- list(X=X.mx[1:(K$W+J$Y),],
                V=V.mx[1:I$Y,][-na.Y,],
                U=U_W.mx[1:K$W,],
                W=W[K$id_W,], 
                Y=Y[IJ$id_Y,][-na.Y,],
                I=I, J=J, IJ=IJ, K=K, 
                na.W=na.W, na.Y=na.Y,
                tax_i=tax_i,
                grd_W.sf=grd_W.sf[K$id_W,])
    if(length(na.W)>0) {
      d.ls$W <- d.ls$W[-na.W,]
      d.ls$X <- d.ls$X[-na.W,]
      d.ls$U <- d.ls$U[-na.W,]
      d.i$W <- d.i$W[-na.W,]
      d.i$X <- d.i$X[-na.W,]
      d.i$U <- d.i$U[-na.W,]
      d.i$grd_W.sf <- d.i$grd_W.sf[-na.W,]
    }
    d.f <- paste0("data/stan_data/", set_type, "_no_pred")
    saveRDS(d.ls, paste0(d.f, "_ls.rds"))
    saveRDS(d.i, paste0(d.f, "_i.rds"))
    rstan::stan_rdump(ls(d.ls),
                      file=paste0(d.f, ".Rdump"),
                      envir=list2env(d.ls))
  }
} else {
  d.ls <- list(K=K$W-length(na.W), K_=K$W_-length(na.W_),
               # part_K_=round(seq(1, K$W_-length(na.W_)+1, length.out=5)),
               J=J$Y, 
               IJ=IJ$Y[-na.Y],
               I=I$Y-length(na.Y),
               S=max(tax_i$sNum), G=max(tax_i$gNum), 
               tax_i=tax_i[,c("sNum", "gNum")], 
               D_prior=tax_i$Dprior, 
               R=dim(X.all)[2], 
               L=dim(V.scale)[2],
               Q=dim(U.all)[2], 
               W=W[K$id_W,],
               W_=W[K$id_W_,],
               Y=Y[IJ$id_Y,][-na.Y,],
               X=X.all[1:(K$W+J$Y),], 
               X_=X.all[K$W+J$Y+(1:(K$W_)),], 
               V=V.scale[1:I$Y,][-na.Y,], 
               U=U.all[1:K$W,], 
               U_=U.all[K$W+(1:K$W_),], 
               h=7.5e-7)
  d.i <- list(X=X.mx[1:(K$W+J$Y),],
              X_=X.mx[K$W+J$Y+(1:(K$W_)),],
              V=V.mx[1:I$Y,][-na.Y,],
              U=U_W.mx[1:K$W,],
              U_=U_W.mx[K$W+(1:K$W_),],
              W=W[K$id_W,], 
              W_=W[K$id_W_,], 
              Y=Y[IJ$id_Y,][-na.Y,],
              I=I, J=J, IJ=IJ, K=K, 
              na.W=na.W, na.W_=na.W_, na.Y=na.Y,
              tax_i=tax_i,
              grd_W.sf=grd_W.sf[c(K$id_W, K$id_W_),])
  if(length(na.W)>0) {
    d.ls$W <- d.ls$W[-na.W,]
    d.ls$X <- d.ls$X[-na.W,]
    d.ls$U <- d.ls$U[-na.W,]
    d.i$W <- d.i$W[-na.W,]
    d.i$X <- d.i$X[-na.W,]
    d.i$U <- d.i$U[-na.W,]
  }
  if(length(na.W_)>0) {
    d.ls$W_ <- d.ls$W_[-na.W_,]
    d.ls$X_ <- d.ls$X_[-na.W_,]
    d.ls$U_ <- d.ls$U_[-na.W_,]
    d.i$W_ <- d.i$W_[-na.W_,]
    d.i$X_ <- d.i$X_[-na.W_,]
    d.i$U_ <- d.i$U_[-na.W_,]
  }
  if(length(na.W)>0 | length(na.W_)>0) {
    d.i$grd_W.sf <- d.i$grd_W.sf[-c(na.W, na.W_),]
  }
  d.f <- paste0("data/stan_data/", set_type)
  saveRDS(d.ls, paste0(d.f, "_ls.rds"))
  saveRDS(d.i, paste0(d.f, "_i.rds"))
  rstan::stan_rdump(ls(d.ls),
                    file=paste0(d.f, ".Rdump"),
                    envir=list2env(d.ls))

} 


