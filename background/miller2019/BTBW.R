# TMS NOTES
# For WinBUGS, univariate CAR models are specified with:
#
# adj[] : A vector listing the ID numbers of the adjacent areas for each area
# (this is a sparse representation of the full adjacency matrix for the study
# region, and can be generated using the Adjacency Tool from the Map menu in
# GeoBUGS.
#
# weights[] : A vector the same length as adj[] giving unnormalised weights
# associated with each pair of areas. For the CAR model described above, taking
# Cij = 1 (equivalently Wij = 1/ ni) if areas i and j are neighbours and 0
# otherwise, gives a vector of 1's for weights[].
#
# num[] : A vector of length N (the total number of areas) giving the number of
# neighbours ni for each area.
#
# tau : A scalar argument representing the precision (inverse variance)
# parameter of the Gaussian CAR prior, or the inverse scale parameter of the
# Laplace prior (for the latter model, the variance = 2 / tau2).


library(R2OpenBUGS)#library(R2WinBUGS)
library(reshape)
library(ggplot2)
load('background/miller2019/DI_Data.Rdata') 
ebird <- sort_df(ebird, "cell")
WINE="/usr/local/Cellar/wine/4.0.2/bin/wine"
WINEPATH="/usr/local/Cellar/wine/4.0.2/bin/winepath"
OpenBUGS.pgm="/Users/tsz/Applications/OpenBUGS323/OpenBUGS.exe"








### ### ### ### ### ##

ni = 30000
nb = 10000
nt = 3
nc = 3


#Bundle data
  DI.data <- list (Y = rowSums(bba[,7:11]),
                    W = ebird$eCount,
                    effort= as.matrix(data.frame(ebird[,2:4])),
                    num = num, 
                    adj = adj,
                    weights = rep(1,length(adj)),
                    ncell = nrow(ebird),
                    cell = bba$cell,
                    nsite = nrow(bba),
                    forest = ebird$Forest1,
                    elev = ebird$Elev1,
                    R = matrix(c(0.02,0,0,0.02),2,2))

  #Set initial values

  DI.inits <- function() {
  list(z = as.numeric(DI.data$Y > 0),
       b.forest = 0,
       b.elev = 0,
       b.effort = rep(0.1,2),
       alpha = runif(2,-1,1),
       omega= matrix(c(.2,0,0,.2),2,2),
       S = matrix(c(rnorm(2*ncell,-0.01,0.01)),2,ncell))
}




  #parameters to monitor
  DI.parameters <- c('psi',"alpha",'lambda','corr',"p",
                     'b.forest','b.effort','b.elev','S')

  #Call winbugs from R
  out.corCARP <- bugs(data =  DI.data,
                      inits = DI.inits, 
                      parameters.to.save = DI.parameters, 
                      model.file = "corCARP.txt", 
                      n.chains = nc, 
                      n.thin = nt, 
                      n.iter = ni, 
                      n.burnin = nb, 
                      OpenBUGS.pgm=OpenBUGS.pgm, 
                      working.directory = paste0(getwd(), "/background/miller2019/"),
                      WINE=WINE, 
                      WINEPATH=WINEPATH,
                      useWINE=T,
                      DIC = FALSE)
  
plot.data <- data.frame(lat = block$Lat1, lon = block$Lon1,psi = out.corCARP$mean$psi)
ggplot(data = plot.data,aes(x = lon, y = lat, fill = psi)) + geom_point(shape = 22, size = 10) + 
  scale_fill_gradientn(colors = magma(100))
plot.data <- data.frame(lat = block$Lat1, lon = block$Lon1,lambda = out.corCARP$mean$lambda)
ggplot(data = plot.data,aes(x = lon, y = lat, fill = lambda)) + geom_point(shape = 22, size = 10) + 
  scale_fill_gradientn(colors = magma(100))

#Bundle data
DI.data <- list (Y = rowSums(bba[,7:11]),
                 W = as.numeric(ebird$eCount>0),
                 effort= as.matrix(data.frame(ebird[,2:4])),
                 num = num, 
                 adj = adj,
                 weights = rep(1,length(adj)),
                 ncell = nrow(ebird),
                 cell = bba$cell,
                 nsite = nrow(bba),
                 forest = ebird$Forest1,
                 elev = ebird$Elev1,
                 R = matrix(c(0.02,0,0,0.02),2,2))

#Set initial values

DI.inits <- function() {
  list(z = as.numeric(DI.data$Y > 0),
       b.forest = 0,
       b.elev = 0,
       b.effort = rep(0.1,2),
       alpha = runif(2,-1,1),
       omega= matrix(c(.2,0,0,.2),2,2),
       S = matrix(c(rnorm(2*ncell,-0.01,0.01)),2,ncell))
}


#parameters to monitor
DI.parameters <- c('psi','delta','corr',"alpha","p",'b.forest','b.effort','b.elev','S')

#Call winbugs from R
out.corCARB <- bugs(data =  DI.data,
                    inits = DI.inits, 
                    parameters.to.save = DI.parameters, 
                    model.file = "corCARB.txt", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni, 
                    n.burnin = nb, 
                    OpenBUGS.pgm=OpenBUGS.pgm, 
                    working.directory = paste0(getwd(), "/background/miller2019/"),
                    WINE=WINE, 
                    WINEPATH=WINEPATH,
                    useWINE=T,
                    debug = F,
                    DIC = FALSE)





plot.data <- data.frame(lat = block$Lat1, lon = block$Lon1,psi = out.corCARB$mean$psi)
ggplot(data = plot.data,aes(x = lon, y = lat, fill = psi)) + geom_point(shape = 22, size = 10) + 
  scale_fill_gradientn(colors = magma(100))
plot.data <- data.frame(lat = block$Lat1, lon = block$Lon1,delta = out.corCARB$mean$delta)
ggplot(data = plot.data,aes(x = lon, y = lat, fill = delta)) + geom_point(shape = 22, size = 10) + 
  scale_fill_gradientn(colors = magma(100))

#Bundle data
DI.data <- list (Y = rowSums(bba[,7:11]),
                 W = ebird$eCount,
                 effort= as.matrix(data.frame(ebird[,2:4])),
                 num = num, 
                 adj = adj,
                 weights = rep(1,length(adj)),
                 ncell = nrow(ebird),
                 cell = bba$cell,
                 nsite = nrow(bba),
                 forest = ebird$Forest1,
                 elev = ebird$Elev1)

#Set initial values

DI.inits <- function() {
  list(z = as.numeric(DI.data$Y > 0),
       b.forest = 0,
       b.elev = 0,
       b.eList = -0.01,
       b.eCount = 0.01,
       alpha.l = 0.5,
       spacesigma = 0.1,
       S = rnorm(ncell,0,0.1))
}


#parameters to monitor
DI.parameters <- c('psi',"alpha","p",'b.forest','b.eList','b.elev','b.eCount','S')

#Call winbugs from R
out.covCAR <- bugs(data =  DI.data,
                   inits = DI.inits, 
                   parameters.to.save = DI.parameters, 
                   model.file = "covCAR.txt", 
                   n.chains = nc, 
                   n.thin = nt, 
                   n.iter = ni, 
                   n.burnin = nb, 
                   OpenBUGS.pgm=OpenBUGS.pgm, 
                   working.directory = paste0(getwd(), "/background/miller2019/"),
                   WINE=WINE, 
                   WINEPATH=WINEPATH,
                   useWINE=T,
                   debug = F,
                   DIC = FALSE)



plot.data <- data.frame(lat = block$Lat1, lon = block$Lon1,psi = out.covCAR$mean$psi)
ggplot(data = plot.data,aes(x = lon, y = lat, fill = psi)) + geom_point(shape = 22, size = 10) + 
  scale_fill_gradientn(colors = magma(100))

