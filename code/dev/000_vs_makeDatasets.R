library(rstan); library(tidyverse)
source("code/00_fn.R")

X_vars <- c("MAT", "AP", "PwarmQ", "TAR",
            "lcH", "Forest", "Edge",
            "npp",
            "aspctN", "aspctE", 
            "bldgPeri", "rdLen",
            "MAT_sq", "AP_sq")[c(1, 8, 9)]
V_vars <- c("SoilTSt", "CnpyOpn", "CnpyMxd", 
            "Pasture", "Crop", 
            "Bare", "Litter",
            "VegTot")[c(1, 2, 5)]
U_vars <- c("rdLen", "bldgPeri")


d.base <- "data/vs/vs_full_30"


subset_vs_data(d.base, "data/vs/vs_test", X_vars, V_vars, U_vars)






#### Set up for rstanarm
library(tidyverse); library(rstanarm); library(bayesplot); library(projpred)

d.base <- "data/vs/vs_full_30"

d.full <- readRDS(paste0(d.base, "_ls.rds"))

Z <- do.call('rbind', map(1:d.full$S, ~cbind(d.full$X[d.full$IJ,-1], d.full$V)))
Z <- cbind(Z, sp=rep(1:d.full$S, each=d.full$I))

d.df <- as.data.frame(cbind(Y=c(d.full$Y), Z))

no_b <- c(1, ncol(Z), which(colnames(Z)=="lcH"), which(colnames(Z)=="npp"))

form.full <- as.formula(paste("Y ~ ",
                              paste(colnames(Z)[-c(1,ncol(Z))],
                                    collapse="+"), "+",
                              paste("(1+",
                                    paste(colnames(Z)[-no_b],
                                          collapse="+"),
                                    "|sp)")))

fit <- readRDS("out/rstanarm_test.rds")
plot(fit, pars="beta")
mcmc_areas(as.matrix(fit), 
           pars=colnames(Z)[colnames(Z) %in% colnames(as.matrix(fit))])
summary(fit, pars="beta")

vs <- varsel(fit, verbose=T, ndraws=1, ndraws_pred=1, nterms_max=2)

