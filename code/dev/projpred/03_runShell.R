# 03_runShell.R
# Tim Szewczyk
#
# This script creates an rstanarm model to be used as a shell for projpred


library(rstanarm); library(tidyverse)
d.f <- paste0("data/stan_data/vs_80")
d.ls <- readRDS(paste0(d.f, "_ls.rds"))
d.i <- readRDS(paste0(d.f, "_i.rds"))

# sampling parameters
nWarm <- 500
nSamp <- 20
nChain <- 12

# # filter species not detected in subset
# detected_spp <- which(colSums(d.ls$Y)>0)
# d.ls$S <- length(detected_spp)
# d.i$tax_i <- d.i$tax_i[detected_spp,]
# d.i$tax_i$sNum <- as.numeric(factor(d.i$tax_i$species))
# d.i$tax_i$gNum <- as.numeric(factor(d.i$tax_i$genus))
# d.ls$Y <- d.ls$Y[,colSums(d.ls$Y)>0]
# d.i$S <- d.ls$S
# d.ls$tax_i <- d.i$tax_i[,3:4]
# d.i$Y <- d.ls$Y
#
# saveRDS(d.ls, paste0(d.f, "_trim_ls.rds"))
# saveRDS(d.i, paste0(d.f, "_trim_i.rds"))
# rstan::stan_rdump(ls(d.ls),
#                   file=paste0(d.f, "_trim.Rdump"),
#                   envir=list2env(d.ls))



# make dataframe
Z <- cbind(d.ls$X[d.ls$IJ,], d.ls$V)
colnames(Z) <- paste0(rep(c("R_", "L_"), c(d.ls$R, d.ls$L)), colnames(Z))
Z_r <- do.call('rbind', map(1:d.ls$S, ~Z))
Y <- c(d.ls$Y)
stan_df <- as.data.frame(cbind(Y, Z_r, sp=rep(1:d.ls$S, each=d.ls$I)))


# run glmm with species random effects
fit <- stan_glmer(
  as.formula(
    paste("Y ~ ", paste(colnames(Z_r)[-1], collapse="+"), "+",
          paste("(1+", paste(colnames(Z_r)[-1], collapse="+"), "|sp)"))),
  data=stan_df,
  family=poisson("log"),
  prior=hs(),
  init_r=0.01,
  chains=nChain, cores=min(nChain, 4), 
  warmup=nWarm, iter=nWarm+nSamp, 
  refresh=10)

# store output
saveRDS(fit, "out/shell_80.rds")
