library(tidyverse)
library(rstanarm)

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
fit <- stan_glmer(form.full, data=d.df, family=poisson("log"), 
                  chains=24, iter=525, refresh=1,
				  prior=hs(), warmup=500)
saveRDS(fit, "out/rstanarm_test.rds")