# Assorted helper functions
# Tim Szewczyk

# Variable selection functions
# Simulation functions
# Output processing functions


#### variable selection functions ##############################################

#' Generate datasets with specified covariates
#' @param d.base Name & path for full dataset
#' @param d.new Name & path for new dataset
#' @param X_vars Vector of X variable names to include
#' @param V_vars Vector of V variable names to include
#' @param U_vars Vector of U variable names to include
#' @return Files created: d.new_ls.rds, d.new_i.rds, d.new.Rdump, d.new_vars.R
subset_vs_data <- function(d.base, d.new, X_vars, V_vars, U_vars) {
  library(rstan); library(tidyverse)
  
  d.ls <- readRDS(paste0(d.base, "_ls.rds"))
  d.i <- readRDS(paste0(d.base, "_i.rds"))
  
  X_cols <- which(colnames(d.ls$X) %in% X_vars)
  V_cols <- which(colnames(d.ls$V) %in% V_vars) 
  U_cols <- which(colnames(d.ls$U) %in% U_vars)
  
  d.ls$R <- length(X_cols) + 1
  d.ls$L <- length(V_cols)
  d.ls$Q <- length(U_cols) + 1
  
  d.ls$X <- d.ls$X[,c(1, X_cols), drop=F]
  d.ls$X_ <- d.ls$X_[,c(1, X_cols), drop=F]
  d.ls$V <- d.ls$V[,V_cols, drop=F]
  d.ls$V_ <- d.ls$V_[,V_cols, drop=F]
  d.ls$U <- d.ls$U[,c(1, U_cols), drop=F]
  
  d.i$X <- d.i$X[,c(1,2, X_cols+1), drop=F]
  d.i$X_ <- d.i$X_[,c(1,2, X_cols+1), drop=F]
  d.i$V <- d.i$V[,c(1,2, V_cols+2), drop=F]
  d.i$V_ <- d.i$V_[,c(1,2, V_cols+2), drop=F]
  d.i$U <- d.i$U[,c(1, U_cols), drop=F]
  
  saveRDS(d.ls, paste0(d.new, "_ls.rds"))
  saveRDS(d.i, paste0(d.new, "_i.rds"))
  rstan::stan_rdump(ls(d.ls),
                    file=paste0(d.new, ".Rdump"),
                    envir=list2env(d.ls))
  
  sink(paste0(d.new, "_vars.R"))
  cat("X_vars.i <-", paste0("c('", paste0(X_vars, collapse="', '"), "')"), "\n")
  cat("V_vars.i <-", paste0("c('", paste0(V_vars, collapse="', '"), "')"), "\n")
  cat("U_vars.i <-", paste0("c('", paste0(U_vars, collapse="', '"), "')"), "\n")
  sink()
  
  return(cat("Created", d.new, "with the following variables:\n",
             "X:", X_vars, "\n",
             "V:", V_vars, "\n",
             "U:", U_vars, "\n"))
}





#' Update rstanarm object with fitted hierarchical model output
#' @param fit_glmer rstanarm::glmer output
#' @param fit_hm stanfit object of hierarchical model output
#' @param d.ls List of data for fit_hm
#' @return Updated rstanarm object
subset_vs_data <- function(fit_glmer, fit_hm, d.ls) {
  library(tidyverse); library(rstan); library(rstanarm)
  
  R <- d.ls$R
  L <- d.ls$L
  S <- d.ls$S
  
  names_glmer <- names(fit_glmer$stanfit@sim$samples[[1]])
  names_hm <- names(fit_hm@sim$samples[[1]])
  name_lookup <- tibble(glmer=c("alpha[1]", 
                                paste0("beta[", 1:(R+L), "]"),
                                paste0("b[", 1:((R+L+1)*S), "]"),
                                "lp__"),
                        hm=c(paste0("beta[", 1:(R+L+1), "]"),
                             paste0(rep(paste0("b[", 1:(R+L+1)), times=S), ",", 
                                    rep(1:S, each=(R+L+1)), "]"),
                             "lp__"))
  
  
  # extract and align predicted values
  llam.summary <- rstan::summary(fit_hm, pars="llambda")$summary
  llam_hm.df <- tibble(llam=llam.summary[,1],
                       par=rownames(llam.summary),
                       site=str_split_fixed(str_remove(par, "llambda\\["), 
                                            ",", 2)[,1],
                       spp=str_split_fixed(str_remove(par, "]"), ",", 2)[,2]) %>%
    mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
    arrange(spp, site)
  
  
  # reference model: take structure of glmer and replace values with fitted hm
  fit_ref <- fit_glmer
  fit_ref$linear.predictors <- llam_hm.df$llam
  fit_ref$fitted.values <- exp(fit_ref$linear.predictors)
  iter_hm <- length(fit_hm@sim$samples[[1]][[1]])
  iter_glmer <- length(fit_glmer$stanfit@sim$samples[[1]][[1]])
  iter_index <- (iter_hm-iter_glmer+1):iter_hm 
  for(i in 1:nrow(name_lookup)) {
    for(j in 1:length(fit_glmer$stanfit@sim$samples)) {
      post_i <- fit_hm@sim$samples[[j]][[name_lookup$hm[i]]]
      fit_ref$stanfit@sim$samples[[j]][[name_lookup$glmer[i]]] <- post_i[iter_index]
    }
  }
  
  return(fit_ref)
}




#### simulation functions ######################################################

#' Simulate regional covariates
#' @param K names list with number of cells in W, W_, Y, Y_
#' @param R number of covariates, including the intercept
#' @param quad logical to indicate whether X[,R] = (X[,R-1])^2
#' @return named list with covariate matrices [K$., R] with elements 'all', 'W',
#'   'Y', "W_', and 'Y_'
make_X <- function(K, R, quad) {
  nTot <- sum(unlist(K))
  X <- cbind(1,
             matrix(rnorm(nTot*(R-1), 0, 1), nrow=nTot, ncol=R-1))
  if(quad) X[,R] <- (X[,R-1])^2
  return(list(all=X, 
              W=X[1:K$W,],
              Y=X[(1:K$Y)+K$W,],
              W_=X[(1:K$W_)+(K$W+K$Y),],
              Y_=X[(1:K$Y_)+(K$W+K$Y+K$W_),]))
}



#' Simulate local covariates
#' @param I names list with number of cells in Y, Y_
#' @param L number of covariates; no intercept for this regression
#' @return named list with covariate matrices [I$., L] with elements 'all', "Y',
#'   and 'Y_'
make_V <- function(I, L) {
  nTot <- sum(unlist(I)) 
  V <- matrix(rnorm(nTot*L, 0, 1), nrow=nTot, ncol=L)
  return(list(all=V,
              Y=V[1:I$Y,],
              Y_=V[(1:I$Y_)+I$Y,]))
}



#' Simulate phylogenetically structured slopes
#' @param agg_true NULL, or vector of true values
#' @param Lam_0 global intercept for Lambda; set to NULL for local env. slopes
#' @param nCov number of covariates, including the intercept
#' @param G number of genera
#' @param S number of species
#' @param tax_i matrix with taxonomic relationships
#' @param sd_sp standard deviation in responses within genera
#' @param quad logical to indicate whether X[,nCov] = (X[,nCov-1])^2
#' @param L_Omega NULL, or cholesky factor for G genera
#' @return named list with slopes for 'agg', 'gen', and 'sp'
make_slopes <- function(agg_true=NULL, Lam_0=NULL, nCov, G, S, tax_i, 
                        sd_sp, quad, L_Omega=NULL) {
  # aggregate: alpha[nCov], beta[nCov]
  if(is.null(agg_true)) {
    if(!is.null(Lam_0)) {
      agg <- cbind(c(Lam_0, rnorm(nCov-1, 0, 1)))
    } else {
      agg <- cbind(rnorm(nCov, 0, 1))
    }
    if(quad) agg[nCov] <- ifelse(agg[nCov] < 0, 2*agg[nCov], -2*agg[nCov])
  } else {
    agg <- cbind(agg_true)
  }
  
  # genus-level: A[G,nCov], B[G,nCov]
  if(is.null(L_Omega)) {
    phy_covMX <- matrix(runif(G^2)*2-1, ncol=G)/2 
    diag(phy_covMX) <- diag(phy_covMX)
    Sigma <- t(phy_covMX) %*% phy_covMX 
  } else {
    Sigma <- L_Omega %*% t(L_Omega)
  }
  gen <- t(apply(agg, 1, function(x) MASS::mvrnorm(1, rep(x,G), Sigma)))
  
  # species-level: a[S,nCov], b[S,nCov]
  sp <- matrix(ncol=S, nrow=nCov)
  for(r in 1:nCov) {
    sp[r,] <- rnorm(S, gen[r,tax_i[,2]], sd_sp)
  }
  return(list(agg=agg, gen=gen, sp=sp))
}






#### output processing functions ###############################################

aggregate_aggSlopes <- function(out.ls, slope.ls, var) {
  gg <- map_dfr(out.ls, ~ggs(., var), .id="model")
  true <- data.frame(value=c(slope.ls$agg),
                     Parameter=paste0(var, "[", 1:length(slope.ls$agg), "]"),
                     R=as.character(1:length(slope.ls$agg)))
  post_mns <- gg %>% group_by(model, Parameter) %>%
    summarise(mean=mean(value), median=median(value))
  return(list(gg=gg, true=true, post_mns=post_mns))
}



aggregate_spSlopes <- function(out.ls, slope.ls, tax_i, var) {
  R <- nrow(slope.ls$sp); S <- ncol(slope.ls$sp)
  gg <- map_dfr(out.ls, ~ggs(., paste0("^", var, "\\[")), .id="model") %>%
    mutate(R=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 3, -1),
           S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2))
  true <- data.frame(true=c(slope.ls$sp), 
                     Parameter=paste0(var, "[", rep(1:R, times=S), 
                                      ",", rep(1:S, each=R), "]"),
                     S=as.character(rep(1:S, each=R)),
                     G=as.character(rep(tax_i[,2], each=R)),
                     R=as.character(rep(1:R, times=S)))
  sum.gg <- gg %>% group_by(model, Parameter) %>%
    summarise(mn=mean(value), 
              med=median(value),
              q025=quantile(value, probs=0.025),
              q975=quantile(value, probs=0.975)) %>%
    full_join(., true, by="Parameter")
  return(list(gg=gg, sum.gg=sum.gg, true=true))
}



aggregate_eta <- function(out.ls, eta) {
  gg <- map_dfr(out.ls, ~ggs(., "^eta"), .id="model")
  true <- data.frame(value=eta,
                       Parameter=paste0("eta[", 1:length(eta), "]"))
  return(list(gg=gg, true=true))
}



aggregate_D <- function(out.ls, D) {
  gg <- map_dfr(out.ls, ~ggs(., "^D\\["), .id="model")
  true <- data.frame(value=D,
                     Parameter=paste0("D[", 1:length(D), "]"))
  return(list(gg=gg, true=true))
}



aggregate_Lambda <- function(out.ls, LAMBDA) {
  Lam_sim <- list(LAMBDA=rbind(LAMBDA$W, LAMBDA$Y),
                  LAMBDA_=rbind(LAMBDA$W_, LAMBDA$Y_))
  gg <- map_dfr(out.ls, ~ggs(., "lLAMBDA"), .id="model")
  true <- map_dfr(Lam_sim, 
                  ~tibble(true=c(.), 
                          K=as.character(rep(1:dim(.)[1], times=dim(.)[2])),
                          S=as.character(rep(1:dim(.)[2], each=dim(.)[1]))),
                  .id="dataset") %>%
    mutate(Parameter=paste0(dataset, "[", K, ",", S, "]"))
  sum.gg <- gg %>% group_by(model, Parameter) %>%
    summarise(mn=mean(exp(value)), 
              med=median(exp(value)),
              q025=quantile(exp(value), probs=0.025, na.rm=T),
              q975=quantile(exp(value), probs=0.975, na.rm=T),
              lmn=mean(value), 
              lmed=median(value),
              lq025=quantile(value, probs=0.025, na.rm=T),
              lq975=quantile(value, probs=0.975, na.rm=T)) %>%
    ungroup %>% mutate(Parameter=str_sub(Parameter, 2L, -1L)) %>%
    full_join(., true, by="Parameter") %>%
    mutate(train=c("train", "test")[grepl("_", dataset)+1])
  return(list(gg=gg, sum.gg=sum.gg, true=true))
}



aggregate_lambda <- function(out.ls, lambda) {
  gg <- map_dfr(out.ls, ~ggs(., "^lambda"), .id="model")
  true <- map_dfr(setNames(lambda, c("lambda", "lambda_")),
                   ~tibble(true=c(.), 
                           K=as.character(rep(1:dim(.)[1], times=dim(.)[2])),
                           S=as.character(rep(1:dim(.)[2], each=dim(.)[1]))),
                   .id="dataset") %>%
    mutate(Parameter=paste0(dataset, "[", K, ",", S, "]"))
  sum.gg <- gg %>% group_by(model, Parameter) %>%
    summarise(mn=mean(value), 
              med=median(value),
              q025=quantile(value, probs=0.025),
              q975=quantile(value, probs=0.975),
              lmn=mean(log(value)), 
              lmed=median(log(value)),
              lq025=quantile(log(value), probs=0.025),
              lq975=quantile(log(value), probs=0.975)) %>%
    full_join(., true, by="Parameter") %>%
    mutate(train=c("train", "test")[grepl("_", dataset)+1])
  return(list(gg=gg, sum.gg=sum.gg, true=true))
}



aggregate_prPres <- function(out.ls, LAMBDA) {
  gg <- map_dfr(out.ls, ~ggs(., "^prPres_"), .id="model")
  true <- map2_dfr(LAMBDA[-1], names(LAMBDA)[-1],
                   ~tibble(true=1-exp(-c(.x)), 
                           K=as.character(rep(1:dim(.x)[1], times=dim(.x)[2])),
                           S=as.character(rep(1:dim(.x)[2], each=dim(.x)[1]))),
                   .id="dataset") %>%
    mutate(Parameter=paste0("prPres_", dataset, "[", K, ",", S, "]"))
  sum.gg <- gg %>% group_by(model, Parameter) %>%
    summarise(mn=mean(value), 
              med=median(value),
              q025=quantile(value, probs=0.025),
              q975=quantile(value, probs=0.975)) %>%
    full_join(., true, by="Parameter") %>%
    mutate(train=c("train", "test")[grepl("_", dataset)+1])
  return(list(gg=gg, sum.gg=sum.gg, true=true))
}