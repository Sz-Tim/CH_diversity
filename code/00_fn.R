# Assorted helper functions
# Tim Szewczyk

# Simulation functions
# Output processing functions


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
#' @param Lam_0 global intercept for Lambda; set to NULL for local env. slopes
#' @param nCov number of covariates, including the intercept
#' @param G number of genera
#' @param S number of species
#' @param tax_i matrix with taxonomic relationships
#' @param sd_sp standard deviation in responses within genera
#' @param quad logical to indicate whether X[,nCov] = (X[,nCov-1])^2
#' @param L_Omega NULL, or cholesky factor for G genera
#' @return named list with slopes for 'agg', 'gen', and 'sp'
make_slopes <- function(Lam_0, nCov, G, S, tax_i, sd_sp, quad, L_Omega=NULL) {
  # aggregate: alpha[nCov], beta[nCov]
  if(!is.null(Lam_0)) {
    agg <- cbind(c(Lam_0, rnorm(nCov-1, 0, 1)))
  } else {
    agg <- cbind(rnorm(nCov, 0, 1))
  }
  if(quad) agg[nCov] <- ifelse(agg[nCov] < 0, 2*agg[nCov], -2*agg[nCov])
  
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
  gg <- map_dfr(out.ls, ~ggs(., "LAMBDA"), .id="model")
  true <- map_dfr(Lam_sim, 
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