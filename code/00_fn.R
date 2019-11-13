# functions for proccessing and aggregating output
# (just for a cleaner script)


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



make_slopes <- function(lam_0, R, G, S, tax_i, sd_b, quad) {
  # beta
  beta <- cbind(c(lam_0, rnorm(R-1, 0, 1)))
  if(quad) beta[R] <- ifelse(beta[R] < 0, 2*beta[R], -2*beta[R])
  # B
  phy_covMX <- matrix(runif(G^2)*2-1, ncol=G)/10 
  diag(phy_covMX) <- diag(phy_covMX)
  Sigma <- t(phy_covMX) %*% phy_covMX
  B <- t(apply(beta, 1, function(x) mvrnorm(1, rep(x,G), Sigma)))
  # b
  b <- matrix(ncol=S, nrow=R)
  for(r in 1:(R)) {
    b[r,] <- rnorm(S, B[r,tax_i[,2]], sd_b)
  }
  return(list(beta=beta, B=B, b=b))
}



aggregate_beta <- function(out.ls, beta.ls) {
  gg <- map_dfr(out.ls, ~ggs(., "beta"), .id="model")
  true <- data.frame(value=c(beta.ls$beta),
                     Parameter=paste0("beta[", 1:length(beta.ls$beta), "]"),
                     R=as.character(1:length(beta.ls$beta)))
  return(list(gg=gg, true=true))
}



aggregate_b <- function(out.ls, beta.ls, tax_i) {
  R <- nrow(beta.ls$b); S <- ncol(beta.ls$b)
  gg <- map_dfr(out.ls, ~ggs(., "^b\\["), .id="model") %>%
    mutate(R=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 3, -1),
           S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2))
  true <- data.frame(true=c(beta.ls$b), 
                     Parameter=paste0("b[", rep(1:R, times=S), 
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
  gg <- map_dfr(out.ls, ~ggs(., "LAMBDA_"), .id="model")
  true <- map2_dfr(LAMBDA[-1], names(LAMBDA)[-1],
                   ~tibble(true=c(.x), 
                           K=as.character(rep(1:dim(.x)[1], times=dim(.x)[2])),
                           S=as.character(rep(1:dim(.x)[2], each=dim(.x)[1]))),
                   .id="dataset") %>%
    mutate(Parameter=paste0("LAMBDA_", dataset, "[", K, ",", S, "]"))
  sum.gg <- gg %>% group_by(model, Parameter) %>%
    summarise(mn=mean(value), 
              med=median(value),
              q025=quantile(value, probs=0.025),
              q975=quantile(value, probs=0.975)) %>%
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