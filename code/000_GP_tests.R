
library(tidyverse)

dgenpois <- function(x, theta, lambda, log=F) {
  # p <- theta * ((theta + lambda*x)^(x-1)) * exp(-theta - lambda*x) / factorial(x)
  
  p <- exp(log(theta - theta*lambda) + 
    (x-1)*log(theta-theta*lambda + x*lambda) - 
    log(factorial(x)) + 
    theta*(lambda-1) - x*lambda)
  
  if(log) {
    return(log(p))
  } else {
    return(p) 
  }
}

rgenpois <- function(n, theta, lambda) {
  r_u <- runif(n)


}

hist(map_dbl(runif(1e6), ~sum(. > cumsum(LaplacesDemon::dgpois(0:10, 2, 0.17)))), freq=F, breaks=(0:15)-0.5)
lines(0:15, LaplacesDemon::dgpois(0:15, 2, 0.17))

x.seq <- 0:10  

dgenpois(x.seq, 1, 0.7)
LaplacesDemon::dgpois(x.seq, 1, 0.7)


theta <- c(1e-30, 0.07, 0.25, 1)
lambda <- c(0, 0.15, 0.18)
lam.cols <- c("#a6611a", "#dfc27d", "#80cdc1", "#018571")

lp <- F

par(mfrow=c(2,2))
for(i in theta) {
  plot(NA, NA, main=paste0("Theta: ", i), xlab="Counts", ylab="density",
       xlim=range(x.seq), 
       ylim=c(min(unlist(map(lambda, ~dgenpois(x.seq, i, .x, lp)))), 
              max(unlist(map(lambda, ~dgenpois(x.seq, i, .x, lp))))))
  
  for(j in seq_along(lambda)) {
    if(lambda[j] < 0) {
      m <- max(x.seq[(i + x.seq*lambda[j]) > 0])
      lam.lim <- max(-1, -theta/m)
    } else {
      m <- max(x.seq)
    }
    lines(x.seq, c(dgenpois(0:m, i, lambda[j], lp), rep(0, length(x.seq)-m-1)), 
          col=lam.cols[j])
  }
  lines(x.seq, dpois(x.seq, i, lp), lty=2)
  
  legend("topright", title="lambda", legend=lambda, col=lam.cols, lty=1, bty="n")
}


# theta > 0
# lambda ≤ 1
# lambda ≥ max(-1, -theta/m)


pmax(-1, -theta/5)


