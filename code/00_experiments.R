# Model experimentation

# Option 1: Citizen science model used to inform priors for structured sampling

# setup
library(rstan); library(tidyverse); library(ggmcmc); theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

opfo.dir <- "~/Documents/unil/opfo_str_sampling/data/"

# Make up parameters
S <- 40  # number of species
G <- 5  # number of genera
K <- 3  # number of regional covariates
J <- 500  # number of grid cells
maxTubes <- 30
sigma_b <- 0.2
sigma_B <- 0.5


tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)))
beta <- cbind(rnorm(K+1))
B <- t(apply(beta, 1, function(x) rnorm(G, x, sigma_B)))
b <- matrix(ncol=S, nrow=K+1)
for(k in 1:(K+1)) {
  b[k,] <- rnorm(S, B[k,tax_i[,2]], sigma_b)
}
b <- t(apply(beta, 1, function(x) rnorm(S, x, .5)))
X <- cbind(1, matrix(rnorm(J*K), nrow=J, ncol=K))  # cell-scale covariates
E <- sample(1:maxTubes, J, TRUE)

# Intensities
LAMBDA <- exp(X %*% b)

# Observations
W <- matrix(NA, nrow=J, ncol=S)
for(j in 1:J) {
  # W[j,] <- rmultinom(1, E[j], LAMBDA[j,]) # probably more realistic...
  W[j,] <- rpois(S, LAMBDA[j,]*E[j])  # but the point process intensity is important
}

# Stan data
stan_data <- list(nCell=J, nSpp=S, nGen=G, K=K, W=W, E=E, X_W=X, tax_i=tax_i)

out <- stan(file="code/mods/00_experiment_W.stan", data=stan_data)
gg.out <- ggs(out, "beta")




b.df <- data.frame(value=c(b), 
                   Parameter=paste0("b[", rep(1:(K+1), times=S), 
                                    ",", rep(1:S, each=(K+1)), "]"),
                   S=as.character(rep(1:S, each=K+1)),
                   G=as.character(rep(tax_i[,2], each=K+1)),
                   K=as.character(rep(1:(K+1), times=S)))
B.df <- data.frame(value=c(B), 
                   Parameter=paste0("B[", rep(1:(K+1), times=G), 
                                    ",", rep(1:G, each=(K+1)), "]"),
                   G=as.character(rep(1:G, each=K+1)),
                   K=as.character(rep(1:(K+1), times=G)))
beta.df <- data.frame(value=c(beta),
                      Parameter=paste0("beta[", 1:length(beta), "]"),
                      K=as.character(1:(K+1)))
ggplot(b.df, aes(K, value)) + geom_point(aes(colour=G), alpha=0.5) + 
  geom_point(data=B.df, aes(colour=G), size=3, shape=1) + 
  geom_point(data=beta.df, colour="red", shape="-", size=10)
ggs_density(gg.out) + 
  geom_vline(data=beta.df, aes(xintercept=value), alpha=0.5) + 
  facet_wrap(~Parameter, scales="free_y")



lam.df <- data.frame(truth=c(LAMBDA),
                     Parameter=paste0("LAMBDA[", rep(1:J, times=S), 
                                      ",", rep(1:S, each=J), "]"),
                     J=as.character(rep(1:J, times=S)),
                     S=as.character(rep(1:S, each=J)))

gg.lam <- ggs(out, "LAMBDA") %>% 
  mutate(J=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 8),
         S=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"))
lam.sum <- gg.lam %>% group_by(Parameter, J, S) %>%
  summarise(mn=mean(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975))

ggplot(lam.sum, aes(S)) + geom_point(aes(y=mn)) + 
  geom_linerange(aes(ymin=q025, ymax=q975)) + 
  facet_wrap(~J, scales="free_y") + coord_flip() + 
  scale_y_log10() +
  geom_point(data=lam.df, aes(y=truth), colour="red", shape=1)




