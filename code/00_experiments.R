# Model experimentation

# Option 1: Citizen science model used to inform priors for structured sampling

# setup
library(rstan); library(tidyverse); library(ggmcmc); library(boot)
theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

opfo.dir <- "~/Documents/unil/opfo_str_sampling/data/"



########
## LAMBDA W
########

# Make up parameters
S <- 40  # number of species
G <- 5  # number of genera
R <- 3  # number of regional covariates
K <- 500  # number of grid cells
maxTubes <- 30
sigma_b <- 0.2
sigma_B <- 0.5


tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)))
beta <- cbind(rnorm(R))
B <- t(apply(beta, 1, function(x) rnorm(G, x, sigma_B)))
b <- matrix(ncol=S, nrow=R)
for(k in 1:(R)) {
  b[k,] <- rnorm(S, B[k,tax_i[,2]], sigma_b)
}
b <- t(apply(beta, 1, function(x) rnorm(S, x, .5)))
X <- cbind(1, matrix(rnorm(K*(R-1)), nrow=K, ncol=(R-1)))  # cell-scale covariates
E <- sample(1:maxTubes, K, TRUE)

# Intensities
LAMBDA <- exp(X %*% b)

# Observations
W <- matrix(NA, nrow=K, ncol=S)
for(k in 1:K) {
  # W[k,] <- rmultinom(1, E[k], LAMBDA[k,]) # probably more realistic...
  W[k,] <- rpois(S, LAMBDA[k,]*E[k])  # but the point process intensity is important
}

# Stan data
stan_data <- list(nCell=K, nSpp=S, nGen=G, R=R, W=W, E=E, X_W=X, tax_i=tax_i)

out <- stan(file="code/mods/00_LAMBDA_W.stan", data=stan_data)
gg.out <- ggs(out, "beta")




b.df <- data.frame(value=c(b), 
                   Parameter=paste0("b[", rep(1:R, times=S), 
                                    ",", rep(1:S, each=(R)), "]"),
                   S=as.character(rep(1:S, each=R)),
                   G=as.character(rep(tax_i[,2], each=R)),
                   R=as.character(rep(1:(R), times=S)))
B.df <- data.frame(value=c(B), 
                   Parameter=paste0("B[", rep(1:R, times=G), 
                                    ",", rep(1:G, each=(R)), "]"),
                   G=as.character(rep(1:G, each=R)),
                   R=as.character(rep(1:R, times=G)))
beta.df <- data.frame(value=c(beta),
                      Parameter=paste0("beta[", 1:length(beta), "]"),
                      R=as.character(1:R))
ggplot(b.df, aes(R, value)) + geom_point(aes(colour=G), alpha=0.5) + 
  geom_point(data=B.df, aes(colour=G), size=3, shape=1) + 
  geom_point(data=beta.df, colour="red", shape="-", size=10)
ggs_density(gg.out) + 
  geom_vline(data=beta.df, aes(xintercept=value), alpha=0.5) + 
  facet_wrap(~Parameter, scales="free_y")



lam.df <- data.frame(truth=c(LAMBDA),
                     Parameter=paste0("LAMBDA[", rep(1:K, times=S), 
                                      ",", rep(1:S, each=K), "]"),
                     K=as.character(rep(1:K, times=S)),
                     S=as.character(rep(1:S, each=K)))

gg.lam <- ggs(out, "LAMBDA") %>% 
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 8),
         S=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"))
lam.sum <- gg.lam %>% group_by(Parameter, K, S) %>%
  summarise(mn=mean(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975))

ggplot(lam.sum, aes(S)) + geom_point(aes(y=mn)) + 
  geom_linerange(aes(ymin=q025, ymax=q975)) + 
  facet_wrap(~K, scales="free_y") + coord_flip() + 
  scale_y_log10() +
  geom_point(data=lam.df, aes(y=truth), colour="red", shape=1)










########
## PSI W
########

# Make up parameters
S <- 40  # number of species
G <- 5  # number of genera
R <- 3  # number of regional covariates
K <- 200  # number of grid cells
maxTubes <- 30
sigma_b <- 0.2
sigma_B <- 0.5


p <- runif(S)
tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)))
beta <- cbind(rnorm(R))
B <- t(apply(beta, 1, function(x) rnorm(G, x, sigma_B)))
b <- matrix(ncol=S, nrow=R)
for(k in 1:(R)) {
  b[k,] <- rnorm(S, B[k,tax_i[,2]], sigma_b)
}
b <- t(apply(beta, 1, function(x) rnorm(S, x, .5)))
X <- cbind(1, matrix(rnorm(K*(R-1)), nrow=K, ncol=(R-1)))  # cell-scale covariates
E <- sample(1:maxTubes, K, TRUE)

# Presences
PSI <- inv.logit(X %*% b)

# Observations
W <- matrix(NA, nrow=K, ncol=S)
for(s in 1:S) {
  for(k in 1:K) {
    W[k,s] <- rbinom(1, 1, PSI[k,s]*(1-(1-p[s])^E[k]))  
  }
}

# Stan data
stan_data <- list(nCell=K, nSpp=S, nGen=G, R=R, W=W, E=E, X_W=X, tax_i=tax_i)

out <- stan(file="code/mods/00_PSI_W.stan", data=stan_data)
gg.out <- ggs(out, "beta")




b.df <- data.frame(value=c(b), 
                   Parameter=paste0("b[", rep(1:R, times=S), 
                                    ",", rep(1:S, each=(R)), "]"),
                   S=as.character(rep(1:S, each=R)),
                   G=as.character(rep(tax_i[,2], each=R)),
                   R=as.character(rep(1:(R), times=S)))
B.df <- data.frame(value=c(B), 
                   Parameter=paste0("B[", rep(1:R, times=G), 
                                    ",", rep(1:G, each=(R)), "]"),
                   G=as.character(rep(1:G, each=R)),
                   R=as.character(rep(1:R, times=G)))
beta.df <- data.frame(value=c(beta),
                      Parameter=paste0("beta[", 1:length(beta), "]"),
                      R=as.character(1:R))
ggplot(b.df, aes(R, value)) + geom_point(aes(colour=G), alpha=0.5) + 
  geom_point(data=B.df, aes(colour=G), size=3, shape=1) + 
  geom_point(data=beta.df, colour="red", shape="-", size=10)
ggs_density(gg.out) + 
  geom_vline(data=beta.df, aes(xintercept=value), alpha=0.5) + 
  facet_wrap(~Parameter, scales="free_y")



psi.df <- data.frame(truth=c(PSI),
                     Parameter=paste0("PSI[", rep(1:K, times=S), 
                                      ",", rep(1:S, each=K), "]"),
                     K=as.character(rep(1:K, times=S)),
                     S=as.character(rep(1:S, each=K)))

gg.psi <- ggs(out, "PSI") %>% 
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 5),
         S=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"))
psi.sum <- gg.psi %>% group_by(Parameter, K, S) %>%
  summarise(mn=mean(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975))

ggplot(psi.sum, aes(S)) + geom_point(aes(y=mn)) + 
  geom_linerange(aes(ymin=q025, ymax=q975)) + 
  facet_wrap(~K, scales="free_y") + coord_flip() + 
  geom_point(data=psi.df, aes(y=truth), colour="red", shape=1)
ggplot(psi.sum, aes(K)) + geom_point(aes(y=mn)) + 
  geom_linerange(aes(ymin=q025, ymax=q975)) + 
  facet_wrap(~S, scales="free_y") + coord_flip() + 
  geom_point(data=psi.df, aes(y=truth), colour="red", shape=1)



p.df <- data.frame(truth=p,
                   Parameter=paste0("p[", 1:S, "]"),
                   S=as.character(1:S))
gg.p <- ggs(out, "p") 
p.sum <- gg.p %>% mutate(S=str_sub(Parameter, 3, -2L)) %>%
  group_by(Parameter, S) %>%
  summarise(mn=mean(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975))
ggs_density(gg.p) + 
  geom_vline(data=p.df, aes(xintercept=truth), alpha=0.5) + 
  facet_wrap(~Parameter, scales="free_y")
ggplot(p.sum, aes(S)) + geom_point(aes(y=mn)) +
  geom_linerange(aes(ymin=q025, ymax=q975)) + coord_flip() + 
  geom_point(data=p.df, aes(y=truth), colour="red", shape=1)




