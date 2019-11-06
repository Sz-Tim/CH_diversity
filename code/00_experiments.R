# Model experimentation

# Option 1: Citizen science model used to inform priors for structured sampling

# setup
library(boot); library(MASS); library(rstan); library(tidyverse); library(ggmcmc)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

opfo.dir <- "~/Documents/unil/opfo_str_sampling/data/"

################################################################################
# LAMBDA as true abundance

# Make up parameters
S <- 16#139  # number of species
G <- 4#28
R <- 4  # number of regional covariates
Q <- 3  # number of W effort covariates
K_W <- 100  # number of grid cells
K_Y <- 44
K_W_ <- 20
K_Y_ <- 10
lambda_0 <- 10.5
effort_0 <- -11
effort_prSoil <- 0.00001885
sigma_b <- 0.4
tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)))
quadratic_R <- TRUE

X <- cbind(1, matrix(rnorm((K_W+K_Y+K_W_+K_Y_)*(R-1)), 
                     nrow=K_W+K_Y+K_W_+K_Y_, ncol=(R-1)))
if(quadratic_R) X[,R] <- (X[,R-1])^2
X_W <- X[1:K_W,]
X_Y <- X[(1:K_Y)+K_W,]
X_W_ <- X[(1:K_W_)+(K_W+K_Y),]
X_Y_ <- X[(1:K_Y_)+(K_W+K_Y+K_W_),]

beta <- cbind(c(lambda_0, rnorm(R-1,0,1)))
if(quadratic_R) beta[R] <- ifelse(beta[R] < 0, 2*beta[R], -2*beta[R])

phy_covMX <- matrix(runif(G^2)*2-1, ncol=G)/10 
diag(phy_covMX) <- diag(phy_covMX)
Sigma <- t(phy_covMX) %*% phy_covMX
B <- t(apply(beta, 1, function(x) mvrnorm(1, rep(x,G), Sigma)))
b <- matrix(ncol=S, nrow=R)
for(r in 1:(R)) {
  b[r,] <- rnorm(S, B[r,tax_i[,2]], sigma_b)
}

LAMBDA_W <- exp(X_W %*% b)  
LAMBDA_Y <- exp(X_Y %*% b)  
LAMBDA_W_ <- exp(X_W_ %*% b)  
LAMBDA_Y_ <- exp(X_Y_ %*% b)  

V <- cbind(1, matrix(rnorm((K_W+K_W_)*(Q-1)), nrow=K_W+K_W_, ncol=Q-1))
V_W <- V[1:K_W,]
V_W_ <- V[(1:K_W_)+K_W,]
eta <- cbind(c(effort_0, runif(Q-1)))
E <- inv.logit(V_W %*% eta)
E_ <- inv.logit(V_W_ %*% eta)
spp_effort <- rbeta(S, 0.8, 1)*2 #runif(S, 0, 2)

W <- matrix(NA, nrow=K_W, ncol=S)
W_ <- matrix(NA, nrow=K_W_, ncol=S)
Y <- matrix(NA, nrow=K_Y, ncol=S)
Y_ <- matrix(NA, nrow=K_Y_, ncol=S)
for(s in 1:S) {
  W[,s] <- rpois(K_W, LAMBDA_W[,s]*E*spp_effort[s])
  W_[,s] <- rpois(K_W_, LAMBDA_W_[,s]*E_*spp_effort[s])
  Y[,s] <- rpois(K_Y, LAMBDA_Y[,s]*effort_prSoil)
  Y_[,s] <- rpois(K_Y_, LAMBDA_Y_[,s]*effort_prSoil)
}

stan_data <- list(nCell_W=K_W, nCell_Y=K_Y, nCell_W_=K_W_, nCell_Y_=K_Y_, 
                  nSpp=S, nGen=G, tax_i=tax_i, R=R, Q=Q, W=W, Y=Y, 
                  X_W=X_W, X_Y=X_Y, X_W_=X_W_, X_Y_=X_Y_, V_W=V_W, V_W_=V_W_,
                  effort_prSoil=effort_prSoil,
                  LAMBDA_W_=LAMBDA_W_, LAMBDA_Y_=LAMBDA_Y_)

out <- stan(file="code/mods/00_LAMBDA_mvPhy_WY_repar.stan",
            data=stan_data, chains=4, thin=5,
            # control=list(adapt_delta=0.99),
            warmup=2000, iter=3000)


## beta ########################################################################
gg.beta <- ggs(out, "beta")
beta.df <- data.frame(value=c(beta),
                      Parameter=paste0("beta[", 1:length(beta), "]"),
                      R=as.character(1:R))
ggs_density(gg.beta) + facet_wrap(~Parameter, scales="free") +
  geom_vline(data=beta.df, aes(xintercept=value), linetype=2)


## eta #########################################################################
gg.eta <- ggs(out, "^eta")
eta.df <- data.frame(value=eta,
                     Parameter=paste0("eta[", 1:Q, "]"))
ggs_density(gg.eta) + facet_wrap(~Parameter, scales="free") +
  geom_vline(data=eta.df, aes(xintercept=value), linetype=2)


## b ###########################################################################
gg.b <- ggs(out, "^b\\[")
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
ggs_density(gg.b) + facet_wrap(~Parameter, scales="free") +
  geom_vline(data=b.df, aes(xintercept=value), linetype=2)
ggplot(b.df, aes(R, value)) + geom_point(aes(colour=G), alpha=0.5) +
  geom_point(data=B.df, aes(colour=G), size=3, shape=1) +
  geom_point(data=beta.df, colour="red", shape="-", size=10)


## spEffort ####################################################################
gg.spEff <- ggs(out, "spp_effort") 
spEff.df <- data.frame(value=spp_effort, 
                       Parameter=paste0("spp_effort[", 1:S, "]"))
ggs_density(gg.spEff) + geom_vline(xintercept=1, size=0.15, colour="gray") +
  geom_vline(data=spEff.df, aes(xintercept=value)) +
  facet_wrap(~Parameter, scales="free_y")


## Facet scale setting #########################################################
sc <- c("fixed", "free", "free_x", "free_y")[2]


## LAMBDA ######################################################################
LAM.df <- rbind(data.frame(truth=c(LAMBDA_W), dataset="W",
                               Parameter=paste0("LAMBDA_W[", rep(1:K_W, times=S), 
                                                ",", rep(1:S, each=K_W), "]"),
                               K=as.character(rep(1:K_W, times=S)),
                               S=as.character(rep(1:S, each=K_W))),
                    data.frame(truth=c(LAMBDA_Y), dataset="Y",
                               Parameter=paste0("LAMBDA_Y[", rep(1:K_Y, times=S), 
                                                ",", rep(1:S, each=K_Y), "]"),
                               K=as.character(rep(1:K_Y, times=S)),
                               S=as.character(rep(1:S, each=K_Y))))
gg.LAM <- ggs(out, "LAMBDA_[A-Z]\\[") %>%
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 11, -1),
         S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
  full_join(., select(LAM.df, "truth", "Parameter", "dataset"), by="Parameter")
sum.LAM <- gg.LAM %>% group_by(Parameter) %>%
  summarise(mn=mean(value), 
            med=median(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975)) %>%
  full_join(., LAM.df, by="Parameter")
ggplot(sum.LAM, aes(x=truth, y=med, ymin=q025, ymax=q975, colour=dataset)) + 
  geom_pointrange(shape=1, alpha=0.7) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~S, scales=sc) + 
  ggtitle(expression(Lambda~vs.~hat(Lambda)~(train)))


## LAMBDA_ #####################################################################
LAM_.df <- rbind(data.frame(truth=c(LAMBDA_W_), dataset="W_",
                               Parameter=paste0("LAMBDA_W_[", rep(1:K_W_, times=S), 
                                                ",", rep(1:S, each=K_W_), "]"),
                               K=as.character(rep(1:K_W_, times=S)),
                               S=as.character(rep(1:S, each=K_W_))),
                    data.frame(truth=c(LAMBDA_Y_), dataset="Y_",
                               Parameter=paste0("LAMBDA_Y_[", rep(1:K_Y_, times=S), 
                                                ",", rep(1:S, each=K_Y_), "]"),
                               K=as.character(rep(1:K_Y_, times=S)),
                               S=as.character(rep(1:S, each=K_Y_))))
gg.LAM_ <- ggs(out, "LAMBDA_[A-Z]_") %>%
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 11, -1),
         S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
  full_join(., select(LAM_.df, "truth", "Parameter", "dataset"), by="Parameter")
sum.LAM_ <- gg.LAM_ %>% group_by(Parameter) %>%
  summarise(mn=mean(value), 
            med=median(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975)) %>%
  full_join(., LAM_.df, by="Parameter")
ggplot(sum.LAM_, aes(x=truth, y=med, ymin=q025, ymax=q975, colour=dataset)) + 
  geom_pointrange(shape=1, alpha=0.7) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~S, scales=sc) + 
  ggtitle(expression(Lambda~vs.~hat(Lambda)~(test)))


## W Y hat #####################################################################
WY.df <- rbind(data.frame(truth=c(W), dataset="What",
                          Parameter=paste0("What[", rep(1:K_W, times=S), 
                                           ",", rep(1:S, each=K_W), "]"),
                          K=as.character(rep(1:K_W, times=S)),
                          S=as.character(rep(1:S, each=K_W))),
               data.frame(truth=c(Y), dataset="Yhat",
                          Parameter=paste0("Yhat[", rep(1:K_Y, times=S), 
                                           ",", rep(1:S, each=K_Y), "]"),
                          K=as.character(rep(1:K_Y, times=S)),
                          S=as.character(rep(1:S, each=K_Y))))

gg.WYhat <- ggs(out, "^[A-Z]hat") %>%
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 6, -1),
         S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
  full_join(., select(WY.df, "truth", "Parameter"), by="Parameter")
sum.WYhat <- gg.WYhat %>% group_by(Parameter) %>%
  summarise(mn=mean(value), 
            med=median(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975)) %>%
  full_join(., WY.df, by="Parameter")
ggplot(sum.WYhat, aes(x=truth, y=med, ymin=q025, ymax=q975, colour=dataset)) + 
  geom_pointrange(shape=1, alpha=0.7) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  facet_wrap(~S, scales=sc) + ggtitle(expression(y~vs.~hat(y)~(train)))


## W Y _hat ####################################################################
WY_.df <- rbind(data.frame(truth=c(W_), dataset="W_hat",
                          Parameter=paste0("W_hat[", rep(1:K_W_, times=S), 
                                           ",", rep(1:S, each=K_W_), "]"),
                          K=as.character(rep(1:K_W_, times=S)),
                          S=as.character(rep(1:S, each=K_W_))),
               data.frame(truth=c(Y_), dataset="Y_hat",
                          Parameter=paste0("Y_hat[", rep(1:K_Y_, times=S), 
                                           ",", rep(1:S, each=K_Y_), "]"),
                          K=as.character(rep(1:K_Y_, times=S)),
                          S=as.character(rep(1:S, each=K_Y_))))

gg.WY_hat <- ggs(out, "^[A-Z]_hat") %>%
  mutate(K=str_sub(str_split_fixed(Parameter, ",", 2)[,1], 7, -1),
         S=str_sub(str_split_fixed(Parameter, ",", 2)[,2], 1, -2)) %>%
  full_join(., select(WY_.df, "truth", "Parameter"), by="Parameter")
sum.WY_hat <- gg.WY_hat %>% group_by(Parameter) %>%
  summarise(mn=mean(value), 
            med=median(value),
            q025=quantile(value, probs=0.025),
            q975=quantile(value, probs=0.975)) %>%
  full_join(., WY_.df, by="Parameter")
ggplot(sum.WY_hat, aes(x=truth, y=med, ymin=q025, ymax=q975, colour=dataset)) + 
  geom_pointrange(shape=1, alpha=0.7) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  facet_wrap(~S, scales=sc) + ggtitle(expression(y~vs.~hat(y)~(test)))
  

## Lambda vs. X ################################################################
lam_sim.df <- as.data.frame(rbind(LAMBDA_W, LAMBDA_Y, LAMBDA_W_, LAMBDA_Y_)) %>%
  mutate(K=row_number())
names(lam_sim.df)[1:S] <- paste0("S_", 1:S)
lam_sim.df <- gather(lam_sim.df, species, Lambda, 1:S)
lam_sim.df$G <- tax_i[as.numeric(str_sub(lam_sim.df$species, 3)),2]
X.df <- as.data.frame(X) %>%
  mutate(K=row_number())
names(X.df)[1:R] <- paste0("R_", 1:R)
sim.df <- full_join(lam_sim.df, X.df, by="K")
ggplot(sim.df, aes(R_2, log(Lambda), group=species, colour=factor(G))) + 
  geom_point(alpha=0.4) + 
  stat_smooth(method="lm", se=F, size=0.3, linetype=3) + facet_wrap(~species)
ggplot(sim.df, aes(R_3, log(Lambda), group=species, colour=factor(G))) + 
  geom_point(alpha=0.4) + 
  stat_smooth(method="lm", se=F, size=0.3, linetype=3) + facet_wrap(~species)
ggplot(sim.df, aes(R_4, log(Lambda), group=species, colour=factor(G))) + 
  geom_point(alpha=0.4) + 
  stat_smooth(method="lm", se=F, size=0.3, linetype=3) + facet_wrap(~species)










################################################################################
# Replicating Miller et al 2019

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







