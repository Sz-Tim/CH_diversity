# Model experimentation

# Option 1: Citizen science model used to inform priors for structured sampling

# setup
library(boot); library(truncnorm); library(MASS); library(rstan); 
library(tidyverse); library(ggmcmc)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("code/00_fn.R")

opfo.dir <- "~/Documents/unil/opfo_str_sampling/data/"

################################################################################
# LAMBDA as true abundance

## make up parameters ##########################################################
S <- 50#139  # number of species
G <- 4#28
R <- 4  # number of regional covariates
Q <- 3  # number of W effort covariates
K <- list(W=100, Y=44, W_=20, Y_=10)  # number of cells
lambda_0 <- 10.5
effort_0 <- -11
h <- 0.00001885  # (25*6*pi*(0.2)^2)/(1000*1000)
sigma_b <- 0.5
tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)), rbinom(S, 1, 0.3)+1)
quadratic_R <- TRUE


## simulate data ###############################################################
# covariates, slopes, and true Lambda
X <- make_X(K, R, quadratic_R)
beta.ls <- make_slopes(lambda_0, R, G, S, tax_i, sigma_b, quadratic_R)
LAMBDA <- map(X, ~exp(. %*% beta.ls$b))

# W effort and bias
U <- cbind(1, matrix(rnorm((K$W+K$W_)*(Q-1)), nrow=K$W+K$W_, ncol=Q-1))
U_W <- U[1:K$W,]; U_W_ <- U[(1:K$W_)+K$W,]
eta <- cbind(c(effort_0, runif(Q-1)))
E <- inv.logit(U_W %*% eta); E_ <- inv.logit(U_W_ %*% eta)
D <- rtruncnorm(S, 0, b=Inf, tax_i[,3]-0.5, 0.5) # species-specific effort

# counts
obs <- map(K, ~matrix(NA, nrow=., ncol=S))
for(s in 1:S) {
  obs$W[,s] <- rpois(K$W, LAMBDA$W[,s]*E*D[s])
  obs$W_[,s] <- rpois(K$W_, LAMBDA$W_[,s]*E_*D[s])
  obs$Y[,s] <- rpois(K$Y, LAMBDA$Y[,s]*h)
  obs$Y_[,s] <- rpois(K$Y_, LAMBDA$Y_[,s]*h)
}


## run stan model ##############################################################
stan_data <- list(K_W=K$W, K_Y=K$Y, K_W_=K$W_, K_Y_=K$Y_, 
                  S=S, G=G, tax_i=tax_i, R=R, Q=Q, W=obs$W, Y=obs$Y, 
                  X_W=X$W, X_Y=X$Y, X_W_=X$W_, X_Y_=X$Y_, 
                  U_W=U_W, U_W_=U_W_, h=h,
                  LAMBDA_W_=LAMBDA$W_, LAMBDA_Y_=LAMBDA$Y_)

out.ls <- list(WY=stan(file="code/mods/00_PPM_mvPhy_WY_repar.stan",
                           data=stan_data, thin=10, warmup=2000, iter=3000),
               W=stan(file="code/mods/00_PPM_mvPhy_W_repar.stan",
                          data=stan_data, thin=10, warmup=2000, iter=3000),
               Y=stan(file="code/mods/00_PPM_mvPhy_Y_repar.stan",
                          data=stan_data, thin=10, warmup=2000, iter=3000))


## plot output: compare W, Y, WY ###############################################
mod_cols <- c(W="#e41a1c", Y="#377eb8", WY="#984ea3")
# beta
beta.out <- aggregate_beta(out.ls, beta.ls)
ggplot(beta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_vline(data=beta.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)

# b
b.out <- aggregate_b(out.ls, beta.ls, tax_i)
ggplot(b.out$gg, aes(x=true, colour=model)) + 
  facet_wrap(~Parameter, scales="free", ncol=S) +
  geom_vline(data=b.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)
ggplot(b.out$sum.gg, aes(mn-true, colour=model)) + 
  geom_vline(xintercept=0, linetype=2) + geom_density() + facet_wrap(~R) +
  scale_colour_manual(values=mod_cols)
ggplot(b.out$sum.gg, aes(true, mn, colour=model)) + 
  geom_abline() + geom_point() + stat_smooth(se=F, linetype=2, method="lm") + 
  facet_wrap(~R, scales="free") + scale_colour_manual(values=mod_cols)

# eta
eta.out <- aggregate_eta(out.ls, eta)
ggplot(eta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") +
  geom_vline(data=eta.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)

# D
D.out <- aggregate_D(out.ls, D)
ggplot(D.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") +
  geom_vline(data=D.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)

# Lambda
Lam.out <- aggregate_Lambda(out.ls, LAMBDA)
ggplot(Lam.out$sum.gg, aes(x=true, y=med, ymin=q025, ymax=q975, colour=model)) +
  geom_pointrange(shape=1, alpha=0.4) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  scale_x_log10() + scale_y_log10() + facet_grid(train~S, scales="free") + 
  ggtitle(expression(Lambda~vs.~hat(Lambda))) +
  scale_colour_manual(values=mod_cols)
ggplot(Lam.out$sum.gg, aes(x=log(mn)-log(true), colour=model)) + 
  geom_density() + geom_vline(xintercept=0, linetype=2) + 
  facet_wrap(~train) +
  scale_colour_manual(values=mod_cols)

# prPres
pP.out <- aggregate_prPres(out.ls, LAMBDA)
ggplot(pP.out$sum.gg, aes(x=true, y=mn, ymin=q025, ymax=q975, colour=model)) + 
  geom_pointrange(shape=1, alpha=0.4) + geom_abline() + 
  stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
  facet_grid(train~S) + 
  ggtitle(expression(prPres~vs.~hat(prPres))) +
  scale_colour_manual(values=mod_cols)
ggplot(pP.out$sum.gg, aes(x=mn-true, colour=model)) + 
  geom_density() + geom_vline(xintercept=0, linetype=2) + 
  facet_wrap(~train, scales="free") +
  scale_colour_manual(values=mod_cols)










## plot output: individual model run ###########################################
## beta ########################################################################
gg.beta <- ggs(out, "beta")
beta.df <- data.frame(value=c(beta.ls$beta),
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
b.df <- data.frame(value=c(beta.ls$b), 
                   Parameter=paste0("b[", rep(1:R, times=S), 
                                    ",", rep(1:S, each=(R)), "]"),
                   S=as.character(rep(1:S, each=R)),
                   G=as.character(rep(tax_i[,2], each=R)),
                   R=as.character(rep(1:(R), times=S)))
B.df <- data.frame(value=c(beta.ls$B), 
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
gg.D <- ggs(out, "^D\\[") 
D.df <- data.frame(value=D, 
                   Parameter=paste0("D[", 1:S, "]"))
ggs_density(gg.D) + geom_vline(xintercept=1, size=0.15, colour="gray") +
  geom_vline(data=D.df, aes(xintercept=value)) +
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







