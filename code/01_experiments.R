# Model experimentation

# Option 1: Citizen science model used to inform priors for structured sampling

# setup
library(boot); library(truncnorm); library(MASS); library(rstan); library(loo);
library(tidyverse); library(ggmcmc)
theme_set(theme_bw() + theme(panel.grid=element_blank()))
options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)
source("code/00_fn.R")

################################################################################
# Simulations based on real data
# LAMBDA, lambda as true abundance
# 1 km2 grid cells and 0.75 m2 soil plots

true.ls <- read_rdump("data/stan_data/test_realData.Rdump")
p.fit <- read_csv("out/test_realWY.csv") %>% rename(Parameter=X1)
tax_i <- read_csv("data/tax_i.csv") %>% filter(., sNum<=true.ls$S) %>% 
  dplyr::select(sNum, gNum, Dprior) %>% as.matrix
tax.df <- tax_i <- read_csv("data/tax_i.csv") %>% filter(., sNum<=true.ls$S) 

# slopes
# beta.ls <- make_slopes(agg_true=p.fit$mean[grepl("beta", p.fit$Parameter)],
#                        nCov=true.ls$R, G=true.ls$G, S=true.ls$S, tax_i=tax_i,
#                        sd_sp=p.fit$mean[grepl("sigma_b", p.fit$Parameter)],
#                        L_Omega=matrix(p.fit$mean[grepl("L_Omega_B",
#                                                        p.fit$Parameter)],
#                                       ncol=true.ls$G, byrow=T))
# alpha.ls <- make_slopes(agg_true=p.fit$mean[grepl("alpha", p.fit$Parameter)],
#                        nCov=true.ls$L, G=true.ls$G, S=true.ls$S, tax_i=tax_i,
#                        sd_sp=p.fit$mean[grepl("sigma_a", p.fit$Parameter)],
#                        L_Omega=matrix(p.fit$mean[grepl("L_Omega_A",
#                                                        p.fit$Parameter)],
#                                       ncol=true.ls$G, byrow=T))
beta.ls <- list(agg=p.fit$mean[grepl("beta", p.fit$Parameter)],
                gen=matrix(p.fit$mean[grepl("^B\\[", p.fit$Parameter)],
                           ncol=true.ls$G, byrow=T),
                sp=matrix(p.fit$mean[grepl("^b\\[", p.fit$Parameter)], 
                          ncol=true.ls$S, byrow=T))
alpha.ls <- list(agg=p.fit$mean[grepl("alpha", p.fit$Parameter)],
                gen=matrix(p.fit$mean[grepl("^A\\[", p.fit$Parameter)],
                           ncol=true.ls$G, byrow=T),
                sp=matrix(p.fit$mean[grepl("^a\\[", p.fit$Parameter)], 
                          ncol=true.ls$S, byrow=T))

# predicted intensities
LAMBDA <- exp(true.ls$X %*% beta.ls$sp)  # rbind(W, Y.km2)
LAMBDA_ <- exp(true.ls$X_ %*% beta.ls$sp)  # rbind(W_, Y_.km2)
lambda <- exp(true.ls$V %*% alpha.ls$sp + 
                log(true.ls$h*LAMBDA[true.ls$IJ + true.ls$K,]))
lambda_ <- exp(true.ls$V_ %*% alpha.ls$sp + 
                log(true.ls$h*LAMBDA_[true.ls$IJ_ + true.ls$K_,]))

# W effort and bias
E <- inv.logit(true.ls$U %*% p.fit$mean[grepl("^eta", p.fit$Parameter)])
E_ <- inv.logit(true.ls$U_ %*% p.fit$mean[grepl("^eta", p.fit$Parameter)])
D <- p.fit$mean[grepl("D", p.fit$Parameter)]

# counts
obs <- list(W=matrix(NA, nrow=true.ls$K, ncol=true.ls$S), 
            W_=matrix(NA, nrow=true.ls$K_, ncol=true.ls$S),
            Yj=matrix(NA, nrow=true.ls$J, ncol=true.ls$S),
            Yj_=matrix(NA, nrow=true.ls$J_, ncol=true.ls$S),
            Y=matrix(NA, nrow=true.ls$I, ncol=true.ls$S), 
            Y_=matrix(NA, nrow=true.ls$I_, ncol=true.ls$S))
for(s in 1:true.ls$S) {
  obs$W[,s] <- rpois(true.ls$K, LAMBDA[1:true.ls$K,s]*E*D[s])
  obs$W_[,s] <- rpois(true.ls$K_, LAMBDA_[1:true.ls$K_,s]*E_*D[s])
  obs$Yj[,s] <- rpois(true.ls$J, LAMBDA[(1:true.ls$J)+true.ls$K,s])
  obs$Yj_[,s] <- rpois(true.ls$J_, LAMBDA_[(1:true.ls$J_)+true.ls$K_,s])
  obs$Y[,s] <- rpois(true.ls$I, lambda[,s])
  obs$Y_[,s] <- rpois(true.ls$I_, lambda_[,s])
}
map(obs, ~summary(c(.)))
map(obs, ~sum(c(.)>0)/length(c(.)))

ShannonH <- prop.table(LAMBDA, 1) %>%
  apply(., 1, function(x) -sum(x * log(x)))
ShannonH_ <- prop.table(LAMBDA_, 1) %>%
  apply(., 1, function(x) -sum(x * log(x)))
ShannonH.obs <- prop.table(rbind(obs$W, obs$Yj), 1) %>%
  apply(., 1, function(x) -sum(x * log(x)))
ShannonH.obs_ <- prop.table(rbind(obs$W_, obs$Yj_), 1) %>%
  apply(., 1, function(x) -sum(x * log(x)))

plot(p.fit$mean[grepl("ShannonH_\\[", p.fit$Parameter)], ShannonH_); abline(0,1)
plot(true.ls$X_[,2], p.fit$mean[grepl("ShannonH_\\[", p.fit$Parameter)])


p.fit %>% filter(grepl("^a\\[", Parameter)) %>% 
  mutate(covariate=str_sub(Parameter, 3,3), 
         species=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"), 
         genus=as.character(tax_i[,2][as.numeric(species)]),
         genName=tax.df$genus[as.numeric(species)]) %>% 
  ggplot(aes(x=species, y=mean, ymin=`2.5%`, ymax=`97.5%`)) + 
  geom_point(alpha=0.5) + geom_linerange(alpha=0.5) + 
  facet_grid(covariate~genName, scales="free", space="free_x") + 
  geom_hline(yintercept=0, linetype=2, size=0.5)
p.fit %>% filter(grepl("^b\\[", Parameter)) %>% 
  mutate(covariate=str_sub(Parameter, 3,3), 
         species=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"), 
         genus=as.character(tax_i[,2][as.numeric(species)]),
         genName=tax.df$genus[as.numeric(species)]) %>% 
  ggplot(aes(x=species, y=mean, ymin=`2.5%`, ymax=`97.5%`)) + 
  geom_point(alpha=0.5) + geom_linerange(alpha=0.5) + 
  facet_grid(covariate~genName, scales="free", space="free_x") + 
  geom_hline(yintercept=0, linetype=2, size=0.5)
p.fit %>% filter(grepl("^a\\[", Parameter)) %>% 
  mutate(covariate=str_sub(Parameter, 3,3), 
         species=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"), 
         genus=as.character(tax_i[,2][as.numeric(species)])) %>% 
  ggplot(aes(x=mean)) + geom_histogram() + facet_wrap(~covariate, scales="free")
p.fit %>% filter(grepl("^b\\[", Parameter)) %>% 
  mutate(covariate=str_sub(Parameter, 3,3), 
         species=str_remove(str_split_fixed(Parameter, ",", 2)[,2], "]"), 
         genus=as.character(tax_i[,2][as.numeric(species)])) %>% 
  ggplot(aes(x=mean)) + geom_histogram() + facet_wrap(~covariate, scales="free")



d.trueSim <- list(K=true.ls$K, J=true.ls$J, K_=true.ls$K_, J_=true.ls$J_, 
                  IJ=true.ls$IJ, IJ_=true.ls$IJ_, I=true.ls$I, I_=true.ls$I_,
                  S=true.ls$S, G=true.ls$G, tax_i=true.ls$tax_i, 
                  D_prior=true.ls$D_prior,
                  R=true.ls$R, L=true.ls$L, Q=true.ls$Q, 
                  W=obs$W, W_=obs$W_, Y=obs$Y, Y_=obs$Y_,
                  X=true.ls$X, X_=true.ls$X_, 
                  V=true.ls$V, V_=true.ls$V_, U=true.ls$U, U_=true.ls$U_, 
                  h=true.ls$h)

out.ls <- list(
  WY=stan(file="code/mods/PPM_mvPhy_WY_IJK.stan", chains=2,
          data=d.trueSim, thin=5, warmup=500, iter=1000, init=0),
  W=stan(file="code/mods/PPM_mvPhy_W_IJK.stan", chains=2,
         data=d.trueSim, thin=5, warmup=500, iter=1000, init=0),
  Y=stan(file="code/mods/PPM_mvPhy_Y_IJK.stan", chains=2,
         data=d.trueSim, thin=5, warmup=500, iter=1000, init=0)
)




################################################################################
# LAMBDA, lambda as true abundance
# 1 km2 grid cells and 0.75 m2 soil plots

## make up parameters ##########################################################
S <- 16  # number of species
R <- 3  # number of regional covariates
L <- 2
Q <- 2  # number of W effort covariates
K <- list(W=200, Y=35, W_=2, Y_=9)  # number of cells (W: 1252, Y: 44)
IJ <- list(Y=rep(1:K$Y, each=25), Y_=rep(1:K$Y_, each=25))
I <- list(Y=length(IJ$Y), Y_=length(IJ$Y_))  # number of soil plots
Lambda_0 <- 6.56#10.5
effort_0 <- -12.16#-11
h <- 7.5e-7  # (6*pi*(0.2)^2)/(1000*1000)
sigma_a <- 0.33
sigma_b <- 0.93
tax_i <- read_csv("data/tax_i.csv") %>% filter(., sNum<=S) %>% 
  dplyr::select(sNum, gNum, Dprior) %>% as.matrix
G <- max(tax_i[,2])
quadratic_R <- TRUE


## simulate data ###############################################################
# covariates, slopes, true Lambda, true lambda
X <- make_X(K, R, quadratic_R)
V <- make_V(I, L)
beta.ls <- make_slopes(NULL, Lambda_0, R, G, S, tax_i, sigma_b, quadratic_R)
alpha.ls <- make_slopes(NULL, NULL, L, G, S, tax_i, sigma_a, FALSE)
LAMBDA <- map(X, ~exp(. %*% beta.ls$sp))
lambda <- map(names(IJ), ~exp(V[[.]] %*% alpha.ls$sp + log(h*LAMBDA[[.]][IJ[[.]],])))
names(lambda) <- names(IJ)

# W effort and bias
U <- cbind(1, matrix(rnorm((K$W+K$W_)*(Q-1)), nrow=K$W+K$W_, ncol=Q-1))
U_W <- U[1:K$W,]; U_W_ <- U[(1:K$W_)+K$W,]
eta <- cbind(c(effort_0, runif(Q-1)))
E <- inv.logit(U_W %*% eta); E_ <- inv.logit(U_W_ %*% eta)
D <- rtruncnorm(S, 0, b=Inf, tax_i[,3]-0.5, 0.5) # species-specific effort

# counts
obs <- list(W=matrix(NA, nrow=K$W, ncol=S), W_=matrix(NA, nrow=K$W_, ncol=S),
            Y=matrix(NA, nrow=I$Y, ncol=S), Y_=matrix(NA, nrow=I$Y_, ncol=S))
for(s in 1:S) {
  obs$W[,s] <- rpois(K$W, LAMBDA$W[,s]*E*D[s])
  obs$W_[,s] <- rpois(K$W_, LAMBDA$W_[,s]*E_*D[s])
  obs$Y[,s] <- rpois(I$Y, lambda$Y[,s])
  obs$Y_[,s] <- rpois(I$Y_, lambda$Y_[,s])
}
map(obs, ~summary(c(.)))

ShannonH <- map(LAMBDA[-1], ~prop.table(., 1)) %>%
  map(~apply(., 1, function(x) -sum(x * log(x))))

## run stan model ##############################################################
stan_data <- list(K=K$W, J=K$Y, K_=K$W_, J_=K$Y_, 
                  IJ=IJ$Y, IJ_=IJ$Y_, I=I$Y, I_=I$Y_,
                  S=S, G=G, tax_i=tax_i[,-3], D_prior=tax_i[,3], R=R, L=L, Q=Q, 
                  W=obs$W, Y=obs$Y, Y_=obs$Y_,
                  X=rbind(X$W, X$Y), X_=rbind(X$W_, X$Y_), 
                  # X_W=X$W, X_Y=X$Y, X_W_=X$W_, X_Y_=X$Y_,
                  V=V$Y, V_=V$Y_, U=U_W, U_=U_W_, h=h)

out.ls <- list(
  WY=stan(file="code/mods/PPM_mvPhy_WY_IJK.stan",
          data=stan_data, thin=5, warmup=1500, iter=2000)#,
  # W=stan(file="code/mods/PPM_mvPhy_W_IJK.stan",
  #        data=stan_data, thin=5, warmup=1500, iter=2000),
  # Y=stan(file="code/mods/PPM_mvPhy_Y_IJK.stan",
  #        data=stan_data, thin=5, warmup=1500, iter=2000)
)


## plot output: compare W, Y, WY ###############################################
mod_cols <- c(W="#e41a1c", Y="#377eb8", WY="#984ea3", "WY_con"="black")

# lpd (Y_ | lambda_)
gg.ll <- map_dfr(out.ls, ~ggs(., "log_lik_lambda_"), .id="model")
gg.ll.lpd <- gg.ll %>% group_by(model, Parameter) %>% 
  summarise(mn=log(mean(exp(value)))) %>% 
  group_by(model) %>% 
  summarise(lpd=-2*sum(mn)/n())
gg.ll.lpd

library(loo)
map(out.ls, ~loo::waic(loo::extract_log_lik(., "log_lik_lambda_"))) %>%
  loo::loo_compare()
map(out.ls, ~loo::loo(loo::extract_log_lik(., "log_lik_lambda_"))) %>%
  loo::loo_compare()



H.out <- ggs(out.ls$W, "ShannonH") %>%
  mutate(site=str_remove(str_split_fixed(Parameter, "\\[", n=2)[,2], "]"),
         site=as.numeric(site),
         par=str_split_fixed(Parameter, "\\[", n=2)[,1]) %>%
  arrange(par, site)

H.out <- rbind(H.out %>% filter(par=="ShannonH") %>%
                 mutate(true=c(ShannonH$W, ShannonH$Y)[site]),
               H.out %>% filter(par=="ShannonH_") %>%
                 mutate(true=c(ShannonH$W_, ShannonH$Y_)[site]))
H.out$diff <- H.out$value - H.out$true

ggs_caterpillar(H.out) + facet_grid(.~par, scales="free_y") + 
  geom_point(aes(x=true, y=Parameter), colour="red", size=0.5, shape=1)

ggplot(H.out, aes(x=diff)) + geom_density() 
ggplot(H.out, aes(true, diff, colour=site)) + geom_point(alpha=0.05) +
  facet_wrap(~par) 

ggplot(H.out, aes(true, value, colour=site)) + geom_point(alpha=0.05) +
  geom_abline(colour="red", linetype=2) + facet_wrap(~par) + ylim(0,3)
H.out %>% group_by(Parameter, par, site) %>% 
  summarise(true=mean(true), postmn=mean(value)) %>%
  ggplot(aes(true, postmn, colour=site)) +  geom_point(alpha=0.5) +
  geom_abline(colour="red", linetype=2) + facet_wrap(~par) + ylim(0,3)



# beta
beta.out <- aggregate_aggSlopes(out.ls, beta.ls, "beta")
ggplot(beta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_rug(data=beta.out$post_mns, aes(x=mean), sides="b") +
  geom_vline(data=beta.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols) 
ggsave("eda/beta_density.pdf", width=8, height=3, units="in")

# alpha
alpha.out <- aggregate_aggSlopes(out.ls, alpha.ls, "alpha")
ggplot(alpha.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_rug(data=alpha.out$post_mns, aes(x=mean), sides="b") +
  geom_vline(data=alpha.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)
ggsave("eda/alpha_density.pdf", width=7, height=3, units="in")

# b
b.out <- aggregate_spSlopes(out.ls, beta.ls, tax_i, "b")
b.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  summarise(mnDiff=sqrt(mean(Diff^2)), seDiff=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(model, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff)) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_pointrange(aes(colour=model)) + 
  labs(title="b", x="Model", y="RMSE among species") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/b_RMSE.pdf", width=8, height=3, units="in")
b.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  ggplot(aes(model, y=Diff)) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_boxplot(aes(colour=model)) + 
  labs(title="b", x="Model", y="Error") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/b_box.pdf", width=8, height=3, units="in")
b.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  ggplot(aes(model, y=abs(Diff))) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_boxplot(aes(colour=model)) + 
  labs(title="b", x="Model", y="abs(Error)") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/b_box_absErr.pdf", width=8, height=3, units="in")

# a
a.out <- aggregate_spSlopes(out.ls, alpha.ls, tax_i, "a")
a.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  summarise(mnDiff=sqrt(mean(Diff^2)), seDiff=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(model, y=mnDiff, ymin=mnDiff-2*seDiff, ymax=mnDiff+2*seDiff)) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_pointrange(aes(colour=model)) + 
  labs(title="a", x="Model", y="RMSE among species") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/a_RMSE.pdf", width=7, height=3, units="in")
a.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  ggplot(aes(model, y=Diff)) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_boxplot(aes(colour=model)) + 
  labs(title="a", x="Model", y="Error") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/a_box.pdf", width=7, height=3, units="in")
a.out$sum.gg %>% group_by(model, R) %>% mutate(Diff=mn-true) %>% 
  ggplot(aes(model, y=abs(Diff))) +
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~R) +
  geom_boxplot(aes(colour=model)) + 
  labs(title="a", x="Model", y="abs(Error)") + 
  scale_colour_manual(values=mod_cols)
ggsave("eda/a_box_absErr.pdf", width=7, height=3, units="in")

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

# Lambda & lambda
Lam.out <- aggregate_Lambda(out.ls, LAMBDA)
lam.out <- aggregate_lambda(out.ls, lambda)

Lam.out$sum.gg %>% mutate(Diff=lmn-log(true)) %>% 
  group_by(model, train) %>%
  summarise(mnDiff=sqrt(mean(Diff^2)), sd=sd(Diff), se=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(x=model, y=mnDiff, ymin=mnDiff-2*se, 
             ymax=mnDiff+2*se, colour=model)) + 
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~train) +
  geom_pointrange() + scale_colour_manual(values=mod_cols) + 
  labs(title=expression(log(Lambda)), x="Model", y="RMSE") 
ggsave("eda/Lambda_RMSE.pdf", width=5, height=4, units="in")

lam.out$sum.gg %>% mutate(Diff=lmn-log(true)) %>% 
  group_by(model, train) %>%
  summarise(mnDiff=sqrt(mean(Diff^2)), sd=sd(Diff), se=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(x=model, y=mnDiff, ymin=mnDiff-2*se, 
             ymax=mnDiff+2*se, colour=model)) + 
  geom_hline(yintercept=0, linetype=2) + facet_wrap(~train) +
  geom_pointrange() + scale_colour_manual(values=mod_cols) + 
  labs(title=expression(log(lambda)), x="Model", y="RMSE") 
ggsave("eda/lambda_RMSE_.pdf", width=5, height=4, units="in")

ggplot(Lam.out$sum.gg, aes(x=model, y=lmn-log(true), fill=model)) +
  geom_hline(yintercept=0, linetype=2) + 
  geom_boxplot(outlier.alpha=0.2, outlier.size=0.25, outlier.shape=1) + 
  scale_fill_manual(values=mod_cols) + facet_wrap(~train) +
  labs(title=expression(Lambda), x="Model", y="Mean error") 
ggsave("eda/Lambda_box.pdf", width=5, height=4, units="in")

ggplot(lam.out$sum.gg, aes(x=model, y=lmn-log(true), fill=model)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_boxplot(outlier.alpha=0.2, outlier.size=0.25, outlier.shape=1) + 
  scale_fill_manual(values=mod_cols) + facet_wrap(~train) +
  labs(title=expression(lambda), x="Model", y="Mean error") 
ggsave("eda/lambda_box_.pdf", width=5, height=4, units="in")

Lam.out$sum.gg %>% 
  mutate(pP.mn=1-exp(-mn), pP.true=1-exp(-true), Diff=pP.mn-pP.true) %>% 
  group_by(model, train) %>%
  summarise(mnDiff=sqrt(mean(Diff^2)), sd=sd(Diff), se=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(x=model, y=mnDiff, ymin=mnDiff-2*se, 
             ymax=mnDiff+2*se, colour=model)) + 
  facet_wrap(~train) +
  geom_pointrange() + scale_colour_manual(values=mod_cols) + 
  labs(y="RMSE in prP")
ggsave("eda/prP_1km_RMSE.pdf", width=5, height=4, units="in")

lam.out$sum.gg %>% 
  mutate(pP.mn=1-exp(-mn), pP.true=1-exp(-true), Diff=pP.mn-pP.true) %>% 
  group_by(model, train) %>%
  summarise(mnDiff=sqrt(mean(Diff^2)), sd=sd(Diff), se=sd(Diff)/sqrt(n())) %>%
  ggplot(aes(x=model, y=mnDiff, ymin=mnDiff-2*se, 
             ymax=mnDiff+2*se, colour=model)) + 
  facet_wrap(~train) +
  geom_pointrange() + scale_colour_manual(values=mod_cols) +
  labs(y="RMSE in prP")
ggsave("eda/prP_075_RMSE.pdf", width=5, height=4, units="in")



Lam.out$sum.gg %>% 
  mutate(pP.mn=1-exp(-mn), pP.true=1-exp(-true), Diff=pP.mn-pP.true) %>% 
  ggplot(aes(x=pP.mn, y=log(true), colour=model)) + 
  facet_grid(model~train) + geom_hline(yintercept=0) + 
  geom_point(alpha=0.2) + scale_colour_manual(values=mod_cols) +
  labs(x="pP", y="log(true Lambda)")
lam.out$sum.gg %>% 
  mutate(pP.mn=1-exp(-mn), pP.true=1-exp(-true), Diff=pP.mn-pP.true) %>% 
  ggplot(aes(x=pP.mn, y=log(true), colour=model)) + 
  facet_grid(model~train) + geom_hline(yintercept=0) + 
  geom_point(alpha=0.2) + scale_colour_manual(values=mod_cols) +
  labs(x="pP", y="log(true  lambda)")







# ggplot(b.out$gg, aes(x=true, colour=model)) + 
#   facet_wrap(~Parameter, scales="free", ncol=S) +
#   geom_vline(data=b.out$true, aes(xintercept=value), linetype=2) +
#   geom_density() + scale_colour_manual(values=mod_cols)
# ggplot(b.out$sum.gg, aes(mn-true, colour=model)) + 
#   geom_vline(xintercept=0, linetype=2) + geom_density() + facet_wrap(~R) +
#   scale_colour_manual(values=mod_cols) + ggtitle("b")
# ggplot(b.out$sum.gg, aes(true, mn, colour=model)) + ggtitle("b") + 
#   geom_abline() + geom_point() + stat_smooth(se=F, linetype=2, method="lm") + 
#   facet_wrap(~R, scales="free") + scale_colour_manual(values=mod_cols)

# ggplot(a.out$gg, aes(x=true, colour=model)) + 
#   facet_wrap(~Parameter, scales="free", ncol=S) +
#   geom_vline(data=a.out$true, aes(xintercept=true), linetype=2) +
#   geom_density() + scale_colour_manual(values=mod_cols)
# ggplot(a.out$sum.gg, aes(mn-true, colour=model)) + 
#   geom_vline(xintercept=0, linetype=2) + geom_density() + facet_wrap(~R) +
#   scale_colour_manual(values=mod_cols) + ggtitle("a")
# ggplot(a.out$sum.gg, aes(true, mn, colour=model)) + ggtitle("a") + 
#   geom_abline() + geom_point() + stat_smooth(se=F, linetype=2, method="lm") + 
#   facet_wrap(~R, scales="free") + scale_colour_manual(values=mod_cols)

# ggplot(Lam.out$sum.gg, aes(x=true, y=med, ymin=q025, ymax=q975, colour=model)) +
#   geom_pointrange(shape=1, alpha=0.4) + geom_abline() + 
#   stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
#   scale_x_log10() + scale_y_log10() + facet_grid(train~S, scales="free") + 
#   ggtitle(expression(Lambda~vs.~hat(Lambda))) +
#   scale_colour_manual(values=mod_cols)
# ggplot(Lam.out$sum.gg, aes(x=log(mn)-log(true), colour=model)) + 
#   geom_density() + geom_vline(xintercept=0, linetype=2) + 
#   facet_wrap(~train) + scale_colour_manual(values=mod_cols)
# ggplot(Lam.out$sum.gg, aes(x=lmn, y=log(true), colour=model)) + 
#   geom_abline() + facet_wrap(~train) +
#   geom_point(alpha=0.5) + stat_smooth(linetype=2, method="lm", se=F) +
#   scale_colour_manual(values=mod_cols) + 
#   labs(title=expression(Lambda), x="posterior mean(log)", y="log(true)") 

# ggplot(lam.out$sum.gg, aes(x=true, y=med, ymin=q025, ymax=q975, colour=model)) +
#   geom_pointrange(shape=1, alpha=0.4) + geom_abline() + 
#   stat_smooth(se=F, method="lm", linetype=2, size=0.5) +
#   scale_x_log10() + scale_y_log10() + facet_grid(train~S, scales="free") + 
#   ggtitle(expression(Lambda~vs.~hat(Lambda))) +
#   scale_colour_manual(values=mod_cols)
# ggplot(lam.out$sum.gg, aes(x=log(mn)-log(true), colour=model)) + 
#   geom_density() + geom_vline(xintercept=0, linetype=2) + 
#   facet_wrap(~train) + scale_colour_manual(values=mod_cols)
# ggplot(lam.out$sum.gg, aes(x=lmn, y=log(true), colour=model)) + 
#   geom_abline() + facet_wrap(~train) +
#   geom_point(alpha=0.25) + stat_smooth(linetype=2, method="lm", se=F) +
#   scale_colour_manual(values=mod_cols) +
#   labs(title=expression(lambda), x="log(posterior mean)", y="log(true)") 
# ggplot(lam.out$sum.gg, aes(x=log(true), y=log(mn)-log(true), colour=model)) + 
#   facet_wrap(~dataset) +
#   geom_point(alpha=0.25) + stat_smooth(linetype=2, method="lm", se=F) +
#   scale_colour_manual(values=mod_cols)

# lam.df <- data.frame(lambda=c(c(lambda$Y), c(lambda$Y_)),
#                      dataset=rep(c("Y", "Y_"), times=c(I$Y*S, I$Y_*S)),
#                      I=c(rep(1:I$Y, times=S), rep(1:I$Y_, times=S)),
#                      S=c(rep(1:S, each=I$Y), rep(1:S, each=I$Y_)),
#                      J=c(rep(J$Y, times=S), rep(J$Y_, times=S)),
#                      LAMBDA_true=c(c(LAMBDA$Y[J$Y,])*h, c(LAMBDA$Y_[J$Y_,])*h))
# lam.df$KS <- paste0(lam.df$J, "_", lam.df$S, "_", lam.df$dataset)
# Lam.out$sum.gg$KS <- paste0(Lam.out$sum.gg$K, "_", 
#                             Lam.out$sum.gg$S, "_", Lam.out$sum.gg$dataset)
# lam.df <- left_join(lam.df, # error because dataset != "Y", "Y_" in Lam.out$sum.gg
#                     Lam.out$sum.gg %>% filter(model=="Y") %>%
#                       select(mn, KS) %>% mutate(LAMBDA_mn=mn*h),
#                     by="KS")
# lam.df$Parameter <- paste0("lambda_", lam.df$dataset, 
#                            "[", lam.df$I, ",", lam.df$S, "]")
# lam.full <- left_join(lam.out$sum.gg, 
#                       select(lam.df, Parameter, J, LAMBDA_true, LAMBDA_mn), 
#                       by="Parameter")












ggplot(Lam.out$sum.gg, aes(x=1-exp(-true), y=1-exp(-mn), colour=model)) + 
  geom_abline() + geom_point(alpha=0.6) + stat_smooth(method="lm", se=F) +
  facet_wrap(~train) + scale_colour_manual(values=mod_cols)
ggplot(Lam.out$sum.gg, aes(x=1-exp(-true), 
                           y=(1-exp(-mn))-(1-exp(-true)), colour=model)) + 
  geom_point() + facet_wrap(~train) + scale_colour_manual(values=mod_cols)
  


ggplot(filter(lam.full, train=="test"), 
       aes(x=log(true), y=log(mn), colour=model)) + 
  geom_abline() + geom_point(alpha=0.5, shape=1) + 
  facet_wrap(~S, scales="free") + scale_colour_manual(values=mod_cols)

ggplot(lam.full, aes(x=1-exp(-true), y=1-exp(-mn), colour=model)) + 
  geom_abline() + geom_point(alpha=0.5, shape=1) + 
  stat_smooth(method="lm", se=F, linetype=2) +
  facet_wrap(~train) + scale_colour_manual(values=mod_cols)

ggplot(lam.full, aes((1-exp(-mn))-(1-exp(-true)), colour=model)) + 
  geom_density() + facet_wrap(~train) + scale_colour_manual(values=mod_cols)





















################################################################################
# LAMBDA as true abundance
# 1 km2 only -- aggregate soil plots Y to 1 km2

## make up parameters ##########################################################
S <- 16#139  # number of species
G <- 4#28
R <- 4  # number of regional covariates
Q <- 3  # number of W effort covariates
K <- list(W=100, Y=44, W_=20, Y_=10)  # number of cells
J <- list(Y=25*K$Y, Y_=25*K$Y_)  # number of soil plots
Lambda_0 <- 10.5
effort_0 <- -11
h <- 0.00001885  # (25*6*pi*(0.2)^2)/(1000*1000)
sigma_b <- 0.5
tax_i <- cbind(1:S, sample(1:G, S, TRUE, runif(G)), rbinom(S, 1, 0.3)+1)
quadratic_R <- TRUE


## simulate data ###############################################################
# covariates, slopes, and true Lambda
X <- make_X(K, R, quadratic_R)
beta.ls <- make_slopes(Lambda_0, R, G, S, tax_i, sigma_b, quadratic_R)
LAMBDA <- map(X, ~exp(. %*% beta.ls$sp))

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
beta.out <- aggregate_aggSlopes(out.ls, beta.ls, "beta")
ggplot(beta.out$gg, aes(x=value, colour=model)) + 
  facet_wrap(~Parameter, scales="free") + 
  geom_vline(data=beta.out$true, aes(xintercept=value), linetype=2) +
  geom_density() + scale_colour_manual(values=mod_cols)

# b
b.out <- aggregate_spSlopes(out.ls, beta.ls, tax_i, "b")
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







