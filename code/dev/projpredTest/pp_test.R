
#### vignette ------------------------------------------------------------------
devtools::install_github('stan-dev/projpred', "develop", build_vignettes=F, 
                         force=TRUE)
# https://github.com/stan-dev/projpred/blob/develop/vignettes/quickstart_glmm.Rmd

library(rstanarm); library(projpred); library(bayesplot); library(tidyverse)
theme_set(theme_bw())

data_pois <- read.table(
  "https://paul-buerkner.github.io/data/data_pois.txt", 
  header = TRUE
)
data_pois$obs <- 1:nrow(data_pois)

fit <- stan_glmer(
  phen_pois ~ cofactor + (1 | phylo) + (1 | obs), data = data_pois,
  family = poisson("log"), chains = 2, iter = 2000,
  control = list(adapt_delta = 0.95)
)
fit2 <- stan_glmer(
  phen_pois ~ cofactor + (1 | phylo), data = data_pois[1:150,],
  family = poisson("log"), chains = 2, iter = 2000,
  control = list(adapt_delta = 0.95)
)

vs <- varsel(fit2, intercept=T)
vs <- varsel(fit, verbose=T)

solution_terms(vs) 
plot(vs, stats = c('elpd', 'rmse'))

proj <- project(vs, nterms = 2, ndraws = 10)
mcmc_areas(as.matrix(proj), pars = solution_terms(vs)[1:2])












data("Orthodont", package = "nlme")
fit <- stan_glmer(distance ~ age * Sex + (age | Subject),
                  chains = 2,  data = Orthodont, seed = 1)
vs <- varsel(fit)
solution_terms(vs) 
plot(vs, stats = c('elpd', 'rmse'))
















#### trying custom model -------------------------------------------------------

# I can fit a model with rstanarm using Y and all of the covariates, then update
# the posteriors and predicted values with the full models. That will give the
# structure that works with projpred, but optimizes for the actual hierarchical
# model (i.e., which include phylogeny and both data sets).
# So: 
# 1) Fit full models for Y, WY, W
# 2) Fit rstanarm model for Y
# 3) Update posteriors and fitted values from 2) with those from 3)
# 4) Select variables to include for each model
# 5) Run optimal models for interpretation

# PROBLEM: how do the fit$stanfit@sim$samples[[1]] b's align with the posterior?
# There are too many of them.... ?? 
# intercept + x1:6 = 7 b's per species
# IT CREATES A NEW SPECIES!! b[x1 sp:_NEW_sp]

# PROBLEM: Is the parameterization an effects or mean? 
# My output will be means: each b is the actual slope
# The glmer output is effects -- I need to subtract 'beta' from each 'b'

library(rstanarm); library(projpred); library(bayesplot); library(tidyverse)
library(rstan)
theme_set(theme_bw())
options(mc.cores=parallel::detectCores())
rstan::rstan_options(auto_write=TRUE)

set.seed(10101)

S <- 10
J <- 20
I <- J*5
IJ <- rep(1:J, each=5)

R <- 3
L <- 2
beta.true <- c(-4, rnorm(R, 0, 1), rnorm(L, 0, 1))
beta.true[sample(2:(R+L), floor((R+L)/3), F)] <- 0
b <- do.call('cbind', map(1:S, ~rnorm(length(beta.true), beta.true, 0.5)))

X <- cbind(1, do.call('cbind', map(1:R, ~rnorm(J)))) # row = site
V <- do.call('cbind', map(1:L, ~rnorm(I))) # row = plot
Z <- cbind(X[IJ,], V) # row = plot

Zt <- do.call('rbind', map(1:S, ~Z))
llambda <- do.call('rbind', map(1:S, ~Zt[(.*I-I+1):(.*I),] %*% b[,.]))
Zt <- cbind(Zt, rep(1:S, each=I))
colnames(Zt) <- c(paste0("x", 1:(ncol(Zt)-1)), "sp")

Y <- rpois(length(llambda), exp(llambda))
# Y <- map2_dbl(runif(length(llambda)), exp(llambda),
#               ~sum(.x > cumsum(LaplacesDemon::dgpois(0:round(.y*2), .y, 0.17))))
fit_glmer <- stan_glmer(
  as.formula(paste("Y ~ ",
                   paste(colnames(Zt)[-c(1, ncol(Zt))],
                         collapse="+"), "+",
                   paste("(1+",
                         paste(colnames(Zt)[-c(1,ncol(Zt))],
                               collapse="+"),
                         "|sp)"))), 
  data=as.data.frame(cbind(Y, Zt)), 
  family=poisson("log"), 
  prior=hs(),
  chains=2, warmup=180, iter=200, refresh=10)
saveRDS(fit_glmer, "out/fit_shell_TEST.rds")

fit_hm <- rstan::stan(
  file="code/mods/full_Y_altRE.stan",
  init=0, 
  # sample_file="out/vs_Y_TEST",
  data=list(K=0, J=J, I=I, IJ=IJ,
            S=S, G=S, R=R+1, L=L,
            tax_i=cbind(1:S, 1:S),
            Y=matrix(Y, nrow=I, ncol=S), 
            X=X, V=V, h=7.5e-7), 
  chains=2, warmup=180, iter=200, refresh=10)
fit_hm_G <- rstan::stan(
  file="code/mods/full_Y_altRE_G.stan",
  init=0, 
  sample_file="out/vs_Y_TEST",
  data=list(K=0, J=J, I=I, IJ=IJ,
            S=S, G=4, R=R+1, L=L,
            tax_i=cbind(1:S, sample(4, S, T)),
            Y=matrix(Y, nrow=I, ncol=S), 
            X=X, V=V, h=7.5e-7), 
  chains=2, warmup=180, iter=200, refresh=10)


# build lookup table to align parameters between glmer and hm
names_glmer <- names(fit_glmer$stanfit@sim$samples[[1]])
names_hm <- names(fit_hm@sim$samples[[1]])
name_lookup <- tibble(glmer=c("alpha[1]", 
                              paste0("beta[", 1:(R+L), "]"),
                              paste0("b[", 1:((R+L+1)*S), "]"),
                              "lp__"),
                      glmer_l=c(
                        rownames(fit_glmer$stan_summary)[1:((R+L+1)*(S+1))],
                        "log-posterior"),
                      hm=c(paste0("beta[", 1:(R+L+1), "]"),
                           paste0(rep(paste0("b[", 1:(R+L+1)), times=S), ",", 
                                  rep(1:S, each=(R+L+1)), "]"),
                           "lp__"), 
                      cov=c(rep(1:(R+L+1), times=S+1), NA))


# extract and align predicted values
llam.summary <- rstan::summary(fit_hm, pars="llambda")$summary
llam_hm.df <- tibble(llam=llam.summary[,1],
                     par=rownames(llam.summary),
                     site=str_split_fixed(str_remove(par, "llambda\\["), 
                                          ",", 2)[,1],
                     spp=str_split_fixed(str_remove(par, "]"), ",", 2)[,2]) %>%
  mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
  arrange(spp, site)

# extract slope means
beta.summary <- rstan::summary(fit_hm, pars="beta")$summary
beta_hm.df <- tibble(beta=beta.summary[,1],
                     beta_md=beta.summary[,6],
                     beta_se=beta.summary[,2],
                     par=rownames(beta.summary)) %>%
  mutate(cov=row_number())
b.summary <- rstan::summary(fit_hm, pars="b")$summary
b_hm.df <- tibble(b=b.summary[,1],
                  b_md=b.summary[,6],
                  b_se=b.summary[,2],
                  par=rownames(b.summary),
                  cov=str_split_fixed(str_remove(par, "b\\["), ",", 2)[,1],
                  spp=str_split_fixed(str_remove(par, "]"), ",", 2)[,2]) %>%
  mutate(cov=as.numeric(cov), spp=as.numeric(spp)) %>%
  arrange(spp, cov) %>%
  left_join(., select(beta_hm.df, beta, beta_md, cov), by='cov') %>%
  mutate(b_eff=b-beta, 
         b_eff_md=b_md-beta_md)


# reference model: take structure of glmer and replace values with fitted hm
fit_ref <- fit_glmer
fit_ref$linear.predictors <- llam_hm.df$llam
fit_ref$fitted.values <- exp(fit_ref$linear.predictors)
fit_ref$coefficients[] <- c(beta_hm.df$beta_md, b_hm.df$b_eff_md)
fit_ref$ses[] <- c(beta_hm.df$beta_se, b_hm.df$b_se)
fit_ref$residuals <- fit_ref$y - fit_ref$fitted.values

iter_hm <- length(fit_hm@sim$samples[[1]][[1]])
iter_glmer <- length(fit_glmer$stanfit@sim$samples[[1]][[1]])
iter_index <- (iter_hm-iter_glmer+1):iter_hm 

sum_hm <- summary(fit_hm, probs=c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))$summary


for(i in 1:nrow(name_lookup)) {
  sum_i <- sum_hm[rownames(sum_hm)==name_lookup$hm[i],]
  # translate to effects parameterization for species effects
  if(grepl('b\\[', name_lookup$hm[i])) {
    sum_i_beta <- sum_hm[rownames(sum_hm)==paste0("beta[", 
                                                  name_lookup$cov[i], "]"),]
    ind <- (1:length(sum_i))[-c(2,length(sum_i)-1, length(sum_i))]
    sum_i[ind] <- sum_i[ind] - sum_i_beta[ind]
  }
  fit_ref$stan_summary[rownames(fit_ref$stan_summary)==name_lookup$glmer_l[i],] <- sum_i
  for(j in 1:length(fit_glmer$stanfit@sim$samples)) {
    # translate to effects parameterization for species effects
    post_i <- fit_hm@sim$samples[[j]][[name_lookup$hm[i]]]
    if(grepl('b\\[', name_lookup$hm[i])) {
      post_i_beta <- fit_hm@sim$samples[[j]][[paste0("beta[", 
                                                     name_lookup$cov[i], "]")]]
      post_i <- post_i - post_i_beta
    }
    fit_ref$stanfit@sim$samples[[j]][[name_lookup$glmer[i]]] <- post_i[iter_index]
  }
}




plot(beta.true[-1], summary(fit_glmer, pars="beta")[,1], 
     xlim=c(-2,2), ylim=c(-2,2), ylab="Estimate", xlab="True")
points(beta.true[-1], rstan::summary(fit_hm, pars="beta")$summary[-1,1], col="red")
abline(a=0,b=1)
legend("topleft", col=c(1, 2), pch=1, c("glmer", "hm"), bty="n")

plot(c(b), 
     as.matrix(summary(fit_glmer, pars="b"))[,1] + 
       summary(fit_glmer, pars=c("alpha","beta"))[,1], 
     xlim=c(-2,2), ylim=c(-2,2), ylab="Estimate", xlab="True")
points(c(t(b)), rstan::summary(fit_hm, pars="b")$summary[,1], col="red")
abline(a=0,b=1)
legend("topleft", col=c(1, 2), pch=1, c("glmer", "hm"), bty="n")


# select variables
vs_glmer <- varsel(fit_glmer, nterms_max=8, verbose=T, validate_search=F)
plot(vs_glmer, stats=c("elpd", "rmse"))
solution_terms(vs_glmer)

vs_ref <- varsel(fit_ref, nterms_max=8, verbose=T, validate_search=F)
plot(vs_ref, stats=c("elpd", "rmse"))
solution_terms(vs_ref)





plot(beta.true[-1], summary(fit_glmer, pars="beta")[,1], 
     xlim=c(-2,2), ylim=c(-2,2), ylab="Estimate", xlab="True")
points(beta.true[-1], rstan::summary(fit_ref, pars="beta")$summary[-1,1], col="red")
abline(a=0,b=1)
legend("topleft", col=c(1, 2), pch=1, c("glmer", "hm"), bty="n")

plot(c(b), 
     as.matrix(summary(fit_glmer, pars="b"))[,1] + 
       summary(fit_glmer, pars=c("alpha","beta"))[,1], 
     xlim=c(-2,2), ylim=c(-2,2), ylab="Estimate", xlab="True")
points(c(t(b)), rstan::summary(fit_hm, pars="b")$summary[,1], col="red")
abline(a=0,b=1)
legend("topleft", col=c(1, 2), pch=1, c("glmer", "hm"), bty="n")














#### stan_glmm with Y only -----------------------------------------------------

library(rstanarm); library(projpred); library(bayesplot); library(tidyverse)
theme_set(theme_bw())
options(mc.cores=parallel::detectCores())

vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
out.dir <- "out/vs/"

d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))

S <- d.ls$S
I <- d.ls$I
Y_vec <- c(d.ls$Y)
Z <- cbind(d.ls$X[(d.ls$K+1):(d.ls$K+d.ls$J),][d.ls$IJ,], d.ls$V)
Zt <- do.call('rbind', map(1:S, ~Z))

YZ <- as.data.frame(cbind(Y_vec, Zt[,-1])) %>%
  mutate(sp=rep(1:S, each=I)) %>% rename(obs=Y_vec) #%>%
  # filter(sp %in% which(colSums(d.i$Y)>30)) %>% 
  # filter(obs>0)

YZ <- rbind(filter(YZ, obs>0),
            filter(YZ, obs==0) %>% sample_n(1000))

var_inc <- colnames(Z)[c(2,14,16,21)]
# form.glmm <- paste("obs ~", paste(var_inc, collapse=" + "), "+",
#                    "(1 | sp) +",
#                    paste("(0 +", paste(var_inc, collapse=" + "), "| sp)"))
form.glmm <- paste("obs ~", paste(var_inc, collapse=" + "), "+",
                   paste("(1 +", paste(var_inc, collapse=" + "), "| sp)"))

fit <- stan_glmer(
  as.formula(form.glmm), data=YZ,
  # obs ~ MAT + SoilTSt + MAT_sq + Crop + (1+MAT+SoilTSt+MAT_sq | sp), data=YZ,
  family=poisson("log"), chains=2, iter=2000
)

plot(fit, pars="beta")

vs <- varsel(fit, validate_search=F)

solution_terms(vs) 
plot(vs, stats = c('elpd', 'rmse'))

proj <- project(vs, nterms=length(solution_terms(vs)), ndraws=10)
mcmc_areas(as.matrix(proj), pars=solution_terms(vs))

p <- mcmc_areas(as.matrix(fit), 
                pars=paste0("b[MAT sp:", unique(YZ$sp), "]"))
ggsave("~/Desktop/test_MAT.pdf", p, width=5, height=20)
p <- mcmc_areas(as.matrix(fit), 
                pars=paste0("b[MAT_sq sp:", unique(YZ$sp), "]"))
ggsave("~/Desktop/test_MAT_sq.pdf", p, width=5, height=20)
p <- mcmc_areas(as.matrix(fit), 
                pars=paste0("b[Litter sp:", unique(YZ$sp), "]"))
ggsave("~/Desktop/test_Litter.pdf", p, width=5, height=20)
p <- mcmc_areas(as.matrix(fit), 
                  pars=paste0("b[SoilTSt sp:", unique(YZ$sp), "]"))
ggsave("~/Desktop/test_SoilTSt.pdf", p, width=5, height=20)
p <- mcmc_areas(as.matrix(fit), 
                  pars=paste0("b[(Intercept) sp:", unique(YZ$sp), "]"))
ggsave("~/Desktop/test_int.pdf", p, width=5, height=20)






#### simulate data with Y structure --------------------------------------------
library(projpred)
S <- 10
J <- 30
I <- J*15
IJ <- rep(1:J, each=15)

R <- 10
L <- 10
beta.true <- c(-4, rnorm(R, 0, 0.5), rnorm(L, 0, 0.5))
beta.true[sample(2:(R+L), floor((R+L)/3), F)] <- 0
b <- do.call('cbind', map(1:S, ~rnorm(length(beta.true), beta.true, 0.1)))

X <- cbind(1, do.call('cbind', map(1:R, ~rnorm(J)))) # row = site
V <- do.call('cbind', map(1:L, ~rnorm(I))) # row = plot
Z <- cbind(X[IJ,], V) # row = plot

Zt <- do.call('rbind', map(1:S, ~Z))

llambda <- do.call('rbind', map(1:S, ~Zt[(.*I-I+1):(.*I),] %*% b[,.]))
Y <- rpois(length(llambda), exp(llambda))
Y <- map2_dbl(runif(length(llambda)), exp(llambda),
             ~sum(.x > cumsum(LaplacesDemon::dgpois(0:round(.y*2), .y, 0.17))))


predfun <- function(zt) {
  out.ls <- vector("list", S)
  for(s in 1:S) {
    sp_col <- grep(paste0(",", s, "\\]"), colnames(b))
    out.ls[[s]] <- b[,sp_col[-1]] %*% t(zt[(s*I-I+1):(s*I),]) + 
      b[,sp_col[1]] #+ log(d.ls$h)
  }
  t(exp(do.call('cbind', out.ls)))
}
fit <- rstan::stan("code/projpredTest/Y_pp.stan", warmup=2000, iter=2500,
                   chains=2,
                   data=list(J=J, I=I, IJ=IJ, S=S, P=length(beta.true), Y=Y, Z=Zt))

draws <- as.matrix(fit)
beta <- draws[,grep('beta', colnames(draws))]
b <- draws[,grep('b', colnames(draws))]

ref <- init_refmodel(Zt[,-1], Y, poisson(), Zt[,-1], predfun,
                     wobs=if_else(Y==0, 1, 1))
vs <- varsel(ref, relax=FALSE, nv_max=R+L)
n_var <- suggest_size(vs, baseline="best", alpha=0.05)
varsel_plot(vs, stats=c('elpd', 'rmse'), deltas=T, baseline="ref")
vs.df <- varsel_stats(vs, stats=c('elpd'), baseline="best", deltas=T) %>%
  mutate(par=rownames(.),
         true=beta.true[-1][vind])

ggplot(vs.df, aes(x=size, y=elpd, label=par, colour=abs(true),
                  ymin=elpd-2*elpd.se, ymax=elpd+2*elpd.se)) +
  geom_point() + geom_linerange() + 
  scale_colour_gradient(low="gray70", high="red3") +
  geom_vline(xintercept=n_var+0.5, linetype=3) + 
  geom_text(size=2.5, nudge_x=0.15, hjust=0)


proj <- project(vs, nv=n_var)
full.beta <- beta[,vs.df$vind[1+(1:n_var)]+1]
colnames(full.beta) <- vs.df$par[1+(1:n_var)]
mcmc_areas(as.matrix(proj)[,-1]) + geom_vline(xintercept=0, linetype=2)
mcmc_areas(full.beta) + geom_vline(xintercept=0, linetype=2)

fit.sum <- summary(fit, pars="llambda")$summary
plot(exp(fit.sum[,1]), exp(llambda)); abline(a=0, b=1, lty=2)

pred <- proj_linpred(vs, xnew=Zt[,-1], ynew=Y, nv=n_var, integrated=T, transform=T)
ggplot() +
  # geom_point(aes(x=pred$pred, y=Y), alpha=0.5) +
  # geom_abline(linetype=3) + stat_smooth(aes(x=pred$pred, y=Y), method="lm") +
  geom_point(aes(x=1-exp(-pred$pred), y=as.numeric(Y>0)), alpha=0.5) +
  stat_smooth(aes(x=1-exp(-pred$pred), y=as.numeric(Y>0)),
              method="glm", size=0.5, method.args=list(family="binomial")) +
  labs(x = 'prediction', y = 'y') 


y1_rep <- proj_predict(vs, xnew=Zt[which.max(Y),-1,drop=F], nv=n_var)
qplot(as.vector(y1_rep), bins=25) +
  geom_vline(xintercept=Y[which.max(Y)], color='red') +
  xlab('y1_rep')







#### try to use real output ----------------------------------------------------
library(rstan); library(tidyverse); library(bayesplot); library(projpred)
theme_set(theme_classic())
# options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
pars.exc <- c("a_std", "A_std", "b_std", "B_std", "f_std", "F_std",
              # "LAMBDA", "LAMBDA_", 
              "prPres", "prPres_",# "prPresL_",
              "lambda", "lambda_",
              "lLAMBDA", "lLAMBDA_", "llambda", #"llambda_",
              "L_Omega_A", "L_Omega_B", "sigma_A", "sigma_B", "L_Omega_F", "sigma_F",
              "p", "p_")

vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
out.dir <- "out/vs/"

d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))
out <- read_stan_csv(dir(out.dir, "^WY_GP_[0-9]", full.names=T))
# out <- stan("code/vs/vs_Y_Pois.stan",
#             data=read_rdump(paste0("data/vs/", vs_d, ".Rdump")), init=0,
#             chains=2, warmup=100, iter=150, pars=pars.exc, include=F)



S <- d.ls$S
I <- d.ls$I
I_ <- d.ls$I_
draws <- as.matrix(out)
b <- draws[,grep('^b\\[', colnames(draws))]
B <- draws[,grep('^B\\[', colnames(draws))]
beta <- draws[,grep('^beta\\[', colnames(draws))]

Y_vec <- c(d.ls$Y)
Y_vec_ <- c(d.ls$Y_)
Z <- cbind(d.ls$X[(d.ls$K+1):(d.ls$K+d.ls$J),][d.ls$IJ,], d.ls$V)
Zt <- do.call('rbind', map(1:S, ~Z))
Z_ <- cbind(d.ls$X_[d.ls$IJ_,], d.ls$V_)
Zt_ <- do.call('rbind', map(1:S, ~Z_))


predfun <- function(zt) {
  out.ls <- vector("list", S)
  for(s in 1:S) {
    sp_col <- grep(paste0(",", s, "\\]"), colnames(b))
    out.ls[[s]] <- b[,sp_col[-1]] %*% t(zt[(s*I-I+1):(s*I),]) + b[,sp_col[1]]
  }
  t(exp(do.call('cbind', out.ls)))
}
ref <- init_refmodel(Zt[,-1], Y_vec, poisson(), NULL, predfun, 
                     wobs=if_else(Y_vec==0, 1, 1))
vs <- varsel(ref, nv_max=d.ls$R+d.ls$L-1)
n_var <- suggest_size(vs, baseline="best", alpha=0.05)
varsel_plot(vs, stats=c('elpd'), deltas=T, baseline="best")
vs.df <- varsel_stats(vs, stats=c('elpd'), baseline="best", deltas=T) %>%
  mutate(par=rownames(.))

ggplot(vs.df, aes(x=size, y=elpd, label=par,
                   ymin=elpd-2*elpd.se, ymax=elpd+2*elpd.se)) +
  geom_point() + geom_linerange() + 
  geom_vline(xintercept=n_var+0.5, linetype=3) + 
  geom_text(size=2.5, nudge_x=0.15, hjust=0)


proj <- project(vs, nv=n_var)
full.beta <- beta[,vs.df$vind[1+(1:n_var)]+1]
colnames(full.beta) <- vs.df$par[1+(1:n_var)]
mcmc_areas(as.matrix(proj)[,-1]) #+ coord_cartesian(c(-2,2))
mcmc_areas(full.beta)

pred <- proj_linpred(vs, xnew=Zt[,-1], ynew=Y_vec, nv=n_var, 
                     integrated=T, transform=F)
ggplot() + xlim(0,1) +
  # geom_point(aes(x=pred$pred, y=Y_vec), alpha=0.5) +
  # geom_abline(linetype=3) + stat_smooth(aes(x=pred$pred, y=Y_vec), method="lm") +
  geom_point(aes(x=1-exp(-pred$pred), y=as.numeric(Y_vec>0)), alpha=0.5) +
  stat_smooth(aes(x=1-exp(-pred$pred), y=as.numeric(Y_vec>0)), fullrange=T,
              method="glm", size=0.5, method.args=list(family="binomial")) +
  labs(x = 'prediction', y = 'y') 


y1_rep <- proj_predict(vs, xnew=Zt[44,-1,drop=F], nv=n_var, seed=7560)
qplot(as.vector(y1_rep), bins=25) +
  geom_vline(xintercept=Y_vec[44], color='red') +
  xlab('y1_rep')

















#### ensemble across species ---------------------------------------------------
library(rstan); library(tidyverse); library(bayesplot); library(projpred)
theme_set(theme_classic())
# options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
pars.exc <- c("a_std", "A_std", "b_std", "B_std", "f_std", "F_std",
              # "LAMBDA", "LAMBDA_", 
              "prPres", "prPres_",# "prPresL_",
              "lambda", "lambda_",
              "lLAMBDA", "lLAMBDA_", "llambda", #"llambda_",
              "L_Omega_A", "L_Omega_B", "sigma_A", "sigma_B", "L_Omega_F", "sigma_F",
              "p", "p_")

vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
out.dir <- "out/vs/"

d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))
out <- read_stan_csv(dir(out.dir, "^WY_Pois_[0-9]", full.names=T))
# out <- stan("code/vs/vs_Y_Pois.stan",
#             data=read_rdump(paste0("data/vs/", vs_d, ".Rdump")), init=0,
#             chains=2, warmup=100, iter=150, pars=pars.exc, include=F)



S <- d.ls$S
I <- d.ls$I
I_ <- d.ls$I_
draws <- as.matrix(out)
b <- draws[,grep('^b\\[', colnames(draws))]
B <- draws[,grep('^B\\[', colnames(draws))]
beta <- draws[,grep('^beta\\[', colnames(draws))]

Z <- cbind(d.ls$X[(d.ls$K+1):(d.ls$K+d.ls$J),][d.ls$IJ,], d.ls$V)
Z_ <- cbind(d.ls$X_[d.ls$IJ_,], d.ls$V_)


predfun <- function(zt) {
  t(exp(b_sp[,-1] %*% t(zt) + b_sp[,1]))
}

ref.ls <- vs.ls <- vs_stats.ls <- proj.ls <- pred.ls <- vector("list", S)
for(s in 1:S) {
  b_sp <- b[,grep(paste0(",", s, "\\]"), colnames(b))]
  ref.ls[[s]] <- init_refmodel(Z[,-1], d.ls$Y[,s], poisson(), NULL, predfun, 
                       wobs=if_else(d.ls$Y[,s]==0, 1, 1))  
  vs.ls[[s]] <- varsel(ref.ls[[s]], nv_max=(d.ls$R+d.ls$L-1))
  vs_stats.ls[[s]] <- varsel_stats(vs.ls[[s]], stats=c('elpd'), deltas=T) %>%
    mutate(par=rownames(.),
           spNum=s, 
           wt_sp=sum(d.ls$Y[,s]>0))
  proj.ls[[s]] <- project(vs.ls[[s]], nv=vs.ls[[s]]$nv_max)
  pred.ls[[s]] <- proj_linpred(vs.ls[[s]], xnew=Z[,-1], ynew=d.ls$Y[,s], 
                               nv=vs.ls[[s]]$nv_max, integrated=T, transform=T)
}

vs.size <- map_dbl(vs.ls, suggest_size) %>% replace_na(0)

vs.best <- map2_dfr(vs_stats.ls, vs.size, ~filter(.x, .x$size <= .y)) %>%
  mutate(spName=d.i$tax_i$species[spNum])
ggplot(filter(vs.best, wt_sp>10), aes(x=size, y=elpd)) + 
  geom_point() + geom_line() + 
  geom_linerange(aes(ymin=elpd-elpd.se, ymax=elpd+elpd.se)) +
  facet_wrap(~spName, scales="free_y")
ggplot(filter(vs.best, par!=""), aes(x=factor(size))) + 
  geom_bar() + facet_wrap(~par)
vs.best %>% filter(par!="") %>% mutate(wt_sp=wt_sp/sum(d.ls$Y>0)) %>%
  group_by(par) %>% 
  summarise(mnOrd=mean(size*wt_sp)) %>% ungroup %>%
  arrange(mnOrd) %>% mutate(par=factor(par, levels=unique(par))) %>%
  ggplot(aes(x=par, y=mnOrd)) + geom_point() 

vs.df <- do.call('rbind', vs_stats.ls) %>%
  mutate(wt_sp=wt_sp/sum(d.ls$Y))

ggplot(vs.df, aes(x=size, y=elpd, group=spNum)) + geom_line(alpha=0.2)
vs.df %>% group_by(size) %>% 
  summarise(elpd=sum(elpd), elpd.se=sum(elpd.se),
            elpd.wt=sum(elpd*wt_sp)) %>%
  ggplot(aes(x=size, y=elpd.wt)) + geom_line() + geom_point() +
  geom_linerange(aes(ymin=elpd-elpd.se, ymax=elpd+elpd.se))
vs.df %>% group_by(par) %>% 
  summarise(mnOrd=mean(size*wt_sp)) %>% ungroup %>%
  arrange(mnOrd) %>% mutate(par=factor(par, levels=unique(par))) %>%
  ggplot(aes(x=par, y=mnOrd)) + geom_point() 
ggplot(filter(vs.df, par!=""), aes(x=size*wt_sp)) + 
  geom_density() + facet_wrap(~par, scales="free")
ggplot(filter(vs.df, par!=""), aes(x=par, y=size*wt_sp)) + 
  geom_boxplot()
vs.df %>% group_by(spNum) %>% mutate(elpd.diff=elpd - lag(elpd)) %>%
  filter(par!="") %>% group_by(par) %>%
  summarise(elpd=mean(elpd.diff)) %>% 
  arrange(-elpd) %>% mutate(par=factor(par, levels=unique(par))) %>%
  ggplot(aes(x=par, y=elpd)) + geom_point()


pred.df <- data.frame(pred=unlist(map(pred.ls, ~.$pred)),
                      obs=c(d.ls$Y),
                      spNum=rep(1:d.ls$S, each=d.ls$I)) %>%
  mutate(spName=d.i$tax_i$species[spNum])
ggplot(pred.df, aes(x=pred, y=obs)) + geom_point(alpha=0.25) + 
  facet_wrap(~spName) + stat_smooth(method="lm", size=0.5)
ggplot(pred.df, aes(x=1-exp(-pred), y=as.numeric(obs>0))) + xlim(0,1) +
  stat_smooth(method="glm", size=0.5, fullrange=T,
              method.args=list(family="binomial")) + 
  geom_point(alpha=0.1) + facet_wrap(~spName)




full.beta <- beta
colnames(full.beta) <- c("int", colnames(Z)[-1])
mcmc_areas(as.matrix(proj.ls[[s]])[,-1]) #+ coord_cartesian(c(-2,2))
mcmc_areas(full.beta)

proj.all <- map(proj.ls, ~matrix(0, nrow=ncol(.$beta), ncol=ncol(Z)-1,
                                 dimnames=list(NULL, colnames(Z)[-1])))
for(i in seq_along(proj.ls)) {
  proj_i.mx <- as.matrix(proj.ls[[i]])[,-1]
  proj.all[[i]][,match(colnames(proj_i.mx), colnames(proj.all[[i]]))] <- proj_i.mx
}

proj.all <- do.call('rbind', proj.all)
mcmc_areas(proj.all)
proj.df <- as.data.frame(proj.all) %>%
  mutate(spNum=rep(1:d.ls$S, each=ncol(proj.ls[[1]]$beta)),
         spName=d.i$tax_i$species[spNum]) %>%
  pivot_longer(1:ncol(proj.all), names_to="par", values_to="value")
proj.mn <- proj.df %>% group_by(par) %>% 
  summarise(mn=mean(value), med=median(value), sd=sd(value))

proj.spMn <- proj.df %>% group_by(spName, par) %>%
  summarise(mn=mean(value)) %>% group_by(par) %>%
  summarise(mnMag=mean(abs(mn)), sd=sd(mn))

beta.df <- as.data.frame(full.beta[,-1]) %>%
  pivot_longer(1:ncol(proj.all), names_to="par", values_to="value") %>%
  group_by(par) %>% summarise(mn=mean(value), med=median(value))

ggplot(proj.df, aes(x=value, group=spName)) + 
  geom_vline(xintercept=0, size=0.5) + 
  geom_density(alpha=0.1, size=0.25, colour=NA, fill="gray30") + 
  geom_vline(data=proj.mn, aes(xintercept=mn), 
             linetype=2, size=0.5, colour="blue") + 
  xlim(-5,5) +
  facet_wrap(~par, scales="free_y") + 
  theme_bw() + theme(panel.grid=element_blank())



b.df <- as.data.frame(b) %>% 
  pivot_longer(1:ncol(.), names_to="p", values_to="value") %>%
  mutate(spNum=as.numeric(str_remove(str_split_fixed(p, ",", 2)[,2], "]")),
         pNum=as.numeric(str_remove(str_split_fixed(p, ",", 2)[,1], "b\\[")),
         spName=d.i$tax_i$species[spNum],
         par=colnames(Z)[pNum])
b.mn <- b.df %>% filter(par!="") %>% group_by(spName, par) %>%
  summarise(mn=mean(value)) %>% group_by(par) %>%
  summarise(mnMag=mean(abs(mn)), sd=sd(mn))








#### ensemble across species ---------------------------------------------------
library(rstan); library(tidyverse); library(bayesplot); library(projpred)
theme_set(theme_classic())
# options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
pars.exc <- c("a_std", "A_std", "b_std", "B_std", "f_std", "F_std",
              # "LAMBDA", "LAMBDA_", 
              "prPres", "prPres_",# "prPresL_",
              "lambda", "lambda_",
              "lLAMBDA", "lLAMBDA_", "llambda", #"llambda_",
              "L_Omega_A", "L_Omega_B", "sigma_A", "sigma_B", "L_Omega_F", "sigma_F",
              "p", "p_")

vs_d <- str_sub(dir("data/vs/", ".Rdump"), 1, -7)[2]
out.dir <- "out/vs/"

d.ls <- readRDS(paste0("data/vs/", vs_d, "_ls.rds"))
d.i <- readRDS(paste0("data/vs/", vs_d, "_i.rds"))
out <- read_stan_csv(dir(out.dir, "^WY_GP_[0-9]", full.names=T))
# out <- stan("code/vs/vs_Y_Pois.stan",
#             data=read_rdump(paste0("data/vs/", vs_d, ".Rdump")), init=0,
#             chains=2, warmup=100, iter=150, pars=pars.exc, include=F)



S <- d.ls$S
I <- d.ls$I
I_ <- d.ls$I_
draws <- as.matrix(out)
b <- draws[,grep('^b\\[', colnames(draws))]
B <- draws[,grep('^B\\[', colnames(draws))]
beta <- draws[,grep('^beta\\[', colnames(draws))]

Y_vec <- c(d.ls$Y)
Z <- cbind(d.ls$X[(d.ls$K+1):(d.ls$K+d.ls$J),][d.ls$IJ,], d.ls$V)
Zt <- do.call('cbind', map(1:S, ~Z))


predfun <- function(zt) {
  t(exp(b %*% t(zt)))
}



ref <- init_refmodel(Zt, Y_vec, poisson(), NULL, predfun, 
                     wobs=if_else(Y_vec==0, 1, 1))
vs <- varsel(ref, nv_max=d.ls$R+d.ls$L+1)
n_var <- suggest_size(vs, baseline="best", alpha=0.05)
varsel_plot(vs, stats=c('elpd'), deltas=T, baseline="best")
vs.df <- varsel_stats(vs, stats=c('elpd'), baseline="best", deltas=T) %>%
  mutate(par=rownames(.))

ggplot(vs.df, aes(x=size, y=elpd, label=par,
                  ymin=elpd-2*elpd.se, ymax=elpd+2*elpd.se)) +
  geom_point() + geom_linerange() + 
  geom_vline(xintercept=n_var+0.5, linetype=3) + 
  geom_text(size=2.5, nudge_x=0.15, hjust=0)




ref.ls <- vs.ls <- vector("list", S)
for(s in 1:S) {
  b_sp <- b[,grep(paste0(",", s, "\\]"), colnames(b))]
  ref.ls[[s]] <- init_refmodel(Z[,-1], d.ls$Y[,s], poisson(), NULL, predfun, 
                               wobs=if_else(d.ls$Y[,s]==0, 1, 1))  
  vs <- varsel(ref.ls[[s]], nv_max=d.ls$R+d.ls$L+1)
  vs.ls[[s]] <- varsel_stats(vs, stats=c('elpd'), baseline="best", deltas=T) %>%
    mutate(par=rownames(.),
           spNum=s)
}

vs.df <- do.call('rbind', vs.ls)

ggplot(vs.df, aes(x=size, y=elpd, group=spNum)) + geom_line(alpha=0.2)
vs.df %>% group_by(size) %>% summarise(elpd=sum(elpd), elpd.se=sum(elpd.se)) %>%
  ggplot(aes(x=size, y=elpd)) + geom_line() + geom_point() +
  geom_linerange(aes(ymin=elpd-elpd.se, ymax=elpd+elpd.se))
vs.df %>% group_by(par) %>% summarise(mnOrd=mean(size)) %>% ungroup %>%
  arrange(mnOrd) %>% mutate(par=factor(par, levels=unique(par))) %>%
  ggplot(aes(x=par, y=mnOrd)) + geom_point() 



vs <- varsel(ref, nv_max=d.ls$R+d.ls$L+1)
n_var <- suggest_size(vs, baseline="best", alpha=0.05)
varsel_plot(vs, stats=c('elpd'), deltas=T, baseline="best")
vs.df <- varsel_stats(vs, stats=c('elpd'), baseline="best", deltas=T) %>%
  mutate(par=rownames(.))

ggplot(vs.df, aes(x=size, y=elpd, label=par,
                  ymin=elpd-2*elpd.se, ymax=elpd+2*elpd.se)) +
  geom_point() + geom_linerange() + 
  geom_vline(xintercept=n_var+0.5, linetype=3) + 
  geom_text(size=2.5, nudge_x=0.15, hjust=0)


proj <- project(vs, nv=n_var)
full.beta <- beta[,vs.df$vind[1+(1:n_var)]+1]
colnames(full.beta) <- vs.df$par[1+(1:n_var)]
mcmc_areas(as.matrix(proj)[,-1]) #+ coord_cartesian(c(-2,2))
mcmc_areas(full.beta)

pred <- proj_linpred(vs, xnew=Z[,-1], ynew=d.ls$Y[,s], nv=n_var, 
                     integrated=T, transform=T)
ggplot() + xlim(0,1) +
  # geom_point(aes(x=pred$pred, y=Y_vec), alpha=0.5) +
  # geom_abline(linetype=3) + stat_smooth(aes(x=pred$pred, y=Y_vec), method="lm") +
  geom_point(aes(x=1-exp(-pred$pred), y=as.numeric(Y_vec>0)), alpha=0.5) +
  stat_smooth(aes(x=1-exp(-pred$pred), y=as.numeric(Y_vec>0)), fullrange=T,
              method="glm", size=0.5, method.args=list(family="binomial")) +
  labs(x = 'prediction', y = 'y') 


y1_rep <- proj_predict(vs, xnew=Zt[44,-1,drop=F], nv=n_var, seed=7560)
qplot(as.vector(y1_rep), bins=25) +
  geom_vline(xintercept=Y_vec[44], color='red') +
  xlab('y1_rep')





