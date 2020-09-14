# 00_fn.R
# Tim Szewczyk
#
# Assorted helper functions



#### variable selection functions ##############################################





#' Make new datasets for forward variable selection
#' @param full_data_base path and base filename for full dataset
#' @param out_base path and base filename for new datasets
#' @param opt_names character vector of already selected variables
#' @return Newly generated data sets and a confirmation message

make_next_datasets <- function(full_data_base, out_base, opt_names) {
  
  # load full dataset
  d.ls <- readRDS(paste0(full_data_base, "_ls.rds"))
  d.i <- readRDS(paste0(full_data_base, "_i.rds"))
  
  # make unambiguous column names
  colnames(d.ls$X) <- paste0("R_", colnames(d.ls$X))
  colnames(d.i$X)[-c(1,2)] <- paste0("R_", colnames(d.i$X)[-c(1,2)])
  colnames(d.ls$X_) <- paste0("R_", colnames(d.ls$X_))
  colnames(d.i$X_)[-c(1,2)] <- paste0("R_", colnames(d.i$X_)[-c(1,2)])
  colnames(d.ls$V) <- paste0("L_", colnames(d.ls$V))
  colnames(d.i$V)[-c(1,2)] <- paste0("L_", colnames(d.i$V)[-c(1,2)])
  colnames(d.ls$V_) <- paste0("L_", colnames(d.ls$V_))
  colnames(d.i$V_)[-c(1,2)] <- paste0("L_", colnames(d.i$V_)[-c(1,2)])
  
  # identify parameters to add
  par_names_all <- c(colnames(d.ls$X), colnames(d.ls$V))
  par_names <- par_names_all[!par_names_all %in% opt_names]
  
  # make datasets
  for(i in seq_along(par_names)) {
    
    d.ls_i <- d.ls
    d.i_i <- d.i
    
    # indexes
    d.ls_i$R <- sum(c(opt_names, par_names[i]) %in% colnames(d.ls$X))
    d.ls_i$L <- sum(c(opt_names, par_names[i]) %in% colnames(d.ls$V))
    
    # X
    d.ls_i$X <- d.ls$X[,which(colnames(d.ls$X) %in% c(opt_names, par_names[i])), drop=F]
    d.i_i$X <- d.i$X[,which(colnames(d.i$X) %in% c("id", "el", par_names[i])), drop=F]
    
    # X_
    d.ls_i$X_ <- d.ls$X_[,which(colnames(d.ls$X_) %in% c(opt_names, par_names[i])), drop=F]
    d.i_i$X_ <- d.i$X_[,which(colnames(d.i$X_) %in% c("id", "el", par_names[i])), drop=F]
    
    # V
    d.ls_i$V <- d.ls$V[,which(colnames(d.ls$V) %in% c(opt_names, par_names[i])), drop=F]
    d.i_i$V <- d.i$V[,which(colnames(d.i$V) %in% c("id", "el", par_names[i])), drop=F]
    
    # V_
    d.ls_i$V_ <- d.ls$V_[,which(colnames(d.ls$V_) %in% c(opt_names, par_names[i])), drop=F]
    d.i_i$V_ <- d.i$V_[,which(colnames(d.i$V_) %in% c("id", "el", par_names[i])), drop=F]
    
    d.f_i <- paste0(out_base, length(opt_names), "__", par_names[i])
    saveRDS(d.ls_i, paste0(d.f_i, "_ls.rds"))
    saveRDS(d.i_i, paste0(d.f_i, "_i.rds"))
    rstan::stan_rdump(ls(d.ls_i),
                      file=paste0(d.f_i, ".Rdump"),
                      envir=list2env(d.ls_i))
    
  }
  
  return(paste("Generated", length(par_names), "new datasets as", 
               paste0(out_base, "_", length(opt_names), "__*")))
  
}





#' Calculate loo metrics for model set
#' @param fit_dir Directory with cmdstan output; loo output will also be saved here if \code{save=T}
#' @param mod Model version; either \code{"Y"} or \code{"WY"}
#' @param mod_size Number of variables in model set
#' @param comp_all If \code{TRUE} (default), any saved loo .rds files for smaller models are also included in the comparison
#' @param save If \code{TRUE} (default), loo metrics are stored as a .rds file and comparison is stored as a .csv
#' @return "loo" with output from loo::loo() and "comp" with output from loo::loo_compare()
compare_models <- function(fit_dir, mod, mod_size, comp_all=T, save=T) {
  
  library(tidyverse)
  
  # load (fit | mod, mod_size) and calculate loo metrics
  fit.f <- unique(str_sub(dir(fit_dir, paste0("^", mod, "_", mod_size, "__")), 1, -7))
  fit.out <- map(fit.f, 
                 ~rstan::read_stan_csv(dir(fit_dir, paste0("^", .x), full.names=T)))
  fit.loo <- map(setNames(fit.out, fit.f), loo::loo)
  
  # load loo metrics for all models of smaller size 
  if(comp_all && mod_size > 1) {
    if(mod_size < 11) {
      smaller.loo <- map(dir(fit_dir, 
                             paste0("loo_", mod, "_[1-", mod_size-1, "].rds"),
                             full.names=T), readRDS)
    } else {
      smaller.loo <- c(
        map(dir(fit_dir, 
                paste0("loo_", mod, "_[1-9].rds"),
                full.names=T), readRDS),
        map(dir(fit_dir, 
                paste0("loo_", mod, "_1[0-", mod_size-11, "].rds"),
                full.names=T), readRDS)
      )
    }
  } else {
    smaller.loo <- NULL
  }
  
  # compare based on elpd
  loo.comp <- loo::loo_compare(c(fit.loo, unlist(smaller.loo, recursive=F)))
  
  if(save) {
    saveRDS(fit.loo, paste0(fit_dir, "/loo_", mod, "_", mod_size, ".rds"))
    write.csv(loo.comp, paste0(fit_dir, "/loo_", mod, "_", mod_size, ".csv"))
  }
  
  return(list(loo=fit.loo, comp=loo.comp))
}





#' Update file containing vector of optimal variables
#' @param out_dir Directory with 'opt_[WY,Y].rds'
#' @param mod Model to update (WY, Y)
#' @param mod_size Model size evaluated
#' @param comp Comparison from loo::loo_compare()
#' @return Status message and automatically updated and saved .rds file
update_opt_vars <- function(out_dir, mod, mod_size, comp) {
  
  library(stringr)
  
  opt.f <- paste0(out_dir, "/opt_", mod, ".rds")
  
  if(grepl(mod_size, rownames(comp)[1])) {
    opt_var <- str_split_fixed(rownames(comp)[1], "__", 2)[,2]
    saveRDS(c(readRDS(opt.f), opt_var), opt.f)
    msg <- paste(opt_var, "  was added to:  ", opt.f, "\n",
                 " updated size: ", length(readRDS(opt.f))-1, "\n",
                 " expected size:", mod_size)
    return(cat(msg))
  } else {
    msg <- paste("No new model is better!", "\n",
                 "Nothing added to", opt.f)
    return(cat(msg))
  }
}






#' Update rstanarm object with fitted hierarchical model output
#' @param fit_glmer rstanarm::glmer output
#' @param fit_hm stanfit object of hierarchical model output
#' @param d.ls List of data for fit_hm
#' @return Updated rstanarm object
update_rstanarm_shell <- function(fit_glmer, fit_hm, d.ls) {
  library(tidyverse); library(rstan); library(rstanarm)
  
  R <- d.ls$R
  L <- d.ls$L
  S <- d.ls$S
  
  
  # build lookup table to align parameters between glmer and hm
  names_glmer <- names(fit_glmer$stanfit@sim$samples[[1]])
  names_hm <- names(fit_hm@sim$samples[[1]])
  par_lu <- tibble(glmer=c("alpha[1]", 
                           paste0("beta[", 1:(R+L-1), "]"),
                           paste0("b[", 1:((R+L)*S), "]"),
                           "lp__"),
                   glmer_l=c(
                     rownames(fit_glmer$stan_summary)[1:((R+L)*(S+1))],
                     "log-posterior"),
                   hm=c(paste0("beta[", 1:(R+L), "]"),
                        paste0(rep(paste0("b[", 1:(R+L)), times=S), ",", 
                               rep(1:S, each=(R+L)), "]"),
                        "lp__"), 
                   hm_post=c(paste0("beta.", 1:(R+L)),
                        paste0(rep(paste0("b.", 1:(R+L)), times=S), ".", 
                               rep(1:S, each=(R+L))),
                        "lp__"), 
                        cov=c(rep(1:(R+L), times=S+1), NA))
  # posterior summaries
  sum_hm <- summary(fit_hm, 
                    probs=c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))$summary
  
  
  # extract and align predicted values
  llam.summary <- sum_hm[grepl("^llambda\\[", rownames(sum_hm)),]
  llam_hm.df <- tibble(llam=llam.summary[,1],
                       par=rownames(llam.summary),
                       site=str_split_fixed(str_remove(par, "llambda\\["), 
                                            ",", 2)[,1],
                       spp=str_split_fixed(str_remove(par, "]"), ",", 2)[,2]) %>%
    mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
    arrange(spp, site)
  
  
  
  
  # extract slope means and translate to effects parameterization for b's
  beta.summary <- sum_hm[grepl("^beta\\[", rownames(sum_hm)),]
  beta_hm.df <- tibble(beta=beta.summary[,1],
                       beta_md=beta.summary[,7],
                       beta_se=beta.summary[,2],
                       par=rownames(beta.summary)) %>%
    mutate(cov=row_number())
  b.summary <- sum_hm[grepl("^b\\[", rownames(sum_hm)),]
  b_hm.df <- tibble(b=b.summary[,1],
                    b_md=b.summary[,7],
                    b_se=b.summary[,2],
                    par=rownames(b.summary),
                    cov=str_split_fixed(str_remove(par, "b\\["), ",", 2)[,1],
                    spp=str_split_fixed(str_remove(par, "]"), ",", 2)[,2]) %>%
    mutate(cov=as.numeric(cov), spp=as.numeric(spp)) %>%
    arrange(spp, cov) %>%
    left_join(., select(beta_hm.df, beta, beta_md, cov), by='cov') %>%
    mutate(b_eff=b-beta, 
           b_eff_md=b_md-beta_md)
  
  
  # iteration indexes
  iter_hm <- length(fit_hm@sim$samples[[1]][[1]])
  iter_glmer <- length(fit_glmer$stanfit@sim$samples[[1]][[1]])
  iter_index <- (iter_hm-iter_glmer+1):iter_hm 
  
  
  # reference model: take structure of glmer and replace values with fitted hm
  fit_ref <- fit_glmer
  
  
  # update point estimates
  fit_ref$linear.predictors <- llam_hm.df$llam
  fit_ref$fitted.values <- exp(fit_ref$linear.predictors)
  fit_ref$coefficients[] <- c(beta_hm.df$beta_md, b_hm.df$b_eff_md)
  fit_ref$ses[] <- c(beta_hm.df$beta_se, b_hm.df$b_se)
  fit_ref$residuals <- fit_ref$y - fit_ref$fitted.values
  
  
  # update posterior summaries and samples
  for(i in 1:nrow(par_lu)) {
    
    # summaries
    sum_i <- sum_hm[rownames(sum_hm)==par_lu$hm[i],]
    if(grepl('b\\[', par_lu$hm[i])) {  # translate to effects parameterization
      ind <- (1:length(sum_i))[-c(2,length(sum_i)-1, length(sum_i))]
      sum_i_beta <- sum_hm[rownames(sum_hm)==paste0("beta[", 
                                                    par_lu$cov[i], "]"),]
      sum_i[ind] <- sum_i[ind] - sum_i_beta[ind]
    }
    fit_ref$stan_summary[rownames(fit_ref$stan_summary)==par_lu$glmer_l[i],] <- sum_i
    
    # samples
    for(j in 1:length(fit_glmer$stanfit@sim$samples)) {  # chains
      post_i <- fit_hm@sim$samples[[j]][[par_lu$hm_post[i]]]
      if(grepl('b\\.', par_lu$hm_post[i])) {  # translate to effects parameterization
        post_i_beta <- fit_hm@sim$samples[[j]][[paste0("beta.", par_lu$cov[i])]]
        post_i <- post_i - post_i_beta
      }
      fit_ref$stanfit@sim$samples[[j]][[par_lu$glmer[i]]] <- post_i[iter_index]
    }
  }
  
  return(fit_ref)
}











#### output processing functions ###############################################

#' Aggregate output from specified hierarchical models for specified parameters
#' @param d.i list of data used in the models
#' @param mods model filename base; must match cmdstan output csv files (e.g., "vs_Y" for "vs_Y_[1-12].csv")
#' @param pars_save Named vector of parameters to aggregate
#' @param out.dir "out"; directory where cmdstan output is stored
#' @return Updated rstanarm object
aggregate_output <- function(d.i, mods, pars_save, out.dir="out") {
  
  # vs for cells with obs only; pred for predictions to new 1km2 cells
  type <- ifelse(any(grepl("vs", mods)), "vs", "pred")
  
  # storage objects
  out.pars <- imap(pars_save, ~vector("list", length(mods)))
  out.stan <- vector("list", length(mods)) %>% setNames(mods)
  
  for(i in seq_along(mods)) {
    cat("\n", format(Sys.time(), "%X"), "-- Beginning model", mods[i])
    cat("\n  Reading", paste0(out.dir, "/", mods[i]))
    out.stan[[i]] <- dir(out.dir, paste0("^", mods[i], "_[0-9]"), full.names=T) %>%
      read_stan_csv()
    cat("\n  Summarizing", paste0(out.dir, "/", mods[i]))
    pars.all <- unique(str_split_fixed(names(out.stan[[i]]), "\\[", 2)[,1])
    pars <- pars.all[pars.all %in% pars_save]
    out <- summary(out.stan[[i]], pars=pars,
                   probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))$summary
    
    out.ls <- data.frame(Parameter=rownames(out),
                         mn=out[,1], se=out[,2],
                         q025=out[,4], q05=out[,5], q25=out[,6],
                         q50=out[,7],
                         q75=out[,8], q95=out[,9], q975=out[,10],
                         model=mods[i], Rhat=out[,12], row.names=NULL) %>%
      mutate(par=str_split_fixed(Parameter, "\\[", n=2)[,1]) %>%
      split(., .$par, drop=T)
    
    ParNames <- c("intercept", 
                  paste0(colnames(d.i$X)[-(1:2)], "_R"),
                  paste0(colnames(d.i$V)[-(1:2)], "_L"))
    
    # lLAMBDA (site)
    out.pars$lLAM[[i]] <- out.ls$lLAMBDA %>%
      mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    if(type=="pred") {
      out.pars$lLAM[[i]] <- rbind(
        out.pars$lLAM[[i]], 
        out.ls$lLAMBDA_ %>%
          mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                      ",", n=2)[,1],
                 spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
    }
    
    # llambda (plot)
    out.pars$llam[[i]] <- out.ls$llambda %>%
      mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(plot, spp) %>%
      mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    # prPres (site)
    if(type=="pred") {
    out.pars$pP_R[[i]] <- rbind(
      out.ls$prPres %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
               Parameter=as.character(Parameter), model=as.character(model)),
      out.ls$prPres_ %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
               Parameter=as.character(Parameter), model=as.character(model))
      )
    }
    
    # prPres (plot)
    out.pars$pP_L[[i]] <- out.ls$prPresL %>%
      mutate(plot=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                  ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>% 
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(plot, spp) %>%
      mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    # beta
    out.pars$beta[[i]] <- out.ls$beta %>%
      mutate(ParName=ParNames) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # B
    if("B" %in% pars) {
    out.pars$B[[i]] <- out.ls$B %>%
      mutate(cov=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                 ",", n=2)[,1],
             gen=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
      mutate(ParName=ParNames[as.numeric(cov)],
             genName=d.i$tax_i$genus[match(gen, d.i$tax_i$gNum)]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    }
    
    # b
    out.pars$b[[i]] <- out.ls$b %>%
      mutate(cov=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                 ",", n=2)[,1],
             spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
      mutate(ParName=ParNames[as.numeric(cov)],
             sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
             genName=str_sub(sppName, 1L, 4L)) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # sigma_b
    out.pars$sig_b[[i]] <- out.ls$sigma_b %>%
      mutate(ParName=ParNames) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # Sigma_B
    if("Sigma_B" %in% pars) {
      out.pars$Sig_B[[i]] <- out.ls$Sigma_B %>%
        mutate(cov=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                   ",", n=3)[,1], 
               gen1=str_split_fixed(Parameter, ",", n=3)[,2],
               gen2=str_remove(str_split_fixed(Parameter, ",", n=3)[,3], "]")) %>%
        mutate(ParName=ParNames[as.numeric(cov)],
               gen1Name=d.i$tax_i$genus[match(gen1, d.i$tax_i$gNum)],
               gen2Name=d.i$tax_i$genus[match(gen2, d.i$tax_i$gNum)]) %>% 
        mutate(Parameter=as.character(Parameter), model=as.character(model))
    }
    
    # dispersion parameter
    out.pars$disp[[i]] <- out.ls$disp_lam %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # D
    if("D" %in% pars) {
      out.pars$D[[i]] <- out.ls$D %>% 
        mutate(spp=str_remove(str_split_fixed(Parameter, "\\[", 2)[,2], "]")) %>%
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
               Parameter=as.character(Parameter), model=as.character(model))
    }
    
    if(type=="pred") {
    out.pars$H[[i]] <- rbind(
      out.ls$ShannonH %>%
        mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                    ",", n=2)[,1],
               spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
        mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
        arrange(site, spp) %>%
        mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
               Parameter=as.character(Parameter), model=as.character(model)), 
        out.ls$ShannonH_ %>%
          mutate(site=str_split_fixed(str_split_fixed(Parameter, "\\[", n=2)[,2],
                                      ",", n=2)[,1],
                 spp=str_sub(str_split_fixed(Parameter, ",", n=2)[,2], 1, -2)) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
    }
    
  }
  return(list(full=out.stan, summaries=map(out.pars, ~do.call('rbind', .))))
}




























#### deprecated functions ######################################################

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