# 00_fn.R
# Tim Szewczyk
#
# Assorted helper functions



#### variable selection functions ##############################################


#' Make new datasets for forward variable selection
#' @param full_data_base path and base filename for full dataset
#' @param out_base path and base filename for new datasets
#' @param opt_names character vector of already selected variables
#' @param type dataset type; "cv" for cross-validation subsets, "vs" for
#'   variable selection, or "pred" for optimal model that predicts whole map
#' @return Newly generated data sets and a confirmation message

make_next_datasets <- function(full_data_base, out_base, opt_names, mod_size, type="cv") {
  
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
  if(type %in% c("cv", "vs")) {  
    colnames(d.ls$V_) <- paste0("L_", colnames(d.ls$V_))
    colnames(d.i$V_)[-c(1,2)] <- paste0("L_", colnames(d.i$V_)[-c(1,2)]) 
  }
  
  if(type %in% c("cv", "vs")) {
    
    # identify parameters to add
    if("L_CnpyOpn" %in% opt_names) opt_names <- c(opt_names, "L_CnpyMxd")
    if("L_Pasture" %in% opt_names) opt_names <- c(opt_names, "L_Crop")
    if("R_grwnDD0" %in% opt_names) opt_names <- c(opt_names, "R_grwnDD0_sq")
    par_names_all <- c(colnames(d.ls$X), colnames(d.ls$V))
    if(is.null(opt_names)) {
      par_names <- "R_"
    } else {
      par_names <- par_names_all[!par_names_all %in% opt_names] 
    } 
    
    # make datasets
    for(i in seq_along(par_names)) {
      
      par_i <- par_names[i]
      if(par_i == "L_CnpyOpn") {
        par_i <- c(par_i, "L_CnpyMxd")
      } else if(par_i == "L_CnpyMxd") {
        next
      } else if(par_i == "L_Pasture") {
        par_i <- c(par_i, "L_Crop")
      } else if(par_i =="L_Crop") {
        next
      } else if(par_i == "R_grwnDD0") {
        par_i <- c(par_i, "R_grwnDD0_sq")
      } else if(par_i =="R_grwnDD0_sq") {
        next
      }
      
      d.ls_i <- d.ls
      d.i_i <- d.i
      
      # indexes
      d.ls_i$R <- sum(c(opt_names, par_i) %in% colnames(d.ls$X))
      d.ls_i$L <- sum(c(opt_names, par_i) %in% colnames(d.ls$V))
      
      # X
      d.ls_i$X <- d.ls$X[,which(colnames(d.ls$X) %in% c(opt_names, par_i)), drop=F]
      d.i_i$X <- d.i$X[,which(colnames(d.i$X) %in% c("id", "el", opt_names, par_i)), drop=F]
      
      # X_
      d.ls_i$X_ <- d.ls$X_[,which(colnames(d.ls$X_) %in% c(opt_names, par_i)), drop=F]
      d.i_i$X_ <- d.i$X_[,which(colnames(d.i$X_) %in% c("id", "el", opt_names, par_i)), drop=F]
      
      # V
      d.ls_i$V <- d.ls$V[,which(colnames(d.ls$V) %in% c(opt_names, par_i)), drop=F]
      d.i_i$V <- d.i$V[,which(colnames(d.i$V) %in% c("Plot_id", "el", opt_names, par_i)), drop=F]
      
      # V_
      d.ls_i$V_ <- d.ls$V_[,which(colnames(d.ls$V_) %in% c(opt_names, par_i)), drop=F]
      d.i_i$V_ <- d.i$V_[,which(colnames(d.i$V_) %in% c("Plot_id", "el", opt_names, par_i)), drop=F]
      
      d.f_i <- paste0(out_base, mod_size, "__", 
                      ifelse(type=="cv", 
                             paste0("k-", str_sub(full_data_base, -1, -1), 
                                    "_")), 
                      par_names[i])
      saveRDS(d.ls_i, paste0(d.f_i, "_ls.rds"))
      saveRDS(d.i_i, paste0(d.f_i, "_i.rds"))
      rstan::stan_rdump(ls(d.ls_i),
                        file=paste0(d.f_i, ".Rdump"),
                        envir=list2env(d.ls_i))
      
    }
    return(paste("Generated", 
                 length(dir("data/fwdSearch", 
                            paste0("^", str_split_fixed(out_base, "/", 3)[1,3], 
                                   mod_size, "__", 
                                   ifelse(type=="cv", 
                                          paste0("k-", str_sub(full_data_base, -1, -1), 
                                                 "_.*Rdump"))))), 
                        "new datasets as", 
                 paste0(out_base, mod_size, "__*")))
    
  } else if(type=="pred") {
    d.ls_i <- d.ls
    d.i_i <- d.i
    
    # indexes
    d.ls_i$R <- sum(opt_names %in% colnames(d.ls$X))
    d.ls_i$L <- sum(opt_names %in% colnames(d.ls$V))
    
    # X
    d.ls_i$X <- d.ls$X[,which(colnames(d.ls$X) %in% opt_names), drop=F]
    d.i_i$X <- d.i$X[,which(colnames(d.i$X) %in% c("id", "el", opt_names)), drop=F]
    
    # X_
    d.ls_i$X_ <- d.ls$X_[,which(colnames(d.ls$X_) %in% opt_names), drop=F]
    d.i_i$X_ <- d.i$X_[,which(colnames(d.i$X_) %in% c("id", "el", opt_names)), drop=F]
    
    # V
    d.ls_i$V <- d.ls$V[,which(colnames(d.ls$V) %in% opt_names), drop=F]
    d.i_i$V <- d.i$V[,which(colnames(d.i$V) %in% c("Plot_id", "el", opt_names)), drop=F]
    
    d.f_i <- paste0(out_base, "__opt_var_set")
    saveRDS(d.ls_i, paste0(d.f_i, "_ls.rds"))
    saveRDS(d.i_i, paste0(d.f_i, "_i.rds"))
    rstan::stan_rdump(ls(d.ls_i),
                      file=paste0(d.f_i, ".Rdump"),
                      envir=list2env(d.ls_i))
    return(paste("Generated dataset with optimal variables as", 
                 paste0(d.f_i, "*")))
  } else {
    return(cat("ERROR: type must be 'vs' or 'pred'."))
  }
  
}





#' Calculate loo metrics for model set
#' @param fit_dir Directory with cmdstan output; loo output will also be saved
#'   here if \code{save=T}
#' @param mod Model version; either \code{"Y"} or \code{"WY"}
#' @param mod_size Number of variables in model set
#' @param comp_all If \code{TRUE}, any saved loo .rds files for smaller models
#'   are also included in the comparison
#' @param save If \code{TRUE} (default), loo metrics are stored as a .rds file
#'   and comparison is stored as a .csv
#' @return "loo" with output from loo::loo() and "comp" with output from
#'   loo::loo_compare()
compare_models <- function(fit_dir, mod, mod_size, comp_all=F, save=T, type="cv") {
  
  library(tidyverse)
  
  # load (fit | mod, mod_size) and calculate loo metrics
  fit.f <- unique(str_sub(dir(fit_dir, paste0("^", mod, "_", mod_size, "__")), 1, -7))
  if(type=="cv") {
    # create lookup table for folds within variables
    fold.lu <- tibble(f=fit.f,
                      k=as.numeric(str_sub(str_split_fixed(fit.f, "k-", 2)[,2], 1, 1)),
                      vars=str_sub(str_split_fixed(fit.f, "k-", 2)[,2], 3, -1))
    # extract log_lik for each fold
    fit.ll <- map(fit.f, 
                  ~rstan::read_stan_csv(dir(fit_dir, paste0("^", .x), full.names=T)) %>%
                    rstan::extract(pars="log_lik")) %>%
      setNames(fit.f) %>%
      map(~matrix(.x$log_lik, nrow=dim(.x$log_lik)[1]))
    
    # concatenate log_lik across folds for each variable set and calculate loo
    fit.loo <- vector("list", n_distinct(fold.lu$vars)) %>% 
      setNames(., unique(fold.lu$vars))
    for(i in seq_along(fit.loo)) {
      ll.index <- match(fold.lu$f[fold.lu$vars==names(fit.loo)[i]], names(fit.ll))
      fit.loo[[i]] <- loo::loo(do.call("cbind", fit.ll[ll.index]))
    }
    names(fit.loo) <- paste0(mod, "_", mod_size, "__", names(fit.loo))
  } else {
    fit.out <- map(fit.f, 
                   ~rstan::read_stan_csv(dir(fit_dir, paste0("^", .x), full.names=T)))
    fit.loo <- map(setNames(fit.out, fit.f), loo::loo)
  }
  
  # load loo metrics for all models of smaller size 
  if(length(dir(fit_dir, paste0("loo_", mod, "_", mod_size)))>0) {
    sub_num <- 1 + length(dir(fit_dir, paste0("loo_", mod, "_", mod_size, ".*rds")))
    same.loo <- map(dir(fit_dir, 
                        paste0("loo_", mod, "_", mod_size, ".*rds"),
                        full.names=T), readRDS)
  } else {
    sub_num <- 1
    same.loo <- NULL
  }
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
  if(mod_size != 0 & mod_size != 14) {
    loo.comp <- loo::loo_compare(c(fit.loo, 
                                   unlist(same.loo, recursive=F),
                                   unlist(smaller.loo, recursive=F)))
  } else {
    loo.comp <- data.frame(elpd_diff=0, se_diff=0, 
                           elpd_loo=fit.loo[[1]]$estimates[1,1],
                           se_elpd_loo=fit.loo[[1]]$estimates[1,2],
                           p_loo=fit.loo[[1]]$estimates[2,1],
                           se_p_loo=fit.loo[[1]]$estimates[2,2],
                           looic=fit.loo[[1]]$estimates[3,1],
                           se_looic=fit.loo[[1]]$estimates[3,2],
                           row.names=ifelse(mod_size==0, 
                                            paste0(mod, "_", 0, "__R_"),
                                            names(fit.loo)))
  }
  
  if(save) {
    saveRDS(fit.loo, paste0(fit_dir, "/loo_", mod, "_", mod_size, "_", sub_num, ".rds"))
    write.csv(loo.comp, paste0(fit_dir, "/loo_", mod, "_", mod_size, "_", sub_num, ".csv"))
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






#### output processing functions ###############################################

#' Aggregate output from specified hierarchical models for specified parameters
#' @param d.f character vector of filenames for data used in the model; should
#'   be same length as \code{mods}, where the elements correspond with one
#'   another
#' @param mods model filename base; must match cmdstan output csv files (e.g.,
#'   "vs_Y" for "vs_Y_[1-12].csv")
#' @param pars_save Named vector of parameters to aggregate
#' @param out.dir "out"; directory where cmdstan output is stored
#' @return List with 1) [["full"]] = list of full stan output for each model,
#'   and 2) [["summaries"]] = list of posterior summaries for specified
#'   parameters, where each element is a dataframe corresponding to each
#'   parameter and includes both models.
aggregate_output <- function(d.f, mods, pars_save, out.dir="out") {
  
  # vs for cells with obs only; pred for predictions to new 1km2 cells
  type <- ifelse(any(grepl("vs", mods), grepl("fwdSearch", out.dir)), "vs", "pred")
  
  # storage objects
  out.pars <- list()
  out.stan <- vector("list", length(mods)) %>% setNames(mods)
  
  for(i in seq_along(mods)) {
    d.i <- readRDS(paste0(d.f[i], "_i.rds"))
    d.ls <- readRDS(paste0(d.f[i], "_ls.rds"))
    
    cat("\n", format(Sys.time(), "%X"), "-- Beginning model", mods[i])
    cat("\n  Reading", paste0(out.dir, "/", mods[i]))
    out.stan[[i]] <- dir(out.dir, paste0("^", mods[i], ".*_[0-9]"), full.names=T) %>%
      read_stan_csv()
    cat("\n  Summarizing", paste0(out.dir, "/", mods[i]))
    pars.all <- unique(str_split_fixed(names(out.stan[[i]]), "\\[", 2)[,1])
    pars <- pars.all[pars.all %in% pars_save]
    pars.trans <- grep("llambda", pars, value=T, ignore.case=T)
    
    ## CALCULATE SUMMARY STATISTICS --------------------------------------------
    # HDIs need to be calculated for lambdas separately on log vs natural scale
    out.mcmc <- c(map(pars, 
                      ~do.call('rbind', 
                               rstan::As.mcmc.list(out.stan[[i]], pars=.x))),
                  map(pars.trans, 
                      ~exp(do.call('rbind', 
                                   rstan::As.mcmc.list(out.stan[[i]], pars=.x))))) %>%
      setNames(c(pars, str_sub(pars.trans, 2, -1)))
    for(p in str_sub(pars.trans, 2, -1)) {
      dimnames(out.mcmc[[p]])[2][[1]] <- str_sub(dimnames(out.mcmc[[p]])[2][[1]], 2, -1)
    }
    out.50 <- map(out.mcmc, ~HDInterval::hdi(.x, 0.5) %>% t %>%
                    as_tibble(rownames="Parameter") %>%
                    rename(L25=lower, L75=upper))
    out.80 <- map(out.mcmc, ~HDInterval::hdi(.x, 0.8) %>% t %>%
                    as_tibble(rownames="Parameter") %>%
                    rename(L10=lower, L90=upper))
    out.90 <- map(out.mcmc, ~HDInterval::hdi(.x, 0.9) %>% t %>%
                    as_tibble(rownames="Parameter") %>%
                    rename(L05=lower, L95=upper))
    out.95 <- map(out.mcmc, ~HDInterval::hdi(.x, 0.95) %>% t %>%
                    as_tibble(rownames="Parameter") %>%
                    rename(L025=lower, L975=upper))
    out.mn <- map(pars, ~summary(out.stan[[i]], pars=.x, probs=0.5)$summary %>%
                    as_tibble(rownames="Parameter") %>%
                    rename(median=`50%`) %>%
                    mutate(Parameter=str_replace_all(Parameter, 
                                                     c("\\[" = ".",
                                                       "," = ".",
                                                       "]" = "")))) %>%
      setNames(pars)
    for(p in seq_along(pars.trans)) {
      p.trans.i <- str_sub(pars.trans[p], 2, -1)
      p.untrans.i <- pars.trans[p]
      out.mn[[p.trans.i]] <- tibble(
        Parameter=out.50[[p.trans.i]]$Parameter,
        mean=apply(out.mcmc[[p.trans.i]], 2, mean),
        sd=apply(out.mcmc[[p.trans.i]], 2, sd),
        median=apply(out.mcmc[[p.trans.i]], 2, median),
        n_eff=out.mn[[p.untrans.i]]$n_eff,
        Rhat=out.mn[[p.untrans.i]]$Rhat,
        se_mean=sd/sqrt(n_eff)
      ) %>% 
        select(Parameter, mean, se_mean, sd, median, n_eff, Rhat)
    }
    
    out.ls <- map(seq_along(c(pars, pars.trans)), 
                  ~full_join(out.50[[.x]], out.80[[.x]], by="Parameter") %>%
                    full_join(., out.90[[.x]], by="Parameter") %>%
                    full_join(., out.95[[.x]], by="Parameter") %>%
                    full_join(., out.mn[[.x]], by="Parameter") %>%
                    mutate(model=mods[i],
                           par=str_split_fixed(Parameter, "\\.", n=2)[,1])) %>%
      setNames(c(pars, str_sub(pars.trans, 2, -1)))
    
    ParNames <- c("intercept", 
                  paste0(colnames(d.ls$X)[-1]),
                  paste0(colnames(d.ls$V)))
    
    
    
    ## CREATE DATAFRAMES WITH SITE/PLOT/COVARIATE INFO -------------------------
    
    #### slopes -----------------------
    
    # beta
    out.pars$beta[[i]] <- out.ls$beta %>%
      mutate(ParName=ParNames) %>%
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # b
    out.pars$b[[i]] <- out.ls$b %>%
      mutate(cov=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(ParName=ParNames[as.numeric(cov)],
             sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
             genName=str_split_fixed(sppName, "_", 2)[,1]) %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # sigma_b
    out.pars$sig_b[[i]] <- out.ls$sigma_b %>%
      mutate(ParName=ParNames) %>%
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    # gamma & zeta
    out.pars$gamma[[i]] <- out.ls$gamma %>%
      mutate(spp=str_split_fixed(Parameter, "\\.", n=3)[,2]) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
             genName=str_split_fixed(sppName, "_", 2)[,1]) %>%
      mutate(Parameter=as.character(Parameter), model=as.character(model))

    out.pars$gamm_sig2[[i]] <- out.ls$gamma_sig2 %>%
      mutate(spp=str_split_fixed(Parameter, "\\.", n=3)[,2]) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
             genName=str_split_fixed(sppName, "_", 2)[,1]) %>%
      mutate(Parameter=as.character(Parameter), model=as.character(model))

    out.pars$gamma_Sigma[[i]] <- out.ls$gamma_Sigma %>%
      mutate(spp1=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp2=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(sppName1=d.i$tax_i$species[match(spp1, d.i$tax_i$sNum)],
             genName1=str_split_fixed(sppName1, "_", 2)[,1],
             sppName2=d.i$tax_i$species[match(spp2, d.i$tax_i$sNum)],
             genName2=str_split_fixed(sppName2, "_", 2)[,1]) %>%
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    out.pars$zeta[[i]] <- out.ls$zeta %>%
      mutate(plot=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
      mutate(plot=as.numeric(plot)) %>%
      arrange(plot) %>%
      mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    
    #### intensities ------------------
    
    # lLAMBDA 
    out.pars$lLAM[[i]] <- out.ls$lLAMBDA %>%
      mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    # LAMBDA
    out.pars$LAM[[i]] <- out.ls$LAMBDA %>%
      mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(site, spp) %>%
      mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    # llambda
    out.pars$llam[[i]] <- out.ls$llambda %>%
      mutate(plot=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(plot, spp) %>%
      mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    # lambda
    out.pars$lam[[i]] <- out.ls$lambda %>%
      mutate(plot=str_split_fixed(Parameter, "\\.", n=3)[,2],
             spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
      mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>%
      mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
      arrange(plot, spp) %>%
      mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
             Parameter=as.character(Parameter), model=as.character(model))
    
    
    #### FULL MODEL ONLY --------------
    
    if(type=="pred") {
      
      #### slopes -----------------------
      
      # B
      out.pars$B[[i]] <- out.ls$B %>%
        mutate(cov=str_split_fixed(Parameter, "\\.", n=3)[,2],
               gen=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
        mutate(ParName=ParNames[as.numeric(cov)],
               genName=d.i$tax_i$genus[match(gen, d.i$tax_i$gNum)]) %>% 
        mutate(Parameter=as.character(Parameter), model=as.character(model))
      
      # Sigma_B
      out.pars$Sig_B[[i]] <- out.ls$Sigma_B %>%
        mutate(cov=str_split_fixed(Parameter, "\\.", n=4)[,2],
               gen1=str_split_fixed(Parameter, "\\.", n=4)[,3],
               gen2=str_split_fixed(Parameter, "\\.", n=4)[,4]) %>%
        mutate(ParName=ParNames[as.numeric(cov)],
               gen1Name=d.i$tax_i$genus[match(gen1, d.i$tax_i$gNum)],
               gen2Name=d.i$tax_i$genus[match(gen2, d.i$tax_i$gNum)]) %>% 
        mutate(Parameter=as.character(Parameter), model=as.character(model))
      
      
      #### intensities ------------------
      
      # lLAMBDA_ 
      out.pars$lLAM[[i]] <- rbind(
        out.pars$lLAM[[i]], 
        out.ls$lLAMBDA_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
                 spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # LAMBDA_
      out.pars$LAM[[i]] <- rbind(
        out.pars$LAM[[i]], 
        out.ls$LAMBDA_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
                 spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # tot_lLAM, tot_lLAM_
      out.pars$tot_lLAM[[i]] <- rbind(
        out.ls$tot_lLAM %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)), 
        out.ls$tot_lLAM_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # tot_LAM, tot_LAM_
      out.pars$tot_LAM[[i]] <- rbind(
        out.ls$tot_LAM %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)), 
        out.ls$tot_LAM_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      ) 
      
      # tot_llam
      out.pars$tot_llam[[i]] <- out.ls$tot_llam %>%
        mutate(plot=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
        mutate(plot=as.numeric(plot)) %>%
        arrange(plot) %>%
        mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               Parameter=as.character(Parameter), model=as.character(model))
      
      # tot_lam
      out.pars$tot_lam[[i]] <- out.ls$tot_lam %>%
        mutate(plot=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
        mutate(plot=as.numeric(plot)) %>%
        arrange(plot) %>%
        mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               Parameter=as.character(Parameter), model=as.character(model))
      
      
      #### presence / absence ---------
      
      # prPres, prPres_
      out.pars$pP_R[[i]] <- rbind(
        out.ls$prPres %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
                 spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)),
        out.ls$prPres_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=3)[,2],
                 spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
          mutate(site=as.numeric(site), spp=as.numeric(spp)) %>%
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(site, spp) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # prPresL
      out.pars$pP_L[[i]] <- out.ls$prPresL %>%
        mutate(plot=str_split_fixed(Parameter, "\\.", n=3)[,2],
               spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
        mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>% 
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
        arrange(plot, spp) %>%
        mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               Parameter=as.character(Parameter), model=as.character(model)) 
      
      
      #### diversity ------------------
      
      # ShannonH, ShannonH_
      out.pars$H[[i]] <- rbind(
        out.ls$ShannonH %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)), 
        out.ls$ShannonH_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # ShannonH_L
      out.pars$H_L[[i]] <- out.ls$ShannonH_L %>%
        mutate(plot=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
        mutate(plot=as.numeric(plot)) %>%
        arrange(plot) %>%
        mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               Parameter=as.character(Parameter), model=as.character(model))
      
      
      #### richness -------------------
      
      # Rich, Rich_
      out.pars$S_R[[i]] <- rbind(
        out.ls$Rich %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X[site,"id"], el=d.i$X[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)), 
        out.ls$Rich_ %>%
          mutate(site=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
          mutate(site=as.numeric(site)) %>%
          arrange(site) %>%
          mutate(id=d.i$X_[site,"id"], el=d.i$X_[site,"el"],
                 Parameter=as.character(Parameter), model=as.character(model))
      )
      
      # RichL
      out.pars$S_L[[i]] <- out.ls$RichL %>%
        mutate(plot=str_split_fixed(Parameter, "\\.", n=2)[,2]) %>%
        mutate(plot=as.numeric(plot)) %>%
        arrange(plot) %>%
        mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
               Parameter=as.character(Parameter), model=as.character(model))
      
    }
    
    
    #### dispersion parameter ---------
    out.pars$disp[[i]] <- out.ls$disp_lam %>% 
      mutate(Parameter=as.character(Parameter), model=as.character(model))
    
    
    #### taxonomic bias ---------------
    if("D" %in% pars) {
      out.pars$D[[i]] <- out.ls$D %>% 
        mutate(spp=str_split_fixed(Parameter, "\\.", 2)[,2]) %>%
        mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)],
               Parameter=as.character(Parameter), model=as.character(model))
    }
    
    
    #### log likelihood ---------------
    if("log_lik" %in% pars) {
      if(type=="pred") {
        out.pars$LL[[i]] <- out.ls$log_lik %>% 
          mutate(Parameter=as.character(Parameter), model=as.character(model))
      } else {
        out.pars$log_lik[[i]] <- out.ls$log_lik %>%
          mutate(plot=str_split_fixed(Parameter, "\\.", n=3)[,2],
                 spp=str_split_fixed(Parameter, "\\.", n=3)[,3]) %>%
          mutate(plot=as.numeric(plot), spp=as.numeric(spp)) %>% 
          mutate(sppName=d.i$tax_i$species[match(spp, d.i$tax_i$sNum)]) %>%
          arrange(plot, spp) %>%
          mutate(id=d.i$V[plot,"Plot_id"], el=d.i$V[plot,"el"],
                 Parameter=as.character(Parameter), model=as.character(model)) 
      }
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
