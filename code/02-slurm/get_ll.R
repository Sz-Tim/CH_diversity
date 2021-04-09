## -----------------------------------------------------------------------------
## Operation Fourmis
## Extract log_lik columns from cmdstan output
## Tim Szewczyk
## -----------------------------------------------------------------------------


mod <- "LV_Y_0__"
fit_dir <- "out/fwdSearch/"
ll_dir <- "out/fwd_ll/"

fit.f <- dir(fit_dir, mod)

library(foreach); library(doParallel)

cl <- makeCluster(24)
registerDoParallel(cl)

foreach(i=seq_along(fit.f)) %dopar% {
  col_names <- scan(paste0(fit_dir, fit.f[i]), sep=",", skip=37, 
                    comment.char="#", nlines=1, what="character")
  ll_cols <- grep("log_lik", col_names)
  column_classes <- vector("list", length(col_names))
  column_classes[ll_cols] <- "character"
  
  full.scan <- scan(paste0(fit_dir, fit.f[i]), sep=",", skip=38, 
                    comment.char="#", what=column_classes)
  ll.df <- as.data.frame(do.call(cbind, full.scan[ll_cols]))
  names(ll.df) <- col_names[ll_cols]
  write.csv(ll.df, paste0(ll_dir, fit.f[i]), row.names=F)
}

stopCluster(cl)





