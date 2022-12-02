
modelfile <- "Stan_scripts/tg_hgf_gamma.stan"
savemodelname <- "M_hgf_gamma.rds"
savedataset <- "Refit_lkj_hgf_gamma" #without .rds 

library(rstan)
library(ggplot2)
library(Stack)
library(tidyr)
library(dplyr)

### load data
data_temp <-data_group

### prepare data
subjList <- unique(data_temp$ID)


removeSubjects<- {} #c(6,19,40)

if (!file.exists("fit_model"))  {
fit_model <- readRDS(paste("Model_results/", savemodelname, sep=""))
}
# choose random predictions from all
extract_parameters <-  rstan::extract(fit_model, pars = c("y_pred"))# , "mu_p", "om_good", "om_bad", "noise", "mu0"))
y_pred_mats <-extract_parameters$y_pred
# mu_p_mats <- extract_parameters$mu_p
# sub_lev_pars <- extract_parameters$om_good
niter <- dim(y_pred_mats)[1]
sample_per_subject <- 5
random_samples <- sample.int(niter, sample_per_subject) # c(3428, 4807, 5115, 6939,7430)


y_pred_mats_subset <- y_pred_mats[random_samples, ,]

# mu_p_mats_subset <- mu_p_mats[random_samples,]
# mu_p_mats_subset
for (sim in 1:sample_per_subject) {
    y_pred_mats_subset_temp <- y_pred_mats_subset[sim,,]

  ankk <- data_temp$ankk[data_temp$trials==1]
  
  numSub <- length(subjList)
  
  
  Tsubj <- as.vector(rep(0, numSub))
  for (ss in 1:numSub) {
    Tsubj[ss] <- max(data_temp$trials[data_temp$s == subjList[ss]]);
  }
  Tsubj_remove <- as.vector(rep(0, numSub))
  maxTrials <- max(Tsubj)
  
  transfer <- array(1,c(numSub, maxTrials))
  backtransfer <- array(1,c(numSub, maxTrials))
  trustee <- array(1,c(numSub, maxTrials))
  
  for (i in 1:numSub) {
    tmp <- subset(data_temp, data_temp$s==subjList[i])
    transfer[i,1:Tsubj[i]] = y_pred_mats_subset_temp[i,]
    backtransfer[i,1:Tsubj[i]] = tmp$backtransfer
    ### now use the prediction instead of the transfer
    trustee[i,1:Tsubj[i]] = tmp$trustee
   
  }
  
  sulpiride <- data_temp$drug[data_temp$trials==1] # sul
  
  serum <- data_temp$serum[data_temp$trials==1] # Nal
  
  
  dataList <- list(N = numSub, T = maxTrials, Tsubj = Tsubj, transfer = transfer, backtransfer=backtransfer,
                   trustee = trustee,  ankk =ankk, sulpiride= sulpiride, serum = serum, simulate = simulate,
                   Tsubj_remove = Tsubj_remove)
  
  
  
  
  # =============================================================================
  #### Running Stan #### 
  # =============================================================================
  rstan_options(auto_write = TRUE)
  options(mc.cores = 4)
  
  
  nIter     <-2000
  nChains   <-4
  nWarmup   <- floor(nIter/4)
  nThin     <- 1
  
  cat("Estimating", modelfile, "model... \n")
  
  startTime = Sys.time(); print(startTime)
  cat("Calling", nChains, "simulations in Stan... \n")
  
  fit_rl <- stan(modelfile, 
                 data    = dataList, 
                 chains  = nChains,
                 iter    = nIter,
                 warmup  = nWarmup,
                 thin    = nThin,
                 init    = "random",#"0", #inits_list,#,
                 seed    = 1450154637,
                 control = list(
                   adapt_delta = 0.99, max_treedepth = 20
                 )
  )
  # 
  # stepsize = 2.0,
  # max_treedepth = 10
  cat("Finishing", modelfile, "model simulation ... \n")
  endTime = Sys.time(); print(endTime)
  cat("It took",as.character.Date(endTime - startTime), "\n")
  cat("Saving in ", savedataset, "... \n")
  saveRDS(fit_rl, file = paste(savedataset, "_sample_no_", sim, ".rds", sep = ""))
  # saveRDS(mu_p_mats_subset, file = paste(savemodelname, "gen_params.rds"))
  # return(fit_rl)
 
  
}




fit_model <- readRDS(paste("Model_results/", savemodelname, sep=""))
par_set_gen <- rstan::extract(fit_model, pars = c("mu_p","sigma", "om_good", "om_bad", "noise", "mu0", "gam"))
all_pars =  matrix( nrow=0 , ncol=2 )
sample_no = 1
# random_samples <- c(3428, 4807, 5115, 6939,7430)
refit_model <- readRDS(paste(savedataset, "_sample_no_", sample_no,".rds", sep = ""))
# par_set <- rstan::extract(refit_model, pars = c("mu_p","sigma", "om_good", "om_bad", "noise", "mu0"))

drf_temp <- tibble(
  om_good = par_set_gen$om_good[random_samples[sample_no],],
  om_bad = par_set_gen$om_bad[random_samples[sample_no],],
  noise = par_set_gen$noise[random_samples[sample_no],],
  mu0 = par_set_gen$mu0[random_samples[sample_no],],
  loggam = log(par_set_gen$gam[random_samples[sample_no],]),
  om_good_rf = get_posterior_mean(refit_model, pars=c('om_good'))[,5],
  om_bad_rf =  get_posterior_mean(refit_model, pars=c('om_bad'))[,5],
  noise_rf =  get_posterior_mean(refit_model, pars=c('noise'))[,5],
  mu0_rf =  get_posterior_mean(refit_model, pars=c('mu0'))[,5],
  loggam_rf =  log(get_posterior_mean(refit_model, pars=c('gam'))[,5]),
  sample_no = sample_no
  )


for (sample_no in 2:5) {
  refit_model <- readRDS(paste(savedataset, "_sample_no_", sample_no,".rds", sep = ""))
 
  
  drf_temp_temp <- tibble(
    om_good = par_set_gen$om_good[random_samples[sample_no],],
    om_bad = par_set_gen$om_bad[random_samples[sample_no],],
    noise = par_set_gen$noise[random_samples[sample_no],],
    mu0 = par_set_gen$mu0[random_samples[sample_no],],
    loggam = log(par_set_gen$gam[random_samples[sample_no],]),
    om_good_rf = get_posterior_mean(refit_model, pars=c('om_good'))[,5],
    om_bad_rf =  get_posterior_mean(refit_model, pars=c('om_bad'))[,5],
    noise_rf =  get_posterior_mean(refit_model, pars=c('noise'))[,5],
    mu0_rf =  get_posterior_mean(refit_model, pars=c('mu0'))[,5],
    loggam_rf =  log(get_posterior_mean(refit_model, pars=c('gam'))[,5]),
    sample_no = sample_no
  )
  
  drf_temp <- Stack(drf_temp, drf_temp_temp)
}
saveRDS(drf_temp, file = paste(savedataset, "_pars_all.rds", sep = ""))


