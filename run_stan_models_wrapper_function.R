# =============================================================================
#### Info #### 
# =============================================================================
# hierarchical model

run_model_fit <- function(modelfile, savemodelname, nIter_set, nWarmup_set=0, simulate = 0) {

   
   if (nIter_set< 3) {
      sampling = "try"
   } else sampling = "sampling"
   
   if (nWarmup_set == 0) {nWarmup_set = ceiling(nIter_set/3)}
      
   
    # =============================================================================
    #### Construct Data #### 
    # =============================================================================
    # clear workspace
    library(rstan)
    library(ggplot2)
    # library(R.matlab)
    library(tidyr)
    library(dplyr)
### load data
   
   
# data_temp <- read.csv ('R_table_beh.csv', header = T, sep="," )
data_temp <- readRDS("Behavioural_data.rds") 

if (simulate ==1) print("simulating from the prior")
 ### prepare data
subjList <- unique(data_temp$ID)

removeSubjects<- {} 

### load data
ankk <- data_temp$ankk[data_temp$Trial==1]

numSub <- length(subjList)
swm_error = data_temp$error_sum_all[data_temp$Trial==1]
swm_error = ave(swm_error, FUN =scale)
mean(swm_error)
Tsubj <- as.vector(rep(0, numSub))
Tsubj_remove <- as.vector(rep(0, numSub))
for (ss in 1:numSub) {
  Tsubj[ss] <- max(data_temp$Trial[data_temp$ID == subjList[ss]]);
}

maxTrials <- max(Tsubj)

transfer <- array(1,c(numSub, maxTrials))
backtransfer <- array(1,c(numSub, maxTrials))
trustee <- array(1,c(numSub, maxTrials))

for (i in 1:numSub) {
  tmp <- subset(data_temp, data_temp$ID==subjList[i])
  transfer[i,1:Tsubj[i]] = tmp$transfer
  backtransfer[i,1:Tsubj[i]] = tmp$backtransfer
  trustee[i,1:Tsubj[i]] = tmp$trustee
  if (subjList[i]%in% removeSubjects) Tsubj_remove[i] = 1;
  
}

sulpiride <- data_temp$drug[data_temp$Trial==1] # sul

serum <- data_temp$serum[data_temp$Trial==1] # ser


dataList <- list(N = numSub, T = maxTrials, Tsubj = Tsubj, transfer = transfer, backtransfer=backtransfer,
                 trustee = trustee,  ankk =ankk, sulpiride= sulpiride, serum = serum, simulate = simulate, 
                 Tsubj_remove = Tsubj_remove, swm_error = swm_error)


if (length(removeSubjects)!=0) {
   print(paste("Removing subjects: "))
   print(removeSubjects)
}

    # =============================================================================
    #### Running Stan #### 
    # =============================================================================
    rstan_options(auto_write = TRUE)
    options(mc.cores = 4)
    
    if (sampling == "try") {
       cat("Trying..  \n")
       nIter     <-2
       nChains   <-1
       nWarmup   <- 1
       nThin     <- 1
    } else {
       nIter     <-nIter_set # 2000
       nChains   <- 4
       nWarmup   <- nWarmup_set
       nThin     <- 1
    }
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
                     adapt_delta = 0.95, max_treedepth = 15
                     )
                   )

    cat("Finishing", modelfile, "model simulation ... \n")
    endTime = Sys.time(); print(endTime)
    cat("It took",as.character.Date(endTime - startTime), "\n")
    cat("Saving in ", savemodelname, "... \n")
    saveRDS(fit_rl, file = savemodelname)
    
}
 

