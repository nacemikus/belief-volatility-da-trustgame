# =============================================================================
# ### Model Summary and Diagnostics ####
# =============================================================================
# source('stan_utility.R')
library(tidyverse)
library(rstan)
library(loo)
library("rstanarm")
library("bayesplot")

source("theme_functions.r")
## compare models 
fit_model_hgf <- readRDS("/stan results/M_full.rds")

fit_model_hgf <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_gamma.rds")
fit_model_hgf_nogamma <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_simple_loo.rds")
fit_model_hgf_ka <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_ka.rds")
fit_model_rw <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_rw_gamma.rds")
fit_model_rw_nogamma <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_rw.rds")


M_hgf_ka <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_ka.rds")
M_hgf_gamma <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_gamma.rds")
M_rw <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_rw.rds")
M_hgf_gamma_wm <- readRDS("~/mnt/p/userdata/mikusn22/data/Sulpride_trustgame/Model_results/M_hgf_gamma_wm.rds")

loo.M_hgf_gamma <- loo::loo(M_hgf_gamma)
loo.M_hgf_ka<- loo::loo(M_hgf_ka)
loo.M_rw <- loo::loo(M_rw)

loo.M_hgf_gamma_wm <- loo::loo(M_hgf_gamma_wm)

cdata <- loo_compare(loo.M_hgf_gamma, loo.M_hgf_ka, loo.M_rw)
cdata%>% print()

M_hgf_gamma_wm %>% View

cdata <- loo_compare(loo.M_hgf_gamma, loo.M_hgf_gamma_wm)
cdata%>% print()

cData <- as.data.frame(cdata)
cData$model <- cData %>% rownames()
cData <- cData%>% mutate(model = c("HGF + \u0263", "HGF","RW"))

cData %>% saveRDS(file = "Model_comparison_single_pt.rds")
# overall model comparison ####
loo.fit_model_hgf<- loo::loo(fit_model_hgf)
loo.fit_model_hgf_nogamma<- loo::loo(fit_model_hgf_nogamma)
loo.fit_model_hgf_ka<- loo::loo(fit_model_hgf_ka)
loo.fit_model_rw<- loo::loo(fit_model_rw)
loo.fit_model_rw_nogamma<- loo::loo(fit_model_rw_nogamma)
loo_list_models = list(loo.fit_model_hgf= loo.fit_model_hgf,
                       loo.fit_model_hgf_nogamma = loo.fit_model_hgf_nogamma,
                       loo.fit_model_hgf_ka = loo.fit_model_hgf_ka,
                       loo.fit_model_rw = loo.fit_model_rw,
                       loo.fit_model_rw_nogamma= loo.fit_model_rw_nogamma
                       )
saveRDS(loo_list_models, file = "loo_list_models.rds")
cdata <- loo_compare(loo.fit_model_hgf, loo.fit_model_hgf_nogamma, loo.fit_model_hgf_ka, loo.fit_model_rw, loo.fit_model_rw_nogamma)
cdata%>% print()


cdata <- loo_compare(loo.fit_model_rw_simple, loo.fit_model_ka, loo.fit_model)

# trial by trial model comparison #####

loo.M_hgf_gamma <- loo::loo(M_hgf_gamma)
loo.M_hgf_ka<- loo::loo(M_hgf_ka)
loo.M_rw <- loo::loo(M_rw)

# loglik_rw_nogamma = extract(fit_model_rw_nogamma, pars = "log_lik")
loglik_rw = extract(M_rw, pars = "log_lik")
loglik_hgf = extract(M_hgf_gamma, pars = "log_lik")
loglik_hgf_ka = extract(M_hgf_ka, pars = "log_lik")

loglik_list = list(loglik_rw= loglik_rw,
                   loglik_hgf = loglik_hgf,
                   loglik_hgf_ka = loglik_hgf_ka)
# looic_hgf_bad = array(NA, dim = c(25,2))
group_trials= 1;
Trial_groups = c(1:24) %>% matrix(nrow = group_trials)

looic_good = array(NA, dim = c(24/group_trials,2,length(loglik_list)))
# looic_hgf_good = array(NA, dim = c(25,2))
looic_bad = array(NA, dim = c(24/group_trials,2,length(loglik_list)))
loo_name = c()

for (t in 1: c(24/group_trials)) {
  
  for (k in 1: length(loglik_list))  {
  loo.fit_model_temp<- loo::loo( loglik_list[[k]]$log_lik[,data_beh$Trial  %in% Trial_groups[,t] & data_beh$Trustee == "Good"] )
  loo_name[k] = names(loglik_list)[[k]]
  # loo.fit_model_hgf <- loo::loo( loglik_hgf$log_lik[,data_beh$Trial %in% Trial_groups[,t] & data_beh$Trustee == "Good"] )
  # loo.fit_model_rw_nogamma <- loo::loo( loglik_rw_nogamma$log_lik[,data_beh$Trial %in% Trial_groups[,t] & data_beh$Trustee == "Good"] )
  # loo.fit_model_hgf_nogamma  <- loo::loo( loglik_hgf_nogamma$log_lik[,data_beh$Trial %in% Trial_groups[,t] & data_beh$Trustee == "Good"] )
  # 
  looic_good[t,1,k] = loo.fit_model_temp$estimates[3,1]
  looic_good[t,2,k] =loo.fit_model_temp$estimates[3,2]
  
  loo.fit_model_temp <- loo::loo( loglik_list[[k]]$log_lik[,data_beh$Trial %in% Trial_groups[,t] & data_beh$Trustee == "Bad"] )
  
  looic_bad[t,1,k] = loo.fit_model_temp$estimates[3,1]
  looic_bad[t,2,k] = loo.fit_model_temp$estimates[3,2]
  }

}
# looic_good %>% glimpse

df = rbind(looic_good[,,1], looic_good[,,2], looic_good[,,3],
           looic_bad[,,1], looic_bad[,,2], looic_bad[,,3])
df %>% dim
df1 <- tibble(looic = df[,1], 
              se = df[,2], 
              type = c(rep(loo_name[1], 2*c(24/group_trials)), rep(loo_name[2], 2*c(24/group_trials)) ,rep(loo_name[3], 2*c(24/group_trials))), 
              trustee = rep(c(rep("good", c(24/group_trials)) ,rep("bad", c(24/group_trials))),3 ), 
              trials = rep(1:c(24/group_trials),6)) 
# df1 %>% saveRDS(file= "Model_comparison_across_trials2.rds")
df1 %>% glimpse
df1 %>% filter(type %in% c("loglik_rw", "loglik_hgf")) %>% ggplot(aes(x=trials, y = looic, group = type)) + 
  geom_line(aes(colour = type)) +
   # geom_smooth(aes(colour = type), se = F)+
  geom_ribbon(aes(ymin = looic -se, ymax = looic + se, fill = type), alpha = 0.2) +
  theme_Publication() + facet_wrap(~trustee)

df1  %>% filter(type %in% c("rw_simple", "hgf_simple"), trials != 25 )%>% ggplot(aes(x=trials, y = looic, group = type, colour = type)) + 
  geom_point(stat = "identity") + 
  geom_smooth()+
  # geom_errorbar(aes(ymin = looic -se, ymax = looic + se)) + 
  theme_Publication() + facet_wrap(~trustee)


# compare_models_df = tibble(looic_rw_good = )
loglik_rw$log_lik[,data_beh$Trial == 1 & data_beh$Trustee == "Good"] %>% glimpse
loo.fit_model_rw <- loo::loo( loglik_rw$log_lik[,data_beh$Trial == 1 & data_beh$Trustee == "Good"] )
loo.fit_model_rw %>% glimpse
loo.fit_model_rw$estimates[3,1] %>% glimpse()
pbma_BB_wts <- pseudobma_weights(cbind(loo.fit_model_hgf$pointwise[,"elpd_loo"],
                                       loo.fit_model_rw_nogamma$pointwise[,"elpd_loo"],
                                       loo.fit_model_hgf_nogamma$pointwise[,"elpd_loo"],
                                       loo.fit_model_hgf_ka$pointwise[,"elpd_loo"])
                                    )

# cdata <- loo_compare(loo.fit_model, loo.fit_model_kool, loo.fit_model_kool2)

cData <- as.data.frame(cdata)
cData$model <- cData %>% rownames()
cData <- cData%>% mutate(model = c("RW + \u0263", "HGF + \u0263", "RW", "HGF", "HGF + ka"), BB_wts = pbma_BB_wts)
saveRDS(cData, file = "Model_comparison.rds")
cData <-readRDS(file = "Model_comparison.rds")
cData %>% glimpse()
g_compare_models <- ggplot(cData, aes(x = model, y =  elpd_diff)) +  
  
  
  geom_errorbar(aes(ymin= elpd_diff - se_diff, ymax = elpd_diff+se_diff), width = 0.2, position = position_dodge(0.9)) + theme_Publication() +
  geom_bar(position = position_dodge(), stat = "identity") #+  scale_x_discrete(name="", limits = c("HGF + \u0263", "HGF", "RW")) + 
  
  ylab("Comparing expected\nlog predictive density")
  g_compare_models
  
  g_compare_models <- ggplot(cData %>% filter(model != ""), aes(x = model, y =  BB_wts)) +  
    
    
    # geom_errorbar(aes(ymin= elpd_diff - se_diff, ymax = elpd_diff+se_diff), width = 0.2, position = position_dodge(0.9)) + theme_Publication() +
    geom_bar(position = position_dodge(), stat = "identity") #+  scale_x_discrete(name="", limits = c("HGF + \u0263", "HGF", "RW")) + 
  
  ylab("Comparing expected\nlog predictive density")
  g_compare_models
## model comparison across drug groups #####
looic_good = array(NA, dim = c(25,2,3,2))
# looic_hgf_good = array(NA, dim = c(25,2))
looic_bad = array(NA, dim = c(25,2,3,2))

t = 1
g_ind = 1
for (t in 1:25) {
  for (g_ind in 1:2) {
    g = c("control", "sulpiride")[g_ind]
  loo.fit_model_rw <- loo::loo( loglik_rw$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Good"& data_beh$Trustee == g] )
  loo.fit_model_hgf <- loo::loo( loglik_hgf$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Good" & data_beh$Trustee == g])
  loo.fit_model_rw_nogamma <- loo::loo( loglik_rw_nogamma$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Good"  & data_beh$Trustee == g])
  
  looic_good[t,1,1,g_ind] = loo.fit_model_hgf$estimates[3,1]
  looic_good[t,2,1,g_ind] = loo.fit_model_hgf$estimates[3,2]
  looic_good[t,1,2,g_ind] = loo.fit_model_rw_nogamma$estimates[3,1]
  looic_good[t,2,2,g_ind] = loo.fit_model_rw_nogamma$estimates[3,2]
  looic_good[t,1,3,g_ind] = loo.fit_model_rw$estimates[3,1]
  looic_good[t,2,3,g_ind] = loo.fit_model_rw$estimates[3,2]
  
  loo.fit_model_rw <- loo::loo( loglik_rw$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Bad"  & data_beh$Trustee == g])
  loo.fit_model_hgf <- loo::loo( loglik_hgf$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Bad"  & data_beh$Trustee == g])
  loo.fit_model_rw_nogamma <- loo::loo( loglik_rw_nogamma$log_lik[,data_beh$Trial == t & data_beh$Trustee == "Bad"  & data_beh$Trustee == g])
  
  looic_bad[t,1,1,g_ind] = loo.fit_model_hgf$estimates[3,1]
  looic_bad[t,2,1,g_ind] = loo.fit_model_hgf$estimates[3,2]
  looic_bad[t,1,2,g_ind] = loo.fit_model_rw_nogamma$estimates[3,1]
  looic_bad[t,2,2,g_ind] = loo.fit_model_rw_nogamma$estimates[3,2]
  looic_bad[t,1,3,g_ind] = loo.fit_model_rw$estimates[3,1]
  looic_bad[t,2,3,g_ind] = loo.fit_model_rw$estimates[3,2]
  }
}

df = rbind(looic_good[,,1,1], looic_good[,,2,1], looic_good[,,3,1], 
           looic_bad[,,1,1], looic_bad[,,2,1], looic_bad[,,3,1],
           looic_good[,,1,2], looic_good[,,2,2], looic_good[,,3,2], 
           looic_bad[,,1,2], looic_bad[,,2,2], looic_bad[,,3,2])
df1 <- tibble(looic = df[,1], 
              se = df[,2], 
              type = c(c(rep("hgf", 50), rep("rw_simple", 50) ,rep("rw", 50)), 
                       c(rep("hgf", 50), rep("rw_simple", 50) ,rep("rw", 50))),
              trustee = rep(c(rep("good", 25) ,rep("bad", 25)),6 ), 
              trials = rep(1:25,12),
              treatment = factor(c(rep(1,150), rep(0,150)) , levels = c(0,1), labels = c("control", "sulpiride"))) 
df1 %>% View

df1 %>% filter(trials != 25) %>% ggplot(aes(x=trials, y = looic, group = type, colour = type)) + 
  geom_point(stat = "identity") + 
  geom_smooth(method = "lm")+
  # geom_errorbar(aes(ymin = looic -se, ymax = looic + se)) + 
  theme_Publication() + facet_wrap(vars(trustee,treatment))

# compare_models_df = tibble(looic_rw_good = )
loglik_rw$log_lik[,data_beh$Trial == 1 & data_beh$Trustee == "Good"] %>% glimpse
loo.fit_model_rw <- loo::loo( loglik_rw$log_lik[,data_beh$Trial == 1 & data_beh$Trustee == "Good"] )
loo.fit_model_rw %>% glimpse
loo.fit_model_rw$estimates[3,1] %>% glimpse()
pbma_BB_wts <- pseudobma_weights(cbind(loo.fit_model_rw_simple$pointwise[,"elpd_loo"],
                                       # loo.fit_model_rw$pointwise[,"elpd_loo"],
                                       loo.fit_model_ka$pointwise[,"elpd_loo"],
                                       loo.fit_model$pointwise[,"elpd_loo"]))

# cdata <- loo_compare(loo.fit_model, loo.fit_model_kool, loo.fit_model_kool2)

## inspect models ####

print(fit_model)
pairs(fit_model, pars = c('c'))#, 'mu_p[6]'))
pairs(fit_model_int, pars = c('mu_p', 'lp__'))#, 'mu_p[6]'))
pairs(fit_model, pars = c('L_Omega'))#, 'mu_p[6]'))
pairs(fit_model, pars = c('mu_p[2]','mu_p[7]', 'mu_p[8]','mu_p[9]','mu_p[10]','lp__'))
divergs <- get_divergent_iterations(fit_model_gen)

# look at rhats

Rhats <- rhat(fit_model)
Rhats %>% hist()
Rhats[is.na(Rhats)]

np_cp <- nuts_params(fit_model)
lp_cp <- log_posterior(fit_model)
head(lp_cp)
head(np_cp)
mcmc_parcoord(fit_model_array, np = np_cp)
ratios_cp <- neff_ratio(fit_model)
mcmc_neff(ratios_cp, size = 2)
rhats <- rhat(fit_model)
rhats %>% range()

mcmc_nuts_divergence(np_cp, lp_cp)
pars_all <- intersect(grep("^beta_",names(fit_model), value = T), grep("mix", names(fit_model), value = T, invert = T))

pars_all <- grep("^beta_",names(fit_model_hgf_ka), value = T)
pars_beta <- grep("^beta_", names(fit_model_hgf_ka), value = T)

traceplot(fit_model, pars = c('mu_p', 'lp__'), inc_warmup = TRUE)
traceplot(fit_model, pars = c('L_Omega'))
traceplot(fit_model, pars = c('sigma'))
traceplot(fit_model, pars = c('c'))
Par_check = rstan::extract(fit_model, pars = c('mu_p[2]', 'sigma[2]', 'mu_p[8]', 'sigma[8]','sess_w', 'w_mean', 'beta_nal', 'beta_ami', 'w', 'sess_w_pr'))
traceplot(fit_model, pars = grep("^beta_", names(fit_model), value = T), inc_warmup = TRUE)


g_m <- stan_plot(M_rw_wm   , pars = grep("^beta_", names(M_rw_wm   ), value = T), outer_level = 0.95) + theme_Publication()+
  geom_vline(xintercept=0) +  labs(subtitle = "80% and 95% CI")#+
# scale_y_discrete(name = "", limits = c("Ami \u03B7","Nal \u03B7","Ami \u0263","Nal \u0263","Ami \u03C9","Nal \u03C9"))+
#scale_x_continuous(name = "", breaks = c(-1,0, 1)) 
g_m





  g1 <- ggplot() + geom_histogram(aes(x=Par_check$lp__[divergs]))  + ylab(paste("parameter",p))
  g2 <- ggplot() + geom_histogram(aes(x=Par_check$lp__[!divergs]))  + ylab(paste("parameter",p))
 plot_grid(g1,g2)


check_hmc_diagnostics(fit_model)

parss<- names(fit_model)
stan_dens(fit_model, pars = parss[11:20])

loo.fit_model<- loo::loo(fit_model_gen)


loo.fit_model_rw <- loo::loo( M_rw_loo)

loo.fit_model_ka<- loo::loo(fit_model_ka)
cdata <- loo_compare(loo.fit_model, loo.fit_model_rw)
cdata%>% summary()
cdata <- loo_compare(loo.fit_model_rw_simple, loo.fit_model_ka, loo.fit_model)

pbma_BB_wts <- pseudobma_weights(cbind(loo.fit_model_rw_simple$pointwise[,"elpd_loo"],
                                       # loo.fit_model_rw$pointwise[,"elpd_loo"],
                                       loo.fit_model_ka$pointwise[,"elpd_loo"],
                                       loo.fit_model$pointwise[,"elpd_loo"]))

# cdata <- loo_compare(loo.fit_model, loo.fit_model_kool, loo.fit_model_kool2)

cData <- as.data.frame(cdata)
cData <- cData%>% mutate(model = c(1,2,3), BB_wts = pbma_BB_wts)
saveRDS(cData, file = "Model_comparison.rds")
cData <-readRDS(file = "Model_comparison.rds")
cData %>% glimpse()
g_compare_models <- ggplot(cData, aes(x = model, y =  elpd_diff)) +  
 
  
  geom_errorbar(aes(ymin= elpd_diff - se_diff, ymax = elpd_diff+se_diff), width = 0.2, position = position_dodge(0.9)) + theme_Publication() +
  geom_bar(position = position_dodge(), stat = "identity") +  scale_x_discrete(name="", limits = c("HGF + \u0263", "HGF", "RW")) + 
  
  ylab("Comparing expected\nlog predictive density")

g_compare_models <- ggplot(cData, aes(x = model, y =  elpd_diff)) +  
  
  
  geom_errorbar(aes(ymin= elpd_diff - se_diff, ymax = elpd_diff+se_diff), width = 0.2, position = position_dodge(0.9)) + theme_Publication() +
  geom_bar(position = position_dodge(), stat = "identity") +  scale_x_discrete(name="", limits = c("HGF + \u0263", "HGF", "RW")) + 
  
  ylab("Comparing expected\nlog predictive density")



