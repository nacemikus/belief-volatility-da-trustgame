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

data_beh$earn = 10 + data_beh$Investment*as.numeric(data_beh$Backtransfer)

data_beh %>% group_by(ID) %>% summarize(earn_sum = sum(earn),
                                        Treatment = Treatment[1],
                                        Genotype = Genotype[1]) %>%
  ggplot(aes(x = Treatment, y = earn_sum, colour = Treatment)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.5))+
  geom_boxplot(position = position_dodge(0.5), alpha = 0.5) + 
  theme_Publication() + facet_wrap(~Genotype)
data_beh$logserum_s = ave(log(data_beh$Serum+1), FUN = scale)
data_beh$Serum01 = data_beh$Serum/max(data_beh$Serum)

fit_model <- readRDS("/stan results/M_full.rds")

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

pars_all <- grep("^beta_",names(fit_model), value = T)
pars_beta <- grep("^beta_", names(fit_model), value = T)

traceplot(fit_model, pars = c('mu_p', 'lp__'), inc_warmup = TRUE)
traceplot(fit_model, pars = c('L_Omega'))
traceplot(fit_model, pars = c('sigma'))
traceplot(fit_model, pars = c('c'))
Par_check = rstan::extract(fit_model, pars = c('mu_p[2]', 'sigma[2]', 'mu_p[8]', 'sigma[8]','sess_w', 'w_mean', 'beta_nal', 'beta_ami', 'w', 'sess_w_pr'))
traceplot(fit_model, pars = grep("^beta_", names(fit_model), value = T), inc_warmup = TRUE)


g_m <- stan_plot(fit_model_int , pars = grep("^beta_", names(pars_all ), value = T), outer_level = 0.95) + theme_Publication()+
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



