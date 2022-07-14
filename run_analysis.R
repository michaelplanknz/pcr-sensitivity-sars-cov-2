##################################################################
######## ---------- PCR SENSITIVITY ANALYSIS --------- ###########
##################################################################

# Code to reproduce results from:
#
# Binny RN, Priest P, French N, Parry M, Lustig A, 
# Hendy SC, Steyn N, Ridings K, Vattiato G, Maclaren O, Plank MJ
# Sensitivity of RT-PCR tests for SARS-CoV-2 through time.
# 
# 23 May 2022

# Run R script to run analysis and save results. Note, model fitting 
# requires test_data_for_analysis_masked.csv and pcr_sensitivity_curve.stan

##################################################################
######## ---------- LOAD PACKAGES AND FUNCTIONS --------- ########
##################################################################

# Load required libraries
library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(rstan)
library(patchwork)
library(ggridges)
library(ggnewscale)
library(gridExtra)
library(binom)


# load incubation period estimates from the literature (function in EpiNow2):
# incubation_periods <- data.table(
#   as_reported = "5.06 (log SD 0.418)",
#   mean = 1.621, mean_sd = 0.0640,
#   sd = 0.418, sd_sd = 0.0691,
#   dist = "lognorm",
#   disease = "SARS-CoV-2",
#   source = "lauer",
#   url = "doi.org/10.7326/M20-0504"
# )

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Plot parameters:
tmax <- 50 
lwdth <- 1
L2 <- "dotted"
col1 <- "steelblue"
col2 <- "darkred"
col3 <- "orange"

#######################################################
######## ---------- ALL TEST DATA --------- ###########
#######################################################

#################
### LOAD DATA ###
#################

# Note, the test dataset has been masked to protect privacy of individuals. Dates of tests and symptom onset dates 
# have been replaced with a single field for the time of test in days since symptom onset. For fitting the model, 
# the day of symptom onset is set to 0 for all individuals. This gives the same results as those presented in the paper. 
# Load masked dataset:
dat0 <- data.table::fread("test_data_for_analysis_masked.csv", header=T) 
dat0$symp_onset_day <- rep(30, nrow(dat0)) # set symptom onset date to day 30 for all individuals. 
dat0$day_of_test <- dat0$symp_onset_day + dat0$days_since_onset

dat <- filter(dat0, included_main_analysis=="Include")

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

#################
# MODEL FITTING #
#################

dfy <- dat

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
nIter <- 8000
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "AllTests_MCMC_results.rds")

res <- rstan::extract(fit)

##########################################
# Parameter posterior summary statistics #
##########################################
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))

write.table(res_tab, "AllTests_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "AllTests_PCR_curve_summary.csv")

#####################################
## Figure 2 Posterior of PCR curve ##
#####################################

pt <- data.frame(top = apply(res$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res$p, 2, quantile, prob = 0.025),
                 y = apply(res$p, 2, median),
                 days = seq(0, tmax, 0.1))

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
res_day <- as.data.table(t(res$T_e)) %>% # T_e estimates from (nIter*nChains x P) to (P x nIter*nChains)
  melt(value.name = "inf_day", variable.name = "iter") # melt to long format

res_day[ , num_id := 1:.N, iter] # re-create num_id variable
res_day <- res_day[, .(iter = 1:.N, inf_day), num_id] # rename iterations 1:N, by num_id

res_day <- merge(dfy, res_day, by = "num_id", allow.cartesian = TRUE) # merge in test data by num_id

res_day[, x := day_of_test - round(inf_day)] # create new variable for 'days since infection' (infection time rounded to nearest day)

res_day <- res_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x) # average (over all individuals) sensitivity for each iteration and each 'day since infection'
][, .(top = quantile(pos, 0.975), # mean and 95% CI of sensitivity (over all iterations) for each 'days since infection'
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0] # keep only positive days since infection

# Plot empirical distribution and posterior distribution
fig2B <- 
  ggplot() +
  geom_ribbon(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, ymin = bottom, ymax = top), fill = add.alpha("grey", 0.5)) +
  ggplot2::geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean), lty = L2, size = lwdth, col="black") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha("steelblue", 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = pt, aes(x = days, y = y), lty = 1, col = "steelblue", size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5)), limits = c(0, tmax)) +
  labs(tag="B") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave(fig2B, filename = "Fig2B.jpg", height = 9, width = 11, units = "cm", dpi=600)

# For Fig S2, change incubation period parameters in pcr_sensitivity_curve.stan and rerun code above.

# ------ EMPIRICAL RESULT ----- #

tSinceOnset <- sort(unique(dfy$days_since_onset))
emp <- data.frame(method="", x=0, n=0, mean=0, lower=0, upper=0)

for (i in 1:length(tSinceOnset)){
  tmp <- filter(dfy, days_since_onset==tSinceOnset[i])
  emp[i,] <- binom.confint(sum(tmp$pcr_result), nrow(tmp), conf.level = 0.95, methods = "exact")
}
emp$tSinceOnset <- tSinceOnset

write.table(emp, "Empirical_sensitivity_summary.csv", row.names=F, sep=",")


fig2A <- ggplot() +
  #coord_cartesian(xlim = c(0, tmax)) +
  geom_ribbon(data = emp, inherit.aes = FALSE, aes(x = tSinceOnset,  ymin = lower, ymax = upper), fill = add.alpha("grey", 0.5)) +
  ggplot2::geom_line(data = emp, inherit.aes = FALSE, aes(x = tSinceOnset, y = mean), lty = L2, size = lwdth, col="black") +
  geom_point(data = emp, inherit.aes = FALSE, aes(x = tSinceOnset, y = mean)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  #ggtitle("PCR positivity over the course of infection") +
  #scale_fill_manual(values = pcols, name = "") +
  #scale_linetype_discrete(name = "") + 
  cowplot::theme_cowplot() + 
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since symptom onset") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  # xlim(0, tmax) +
  ggplot2::scale_x_continuous(breaks = seq(-10, 50, 5), limits = c(-10, 40)) +
  labs(tag="A") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave(fig2A, filename = "Fig2A.jpg", height = 9, width = 11, units = "cm", dpi=600)


#######################################################
######## ---------- SUBGROUP ANALYSIS --------- #######
#######################################################

####################################
### SUBSET BY IMMUNITY FOR DELTA ###
####################################

#### ---- DELTA, VACCINATED ----####
rm(dat, fit, res, res_day, pt, dfy)
dat <- filter(dat0, included_main_analysis=="Include", time_period=="Third", immunised=="Yes") # after 1 July 2021 and vaccinated

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)


### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "VaccinatedDelta_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "VaccinatedDelta_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "VaccinatedDelta_PCR_curve_summary.csv")


#### ---- DELTA, UNVACCINATED ----####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", time_period=="Third", immunised=="No") # after 1 July 2021 and unvaccinated

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "UnvaccinatedDelta_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "UnvaccinatedDelta_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "UnvaccinatedDelta_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###

# vaccinated, delta
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "Vaccinated"
# unvaccinated, delta
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "Unvaccinated"
ptComb <- rbind(pt, pt2)

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# vaccinated, delta
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "Vaccinated"

# unvaccinated, delta
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "Unvaccinated"

res_day <- rbind(res1_day, res2_day)

fig3A <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Vaccination (Delta)") +
  scale_colour_manual(values=c(col1, col2), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="A") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.6, 0.8), plot.margin = margin(7,9,7,9)) 


ggsave(fig3A, filename = "Fig3A.jpg", height = 9, width = 11, units = "cm", dpi=600)

####################################
######## SUBSET BY SOURCE ##########
####################################

#### ---- DOMESTIC ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, included_main_analysis=="Include", overseas=="No") # domestic

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Domestic_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Domestic_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Domestic_PCR_curve_summary.csv")


#### ---- OVERSEAS ----####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", overseas=="Yes") # overseas

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Overseas_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Overseas_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Overseas_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###

# domestic
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "Domestic"
# overseas
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "Overseas"
ptComb <- rbind(pt, pt2)

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# domestic
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "Domestic"

# overseas
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "Overseas"

res_day <- rbind(res1_day, res2_day)

fig3B <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Source") +
  scale_colour_manual(values=c(col1, col2), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="B") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.8), plot.margin = margin(7,12,7,9))

ggsave(fig3B, filename = "Fig3B.jpg", height = 9, width = 11, units = "cm", dpi=600)

####################################
######## SUBSET BY GENDER ##########
####################################

#### ---- FEMALE ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, included_main_analysis=="Include", gender=="Female") # female

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Female_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Female_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Female_PCR_curve_summary.csv")


#### ---- MALE ----####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", gender=="Male") # male

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Male_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Male_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Male_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###

# female
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "Female"
# male
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "Male"
ptComb <- rbind(pt, pt2)

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# female
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "Female"

# male
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "Male"

res_day <- rbind(res1_day, res2_day)

fig3C <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Gender") +
  scale_colour_manual(values=c(col1, col2), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="C") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.6, 0.8), plot.margin = margin(7,9,7,9))

ggsave(fig3C, filename = "Fig3C.jpg", height = 9, width = 11, units = "cm", dpi=600)

####################################
######## SUBSET BY AGE #############
####################################

#### ---- over 40s ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, included_main_analysis=="Include", age_group=="over40") # over 40

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Over40_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Over40_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Over40_PCR_curve_summary.csv")


#### ---- under 40s ----####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", age_group=="under40") # under 40s

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Under40_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Under40_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Under40_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###
# over40
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "> 40"
# under40
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "<= 40"
ptComb <- rbind(pt, pt2)

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# over40
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "> 40"

# under40
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "<= 40"

res_day <- rbind(res1_day, res2_day)

fig3D <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Age") +
  scale_colour_manual(values=c(col1, col2), labels=c(expression(phantom()<= 40), "> 40"), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="D") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.8), plot.margin = margin(7,12,7,9)) 

ggsave(fig3D, filename = "Fig3D.jpg", height = 9, width = 11, units = "cm", dpi=600)

########################################
######## SUBSET BY COMORBIDITIES #######
########################################

#### ---- WITH COMORBIDITIES ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, included_main_analysis=="Include", comorbidities==1) # With comorbidities

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Comorb_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Comorb_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Comorb_PCR_curve_summary.csv")


#### ---- NO COMORBIDITIES ---- ####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", comorbidities==0) # no comorbs

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "NoComorb_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "NoComorb_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "NoComorb_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###
# comorbidities
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "At least one"
# no comorbidities
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "None"
ptComb <- rbind(pt, pt2)

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# comorbidities
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "At least one"

# no comorbidities
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "None"

res_day <- rbind(res1_day, res2_day)

fig3E <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Comorbidities") +
  scale_colour_manual(values=c(col1, col2), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="E") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.6, 0.8), plot.margin = margin(7,9,7,9)) 

ggsave(fig3E, filename = "Fig3E.jpg", height = 9, width = 11, units = "cm", dpi=600)

########################################
########## SUBSET BY ETHNICITY #########
########################################

#### ---- MAORI ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, included_main_analysis=="Include", ethnicity=="Maori") # Maori

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Maori_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Maori_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Maori_PCR_curve_summary.csv")


#### ---- PACIFIC PEOPLES ---- ####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", ethnicity=="Pacific") # Pacific

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Pacific_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Pacific_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Pacific_PCR_curve_summary.csv")

#### ---- OTHER ETHNICITIES ---- ####
rm(dat, fit, res)
dat <- filter(dat0, included_main_analysis=="Include", ethnicity=="Other") # other

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy3 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "Other_MCMC_results.rds")

res <- rstan::extract(fit)
res3 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "Other_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "Other_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###
# Maori
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "Maori"
# Pacific
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "Pacific"
# nonMaoriPacific
pt3 <- data.frame(top = apply(res3$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res3$p, 2, quantile, prob = 0.025),
                  y = apply(res3$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt3$group <- "non-Maori & non-Pacific"
ptComb <- rbind(pt, pt2, pt3)
ptComb$group <- factor(ptComb$group, levels = c("Maori","Pacific","non-Maori & non-Pacific"))

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# Maori
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "Maori"

# Pacific
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "Pacific"

# non-Maori/non-Pacific (Other)
res3_day <- as.data.table(t(res3$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res3_day[ , num_id := 1:.N, iter]
res3_day <- res3_day[, .(iter = 1:.N, inf_day), num_id]

res3_day <- merge(dfy3, res3_day, by = "num_id", allow.cartesian = TRUE)

res3_day[, x := day_of_test - round(inf_day)]

res3_day <- res3_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res3_day$group <- "non-Maori & non-Pacific"

res_day <- rbind(res1_day, res2_day, res3_day)

fig3F <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Ethnicity") +
  scale_colour_manual(values=c(col3, col2, col1), labels=c("M\u101ori", "Pacific", "Other"), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col3, 0.5)) + 
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt3, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  labs(tag="F") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.65, 0.8), plot.margin = margin(7,12,7,9)) 

ggsave(fig3F, filename = "Fig3F.jpg", height = 9, width = 11, units = "cm", dpi=600)

############################################
############# FINAL PLOTS ##################
############################################

cairo_pdf(filename = "PCRcurve_Subgroups.pdf", height=12, width=9)
grid.arrange(fig3A, fig3B, fig3C, fig3D, fig3E, fig3F, nrow = 3, ncol=2, heights=unit(c(8, 8, 8), "cm"), widths=unit(c(10, 10), "cm"), as.table=T)
dev.off()

#######################################################################
############# SUPPLEMENTARY MATERIAL: TIME PERIODS ANALYSIS ###########
#######################################################################

##########################################
########## SUBSET BY TIME PERIOD #########
##########################################

#### ---- Pre June 2020 ---- ####
rm(dat, fit, res, res1, res2, res1_day, res2_day, pt, pt2, dfy, dfy1, dfy2)
dat <- filter(dat0, time_period=="First") # Pre 15 June 2020

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

### ---- Model Fitting ---- ###

dfy <- dat
dfy1 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "PreJune20_MCMC_results.rds")

res <- rstan::extract(fit)
res1 <- res


### ---- Parameter posterior summary statistics ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "PreJune20_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]


# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "PreJune20_PCR_curve_summary.csv")


#### ---- JUNE 2020 TO JUNE 2021---- ####
rm(dat, fit, res)
dat <- filter(dat0, time_period=="Second") # 15 June 2020 - 30 June 2021

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy2 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "June20-June21_MCMC_results.rds")

res <- rstan::extract(fit)
res2 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "June20-June21_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "June20-June21_PCR_curve_summary.csv")

#### ---- POST JULY 2021 ---- ####
rm(dat, fit, res)
dat <- filter(dat0, time_period=="Third") # Post 1 July 2021

# Generate new IDs:
newIDs <- data.frame(num_id=unique(dat$num_id), newID=seq(1,length(unique(dat$num_id))))
dat <- left_join(dat, newIDs, by="num_id") %>% select(., -num_id)
names(dat)[names(dat)=="newID"] <- "num_id"
setkey(dat, num_id)

patient_onsets <- unique(select(dat, num_id, symp_onset_day)) # day of symptom onset for each individual

setkey(patient_onsets, num_id) # sort by num_id
setkey(dat, num_id)

## ---- MODEL FITTING ---- ###

dfy <- dat
dfy3 <- dfy

# Generate list of data for stan
dat <- list()
dat$P <- patient_onsets[, .N]
dat$N <- dfy[, .N]
dat$day_of_test <- dfy$day_of_test
dat$test_result <- dfy$pcr_result %>% as.numeric()
dat$patient_ID <- dfy$num_id
dat$time_first_symptom <- patient_onsets$symp_onset_day

# Upper bound on time of infection, infection must occur before
# symptom onset or first positive PCR, whichever
# is first
dat$te_upper_bound <- dfy[, te_upper_bound := ifelse(
  any(day_of_test[pcr_result == TRUE] < symp_onset_day[pcr_result == TRUE]),
  min(day_of_test[pcr_result == TRUE & day_of_test < symp_onset_day]),
  symp_onset_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
#options(mc.cores = parallel::detectCores())
#mod <- rstan::stan_model("pcr_sensitivity_curve.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       cores = 4, 
                       iter = nIter,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      max_treedepth = 15))

saveRDS(fit, "PostJuly21_MCMC_results.rds")

res <- rstan::extract(fit)
res3 <- res

### ---- Parameter posterior summary statistics  ---- ###
# Table:
res_tab <- data.frame(Parameter=c("Cutpoint","Beta1","Beta2","Beta3","Inc1"), Posterior.median=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=median),
                      Posterior.mean=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=mean),
                      CI95Low=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.025),
                      CI95Up=apply(cbind(res$cutpoint,res$beta1,res$beta2,res$beta3,res$Inc[,1]), MARGIN=2, FUN=quantile, probs=0.975))
write.table(res_tab, "PostJuly21_parameter_summary.csv", row.names=F, sep=",")


# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, tmax, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(((nIter-1000)*4), length(p_vals)))
p_tab$iter <- rep(1:((nIter-1000)*4), length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value),
                 lower_95 = quantile(value, 0.025),
                 upper_95 = quantile(value, 0.975)),
             by = list(days_since_infection = diff)],
       file = "PostJuly21_PCR_curve_summary.csv")

### ---- Figure 3 Posterior of PCR curve ---- ###

# pre-June 2020 
pt <- data.frame(top = apply(res1$p, 2, quantile, prob = 0.975), 
                 bottom = apply(res1$p, 2, quantile, prob = 0.025),
                 y = apply(res1$p, 2, median),
                 days = seq(0, tmax, 0.1))
pt$group <- "-15/6/20"
# June 20 to June 21
pt2 <- data.frame(top = apply(res2$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res2$p, 2, quantile, prob = 0.025),
                  y = apply(res2$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt2$group <- "15/6/20-30/6/21"
# July 21 onwards
pt3 <- data.frame(top = apply(res3$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res3$p, 2, quantile, prob = 0.025),
                  y = apply(res3$p, 2, median),
                  days = seq(0, tmax, 0.1))
pt3$group <- "1/7/21-"
ptComb <- rbind(pt, pt2, pt3)
ptComb$group <- factor(ptComb$group, levels = c("-15/6/20","15/6/20-30/6/21","1/7/21-"))

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times
# Pre-June 2020
res1_day <- as.data.table(t(res1$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res1_day[ , num_id := 1:.N, iter]
res1_day <- res1_day[, .(iter = 1:.N, inf_day), num_id]

res1_day <- merge(dfy1, res1_day, by = "num_id", allow.cartesian = TRUE)

res1_day[, x := day_of_test - round(inf_day)]

res1_day <- res1_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res1_day$group <- "-15/6/20"

# June 20 to June 2021
res2_day <- as.data.table(t(res2$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res2_day[ , num_id := 1:.N, iter]
res2_day <- res2_day[, .(iter = 1:.N, inf_day), num_id]

res2_day <- merge(dfy2, res2_day, by = "num_id", allow.cartesian = TRUE)

res2_day[, x := day_of_test - round(inf_day)]

res2_day <- res2_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res2_day$group <- "15/6/20-30/6/21"

# July 21 onwards:
res3_day <- as.data.table(t(res3$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res3_day[ , num_id := 1:.N, iter]
res3_day <- res3_day[, .(iter = 1:.N, inf_day), num_id]

res3_day <- merge(dfy3, res3_day, by = "num_id", allow.cartesian = TRUE)

res3_day[, x := day_of_test - round(inf_day)]

res3_day <- res3_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]
res3_day$group <- "1/7/21-"

res_day <- rbind(res1_day, res2_day, res3_day)

figS6 <- 
  ggplot() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  ggtitle("Testing Period") +
  scale_colour_manual(values=c(col3, col2, col1), aesthetics = c("colour", "fill"), name="") +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt2, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col2, 0.5)) +
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col3, 0.5)) + 
  ggplot2::geom_ribbon(inherit.aes = FALSE, data = pt3, ggplot2::aes(x = days, y = y,  ymin = bottom, ymax = top), fill = add.alpha(col1, 0.5)) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(inherit.aes = FALSE, data = ptComb, aes(x = days, y = y, col=group), lty = 1, size = lwdth) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, color=group), lty = L2, size = lwdth) +
  ggplot2::labs(y = "RT-PCR sensitivity (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10)), limits = c(0,1)) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, tmax, 5))) +
  coord_cartesian(xlim = c(0, tmax)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.57, 0.8), plot.margin = margin(7,12,7,9)) 

ggsave(figS6, filename = "FigS6.jpg", height = 9, width = 11, units = "cm", dpi=600)