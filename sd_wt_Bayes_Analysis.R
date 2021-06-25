#' ---
#' output: github_document
#' ---

#' title: Individual Project Analysis, Data Science Fall 2020
#' author: Peter Innes
#' date: Dec 12

#+ results=FALSE, message=FALSE, warning=FALSE
library(dplyr)
library(magrittr)
library(tidyr)
library(gridExtra)
library(lme4)
library(arm) #for se.ranef()
library(rstan) #for extract()
library(AICcmodavg)
library(ggplot2)
library(rstanarm)
options(mc.cores = parallel::detectCores())
options(max.print=1000000)
theme_set(theme_grey())


#' Read-in collection data for population-level predictors
env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>%
  dplyr::select(source,population,Lat,Long,Elev_m)
env_data$source %<>% as.factor
env_data$population %<>% as.factor
# scale the predictors
env_data$Lat_s <- scale(env_data$Lat)
env_data$Long_s <- scale(env_data$Long)
env_data$Elev_m_s <- scale(env_data$Elev_m)

#' Read-in seed weight data. Seed weight will be the only trait I analyze for the final project.
sd_wt_data <- read.csv("data/cleaned_LILE_yield_data_2013_seed_wt.csv", header = T)
sd_wt_data$source %<>% as.factor
sd_wt_data$block %<>% as.factor
sd_wt_data <- filter(sd_wt_data, is.na(notes)) %>% #filter out observations with fewer than 50 seeds (described in notes column). All other obs are weights of 50 seeds. 
  dplyr::select(!notes) %>% #don't need notes column anymore
  filter(!source %in% c(2,5,32,38)) %>% # exclude these sources bc they were found to be mostly 'Appar', which is already represented (source 41)
  left_join(dplyr::select(env_data,source,population,Lat_s, Long_s, Elev_m_s))



#' 2. MODEL FITS \n
#' Fit models with all combinations of geographic predictors for model selection. For now will do AICc insted of LOOIC because of computation costs.
# Three-factor models
fit_LaxLoxE <- lmer(sd_wt_50_ct ~ Lat_s*Long_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE) #This model doesn't converge. Not sure why, if the previous, more complex one did.

fit_LaLoE_LaxLo_LaxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + 
                               Lat_s:Long_s + Lat_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                               Lat_s:Long_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                              Lat_s:Elev_m_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LaxE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                                    Lat_s:Long_s + Lat_s:Elev_m_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                          Lat_s:Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                         Lat_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                         Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# Two-factor models
fit_LaxLo <- lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaxE <- lmer(sd_wt_50_ct ~ Lat_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LoxE <-lmer(sd_wt_50_ct ~ Long_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLo <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaE <-lmer(sd_wt_50_ct ~ Lat_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LoE <- lmer(sd_wt_50_ct ~ Long_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# Single-factor models
fit_Lat <- lmer(sd_wt_50_ct ~ Lat_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_Long <- lmer(sd_wt_50_ct ~ Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_Elev <- lmer(sd_wt_50_ct ~ Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# No predictors. Not sure if there is a point to including this. 
fit_no_pred <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

#' AICc
models <- list(fit_LaxLoxE, fit_LaLoE, fit_LaLoE_LaxLo_LaxE, fit_LaLoE_LaxLo_LoxE, fit_LaLoE_LaxE_LoxE, fit_LaLoE_LaxLo_LaxE_LoxE, fit_LaLoE_LaxLo, fit_LaLoE_LaxE, fit_LaLoE_LoxE, fit_LaxLo, fit_LaxE, fit_LoxE, fit_LaLo, fit_LaE, fit_LoE, fit_Lat, fit_Long, fit_Elev, fit_no_pred)
model_names <- c('LaxLoxE', 'LaLoE', 'LaLoE_LaxLo_LaxE', 'LaLoE_LaxLo_LoxE', 'LaLoE_LaxE_LoxE', 'LaLoE_LaxLo_LaxE_LoxE', 'LaLoE_LaxLo', 'LaLoE_LaxE', 'LaLoE_LoxE', 'LaxLo', 'LaxE', 'LoxE', 'LaLo', 'LaE', 'LoE', 'Lat', 'Long', 'Elev', 'No_predictor')
aictab(cand.set = models, modnames = model_names)

#' Fit a Bayesian version of the model with lowest AICc score, the Lat*Long model:
# bayes_LaxLo <- stan_lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, iter=8000)
# save(bayes_LaxLo, file="data/12_13_bayes_sw_LaxLo.RData")
load(file="data/12_13_bayes_sw_LaxLo.RData")
#' 3. MODEL CHECKS
#' Diagnostics for max-likelihood model 
summary(fit_LaxLo)
# ci_LaxLo <- confint(fit_LaxLo)

# check for homoskedasticity 
plot(fit_LaxLo)
# qqplot
qqnorm(resid(fit_LaxLo)) 
# histogram of residuals
sw_resid <- data.frame(resid=resid(fit_LaxLo))
ggplot(data=sw_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 100) 
#' The above diagnostics look reasonable, barring a couple outliers. It shouldn't be much of an issue given how much data we have. Overall, assumptions of equal variance and normality of residuals are met in the frequentist mixed model.
#' Check parameter estimates. The predictor slopes and population-level estimates of seed weight are of key interest
# ranef(fit_LaxLo)$population
coef(fit_LaxLo)$population

#' These population level estimates seem reasonable, however, they are substantially different than the results from a model fit with no predictors. Why does this happen? Does 'shrinkage' influence values if there are relatively fewer observations from one end of the latitude (or longitude) range versus the other?


#' Diagnostics for bayesian model.
prior_summary(bayes_LaxLo)
# population-level ranefs look reasonable—very close to the max-likelihood model. Rhat are all close to 1, n_eff is high across the board.
print(summary(bayes_LaxLo)[1:4,c("mean","sd","n_eff","Rhat")],digits=3)
print(summary(bayes_LaxLo)[284:333,c("mean","sd","n_eff","Rhat")],digits=3)

results <- print(summary(bayes_LaxLo)[,c("mean","sd","n_eff","Rhat")],digits=3)
samples <- extract(bayes_LaxLo$stanfit)

dim(results[284:319,]) #36 population 'b' params
dim(results[5:283,]) #279 population:block 'b' params
dim(results[320:327,]) #8 block 'b' params 
dim(results[5:327,]) #323 total 'b' params
str(samples$b) #326 params. Why 3 more than the in the results???? It's because the samples$b come from the extract() fxn, which is adding an additional 'b' param for each grouping scale e.g. 'b[(Intercept) population:block:_NEW_population:block]' 



# launch_shinystan(bayes_LaxLo) #again, no red flags in diagnostics of the bayesian model. Chains appear well mixed.
pp_check(bayes_LaxLo) # this looks good—replications match the data very closely.

#' 4. MODEL INFERENCES. For parameters or predictions that answer the scientific question, preferably presented graphically.

# Make point estimates of population means with uncertainty.
# First need a simple data frame of population coordinates to be used later on in model equation.
coords <- dplyr::select(env_data, source, population, Lat, Long, Lat_s, Long_s) %>%
  na.omit() %>%
  filter(!source %in% c(2,5,32,38)) %>%
  dplyr::select(population, Lat, Lat_s, Long, Long_s) %>%
  arrange(population) #Must sort alphabetically because that is the order in the bayes fit samples!

# Next derive posterior samples for population means
samples <- extract(bayes_LaxLo$stanfit)
pops <- samples$b[,281:316] #Samples of population deviations. Have to pick out just the population-level samples because block and population:block samples are also included in the samples$b vector. The correct subsetting here was tricky to figure out bc the samples$b vector had an extra param for each grouping scale.
pop_samples <- pops * NA

for ( i in 1:36 ) {
  pop_samples[,i] <- samples$alpha + pops[,i] + samples$beta[,1]*coords$Lat_s[i] + samples$beta[,2]*coords$Long_s[i] + samples$beta[,3]*coords$Lat_s[i]*coords$Long_s[i] #I'm not sure if this is the correct equation if I am interested in population-level means. Why not just do alpha plus population deviations? Wouldn't this match the coefs() of ML fit?
  #pop_samples[,i] <- samples$alpha + pops[,i]
}
# Now calculate mean and standard deviation of the posterior distributions for the county means.
pop_postmns <- rep(NA,36)
pop_postses <- rep(NA,36)
for ( i in 1:36 ) {
  pop_postmns[i] <- mean(pop_samples[,i])
  pop_postses[i] <- sd(pop_samples[,i])
}

#' Plot of posterior means and standard errors
bayes_pred_df <- data.frame(pop_mn=pop_postmns, pop_se=pop_postses)
bayes_pred_df <- cbind(bayes_pred_df, coords)

plot_sw_bayes <- ggplot(data=bayes_pred_df) +
  geom_abline(intercept=mean(samples$alpha),
              slope=mean(samples$beta[,1]),
              col="blue",
              lty=2) +
  geom_point(mapping=aes(x=Lat_s, y=pop_mn)) +
  geom_linerange(mapping=aes(x=Lat_s, ymin=pop_mn-pop_se, ymax=pop_mn+pop_se)) +
  labs(x="Latitude (scaled)", y="Estimated intercept in population l", title="Bayesian")

# extract population-level means from ML fit using coef()
ml_pred_df <- data.frame(coef(fit_LaxLo)$population,
                         se.ranef(fit_LaxLo)$population[,1])
names(ml_pred_df) <- c("pop_b0", "b1", "b2", "b3", "pop_b0_se")
ml_pred_df <- ml_pred_df %>%
  tibble::rownames_to_column("population") %>%
  inner_join(coords)

# Calculate intercepts for each population 
ml_pred_df$pop_b0 <- ml_pred_df$pop_b0 + ml_pred_df$b1*ml_pred_df$Lat_s + ml_pred_df$b2*ml_pred_df$Long_s + ml_pred_df$b3*ml_pred_df$Lat_s*ml_pred_df$Long_s
# Compare bayesian vs ML 
cbind(ml_pred_df$pop_b0, bayes_pred_df$pop_mn)
#' The calculated intercepts match very well. However I'm still unclear on what the intercepts I've calculated here represent for this model. Again, how is it different than just the population-level means i.e. the coefs() from lmer()?

# compare also to the model with no predictors, just for fun
nopred_mns <- coef(fit_no_pred)$population[1] %>% tibble::rownames_to_column("population")
mns <- full_join(nopred_mns, dplyr::select(ml_pred_df, population, pop_b0)) %>% full_join(dplyr::select(bayes_pred_df, population, pop_mn)) 

# take mean 'b' parameter values from model results and add alpha to compare to the ML coefs(). 
pop_dev_mns <- summary(bayes_LaxLo)[284:319, c("mean")]
alpha <- summary(bayes_LaxLo)[1,c("mean")]
pop_mns <- pop_dev_mns*NA
for (i in 1:length(pop_dev_mns) ) {
  pop_mns[i] <- pop_dev_mns[i] + alpha
}
cbind(pop_mns, coef(fit_LaxLo)$population[1]) #these also match well

#' Plot the comparison of intercepts
plot_sw_ML <- ggplot(data=ml_pred_df) +
  geom_abline(intercept=fixef(fit_LaxLo)[1], slope=fixef(fit_LaxLo)[2], col="blue", lty=2) +
  geom_point(mapping=aes(x=Lat_s, y=pop_b0)) +
  geom_linerange(mapping=aes(x=Lat_s, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se)) +
  labs(x="Latitude (scaled)", y="Estimated intercept in population l", title="Maximum likelihood")

grid.arrange(plot_sw_ML, plot_sw_bayes, nrow = 1)
#' They look nearly identical. Even uncertainty is hardly different. Would it make sense to have uncertainty shown for the trend line? \n

#' I'm not sure how to best visualize results for the interaction b/w two continuous predictors (e.g. Lat and Long). Perhaps an effect size plot for slope parameters (lat, long, and interaction)? \n

#' 5. CONCLUSIONS. Lewis flax seed weight has substantial biogeographic variation across its range in the Intermountain West and this variation is fit best by a model with an interaction of latitude and longitude. Seed weight tends to be greater in more southern and south-western areas of its range. These data (from a common garden experiment) suggest that seed weight has some degree of genetic basis and is a potentially adaptive trait. Future studies should test the hypothesis that climate variables such as temperature and precipitation, which also follow latitudinal and longitudinal clines, act as selective pressures on seed weight in Lewis flax. 
