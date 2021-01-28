#' ---
#' output: github_document
#' ---


# LILE study: Trait analyses
# Peter Innes
# 9.29.20
install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos", type="binary")

#+ results=FALSE, message=FALSE, warning=FALSE
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(modelsummary)
library(magrittr)
library(reshape2)
library(rstan)
library(rstanarm)
library(arm) #for se.ranef()

#' Collection/environmental data
env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population), Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) #scale predictors

#' Seed weight data
sd_wt_data <- read.csv("data/sd_wt_data.csv", header = TRUE) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' Fruit fill data
ff_data <- read.csv("data/ff_data.csv", header = TRUE) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' Stem data
stem_data <- read.csv("data/stem_data.csv", header = T) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' Read-in yield data (composite of all the other traits, data frame was created in 01_traits_EDA.R script)
yield_df <- read.csv("data/yield_df.csv", header=T) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' #### FITTING LINEAR MODELS

# see: https://stackoverflow.com/questions/45788123/general-linear-mixed-effect-glmer-heteroscedasticity-modelling regarding modeling heteroscedasticity 

#' 1. Seed weight
# Fit a fixed-effect model. Population:block interaction accounts for subsampling at the 'plot' level (i.e. multiple plants of the same population grown together in the same plot). Sample sizes are mostly consistent (balanced design), so shrinkage would be minimal anyways if we fitted population as a random effect here.
fit_sd_wt <- lmer(sd_wt_50_ct ~ -1 + population + (1|block) + (1|population:block), data = sd_wt_data)# [-(252:253),]) 

fit_sd_wt2 <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data = sd_wt_data)

# Fit a model with geographic predictors. We will actually come back to this, in the trait-env script.
fit_sd_wt3 <- lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data = sd_wt_data, REML = FALSE) #Lat and Long as population-level predictors.

# Bayesian version of above model
#bayes_sd_wt <- stan_lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, iter=8000)

summary(fit_sd_wt) #Notice right-skew of residuals.
# compare pop means fixed vs random. Essentially the same
sw_mod_comps <- cbind(fixef(fit_sd_wt), coef(fit_sd_wt2)$population) %>%
  rename("coef(fit_sd_wt2)"="(Intercept)") %>%
  tibble::rownames_to_column("population")
sw_mod_comps$population <- gsub("population", "", sw_mod_comps$population)

#' 2. Fruit fill mod. Normal distro could suffice here.  Or possion? Might need bayesian to get it to fit
fit_ff <- lmer(good_fill ~ -1 + population + (1|block) + (1|population:block), data = ff_data)
ff_fit_summary <- summary(fit_ff)

ff_data$obsv <- 1:nrow(ff_data)
fit_ff2 <- glmer(good_sds ~ -1 + population + (1|block) + (1|population:block) + (1|obsv), family="poisson", data = ff_data) #singular fit. use Bayesian instead?
ff_ff2_summary <- summary(fit_ff2)

#' 3. Number of stems per plant mod. Normal distro
stems <- stem_data %>% dplyr::select(population,trt,block,row,plot,plant,num_of_stems) %>% unique()
fit_num_stems <- lmer(num_of_stems ~ -1 + population + (1|block) + (1|population:block), data=stems)
ns_fit_summary <- summary(fit_num_stems)

#' 4. Fruit per stem mod. For this trait and subsequent per-stem traits, we need an additional mofel term to account for multiple measurements taken from same plant
fit_fruit <- lmer(fruits ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) 
summary(fit_fruit)

fit_fruit_nb <- glmmTMB(fruits ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2") #try negative binomial instead

#' 5. Buds/flowers per stem mod
fit_bf <- lmer(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.
summary(fit_bf)

fit_bf_nb <- glmmTMB(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2")

fit_bf_p <- glmmTMB(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "poisson")
  
#' 6. Forks per stem mod
fit_forks <- lmer(forks ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.

fit_forks_g <- glmmTMB(forks ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "gaussian")

#' 7. Stem diameter mod
fit_stem_diam <- lmer(log(diam_stem) ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) #log transform to account for right skew
summary(fit_stem_diam)

#' 8. Capsule diameter mod
fit_caps_diam <- lmer(diam_caps ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.
summary(fit_caps_diam)

#' 9. Estimated yield mod. (see 01_traits_EDA.R for how we calculate yeild, plus exploratory analysis and fit assessments)
fit_log_yield <- lmer(log(EST_YIELD) ~ -1 + population + (1|block) + (1|population:block), data=yield_df) 
ly_fit_summary <- summary(fit_log_yield)

fit_log_yield2 <- lmer(log(EST_YIELD) ~ (1|population) + (1|block) + (1|population:block), data=yield_df) 

fit_yield <- lmer(EST_YIELD ~ -1 + population + (1|block) + (1|population:block), data=yield_df)
summary(fit_yield)
yield_fit_summary <- summary(fit_yield)

fit_yield2 <- lmer(EST_YIELD ~ (1|population) + (1|block) + (1|population:block), data=yield_df)
summary(fit_yield2)

fit_yield_glmer <- glmer(EST_YIELD ~ (1|population) + (1|block) + (1|population:block), family=Gamma(link="log"), data=yield_df) #fails to converge when population is a fixed effect.

summary(fit_yield_glmer)


#' #### Model DIAGNOSTICS

#' 1. Seed weight diagnostics
plot(fit_sd_wt) #looks okay except for a couple outliers in the mid-range.
qqnorm(resid(fit_sd_wt), main="Seed weight")
title(main="seed_weight") #looks okay, just a couple outliers 

# Histogram of residuals
sw_fit_summary <- summary(fit_sd_wt)
sw_resid <- data.frame(resid(fit_sd_wt))
sw_scaled_resid <- data.frame(sw_fit_summary$residuals) 
sw_fitted <- data.frame(fitted(fit_sd_wt))
ggplot(data=sw_resid, aes(x=resid.fit_sd_wt., y=stat(density))) +
  geom_histogram(bins = 100) #unscaled resids
ggplot(data=sw_scaled_resid, aes(x=sw_fit_summary.residuals, y=stat(density))) +
  geom_histogram(bins = 100) #how are residuals scaled? 

# Outliers
which(sw_scaled_resid>=5 | sw_scaled_resid <=-5) #obs 252 and 253 are the culprits. these observations are part of source:block groups with only 3 technical replicates instead of 4. 
which(sw_resid>=.005 | sw_resid<=-.005)

#' 2. Fruit fill diagnostics
plot(fit_ff) #doesn't look terrible...residuals tend to be larger at lower fruit fill values
qqnorm(resid(fit_ff), main = "Fruit fill") #fine?
ff_resid <- data.frame(resid=resid(fit_ff))
ggplot(data=ff_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks normal-ish. slightly left skewed.

plot(fit_ff2) #doesnt look much better than the first
qqnorm(resid(fit_ff2)) #about the same as first fit
ff2_resid <- data.frame(resid=resid(fit_ff2))
ggplot(data=ff2_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 100) #resids skewed left

#' 3. Num stems diagnostics
plot(fit_num_stems) 
qqnorm(resid(fit_num_stems)) #looks pretty good actually
ns_resid <- data.frame(resid=resid(fit_num_stems))
ggplot(data=ns_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks pretty good. a single positive outlier

#' 4. Fruit per stem diagnostics
plot(fit_fruit)
qqnorm(resid(fit_fruit)) # a bit curvy. might be okay. 
plot(fit_fruit2) #doesn't really change.
qqnorm(resid(fit_fruit2))

plot(predict(fit_fruit), na.omit(stem_data$fruits))
abline(a=0, b=1)

frt_resid <- data.frame(resid=resid(fit_fruit))
ggplot(data=frt_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks good? a few outliers?

fruit_nb_simres <- simulateResiduals(fit_fruit_nb)
plot(fruit_nb_simres)
plotQQunif(fruit_nb_simres)

#' 5. Buds/flowers per stem diagnostics
plot(fit_bf) #normal probably not best fit
qplot(fitted(fit_bf), resid(fit_bf)) #ggplot alternative to plot()

qqnorm(resid(fit_bf)) #curvy
# qqplot using broom.mixed and ggplot2
library(broom.mixed)
ggplot(broom.mixed::augment(fit_bf), aes(sample=.resid/sd(.resid))) + #scale to variance=1 
  stat_qq() +
  labs(title="buds and flowers per stem")

# histogram of residuals
bf_resid <- data.frame(resid=resid(fit_bf))
ggplot(data=bf_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) 

# diagnostics for nbinom and poisson fits
summary(fit_bf_nb)
bf_nb_simres <- simulateResiduals(fit_bf_nb)
plot(bf_nb_simres) #not sure how to interpret the DHARMa diagnostics
plotQQunif(bf_nb_simres)
plotResiduals(bf_nb_simres)

bf_pois_simres <- simulateResiduals(fit_bf_p)
plot(bf_pois_simres)


#' 6. Forks per stem diagnostics.
summary(fit_forks)
plot(fit_forks)
qqnorm(resid(fit_forks)) #not bad. normal might be okay.
forks_resid <- data.frame(resid=resid(fit_forks))
ggplot(data=forks_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks good. a few outliers?

forks_simres <- simulateResiduals(fit_forks_g)
plot(forks_simres)
plotResiduals(forks_simres)
plotQQunif(forks_simres)

#' 7. Stem diam diagnostics. Stem diam was log-transformed. Looks okay, a few right outliers.
plot(fit_stem_diam) 
qqnorm(resid(fit_stem_diam)) 
stem_diam_resid <- data.frame(resid=resid(fit_stem_diam))
ggplot(data=stem_diam_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50)

#' 8. Caps diam show a left skew
plot(fit_caps_diam)
qqnorm(resid(fit_caps_diam)) 
caps_resid <- data.frame(resid=resid(fit_caps_diam))
ggplot(data=caps_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) 
which(caps_resid<=-1) #is it reasonable to exclude these outliers? 

temp <- stem_data %>% dplyr::select(population, block, row, plant, diam_caps) %>% 
  na.omit()
rownames(temp) <- 1:nrow(temp)

#' 9. Yield diagnostics
# compare yield estimates from dif models
yield_mod_comps <- cbind(fixef(fit_yield), coef(fit_yield2)$population, exp(fixef(fit_log_yield)), exp(coef(fit_log_yield2)$population), exp(coef(fit_yield_glmer)$population)) 
colnames(yield_mod_comps) <- c("fy","fy2", "fly", "fly2","flg")
yield_mod_comps

# diagnostics
plot(fit_log_yield)
plot(fit_yield) #megaphone effect
 #log transform looks better
plot(fit_log_yield2) 
plot(fit_yield_glmer) #gamma distro also looks okay

qqnorm(resid(fit_log_yield)) #bout the same as untransformed
qqnorm(resid(fit_yield)) #actually looks fine
qqnorm(resid(fit_yield_glmer)) #not as good

yield_resid <- data.frame(resid=resid(fit_yield))
ggplot(data=yield_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks pretty good, slight right skew

ly_resid <- data.frame(resid=resid(fit_log_yield))
ggplot(data=ly_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #gets rid of the right-skew.

#' #### gather and summarise model diagnostics
quartz()
par(mfrow=c(3,3))
qqnorm(resid(fit_ff), main = "Fruit fill")
qqnorm(resid(fit_sd_wt), main="Seed weight")
qqnorm(resid(fit_num_stems), main = "Number of stems")
qqnorm(resid(fit_fruit), main = "Fruit per stem")
qqnorm(resid(fit_bf), main = "Buds and flowers per stem") #curvy
qqnorm(resid(fit_forks), main = "Forks per stem")
qqnorm(resid(fit_stem_diam), main = "log Stem diameter") 
qqnorm(resid(fit_caps_diam), main = "Capsule diameter") 
qqnorm(resid(fit_log_yield), main = "log Yield")


#' #### Gather ls-means of all traits together
# 1. Seed weight
sw_means <- as.data.frame(fixef(fit_sd_wt)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("seed_weight"="fixef(fit_sd_wt)") 
sw_means$population <- gsub("population", "",sw_means$population)

# sw confidence intervals
confint(fit_sd_wt)
sd_wt_CI <- as.data.frame(confint(fit_sd_wt)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("lwr"="2.5 %", "upr"="97.5 %")
sd_wt_CI$population <- gsub("population", "",sd_wt_CI$population)

# sw effect size plot
sw_results <- inner_join(sw_means, sd_wt_CI)
ggplot(data=sw_results, aes(x=seed_weight, y=reorder(population, seed_weight), xmin=lwr,xmax=upr)) +
  geom_point() +
  geom_errorbar() +
  ylab("Population") +
  xlab("Mean 50-count seed weight ± 95%CI")

# 2. Fruit fill
ff_means <- as.data.frame(fixef(fit_ff)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("fruit_fill"="fixef(fit_ff)") 
ff_means$population <- gsub("population", "", ff_means$population)

# 3. Number of stems
ns_means <- as.data.frame(fixef(fit_num_stems)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("number_of_stems"="fixef(fit_num_stems)")
ns_means$population <- gsub("population", "",ns_means$population)

# 4. Fruits per stem
fruit_means <- as.data.frame(fixef(fit_fruit)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("fruits_per_stem"="fixef(fit_fruit)")
fruit_means$population <- gsub("population", "",fruit_means$population)

# 5. Buds/flowers per stem
bf_means <- as.data.frame(fixef(fit_bf)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("buds_and_flowers_per_stem"="fixef(fit_bf)")
bf_means$population <- gsub("population", "",bf_means$population)

# 6. Forks per stem
fork_means <- as.data.frame(fixef(fit_forks)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("forks_per_stem"="fixef(fit_forks)")
fork_means$population <- gsub("population", "",fork_means$population)

# 7. Stem diameter
stem_diam_means <- as.data.frame(fixef(fit_stem_diam)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("stem_diameter"="fixef(fit_stem_diam)")
stem_diam_means$population <- gsub("population", "", stem_diam_means$population)

# 8. Capsule diameter
caps_diam_means <- as.data.frame(fixef(fit_caps_diam)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("capsule_diameter"="fixef(fit_caps_diam)")
caps_diam_means$population <- gsub("population", "",caps_diam_means$population)

# 9. Yield
yield_means <- as.data.frame(fixef(fit_log_yield)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("log_est_yield"="fixef(fit_log_yield)")
yield_means$population <- gsub("population", "", yield_means$population)
yield_CI <- as.data.frame(confint(fit_log_yield)) %>% #confidence intervals
  tibble::rownames_to_column(c("population"))
yield_CI$population <- gsub("population", "", yield_CI$population) #need to get rid of 'source' in front of all the source numbers
yield_results <- inner_join(yield_means, yield_CI) %>% 
  `colnames<-`(c("population","log_EST_YIELD","lwr","upr")) #%>% 
#mutate(exp_yield=exp(EST_YIELD), exp_lwr=exp(lwr), exp_upr=exp(upr))
ggplot(data=yield_results, aes(x=log_EST_YIELD, y=reorder(population, log_EST_YIELD), xmin=lwr,xmax=upr)) +
  geom_point() +
  geom_errorbar() +
  ylab("Population") +
  xlab("Mean log Estimated Yield (g seed/plant) ± 95%CI")

pop_trait_means <- full_join(sw_means, ff_means)
pop_trait_means <- full_join(pop_trait_means, ns_means)
pop_trait_means <- full_join(pop_trait_means, fruit_means)
pop_trait_means <- full_join(pop_trait_means, fork_means)
pop_trait_means <- full_join(pop_trait_means, bf_means)
pop_trait_means <- full_join(pop_trait_means, stem_diam_means)
pop_trait_means <- full_join(pop_trait_means, caps_diam_means)
pop_trait_means <- full_join(pop_trait_means, yield_means)
write.csv(pop_trait_means, file="data/pop_trait_means.csv", row.names = FALSE)


#' #### POST-HOC comparisons
# Tukey all pairwise comparisons for Yield
post_hoc <- emmeans(fit_yield, list(pairwise ~ population), adjust = "tukey", lmer.df = "satterthwaite")
contrasts <- post_hoc$`pairwise differences of population` %<>% 
  as.data.frame %>%
  filter(p.value <= 0.05) #keep only significant comparisons
dim(contrasts)

# Dunnett post-hoc test to compare everything versus just Maple Grove https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
post_hoc_MG <- emmeans(fit_yield, specs = trt.vs.ctrlk ~ population, ref = c("MAPLE_GROVE"), lmer.df = "satterthwaite") 
mg_contrasts <- post_hoc_MG$contrasts %<>% 
  as.data.frame %>%
  filter(p.value <=0.05) 

#' #### Trait pairwise correlations, performed with population-level trait means
# using the ggcorrplot package
library(ggcorrplot)

trait_corr_df <- pop_trait_means %>%
  mutate(ratio_bf_to_fruit=buds_and_flowers_per_stem/fruits_per_stem) %>%
  filter(!population=='APPAR')
corr_mat <- round(cor(trait_corr_df[,2:10], method=c("pearson"), use = "complete.obs"),4)
p_mat <- cor_pmat(trait_corr_df[,2:10])
head(p_mat)
quartz()
corr_plot <- ggcorrplot(corr_mat, hc.order = TRUE,type = "lower", lab = TRUE, p.mat=p_mat, insig = "blank", lab_size = 3, tl.cex = 7, show.legend = FALSE, title = "Correlations of least-square means")
#png("lsmeans_corrplot.png", width=8, height=7, res=300, units="in")
#corr_plot
#dev.off()

# Without the ggcorrplot package
get_upper_tri <- function(corr_mat){ # Get upper triangle of the correlation matrix
  corr_mat[lower.tri(corr_mat)]<- NA
  return(corr_mat)
}
upper_tri <- get_upper_tri(corr_mat)
melted_corr_mat <- melt(upper_tri)
ggplot(data = melted_corr_mat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  ylab("") +
  xlab("")
coord_fixed()

#' #### PCA of trait means
library(vegan)
library(ggvegan)
pop_trait_means <- full_join(pop_trait_means, dplyr::select(env_data, source, population))
temp <- pop_trait_means[,-1]
rownames(temp) <- temp[,10]
temp <- temp[,-10]
pop_trait_means <- temp

my_trait_rda <- rda(na.omit(pop_trait_means), scale = T) #scale everything bc they are vastly different units. Also don't use EST_YIELD (14th col) since that was a composite index of the other values
summary(my_trait_rda)
# Base R biplot with labels
biplot(my_trait_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_trait_rda, dis="sites", cex=0.5)

# ggvegan version
autoplot(my_trait_rda, arrows = TRUE, geom = "text", legend = "none") #basic version

# Extract df of PC scores to manually build a plot
pops_rda <- fortify(my_trait_rda, display='sites')
traits_rda <- fortify(my_trait_rda, display='species')
quartz()
ggplot() +
  geom_point(data=pops_rda, aes(x = PC1, y = PC2)) +
  geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_segment(data=traits_rda, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=traits_rda, 
            aes(x=PC1,y=PC2,label=Label,
                hjust="inward",vjust=0.75*(1-sign(PC2))), 
            color="red", size=4)
