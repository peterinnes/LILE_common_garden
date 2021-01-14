# LILE study: yield trait analyses
# Peter Innes
# 9.29.20

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(magrittr)
library(reshape2)
library(rstan)
library(rstanarm)
library(arm) #for se.ranef()

# collection/environmental data
env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population), Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) #scale predictors

# seed weight data
sd_wt_data <- read.csv("data/cleaned_LILE_yield_data_2013_seed_wt.csv", header = T) %>%
  filter(is.na(notes)) %>% #filter out rows with fewer than 50 seeds (described in notes column in spreadsheet, obs with standard 50 seeds have 'NA' in notes column)
  mutate(source=as.factor(source), block=as.factor(block)) %>%
  left_join(dplyr::select(env_data, source, population, Lat, Lat_s, Long, Long_s, Elev_m, Elev_m_s))

# alternatively, for rows with fewer than 50 seeds, scale the seed weight measurements up to approx 50-count value
#for(i in 1:length(sd_wt_data$sd_wt_50_ct)){
#  if(!is.na(sd_wt_data$num_seeds[i])){
#    sd_wt_data$sd_wt_50_ct[i] <- 1/sd_wt_data$num_seeds[i]*50*sd_wt_data$sd_wt_50_ct[i]
#  }
#}

sd_wt_data <- sd_wt_data %>% 
  dplyr::select(!notes) %>%
  filter(!source %in% c(2,5,22,32,38)) # exclude these sources bc they were found to be mostly 'Appar', which is already represented (source 41). Source 22 should be excluded as well—6 of  8 source 22 plants are Appar.

# fruit fill data
ff_data <- read.csv("data/cleaned_LILE_yield_data_2013_fruit_fill.csv", header=T) %>%
  mutate(source=as.factor(source), block=as.factor(block)) %>%
  filter(trt=="B") %>% #exclude 'trt A' (non-study/non-harvested plants)
  left_join(dplyr::select(env_data,source,population,Lat,Lat_s,Long,Long_s,Elev_m,Elev_m_s)) 

Appar <- ff_data %>% # list of individual Appar plants as noted in Frt Fill data. 6 plants from source 22 are listed as Appar. Brent says go ahead and exclude this entire source.
  filter(notes==c("Appar")) %>%
  filter(trt=="B") %>%
  dplyr::select(source,trt,block,row,plot,plant)

ff_data <- ff_data %>% 
  filter(!source %in% c(2,5,22,32,38)) #exclude mistaken Appar sources 
  
# stem and fruit data
stem_data <- read.csv("data/cleaned_LILE_yield_data_2013_stem_and_fruit.csv", header=T) %>%
  mutate(source=as.factor(source), block=as.factor(block)) %>%
  filter(trt=="B") %>%
  left_join(dplyr::select(env_data,source,population,Lat,Lat_s,Long,Long_s,Elev_m,Elev_m_s)) 

stem_data <- stem_data %>% 
  filter(!source %in% c(2,5,22,32,38)) #exclude mistaken Appar sources
  
# read-in yield data (composite of all the other traits, data frame was created in 01_traits_EDA.R script)
yield_df <- read.csv("data/yield_data.csv", header=T)

####
#### linear models
####
# see: https://stackoverflow.com/questions/45788123/general-linear-mixed-effect-glmer-heteroscedasticity-modelling regarding modeling heteroscedasticity 

#### seed weight
# fit a fixed-effect model, sample sizes are mostly consistent (balanced design), so shrinkage should be minimal anyways. 
fit_sd_wt <- lmer(sd_wt_50_ct ~ -1 + population + (1|block) + (1|population:block), data = sd_wt_data)# [-(252:253),]) 

fit_sd_wt2 <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data = sd_wt_data)

fit_sd_wt3 <- lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data = sd_wt_data, REML = FALSE) #Lat and Long as population-level predictors.

# bayesian 
bayes_sd_wt <- stan_lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, iter=8000) #Bayesian version of fit_sd_wt3

summary(fit_sd_wt) #notice right-skew of residuals.
# compare pop means fixed vs random
sw_mod_comps <- cbind(fixef(fit_sd_wt), coef(fit_sd_wt2)$population) %>%
  rename("coef(fit_sd_wt2)"="(Intercept)") %>%
  tibble::rownames_to_column("population")
sw_mod_comps$population <- gsub("population", "", sw_mod_comps$population)

coef_fit_sd_wt3 <- coef(fit_sd_wt3)$population[1] %>%
  tibble::rownames_to_column("population")
sw_mod_comps <- full_join(sw_mod_comps, coef_fit_sd_wt3)

# sd wt diagnostic plots
plot(fit_sd_wt) #looks okay except for a couple outliers in the mid-range.
qqnorm(residuals(fit_sd_wt)) #doesn't look perfect...but pretty good? 

#histogram of residuals
sw_fit_summary <- summary(fit_sd_wt)
sw_resid <- data.frame(resid(fit_sd_wt))
sw_scaled_resid <- data.frame(sw_fit_summary$residuals)
sw_fitted <- data.frame(fitted(fit_sd_wt))

ggplot(data=sw_resid, aes(x=resid.fit_sd_wt., y=stat(density))) +
  geom_histogram(bins = 100) #unsacled resids
ggplot(data=sw_scaled_resid, aes(x=sw_fit_summary.residuals, y=stat(density))) +
  geom_histogram(bins = 100) #how are residuals scaled? 

#outliers
which(sw_scaled_resid>=5 | sw_scaled_resid <=-5) #obs 252 and 253 are the culprits. these observations are part of source:block groups with only 3 technical replicates instead of 4. 
which(sw_resid>=.005 | sw_resid<=-.005)

# pull out fixed effect estimates (mean seed weights for each population)
sd_wt_means <- as.data.frame(fixef(fit_sd_wt)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("seed_weight"="fixef(fit_sd_wt)") 
sd_wt_means$population <- gsub("population", "",sd_wt_means$population)
# confidence intervals
sd_wt_CI <- as.data.frame(confint(fit_sd_wt)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("lwr"="2.5 %", "upr"="97.5 %")
sd_wt_CI$population <- gsub("population", "",sd_wt_CI$population)

# gather means and CIs for effect size plot
sd_wt_df <- inner_join(sd_wt_means, sd_wt_CI)
ggplot(data=sd_wt_df, aes(x=mean_sd_wt_50_ct, y=reorder(population, mean_sd_wt_50_ct), xmin=lwr,xmax=upr)) +
  geom_point() +
  geom_errorbar() +
  ylab("population") +
  xlab("mean 50-count seed weight ± 95%CI")

# note about the rest of the traits. up to 8 plants per accession we measured for these traits, sometimes fewer, and sometimes there were multiple plants from the same block. We may have to choose only one plant per block? Also, for now, only looking at 'study plants' i.e. 'trt B' and thus our sample is non-random, and we have source as a fixed effect. Also I think we should maybe use a Poisson distribution for fruit fill and fruit per stem, stem per plant. 

#### Fruit fill. Normal distro could suffice here.  Or possion? Might need bayesian to get it to fit
fit_ff <- lmer(good_fill ~ -1 + population + (1|block) + (1|population:block), data = ff_data)
ff_fit_summary <- summary(fit_ff)

ff_data$obsv <- 1:nrow(ff_data)
fit_ff2 <- glmer(good_sds ~ -1 + population + (1|block) + (1|population:block) + (1|obsv), family="poisson", data = ff_data) #singular fit. use Bayesian instead?
ff_ff2_summary <- summary(fit_ff2)

# diagnostics
plot(fit_ff) #doesn't look terrible...residuals tend to be larger at lower fruit fill values
qqnorm(resid(fit_ff)) #fine?
ff_resid <- data.frame(resid=resid(fit_ff))
ggplot(data=ff_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks normal-ish. slightly left skewed.

plot(fit_ff2) #doesnt look much better than the first
qqnorm(resid(fit_ff2)) #about the same as first fit
ff2_resid <- data.frame(resid=resid(fit_ff2))
ggplot(data=ff2_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 100) #resids skewed left
# gather means
ff_means <- as.data.frame(fixef(fit_ff)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("fruit_fill"="fixef(fit_ff)") 
ff_means$population <- gsub("population", "", ff_means$population)

#### Number of stems per plant. Again normal distro probably will suffice
stems <- stem_data %>% dplyr::select(population,trt,block,row,plot,plant,num_of_stems, Lat_s) %>% unique()
fit_num_stems <- lmer(num_of_stems ~ -1 + population + (1|block) + (1|population:block), data=stems)
ns_fit_summary <- summary(fit_num_stems)

# diagnostics
plot(fit_num_stems) 
qqnorm(resid(fit_num_stems)) #looks pretty good actually
ns_resid <- data.frame(resid=resid(fit_num_stems))
ggplot(data=ns_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks pretty good. a single positive outlier

# fit with glmer.
stems$plant <- 1:nrow(stems) #hack that Brent showed us to deal with overdispersion in glmer
fit_num_stems2 <- glmer(num_of_stems ~ (1|population) + (1|block) + (1|population:block) + (1|plant), family = poisson(link="log"), data=stems) #doesn't converge. try Bayesian?
plot(fit_num_stems2) #reverse megaphone-ish. this actually looks worse than with normal distro
summary(fit_num_stems2) 

# gather means
ns_means <- as.data.frame(fixef(fit_num_stems)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("number_of_stems"="fixef(fit_num_stems)")
ns_means$population <- gsub("population", "",ns_means$population)
# confidence intervals
ns_CI <- as.data.frame(confint(fit_num_stems)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("lwr"="2.5 %", "upr"="97.5 %")
ns_CI$population <- gsub("population", "",ns_CI$population)

# gather means and CIs for effect size plot
ns_df <- inner_join(ns_means, ns_CI)
ggplot(data=ns_df, aes(x=mean_stem_number, y=reorder(population, mean_stem_number), xmin=lwr,xmax=upr)) +
  geom_point() +
  geom_errorbar() +
  ylab("population") +
  xlab("mean stems per plant ± 95%CI")

#### Fruit per stem. 
fit_fruit <- lmer(fruits ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) #nested structure to account for multiple measurements taken from same plant, and in come cases, sub-sampling of populations in blocks.
summary(fit_fruit)

#stem_data$obsv <- 1:nrow(stem_data)
fit_fruit2 <- glmer(fruits ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), family=poisson(link="log"), data=stem_data) #try poisson instead...singular fit.

# diagnostics
plot(fit_fruit)
qqnorm(resid(fit_fruit)) # a bit curvy
plot(fit_fruit2) #doesn't really change.
qqnorm(resid(fit_fruit2))

frt_resid <- data.frame(resid=resid(fit_fruit))
ggplot(data=frt_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks good? a few outliers?

# get ls means
fruit_means <- as.data.frame(fixef(fit_fruit)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("fruits_per_stem"="fixef(fit_fruit)")
fruit_means$population <- gsub("population", "",fruit_means$population)

#### Buds/flowers per stem
fit_bf <- lmer(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data)
summary(fit_bf)
plot(fit_bf) #normal distro might not be best fit

qqnorm(resid(fit_bf)) #curvy
bf_resid <- data.frame(resid=resid(fit_bf))
ggplot(data=bf_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) 

# get ls means
bf_means <- as.data.frame(fixef(fit_bf)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("buds_and_flowers_per_stem"="fixef(fit_bf)")
bf_means$population <- gsub("population", "",bf_means$population)

#### Forks per stem
fit_forks <- lmer(forks ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data)
fit_forks2 <- lmer(forks ~ Lat_s + (1|population) + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data)
fit_forks3 <- lmer(forks ~ (1|population) + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data)

# diagnostics
summary(fit_forks)
plot(fit_forks)
qqnorm(resid(fit_forks))
forks_resid <- data.frame(resid=resid(fit_forks))
ggplot(data=forks_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks good. a few outliers?

# gather means
fork_means <- as.data.frame(fixef(fit_forks)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("forks_per_stem"="fixef(fit_forks)")
fork_means$population <- gsub("population", "",fork_means$population)
fork_means <- inner_join(fork_means, dplyr::select(env_data, population, Lat))

ranef_forks <- coef(fit_forks2)$population %>%
  tibble::rownames_to_column(var=c("population")) %>%
  rename("ranef"="(Intercept)") %>%
  select(population, ranef) %>%
  inner_join(select(env_data, population, Lat))
ranef_forks3 <- coef(fit_forks3)$population %>%
  tibble::rownames_to_column(var=c("population")) %>%
  rename("ranef3"="(Intercept)") %>%
  select(population, ranef3) %>%
  inner_join(select(env_data, population, Lat))

# plot just the fixed effects vs lat
# quick plots
ggplot(data=filter(fork_means, !population %in% c('MAPLE_GROVE')), aes(x=Lat, y=fixef)) +
  geom_point() #

# compare fixed vs random effects
fork_means <- full_join(ranef_forks, fork_means) %>%
  full_join(ranef_forks3) %>%
  pivot_longer(c(ranef,ranef3,fixef), names_to=c("effect_type"), values_to ="forks_per_stem")
ggplot(data=filter(fork_means, !population %in% c('MAPLE_GROVE')), aes(x=Lat, y=forks_per_stem, color=effect_type)) +
  geom_point()
cor.test(ranef_forks3$ranef, ranef_forks3$Lat)

#### stem diameter
fit_stem_diam <- lmer(diam_stem ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) #not picking up plant grouping??/
summary(fit_stem_diam)

# diagnostics
plot(fit_stem_diam) #looks good, maybe an outlier or two.
qqnorm(resid(fit_stem_diam))
stem_diam_resid <- data.frame(resid=resid(fit_stem_diam))
ggplot(data=stem_diam_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks okay, a number of right-skewed outliers

# gather lsmeans
stem_diam_means <- as.data.frame(fixef(fit_stem_diam)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("stem_diameter"="fixef(fit_stem_diam)")
stem_diam_means$population <- gsub("population", "", stem_diam_means$population)

#### capsule diameter
fit_caps_diam <- lmer(diam_caps ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term, leaving this term out fixes issue
summary(fit_caps_diam)

plot(fit_caps_diam)
qqnorm(resid(fit_caps_diam))
caps_resid <- data.frame(resid=resid(fit_caps_diam))
ggplot(data=caps_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) 
which(caps_resid<=-1)

# gather lsmeans
caps_diam_means <- as.data.frame(fixef(fit_caps_diam)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("capsule_diameter"="fixef(fit_caps_diam)")
caps_diam_means$population <- gsub("population", "",caps_diam_means$population)

#### gather lsmeans of all traits together
pop_trait_means <- full_join(sd_wt_means, ff_means)
pop_trait_means <- full_join(pop_trait_means, ns_means)
pop_trait_means <- full_join(pop_trait_means, fruit_means)
pop_trait_means <- full_join(pop_trait_means, fork_means)
pop_trait_means <- full_join(pop_trait_means, bf_means)
pop_trait_means <- full_join(pop_trait_means, stem_diam_means)
pop_trait_means <- full_join(pop_trait_means, caps_diam_means)
write.csv(pop_trait_means, file="data/pop_trait_means.csv", row.names = FALSE)

#### analyzing yield. (see 01_traits_EDA.R for how we calculate yeild, plus exploratory analysis and fit assessments

# yield linear model 
fit_yield <- lmer(EST_YIELD ~ -1 + population + (1|block) + (1|population:block), data=yield_df)
summary(fit_yield)
yield_fit_summary <- summary(fit_yield)
fit_yield2 <- lmer(EST_YIELD ~ (1|population) + (1|block) + (1|population:block), data=yield_df)
summary(fit_yield2)

fit_log_yield <- lmer(log(EST_YIELD) ~ -1 + population + (1|block) + (1|population:block), data=yield_df) 
ly_fit_summary <- summary(fit_log_yield)

fit_log_yield2 <- lmer(log(EST_YIELD) ~ (1|population) + (1|block) + (1|population:block), data=yield_df) 

fit_yield_glmer <- glmer(EST_YIELD ~ (1|population) + (1|block) + (1|population:block), family=Gamma(link="log"), data=yield_df) #fails to converge when population is a fixed effect.
summary(fit_yield_glmer)

# compare yield estimates from dif models
cbind(exp(coef(fit_yield_glmer)$population), exp(coef(fit_log_yield2)$population), exp(fixef(fit_log_yield)), fixef(fit_yield), coef(fit_yield2)$population)

# diagnostics
plot(fit_yield) #slight megaphone effect. but maybe okay?
plot(fit_log_yield) #log transform looks a little better
plot(fit_log_yield2) 
plot(fit_yield_glmer) #gamma distro also looks good

qqnorm(resid(fit_yield)) #actually looks fine
qqnorm(resid(fit_log_yield)) #bout the same as untransformed
qqnorm(resid(fit_yield_glmer))

yield_resid <- data.frame(resid=resid(fit_yield))
ggplot(data=yield_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #looks pretty good, slight right skew

ly_resid <- data.frame(resid=resid(fit_log_yield))
ggplot(data=ly_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 50) #gets rid of the right-skew.

# gather lsmeans and confidence intervals
yield_means <- as.data.frame(fixef(fit_yield)) %>%
  tibble::rownames_to_column(c("population")) %>%
  rename("est_yield"="fixef(fit_yield)")
yield_means$population <- gsub("population", "", yield_means$population)
# confidence intervals
yield_CI <- as.data.frame(confint(fit_yield)) %>%
  tibble::rownames_to_column(c("population"))
yield_CI$population <- gsub("population", "", yield_CI$population)#need to get rid of 'source' in front of all the source numbers

yield_df <- inner_join(yield_means, yield_CI) %>% 
  `colnames<-`(c("population","EST_YIELD","lwr","upr")) #%>% 
  #mutate(exp_yield=exp(EST_YIELD), exp_lwr=exp(lwr), exp_upr=exp(upr))

quartz()
ggplot(data=yield_df, aes(x=EST_YIELD, y=reorder(population, EST_YIELD), xmin=lwr,xmax=upr)) +
  geom_point() +
  geom_errorbar() +
  ylab("population") +
  xlab("mean estimated seed yield (g/plant) ± 95%CI")

#### Tukey post-hoc test all pairwise comparisons
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

#### Trait correlations. Brent wants me to calculate correlations based on the population means (right now they are calculated at the individual plant (ish) level. Good to aggregate into means, or no?
install.packages("ggcorrplot")

trait_means_df <- inner_join(ns_means, fruit_means) %>%
  inner_join(ff_means) %>%
  inner_join(sd_wt_means) %>%
  inner_join(stem_diam_means) %>%
  inner_join(caps_diam_means) %>%
  inner_join(bf_means) %>%
  inner_join(fruit_means) %>%
  inner_join(fork_means) %>%
  inner_join(yield_means) %>% 
  mutate(ratio_bf_to_fruit=buds_and_flowers_per_stem/fruits_per_stem) %>%
  inner_join(select(env_data, population, Lat, Long, Elev_m)) 

#correlations @plant 
corr_mat <- round(cor(yield_df[,6:14], method=c("pearson"), use = "complete.obs"),4)
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr_mat){
  corr_mat[lower.tri(corr_mat)]<- NA
  return(corr_mat)
}
upper_tri <- get_upper_tri(corr_mat)
melted_corr_mat <- melt(upper_tri)
quartz()
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

# correlations @pop means
corr_mat <- round(cor(trait_means_df[,2:13], method=c("pearson"), use = "complete.obs"),4)
# Get upper triangle of the correlation matrix
get_upper_tri <- function(corr_mat){
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

# using the ggcorrplot package
library(ggcorrplot)
trait_means_df <- filter(trait_means_df, !population=='APPAR')
corr_mat <- round(cor(trait_means_df[,2:13], method=c("pearson"), use = "complete.obs"),4)
p_mat <- cor_pmat(trait_means_df[,2:13])
head(p_mat)
quartz()
corr_plot <- ggcorrplot(corr_mat, hc.order = TRUE,type = "lower", lab = TRUE, p.mat=p_mat, insig = "blank", lab_size = 3, tl.cex = 7, show.legend = FALSE, title = "Correlations of least-square means")
png("lsmeans_corrplot.png", width=8, height=7, res=300, units="in")
corr_plot
dev.off()


