#' ---
#' output: github_document
#' ---


# LILE study: Trait analyses
# Peter Innes
# 9.29.20

#+ results=FALSE, message=FALSE, warning=FALSE
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)
library(multcomp)
library(multcompView)
library(broom.mixed)
#library(glmmTMB)
#library(DHARMa)
#library(modelsummary)
library(magrittr)
library(reshape2)
library(sjstats)
#library(rstan)
#library(rstanarm)
#library(arm) #for se.ranef()
library(vegan)
library(ggvegan)
library(climatedata)
library(raster)
library(sp)
library(rgdal)

options(contrasts = c("contr.sum","contr.poly"))
#### Read in the data ####
# For trait data, we'll use the cleaned/filtered data frames created and written to files in the traits_EDA.R script 

# ENVIRONMENT/CLIMATE DATA
env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
  filter(!is.na(Lat) | !is.na(Long)) %>% #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)
  mutate(Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) # scale predictors

#' Seed weight data
sd_wt_data <- read.csv("data/sd_wt_data.csv", header = TRUE) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' Fruit fill data
ff_data <- read.csv("data/ff_data.csv", header = TRUE) %>%
  mutate(source=as.factor(source), block=as.factor(block))

#' Stem data
stem_data <- read.csv("data/stem_data.csv", header = T) %>%
  mutate(source=as.factor(source), block=as.factor(block))
stems <- stem_data %>% dplyr::select(source,population,trt,block,row,plot,plant,num_of_stems) %>% unique()

#' Read-in yield data (composite of all the other traits, data frame was created in 01_traits_EDA.R script)
yield_df <- read.csv("data/yield_df.csv", header=T) %>%
  mutate(source=as.factor(source), block=as.factor(block))

# height data
eph_ht_rust <- read.csv("data/StanHtCrownRustJune3_4_13.csv", header=T) %>%
  rename(source=Species_source) %>%
  mutate(source=as.factor(source)) %>%
  full_join(dplyr::select(env_data, source, population)) %>%
  filter(!source %in% c(2,5,22,32,38))

#### DOWNLOAD climate data from the CHELSA database ####
# CHELSA has the bioclim data at high resolution (30 arc sec, ~1km()
BioClim_codes <- read.csv("BioClim_codes.csv") #this file matches the vague bioclim codes (e.g. bio01) with actual descriptions (e.g. Mean Annual Temp)
chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

coords <- data.frame(Long=geo_data$Long, Lat=geo_data$Lat,
                     row.names = geo_data$source) %>% na.omit()

points <- SpatialPoints(coords, proj4string = chelsa@crs)

values <- raster::extract(chelsa,points) #previously raster::extract(r,points)

clim_data <- cbind.data.frame(coordinates(points),values) %>%
  tibble::rownames_to_column("source")
colnames(clim_data)[4:22] <- lapply(colnames(clim_data)[4:22], gsub, pattern = "CHELSA_bio10_", replacement = "bio") #simplify column names

# Combine climate and geography predictors
geo_clim_df <- inner_join(geo_data[1:5], clim_data) 
rownames(geo_clim_df) <- geo_clim_df[,2] #set site to the rownames, otherwise we lose these labels in the PCA below.
head(geo_clim_df)


#### FITTING LINEAR MODELS for each trait ####
# see: https://stackoverflow.com/questions/45788123/general-linear-mixed-effect-glmer-heteroscedasticity-modelling regarding modeling heteroscedasticity 

#' 1. Seed weight
# Fit a fixed-effect model. Population:block interaction accounts for sub-sampling at the 'plot' level (i.e. multiple plants of the same population grown together in the same plot). Sample sizes are mostly consistent (balanced design) across populations and blocks, so shrinkage would be minimal anyways if we fitted population as a random effect.
sd_wt_data$sd_wt_50_ct <- sd_wt_data$sd_wt_50_ct*1000 #convert to mg (per 50 seeds)
fit_sd_wt <- lmer(sd_wt_50_ct ~ population + (1|block) + (1|population:block), data = sd_wt_data)# [-(252:253),])
#summary(fit_sd_wt) #right-skew of residuals?
#fit_sd_wt2 <- lmer(sd_wt_50_ct ~ population + (1|block), data = sd_wt_data) #same model but without plot RE, for simpler estimation of the coefficient of variation.

#fit_sd_wt2 <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data = sd_wt_data)

# Compare pop means fixed vs random. Essentially the same
#sw_mod_comps <- cbind(fixef(fit_sd_wt), coef(fit_sd_wt2)$population) %>%
  #rename("coef(fit_sd_wt2)"="(Intercept)") %>%
  #tibble::rownames_to_column("population")
#sw_mod_comps$population <- gsub("population", "", sw_mod_comps$population)

#' 2. Fruit (capsule) fill. Normal distro probably suffices here.
fit_ff <- lmer(good_fill ~ population + (1|block) + (1|population:block), data = ff_data)
#fit_ff2 <- lmer(good_fill ~ (1|population) + (1|block) + (1|population:block), data = ff_data)
#fit_ff2 <- lmer(good_fill ~ population + (1|block), data = ff_data)

#ff_data$obsv <- 1:nrow(ff_data)
#fit_ff_glm <- glmer(good_sds ~ -1 + population + (1|block) + (1|population:block), family="poisson", data = ff_data) #Trying Poisson distro, but getting singular fit. Diagnostic plots don't look great either. Normal distro probably fine

#' 3. Number of stems per plant. Normal distro.
stems <- stem_data %>% dplyr::select(source,population,trt,block,row,plot,plant,num_of_stems) %>% unique() #Need to filter stem_data to just get rows with number of stems, otherwise this info is duplicated bc there is also fruit, fork count etc for each stem.
fit_num_stems <- lmer(num_of_stems ~ population + (1|block) + (1|population:block), data=stems)
#fit_num_stems2 <- lmer(num_of_stems ~ (1|population) + (1|block) + (1|population:block), data=stems)
#fit_num_stems2 <- lmer(num_of_stems ~ population + (1|block), data=stems)

#' 4. Est. ttl capsules per stem (capsules + buds/flowers). or this trait and subsequent per-stem traits, we need an additional model term to account for multiple measurements taken from same plant
stem_data$EST_ttl_caps <- stem_data$fruits + stem_data$bds_flow
fit_log_EST_ttl_caps <- lmer(log(EST_ttl_caps) ~ population + (1|block) + (1|population:block:plant), data=stem_data) #(1|population:block) explains 0 variance
summary(fit_log_EST_ttl_caps)
qqnorm(resid(fit_log_EST_ttl_caps))
#fit_log_EST_ttl_caps2 <- lmer(log(EST_ttl_caps) ~ population + (1|block), data=stem_data)
#fit_log_EST_ttl_caps2 <- lmer(log(EST_ttl_caps) ~ (1|population) + (1|block) + (1|population:block:plant), data=stem_data)

#' 5. Indeterminancy index (ratio of buds and flowers per plant to capsules per plant). per-plant=per 20 stems. square-root transform. 
stem_data_DI <- stem_data %>%
  group_by(population, block, row, plant) %>% 
  na.omit() %>%
  summarise(ttl_caps=sum(fruits), ttl_bds_flow=sum(bds_flow)) %>% #sum fruits/flowers/buds for each plant
  mutate(DI=(ttl_bds_flow/ttl_caps))

fit_sqr_DI <- lmer(sqrt(DI) ~ population + (1|population:block), data=stem_data_DI)
plot(fit_sqr_DI)
qqnorm(resid(fit_sqr_DI))
#fit_sqr_DI2 <- lmer(sqrt(DI) ~ population + (1|block), data=stem_data_DI)

#fit_sqr_DI2 <- lmer(sqrt(DI) ~ (1|population) + (1|population:block), data=stem_data_DI)

#' 6. Forks per stem. square root transform.
fit_sqr_forks <- lmer(sqrt(forks) ~ population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. ~0 variance due to block, leaving this term out. 
#fit_sqr_forks2 <- lmer(sqrt(forks) ~ (1|population) + (1|block) + (1|population:block:plant), data=stem_data)
#fit_sqr_forks2 <- lmer(sqrt(forks) ~ population + (1|block), data=stem_data)

#fit_forks_nb <- glmmTMB(forks ~ (1|population) + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2") #trying negative binom with glmmTMB. Doesn't look great. will stick with square root transform I think.

#' 7. Stem diameter. log transform
fit_stemd <- lmer(diam_stem ~ population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) #log transform to account for right skew
#fit_stemd2 <- lmer(diam_stem ~ (1|population) + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data)
#fit_stemd2 <- lmer(diam_stem ~ population + (1|block), data=stem_data)


#' 8. Capsule diameter
fit_capsd <- lmer(diam_caps ~ population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term so leave it out.
#fit_capsd2 <- lmer(diam_caps ~ population + (1|block), data=stem_data)

#' 9. Estimated yield. (see 01_traits_EDA.R for how we estimate yield)
fit_log_yield <- lmer(log(EST_YIELD) ~ population + (1|block) + (1|population:block), data=yield_df)
#fit_log_yield2 <- lmer(log(EST_YIELD) ~ population + (1|block), data=yield_df)

#' 10. Estimated fecundity (seeds per plant): Fruit fill x (Fruits per stem + Buds_flowers per stem) x Stems per plant
yield_df$EST_fecundity <- yield_df$num_of_stems *
  (yield_df$fruit_per_stem + yield_df$bds_flws_per_stem) *
  yield_df$good_fill

fit_log_fecundity <- lmer(log(EST_fecundity) ~ population + (1|block) + (1|population:block), data=yield_df)
summary(fit_log_fecundity) # a lot of variation, mostly at the individual plant level (residual)
plot(fit_log_fecundity)
qqnorm(resid(fit_log_fecundity))

fit_log_fecundity2 <- lmer(log(EST_fecundity) ~ population + (1|block), data=yield_df)

#' 11. Height
fit_eph_height <- lmer(Height ~ population + (1|Block) + (1|population:Block), data=eph_ht_rust)
#fit_eph_height2 <- lmer(Height ~ population + (1|Block), data=eph_ht_rust)

#' 12. Rust presence
fit_rust <- lmer(Rust_17th ~ population + (1|Block) + (1|population:Block), data=eph_ht_rust)
#fit_rust2 <- lmer(Rust_17th ~ population + (1|Block), data=eph_ht_rust)

#### Model DIAGNOSTICS ####
eph_fit_list <- c(fit_sd_wt, fit_capsd, fit_ff, fit_num_stems, fit_stemd, fit_log_EST_ttl_caps, fit_sqr_DI, fit_sqr_forks, fit_log_yield, fit_log_fecundity, fit_eph_height, fit_rust)
eph_trait_list <- c("Seed_mass", "Capsule_diam", "Capsule_fill", "Stems_per_plant","Stem_diam", "Est_Capsules_per_stem", "Indeterminancy_index", "Forks_per_stem", "Est_yield", "Est_fecundity", "Height", "Rust_presence")
fvr_list <- list()
qqp_list <- list()
hist_list <- list()
for ( i in 1:length(eph_fit_list) ){
  
  fit <- eph_fit_list[[i]]
  df <- data.frame(fitted=fitted(fit), resid=resid(fit))
  # fitted vs residuals
  fvr <- ggplot(df, aes(x=scale(fitted), y=scale(resid))) +
    geom_point(alpha=0.5, shape=1) +
    labs(title=eph_trait_list[[i]], x="Fitted", y="Resid") 
  fvr_list[[i]] <- fvr
  
  # qqplot
  qqp <- ggplot(broom.mixed::augment(fit), aes(sample=.resid/sd(.resid))) + #scale to variance=1 
    stat_qq(alpha=0.5, shape=1) +
    labs(title=eph_trait_list[[i]])
  qqp_list[[i]] <- qqp
  # histogram of residuals
  
  hist <- ggplot(data=df, aes(x=scale(resid))) +
    geom_histogram() +
    labs(title=eph_trait_list[[i]], x="Resid")
  hist_list[[i]] <- hist
}


fvr_grid <- cowplot::plot_grid(plotlist = fvr_list, ncol = 3)
qqp_grid <- cowplot::plot_grid(plotlist = qqp_list, ncol = 3)
hist_grid <- cowplot::plot_grid(plotlist = hist_list, ncol = 3)

jpeg("plots/eph_mod_diagnostics_resid_v_fitted.jpg", width=17, height=23, res=600, units="cm")
fvr_grid
dev.off()

jpeg("plots/eph_mod_diagnosticss_qqplots.jpg", width=17, height=23, res=600, units="cm")
qqp_grid
dev.off()

jpeg("plots/eph_mod_diagnostics_resid_histograms.jpg", width=17, height=23, res=600, units="cm")
hist_grid
dev.off()


##' misc model diagnostics
##' 1. Seed weight
## sw outliers
#sw_resids <- data.frame(fitted=scale(fitted(fit_sd_wt)), resid=scale(resid(fit_sd_wt)))
#which(sw_resids$resid>=5 | sw_resids$resid <=-5) #obs 252 and 253 are the culprits. these #observations are part of source:block groups with only 3 technical replicates instead of 4. 
#
##' 2. Fruit fill
#plot(fit_ff2) #doesnt look much better than the first
#qqnorm(resid(fit_ff2)) #about the same as first fit
#ff2_resid <- data.frame(resid=resid(fit_ff2))
#ggplot(data=ff2_resid, aes(x=resid, y=stat(density))) +
#  geom_histogram(bins = 100) #resids skewed left
#
##' 3. Num stems
#
##' 6. Forks per stem
##forks_simres <- simulateResiduals(fit_forks_nb)
##plotResiduals(forks_simres)
##plotQQunif(forks_simres)
#
##' 7. Stem diam diagnostics. Stem diam was log-transformed. Looks okay, a few right outliers.
#
##' 8. Caps diam show a slight left skew in resids
#
##' 9. Yield diagnostics
## compare yield estimates from dif models
#yield_mod_comps <- cbind(fixef(fit_yield), coef(fit_yield2)$population, exp(fixef(fit_log_yield)), #exp(coef(fit_log_yield2)$population), exp(coef(fit_yield_glmer)$population)) 
#colnames(yield_mod_comps) <- c("fy","fy2", "fly", "fly2","fyg")
#yield_mod_comps
#
#plot(fit_log_yield)
#plot(fit_yield) #megaphone effect
# #log transform looks better
#plot(fit_log_yield2) 
#plot(fit_yield_glmer) #gamma distro also looks okay
#
#qqnorm(resid(fit_log_yield)) #bout the same as untransformed
##qqnorm(resid(fit_yield)) #actually looks fine
##qqnorm(resid(fit_yield_glmer)) #not as good
#
#yield_resid <- data.frame(resid=resid(fit_yield))
#ggplot(data=yield_resid, aes(x=resid, y=stat(density))) +
#  geom_histogram(bins = 50) #looks pretty good, slight right skew
#
#ly_resid <- data.frame(resid=resid(fit_log_yield))
#ggplot(data=ly_resid, aes(x=resid, y=stat(density))) +
#  geom_histogram(bins = 50) #gets rid of the right-skew.

#### Summarize ls-means of all (non-oil) traits  ####
eph_fit_list <- c(fit_sd_wt, fit_capsd, fit_ff, fit_num_stems, fit_stemd, fit_log_EST_ttl_caps, fit_sqr_DI, fit_sqr_forks, fit_log_yield, fit_log_fecundity, fit_eph_height, fit_rust)
eph_trait_list <- c("Seed_mass", "Capsule_diam", "Capsule_fill", "Stems_per_plant","Stem_diam", "Est_Capsules_per_stem", "Indeterminancy_index", "Forks_per_stem", "Est_yield", "Est_fecundity", "Height", "Rust_presence")
eph_results <- list() #list to store means and confidence intervals
eph_results_bt <- list()
eph_emm_list <- list() #list to store comparison plots
#eph_esp_list <- list() #list to store effect size plots 
emm_options(pbkrtest.limit = 5565) #need to increase mem usage limit.

# This will take a while
for (i in 1:length(eph_fit_list) ){
 
  fit <- eph_fit_list[[i]]
  lsmeans <- as.data.frame(lsmeans(fit, "population")) #as data frame
  emm1 <- lsmeans(fit, "population") #same as lsmeans above but dont convert to df
  emm2 <- as.data.frame(lsmeans(fit, "population", type="response")) #separate object for backtransformed lsmeans
  
  # Plotting means and CIs with emmeans:::plot.emmGrid. in order to reorder the populations on y-axis, we need to edit the .plot.srg function in emmeans package:
  # trace(emmeans:::.plot.srg, edit=T). Edit lines 239-240, change aes_() to aes() and delete tilde from in front of x and y variables. Then use reorder() on the y variable as such: y = reorder(pri.fac, the.emmean).
  emm_plot <- plot(emm1, type = "response", comparisons = T, colors = c("salmon", "blue", "black")) +
    theme_minimal() +
    xlab(eph_trait_list[i]) +
    ylab("")
  eph_emm_list[[i]] <- emm_plot #store plot
  
  # Compact letter display
  contrasts <- emmeans::emmeans(object=fit, type="response", pairwise ~ "population", adjust="tukey")
  cld <- emmeans:::cld.emmGrid(object=contrasts$emmeans, Letters=letters, sort=F)
  cld_df <- data.frame(cld)

  # Renaming columns and sorting results  
  names(lsmeans)[2] <- eph_trait_list[[i]] #Change to actual trait name before storing in results
  eph_results[[i]] <- lsmeans #store means and confidence intervals
  lsmeans <- lsmeans %>% arrange(-lsmeans[2]) #sort by descending trait value to make more readable
  # Same thing but with backtransformed emms and clds. This will be results table 2(?) in manuscript.
  names(cld_df)[2] <- eph_trait_list[[i]] 
  eph_results_bt[[i]] <- cld_df
  #emm2 <- emm2 %>% arrange(-emm2[2])
}
names(eph_results) <- eph_trait_list
names(eph_results_bt) <- eph_trait_list

# store emms in one dataframe with population as rowname. transformed variables are not back-transformed here.
eph_means_df <- data.frame(matrix(ncol = length(eph_trait_list), nrow = length(unique(sd_wt_data$population))))
names(eph_means_df) <- eph_trait_list
rownames(eph_means_df) <- eph_results[[1]]$population
for (i in 1:length(eph_trait_list) ){
  eph_means_df[i]  <- eph_results[[i]][2]
}
#eph_means_df <- eph_means_df %>% arrange(desc(Seed_mass)) %>%
  #tibble::rownames_to_column("Accession") %>%
  #relocate(Accession, .before = Seed_mass)
#write.csv(eph_means_df, file="plots/millville_trait_means_table.csv", row.names = F)

# dataframe for backtransformed ls-means without clds. for PCA/RDA
eph_means_df2 <- data.frame(matrix(ncol = length(eph_trait_list), nrow = length(unique(sd_wt_data$population))))
names(eph_means_df2) <- eph_trait_list
rownames(eph_means_df2) <- eph_results_bt[[1]]$population
for (i in 1:length(eph_trait_list) ){
  eph_means_df2[i] <- eph_results_bt[[i]][2]
}
# join with the oil emms
eph_means_df2 <- cbind(eph_means_df2, oil_means_df)

# Store EMMs with CLDs in dataframe with column for accession/population. This will be for publication, table 4b. 
eph_means_df3 <- data.frame(matrix(ncol = length(eph_trait_list), nrow = length(unique(sd_wt_data$population))))
names(eph_means_df3) <- eph_trait_list
rownames(eph_means_df3) <- eph_results[[1]]$population
eph_means_df3 <- eph_means_df3 %>%
  tibble::rownames_to_column("Accession")
for (i in 1:length(eph_trait_list) ){
  
  emm_sf <- data.frame(apply(eph_results_bt[[i]][c(2:6)], 1:2,
                             function(x) signif(x, 3))) %>% #change sig figs
    mutate(letter_group=eph_results_bt[[i]][7])
  if ( names(eph_results_bt[[i]][2])=="Rust_presence" ){
    emm_sf[1] <- round(emm_sf[1], digits=3)
  }
  eph_means_df3[i+1] <- apply(emm_sf[c(1,6)], 1, paste, collapse="") #combine emmeans and letters into single column
}

eph_means_df3 <- eph_means_df3 %>% arrange(desc(Seed_mass))
#names(eph_means_df3)[2:8] <- c("Capsules per plot", "Capsules per stem", "Stems per plot", "2013 Biomass per plot (g)", "2014 Biomass per plot (g)", "Plant height (cm)", "Plant diameter (cm)")
write.csv(eph_means_df3, "plots/Ephraim_trait_means_table.csv", row.names = F) #Table 4b.
#test <- round(eph_means_df2$Rust_presence, digits = 3)

# Join emm plots together. This figure didn't end up in the man. 
eph_emm_list[[1]] <- eph_emm_list[[1]] + ylab("Population") #set y axis labels for these two plots, which are the left most of the two rows, so that we don't have redundant axes. 
eph_emm_list[[6]] <- eph_emm_list[[6]] + ylab("Population")

# Join all the emm plots together
library(patchwork) #patchwork not working but cowplot is
#layout <- "
#ABCDE
#FGHIJ
#"

#emm_patchwork <- eph_emm_list[[1]] + eph_emm_list[[2]] + eph_emm_list[[3]] + eph_emm_list[[4]] + eph_emm_list[[5]] + eph_emm_list[[6]] + eph_emm_list[[7]] + eph_emm_list[[8]] + eph_emm_list[[9]] + eph_emm_list[[10]] +
  #plot_layout(design = layout)

eph_emm_grid <- cowplot::plot_grid(plotlist = eph_emm_list, ncol = 3)

#jpg("plots/Ephraim_traits_emm_plots.png", width=17, height=23, res=300, units="cm")
#eph_emm_grid
#dev.off()

#### Cet coefficients of variation ####Or just use sjstats?

# Use same models as above
eph_cvs <- c()
# Use the same models as above 
for (fit in eph_fit_list) {
  #resid_sd <- data.frame(VarCorr(fit))$sdcor[2]
  rmse <- performance::rmse(fit)
  grand_mean <- fixef(fit)[1]
  cv <- rmse/grand_mean
  eph_cvs <- append(eph_cvs, cv)
}
names(eph_cvs) <- eph_trait_list



# Gather BLUPs/conditional modes (for PCA/RDA?). incl Oil content, ALA, Linoleic, Oleic, Palmitic, Stearic 
#eph_traits_blups <- c("Seed_mass", "Oil_content", "ALA", "Linoleic", #"Capsule_fill", "Stems_per_plant", "Est_Capsules_per_stem", #"Indeterminancy_index", "Forks_per_stem", "Stem_diam", "Capsule_diam", #"Est_yield", "Est_fecundity")
#
#eph_fit_list_blups <- c(fit_sd_wt2, fit_oil_cont2, fit_ala2, fit_linoleic2, #fit_ff2, fit_num_stems2, fit_log_EST_ttl_caps2, fit_sqr_DI2, fit_sqr_forks2, #fit_stemd2, fit_capsd2, fit_log_yield2, fit_log_fecundity2)
#
#eph_blup_df <- data.frame(matrix(ncol = length(eph_traits), nrow = 37)) #we #have 37 accessions in ephraim garden. 
#names(eph_blup_df) <- eph_traits_blups
#
#for (i in 1:length(eph_fit_list2) ){
#  fit <- eph_fit_list2[[i]]
#  blups <- coef(fit)$population
#  eph_blup_df[i] <- blups
#}
#rownames(eph_blup_df) <- rownames(blups)

#### Trait PCA of accession level means (emmeans) for ALL traits incl oil content and fatty acid compositions ####

#eph_trait_pca <- rda(eph_blup_df, scale = T) #scale everything to unit variance bc different units.
#eph_trait_pca_noAppar <- rda(eph_blup_df[-2,], scale = T)
eph_trait_pca_noAppar2 <- rda(eph_means_df2[-2,], scale = T) #PCA with ALL response-scale emmeans, instead of blups. includes all oil traits. Exclude appar (row 2) and Rust level (column 12)

summary(eph_trait_pca_noAppar2)
biplot(eph_trait_pca_noAppar2)

# Species (traits) loadings.
eph_trait_PC_loadings <- round(data.frame(scores(eph_trait_pca_noAppar2, choices=1:3, display = "species", scaling = 0)), digits=3) %>%
  arrange(desc(abs(PC1)))
eph_trait_PC_loading_cutoff <- sqrt(1/ncol(eph_means_df2)) #loading of a single variable if each variable contributed equally; sum of squares of all loadings for an individual principal components must sum to 1.
eph_PCA_eigenvals <- round(summary(eigenvals(eph_trait_pca_noAppar2))[,1:3], digits = 3)
write.csv(rbind(eph_trait_PC_loadings, eph_PCA_eigenvals), "plots/Ephraim_trait_PC_loadings_eigenvals.csv") #Table 2b


# Eigenvalue plot
data.frame(summary(eigenvals(eph_trait_pca_noAppar2)))[2,1:10] %>% 
  pivot_longer(1:10, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()

# Extract df of PC scores to manually build a plot
library(ggrepel)
eph_fort_pops <- fortify(eph_trait_pca_noAppar2, display='sites', scaling=0)
eph_fort_traits <- fortify(eph_trait_pca_noAppar2, display='species', scaling=0)
ephraim_trait_pca_plot <- ggplot() +
  geom_point(data=eph_fort_pops, aes(x = PC1, y = PC2), size=2, alpha=0.5) +
  geom_text_repel(data=eph_fort_pops, aes(x = PC1, y = PC2, label=Label),size=4, alpha=.5) +
  #geom_text_repel(data=pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=eph_fort_traits, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.02,"npc"))) +
  geom_text_repel(data=eph_fort_traits, 
                  aes(x=PC1,y=PC2,label=Label), 
                  color="red", size=4) +
  theme_bw() +
  labs(x=paste0("PC1 ","(",100*eph_PCA_eigenvals[2,1],"%)"), y=paste0("PC2 ", "(",100*eph_PCA_eigenvals[2,2],"%)"), title="Ephraim") +
  theme(text = element_text(size=14))

ephraim_trait_pca_plot2 <- ggplot() +
  geom_point(data=eph_fort_pops, aes(x = PC2, y = PC3), size=2, alpha=0.5) +
  geom_text_repel(data=eph_fort_pops, aes(x = PC2, y = PC3, label=Label),size=4, alpha=.5) +
  #geom_text_repel(data=pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=eph_fort_traits, aes(x=0, xend=PC2, y=0, yend=PC3), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.02,"npc"))) +
  geom_text_repel(data=eph_fort_traits, 
                  aes(x=PC2,y=PC3,label=Label), 
                  color="red", size=4) +
  theme_bw() +
  labs(x=paste0("PC2 ","(",100*eph_PCA_eigenvals[2,2],"%)"), y=paste0("PC3 ", "(",100*eph_PCA_eigenvals[2,3],"%)"), title="Ephraim") +
  theme(text = element_text(size=14))

# ggvegan plot
autoplot(eph_trait_pca, arrows = TRUE, geom = "text", legend = "none") #basic version

#png("plots/Ephraim_traits_PCA.png", width=9, height=9, res=300, units="in")
#ephraim_trait_pca_plot
#dev.off()

# Join Millville and Ephraim plots together. FIG 2
fig2 <- plot_grid(millville_trait_pca_plot, ephraim_trait_pca_plot, labels=c("a)","b)"), ncol=1, nrow=2)
png("plots/Fig2.jpg", width=17, height = 23, res=600, units = "cm")
fig2
dev.off()

#### Full RDA (all env predictors and all traits, incl oil content, ALA, linoleic) of accession-level data ####
eph_means_df2$population <- rownames(eph_means_df2)
#eph_rda_means <- eph_means_df2 %>%
  #mutate(Saturated_FA = Stearic + Palmitic, Unsaturated_FA = Alphalinolenic + Linoleic + Oleic) %>%
  #dplyr::select(-c(Alphalinolenic, Linoleic, Oleic, Palmitic, Stearic))
eph_fullRDA_df <- inner_join(eph_means_df2, geo_clim_df) 
rownames(eph_fullRDA_df) <- eph_fullRDA_df$population
eph_fullRDA_df <- eph_fullRDA_df %>% dplyr::select(-c(population, source))
eph_RDA_traits <- eph_fullRDA_df[1:18]
eph_RDA_preds <- eph_fullRDA_df[19:40]

eph_full_rda <- rda(eph_RDA_traits ~ ., data=eph_RDA_preds, scale = T)
#eph_full_rda <- rda(eph_RDA_traits[-12] ~ ., data=eph_RDA_preds, scale = T) # exclude Rust_presence trait (column 12)

# trait loadings
eph_rda_trait_loadings <- round(data.frame(scores(eph_full_rda, choices=1:4, display = "species", scaling = 0)), digits=3) %>%
  arrange(desc(abs(RDA1)))
eph_rda_trait_loading_cutoff <- sqrt(1/ncol(eph_RDA_traits)) # 

# get env_loadings
eph_rda_env_loadings <- round(data.frame(scores(eph_full_rda, choices=1:4, display = "bp", scaling = 0)), digits=3) %>% 
  arrange(desc(abs(RDA1)))
write.csv(eph_rda_env_loadings, "plots/Ephraim_rda_env_loadings.csv")

eph_rda_eigenvals <- round(summary(eigenvals(eph_full_rda, model = "constrained"))[,1:4], digits = 3) #constrained by climate
eph_rda_eigenvals_adj <- round(rbind(eph_rda_eigenvals["Eigenvalue",], data.frame(eph_rda_eigenvals[2:3,]) * RsquareAdj(eph_full_rda)[2]), digits = 3) 
rownames(eph_rda_eigenvals_adj)[1] <- "Eigenvalue"

# bind trait loadings and eigenvals for TABLE 3b.
write.csv(rbind(eph_rda_trait_loadings, eph_rda_eigenvals_adj), "plots/Ephraim_rda_loadings_eigenvals.csv")

#eph_rda.sp_sc0 <- scores(eph_full_rda, choices = 1:2, scaling=0, display="sp") #scaling 0
#eph_rda.sp_sc1 <- scores(eph_full_rda, choices = 1:2, scaling=1, display="sp") #scaling 1
eph_rda.sp_sc <- data.frame(scores(eph_full_rda, choices = 1:2, scaling = 0, display="sp")) #scaling 0. (scaling 2 is default)
eph_rda.env_sc <- data.frame(scores(eph_full_rda, choices = 1:2, scaling = 0, display = "bp"))
#eph_mul <- ordiArrowMul(scores(eph_full_rda, choices = 1:2, scaling = 2, display = "bp")) #multiplier for the coordinates of the head of the env vectors such that they fill set proportion of the plot region. This function is used in the default plot() function for rda objects. value is 3.297
#eph_rda.env_sc <- eph_rda.env_sc*eph_mul

# site scores
eph_rda.site_sc <- data.frame(scores(eph_full_rda, choices = 1:2, scaling = 0, display = "wa"))

# RDA PLOTS
eph_rda_triplotgg <- ggplot() +
  geom_point(data=eph_rda.site_sc, aes(x = RDA1, y = RDA2), size=2, alpha=0.5) +
  #geom_text(data=eph_rda.site_sc, aes(x = RDA1, y = RDA2, label=rownames(eph_rda.site_sc)), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=eph_full_rda, aes(x = RDA1, y = RDA2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=eph_rda.sp_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(.02, "npc"))) +
  geom_text_repel(data=eph_rda.sp_sc, 
                  aes(x=RDA1,y=RDA2,label=rownames(eph_rda.sp_sc)), 
                  color="red", size=4) +
  geom_segment(data=eph_rda.env_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2),
               alpha=0.5, color="blue", arrow=arrow(length=unit(.02,"npc"))) +
  #geom_text_repel(data=eph_rda.env_sc, 
                  #aes(x=RDA1,y=RDA2,label=rownames(eph_rda.env_sc)), 
                  #color="blue", size=4) +
  annotate("text", x = -.68, y = -.245, label = "Lat", color='blue') +
  annotate("text", x = .325, y = .45, label = "Temperature", color='blue') +
  annotate("text", x = -.13, y = -.425, label = "Precipitation", color='blue') +
  annotate("text", x = .6, y = -.17, label = "MDR", color='blue') +
  annotate("text", x = -.17, y = .46, label = "bio08", color='blue') +
  annotate("text", x = -.03, y = .46, label = "bio04", color='blue') +
  annotate("text", x = -.18, y = -.25, label = "Long", color='blue') +
  labs(x=paste0("RDA1 ","(",100*eph_rda_eigenvals_adj[2,1],"%)"), y=paste0("RDA2 ", "(",100*eph_rda_eigenvals_adj[2,2],"%)"), title="Ephraim") +
  theme(axis.text=element_text(size=12),axis.title = element_text(size=16)) +
  #geom_hline(yintercept = 0, lty=2, alpha=0.5) +
  #geom_vline(xintercept = 0, lty=2, alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size = 14))

eph_rda_triplotgg

# combine with Milville RDA (code for Milville plot is in separate file) FIG 3
fig3 <- plot_grid(mv_rda_triplotgg, eph_rda_triplotgg, labels=c("a)","b)"), ncol=1, nrow=2)
jpeg("plots/fig3.jpg", width=17, height=23, res=600, units="cm")
fig3
dev.off()

# Plot for SUPP mat with labels for every predictor arrow (e.g. bio01, bio02, etc)
eph_rda_triplotgg_SUPP <- ggplot() +
  geom_point(data=eph_rda.site_sc, aes(x = RDA1, y = RDA2), size=2, alpha=0.5) +
  #geom_text(data=eph_rda.site_sc, aes(x = RDA1, y = RDA2, label=rownames(eph_rda.site_sc)), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=eph_full_rda, aes(x = RDA1, y = RDA2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=eph_rda.sp_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(.02, "npc"))) +
  geom_text_repel(data=eph_rda.sp_sc, 
                  aes(x=RDA1,y=RDA2,label=rownames(eph_rda.sp_sc)), 
                  color="red", size=4) +
  geom_segment(data=eph_rda.env_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2),
               alpha=0.5, color="blue", arrow=arrow(length=unit(.02,"npc"))) +
  geom_text_repel(data=eph_rda.env_sc, 
  aes(x=RDA1,y=RDA2,label=rownames(eph_rda.env_sc)), 
  color="blue", size=4, alpha=0.5) +
  #annotate("text", x = -.6, y = .3, label = "Lat", color='blue') +
  #annotate("text", x = .25, y = -.5, label = "Temperature", color='blue') +
  #annotate("text", x = -.05, y = .55, label = "Precip", color='blue') +
  #annotate("text", x = .6, y = .05, label = "MDR", color='blue') +
  labs(x=paste0("RDA1 ","(",100*eph_rda_eigenvals_adj[2,1],"%)"), y=paste0("RDA2 ", "(",100*eph_rda_eigenvals_adj[2,2],"%)"), title="Ephraim") +
  theme(axis.text=element_text(size=12),axis.title = element_text(size=16)) +
  #geom_hline(yintercept = 0, lty=2, alpha=0.5) +
  #geom_vline(xintercept = 0, lty=2, alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size = 14))

# SUPP FIG 3. 
figS3 <- plot_grid(mv_rda_triplotgg_SUPP, eph_rda_triplotgg_SUPP, labels=c("a)","b)"), ncol=1, nrow=2)
jpeg("plots/figS3.jpg", width=17, height=23, res=600, units="cm")
figS3
dev.off()
# Alt/base RDA plots for comparison. Custom plot above should have same arrow positions as the base plot() below. The autoplot() from ggvegan is slighty different--different scaling I think?
plot(eph_full_rda, display = c("sp", "wa", "bp"))
arrows(0,0,eph_rda.sp_sc[,1], eph_rda.sp_sc[,2], length=0, col="red")
autoplot(eph_full_rda, arrows=FALSE, geom="text", legend= "none", scaling=2)

# Tests of significance, variance partition, and variable selection--something we should do? Not doing variable selection.
set.seed(9)
# Global test of RDA result
anova.cca(eph_full_rda, permutations=1000)
# Marginal effects of the terms (each marginal term analysed in a model with all other variables). 
anova.cca(eph_full_rda, by="term", permutations=1000)
anova.cca(eph_full_rda, by="margin", permutations=1000) #or should it be by="margin"
# Test significance of axes
anova.cca(eph_full_rda, by="axis", permutations=1000)

## Variance partition? Didn't end up using this. 
#eph_RDA_geog <- eph_RDA_preds[1:3]
#eph_RDA_clim <- eph_RDA_preds[4:22]
#RsquareAdj(eph_full_rda)
#varpart(eph_RDA_traits, eph_RDA_geog, eph_RDA_clim, scale = T, permutations = 1000)
#showvarparts(2)

## Variable selection with vegan's ordistep()?
#eph_mod0 <- rda(eph_RDA_traits ~ 1, data=eph_RDA_preds, scale = T)
#eph_mod1 <- rda(eph_RDA_traits ~ ., data=eph_RDA_preds, scale = T)
#step.forward <- ordistep(eph_mod0, scope=formula(eph_mod1), direction="forward", pstep=1000)
#eph_R2step <- ordiR2step(eph_mod0, scope=formula(eph_mod1), direction="forward", pstep=1000)

## Variable selection with forward.sel from packfor
#library(packfor)
#eph_R2a.full <- RsquareAdj(eph_full_rda)$adj.r.squared
#forward.sel(eph_RDA_traits, eph_RDA_preds, adjR2thresh =  eph_R2a.full) #bio02,bio08.

#eph_pars_rda <- rda(eph_RDA_traits ~ Lat + bio08 + bio02 + bio16 + bio12 + bio13, data = eph_RDA_preds, scale = T)
#autoplot(eph_pars_rda, arrows=FALSE, geom="text", legend= "none", scaling=2)


#### Enrionment transfer distance test for local adaptation ####
# *moved to misc_Env_Analyses.R script*


#### Trait pairwise Pearson correlations ####
# performed with population trait means from above, scaled and centered
library(ggcorrplot)
eph_trait_corr_df <- eph_means_df2[-2,] #exclude Appar bc different species

scaled_eph_trait_corr_df <- scale(eph_trait_corr_df, center = T, scale = T)

eph_corr_mat <- round(cor(scaled_eph_trait_corr_df, method=c("pearson"), use = "complete.obs"),4)

eph_p_mat <- cor_pmat(scaled_eph_trait_corr_df)
head(eph_p_mat)
eph_corr_plot <- ggcorrplot(eph_corr_mat, hc.order = TRUE,type = "lower", lab = TRUE, p.mat=eph_p_mat, insig = "blank", lab_size = 4, tl.cex = 10, show.legend = FALSE, title = "Ephraim") + #Using default sig level of .05
  theme(text = element_text(size = 14))

# FIG S4b
jpeg("plots/ephraim_trait_corrplot.jpg", width=17, height=23, res=600, units="cm")
eph_corr_plot
dev.off()

both_corr_plots <- plot_grid(mv_corr_plot, eph_corr_plot, labels=c("a)","b)"), ncol=2, nrow=1)

jpeg("plots/corr_plots.jpg", width=17, height=17, res=600, units="cm")
both_corr_plots
dev.off()

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

#### Fecundity vs. seed mass ####
# FIG S5
S5_df <- dplyr::select(eph_means_df2, Seed_mass, Est_fecundity) %>%
  mutate(species="L. lewisii")
S5_df[2,3] <- "L. perenne"

figS5 <- ggplot(data=S5_df, aes(x=Seed_mass, y=Est_fecundity, color=species, shape=species)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(data=S5_df[-2,], aes(x=Seed_mass, y=Est_fecundity), method=lm) +
  labs(x="Seed mass (mg per 50 seeds)", y="Est. fecundity (seeds per plant)") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.7)) +
  theme(text = element_text(size = 14))

jpeg(file="plots/figS5.jpg",
     width=17, height=12, res=600, units="cm")
figS5
dev.off()

cor.test(S5_df$Est_fecundity[-2], S5_df$Seed_mass[-2])
  

