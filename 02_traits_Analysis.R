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
site_info

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

#' #### FITTING LINEAR MODELS fir each trait
# see: https://stackoverflow.com/questions/45788123/general-linear-mixed-effect-glmer-heteroscedasticity-modelling regarding modeling heteroscedasticity 

#' 1. Seed weight
# Fit a fixed-effect model. Population:block interaction accounts for subsampling at the 'plot' level (i.e. multiple plants of the same population grown together in the same plot). Sample sizes are mostly consistent (balanced design), so shrinkage would be minimal anyways if we fitted population as a random effect here.
fit_sd_wt <- lmer(sd_wt_50_ct ~ -1 + population + (1|block) + (1|population:block), data = sd_wt_data)# [-(252:253),]) 
summary(fit_sd_wt) #Notice right-skew of residuals.

fit_sd_wt2 <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data = sd_wt_data)

# Compare pop means fixed vs random. Essentially the same
sw_mod_comps <- cbind(fixef(fit_sd_wt), coef(fit_sd_wt2)$population) %>%
  rename("coef(fit_sd_wt2)"="(Intercept)") %>%
  tibble::rownames_to_column("population")
sw_mod_comps$population <- gsub("population", "", sw_mod_comps$population)

#' 2. Fruit fill mod. Normal distro could suffice here.
fit_ff <- lmer(good_fill ~ -1 + population + (1|block) + (1|population:block), data = ff_data)
ff_fit_summary <- summary(fit_ff)

ff_data$obsv <- 1:nrow(ff_data)
fit_ff2 <- glmer(good_sds ~ -1 + population + (1|block) + (1|population:block) + (1|obsv), family="poisson", data = ff_data) #singular fit, Normal distro probably fine
ff_ff2_summary <- summary(fit_ff2)

#' 3. Number of stems per plant. Normal distro
stems <- stem_data %>% dplyr::select(population,trt,block,row,plot,plant,num_of_stems) %>% unique() #Need to filter stem_data to just get rows with number of stems, otherwise this info is duplicated bc there is also fruit, fork count etc for each stem.
fit_num_stems <- lmer(num_of_stems ~ -1 + population + (1|block) + (1|population:block), data=stems)
ns_fit_summary <- summary(fit_num_stems)

#' 4. Fruit per stem mod. For this trait and subsequent per-stem traits, we need an additional model term to account for multiple measurements taken from same plant
fit_fruit <- lmer(fruits ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) 
summary(fit_fruit)

# try negative binomial instead
fit_fruit_nb <- glmmTMB(fruits ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2") 
#' 5. Buds/flowers per stem mod
fit_bf <- lmer(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.
summary(fit_bf)

# Trying other distros. neither seem that promising
fit_bf_nb <- glmmTMB(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2")
fit_bf_p <- glmmTMB(bds_flow ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data, family = "poisson")
  
#' 6. Forks per stem
fit_forks <- lmer(forks ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.

fit_forks_nb <- glmmTMB(forks ~ (1|population) + (1|population:block) + (1|population:block:plant), data=stem_data, family = "nbinom2") #trying negative binom with glmmTMB. Doesn't seem much better

#' 7. Stem diameter
fit_log_stemd <- lmer(log(diam_stem) ~ -1 + population + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data) #log transform to account for right skew
summary(fit_log_stemd)

#' 8. Capsule diameter
fit_capsd <- lmer(diam_caps ~ -1 + population + (1|population:block) + (1|population:block:plant), data=stem_data) #singular fit with (1|block) term. Leaving this term out fixes issue.
summary(fit_capsd)

#' 9. Estimated yield. (see 01_traits_EDA.R for how we estimate yield)
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
fit_list <- c(fit_sd_wt, fit_ff, fit_num_stems, fit_fruit, fit_bf, fit_forks, fit_log_stemd, fit_capsd, fit_log_yield)
trait_list <- c("Seed_weight", "Fruit_fill", "Stems", "Fruits", "Buds_flowers", "Forks", "log_Stem_dia", "Capsule_dia", "log_Est_yield")
fvr_list <- list()
qqp_list <- list()
hist_list <- list()
for ( i in 1:length(fit_list) ){
  
  fit <- fit_list[[i]]
  df <- data.frame(fitted=fitted(fit), resid=resid(fit))
  # fitted vs residuals
  fvr <- ggplot(df, aes(x=scale(fitted), y=scale(resid))) +
    geom_point(alpha=0.5, shape=1) +
    labs(title=trait_list[[i]], x="Fitted", y="Resid") 
  fvr_list[[i]] <- fvr
  
  # qqplot
  qqp <- ggplot(broom.mixed::augment(fit), aes(sample=.resid/sd(.resid))) + #scale to variance=1 
    stat_qq(alpha=0.5, shape=1) +
    labs(title=trait_list[[i]])
  qqp_list[[i]] <- qqp
  # histogram of residuals
  
  hist <- ggplot(data=df, aes(x=scale(resid))) +
    geom_histogram() +
    labs(title=trait_list[[i]], x="Resid")
  hist_list[[i]] <- hist
}


fvr_grid <- cowplot::plot_grid(plotlist = fvr_list, ncol = 3)
qqp_grid <- cowplot::plot_grid(plotlist = qqp_list, ncol = 3)
hist_grid <- cowplot::plot_grid(plotlist = hist_list, ncol = 3)

png("plots/rvf_plots.png", width=11, height=9, res=300, units="in")
fvr_grid
dev.off()

png("plots/qq_plots.png", width=11, height=9, res=300, units="in")
qqp_grid
dev.off()

png("plots/hist_resids.png", width=11, height=9, res=300, units="in")
hist_grid
dev.off()



#' misc additional diagnostics and alternative model diagnostics
#' 1. Seed weight
# sw outliers
sw_resids <- data.frame(fitted=scale(fitted(fit_sd_wt)), resid=scale(resid(fit_sd_wt)))
which(sw_resids$resid>=5 | sw_resids$resid <=-5) #obs 252 and 253 are the culprits. these observations are part of source:block groups with only 3 technical replicates instead of 4. 

#' 2. Fruit fill
plot(fit_ff2) #doesnt look much better than the first
qqnorm(resid(fit_ff2)) #about the same as first fit
ff2_resid <- data.frame(resid=resid(fit_ff2))
ggplot(data=ff2_resid, aes(x=resid, y=stat(density))) +
  geom_histogram(bins = 100) #resids skewed left

#' 3. Num stems

#' 4. Fruit per stem
fruit_nb_simres <- simulateResiduals(fit_fruit_nb)
plot(fruit_nb_simres) #not any better, really
plotQQunif(fruit_nb_simres)

#' 5. Buds/flowers per stem
# diagnostics for nbinom and poisson fits
summary(fit_bf_nb)
bf_nb_simres <- simulateResiduals(fit_bf_nb)
# not entirely sure how to interpret the DHARMa diagnostics
plotQQunif(bf_nb_simres) # all the model deviations are significant
plotResiduals(bf_nb_simres)

bf_pois_simres <- simulateResiduals(fit_bf_p)
plotQQunif(bf_pois_simres) # all the model deviations are significant
plotResiduals(bf_pois_simres)


#' 6. Forks per stem
forks_simres <- simulateResiduals(fit_forks_nb)
plotResiduals(forks_simres)
plotQQunif(forks_simres)

#' 7. Stem diam diagnostics. Stem diam was log-transformed. Looks okay, a few right outliers.

#' 8. Caps diam show a slight left skew in resids

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


#' #### Summarize ls-means of all traits
fit_list <- c(fit_sd_wt, fit_ff, fit_num_stems, fit_fruit, fit_bf, fit_forks, fit_log_stemd, fit_capsd, fit_log_yield)
trait_list <- c("Seed_weight", "Fruit_fill", "Stems", "Fruits", "Buds_flowers", "Forks", "log_Stem_dia", "Capsule_dia", "log_Est_yield")
results <- list() #list to store means and confidence intervals
esp_list <- list() #list to store effect size plots
for (i in 1:length(fit_list) ){
  fit <- fit_list[[i]]
  
  means <- as.data.frame(fixef(fit)) %>% #Get means
    tibble::rownames_to_column(c("population"))
  names(means)[2] <- "trait_value"
  
  ci <- as.data.frame(confint(fit)) %>% #Get confidence intervals
    tibble::rownames_to_column(c("population"))
    names(ci)[2:3] <- c("lwr", "upr")

  means_ci <- inner_join(means, ci)
  means_ci$population <- gsub("population", "", means_ci$population)
  means_ci <- dplyr::select(env_data, population, source, site) %>% inner_join(means_ci)
  # Create effect size plot
  esp <- ggplot(data=means_ci, aes(x=trait_value, y=reorder(site, trait_value), xmin=lwr, xmax=upr)) +
    geom_point() +
    geom_errorbar() +
    ylab("Population") +
    xlab(trait_list[[i]]) +
    theme(axis.text.y = element_text(size = 6))
  
  esp_list[[i]] <- esp #Store plot
  
  names(means_ci)[4] <- trait_list[[i]] #Change to actual trait name before storing in results
  results[[i]] <- means_ci #store means and confidence intervals
  
  # Tweak dataframe and write to csv for summary to send to Scott J et al
  means_ci <- means_ci %>% arrange(-means_ci[4]) #sort descending trait value to make more readable
  write.csv(means_ci, file=paste0("results_summaries/", names(means_ci)[4], "_summary.csv"))
}

# Join all the effect size plots together
esp_grid <- cowplot::plot_grid(plotlist = esp_list, ncol = 3)
png("plots/traits_esp.png", width=12, height=9, res=300, units="in")
esp_grid
dev.off()

#' #### Trait pairwise pearson correlations, performed with population trait means from above
# Gather just the trait ls-means together (no conf intervals)
#pop_trait_means <- data.frame(population=results[[1]]$population)
pop_trait_means <- data.frame(site=results[[1]]$site)
for ( i in 1:length(results) ){
  pop_trait_means <- cbind(pop_trait_means, results[[i]][4]) #add trait means to growing df
}
write.csv(pop_trait_means, file="data/pop_trait_means.csv", row.names = FALSE) #save df
# using the ggcorrplot package
library(ggcorrplot)

trait_corr_df <- pop_trait_means %>%
  mutate(ratio_bf_to_fruit=buds_and_flowers_per_stem/fruits_per_stem) %>%
  filter(!population=='APPAR') #exclude Appar bc different species
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

#' #### PCA of population trait means
library(vegan)
library(ggvegan)
# with source as rownames/labels
pop_trait_means <- full_join(pop_trait_means, dplyr::select(env_data, source, population))
temp <- pop_trait_means[,-1]
rownames(temp) <- temp[,10]
temp <- temp[,-10]
pop_trait_means <- temp
# with population id as rownames/labels
temp <- pop_trait_means[,-1]
rownames(temp) <- pop_trait_means[,1]
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
library(ggrepel)
pops_rda <- fortify(my_trait_rda, display='sites')
traits_rda <- fortify(my_trait_rda, display='species')
quartz()
trait_pca <- ggplot() +
  geom_point(data=pops_rda, aes(x = PC1, y = PC2), shape=1, size=2, alpha=0.75) +
  #geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_text_repel(data=pops_rda, aes(x = PC1, y = PC2, label=Label), size=3) +
  geom_segment(data=traits_rda, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.5, arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=traits_rda, 
            aes(x=PC1,y=PC2,label=Label,
                hjust="inward",vjust=0.75*(1-sign(PC2))), 
            color="red", size=3, alpha=0.5) +
  theme_minimal()
png("plots/trait_pca2.png", width=9, height=9, res=300, units="in")
trait_pca
dev.off()


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


