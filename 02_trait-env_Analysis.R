#' ---
#' output: github_document
#' ---

#' trait-environment correlations
#' Peter Innes
#' 10.5.20 / last updated 1.15.21

#+ results=FALSE, message=FALSE, warning=FALSE

#devtools::install_github("jimhester/archive") #for CHELSA climate data
#remotes::install_github("MirzaCengic/climatedata") #for CHELSA climate data
library(dplyr)
library(ggplot2)
#library(interactions) #for interaction plots
library(rgdal)
library(raster)
library(sp)
library(rgdal)
library(remotes)
library(climatedata)
library(AICcmodavg)
library(magrittr)
library(lme4)
library(arm) #for se.ranef()
library(vegan)
library(ggvegan)

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population), Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
  filter(!is.na(Lat) | !is.na(Long)) #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)

#### Get climate data from the CHELSA database ####
# CHELSA has the bioclim data at high resolution (30 arc sec, ~1km()
chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

coords <- data.frame(Long=geo_data$Long, Lat=geo_data$Lat,
                     row.names = geo_data$population) %>% na.omit()
points <- SpatialPoints(coords, proj4string = chelsa@crs)

values <- raster::extract(chelsa,points) #previously raster::extract(r,points)

clim_df <- cbind.data.frame(coordinates(points),values) %>%
  tibble::rownames_to_column("population")
colnames(clim_df)[4:22] <- lapply(colnames(clim_df)[4:22], gsub, pattern = "CHELSA_bio10_", replacement = "bio") #simplify column names

geo_clim_df <- inner_join(geo_data, clim_df) 

#' #### Checking for association b/w environmental variables ####
# Check correlations b/w geographic predictors. No significant correlations, that's good.
plot(geo_data$Elev_m_s ~ geo_data$Long_s)
plot(geo_data$Elev_m_s ~ geo_data$Lat_s)
plot(geo_data$Lat ~ geo_data$Long)
cor.test(geo_data$Long, geo_data$Elev_m)
cor.test(geo_data$Lat, geo_data$Elev_m)
cor.test(geo_data$Long, geo_data$Lat)

# PCA of environmental variables: which ones are redundant?
# First plot everything, then we will check for colinearity.
rownames(geo_clim_df) <- geo_clim_df[,1]
my_env_rda <- rda(geo_clim_df[3:24], scale = T)
biplot(my_env_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_env_rda, dis="sites", cex=0.5)
write.csv(summary(eigenvals(my_env_rda)), file = "results_summaries/env_pca_full.eigenvals.csv")

autoplot(my_env_rda, arrows = TRUE, geom = "text", legend = "none")
png("plots/env_pca.png", width=9, height=9, res=300, units="in")
my_env_rda
dev.off()

# using base R
my_env_pca <- prcomp(geo_clim_df[3:24], scale = T) #same as vegan RDA above
summary(my_env_pca)

# pca of just climate vars
my_clim_rda <- rda(geo_clim_df[6:24], scale = T)
biplot(my_clim_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
summary(eigenvals(my_clim_rda))


# check colinearity of climate and geo predictors. Which predictors should we remove?
geo_clim_corr <- round(cor(geo_clim_df[,3:24], method=c("pearson"), use = "complete.obs"),4)
redun_preds <- caret::findCorrelation(geo_clim_corr, names = T)

clim_corr <- round(cor(geo_clim_df[,6:24], method=c("pearson"), use = "complete.obs"),4) #just climate variables, no geo vars
redun_clim_preds <- caret::findCorrelation(clim_corr, names = T)

geo_clim_df_sub <- dplyr::select(geo_clim_df, -all_of(redun_preds))
my_env_sub_rda <- rda(geo_clim_df_sub[3:16], scale = T)
biplot(my_env_sub_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_env_sub_rda, dis="sites", cex=0.5)
summary(my_env_sub_rda)
summary(eigenvals(my_env_sub_rda))

clim_df_sub <- dplyr::select(geo_clim_df[,6:24], -all_of(redun_clim_preds))
my_clim_sub_rda <- rda(clim_df_sub, scale = T)
biplot(my_clim_sub_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))

summary(eigenvals(my_clim_sub_rda))

#' #### Hypothesis testing for latitudinal clines ####
# make_pred_df is a function that takes a LMM as its argument and returns a data frame with estimated group means (intercepts). I use it to find the population trait means of trait~latitude models, which have population as a random effect.
make_pred_df <- function(fit){
  pred_df <- data.frame(coef(fit)$population,
                               se.ranef(fit)$population[,1])
  names(pred_df) <- c("pop_b0", "b1", "pop_b0_se")
  pred_df <- pred_df %>%
    tibble::rownames_to_column("population") %>%
    inner_join(dplyr::select(env_data, population, Lat))
  # Calculate means (intercepts?) for each population 
  pred_df$pop_b0 <- pred_df$pop_b0 + pred_df$b1*pred_df$Lat
  return(pred_df)
}

# Function for making predictions from LMMs, to obtain confidence intervals subsequently with bootMer() and confint()
predict_fun <- function(fit) {
  predict(fit, newd, re.form=NA)   # This is predict.merMod 
}

# Fit models for each trait vs Latitude. Would prefer to do this with a for() loop but different models come from different data frames, or have slightly different parameterizations due to singularity issues
sd_wt_data <- sd_wt_data %>% inner_join(dplyr::select(env_data,source,population, Lat, Lat_s, Long, Long_s, Elev_m, Elev_m_s)) 

fit_sw_Lat <- lmer(sd_wt_50_ct ~ Lat + (1|population) + (1|population:block),
                   data=sd_wt_data) #leave out (1|block), variance is essentially zero 

ff_data <- ff_data %>% inner_join(dplyr::select(env_data,source,population, Lat))
fit_ff_Lat <- lmer(good_fill ~ Lat + (1|population) + (1|block) + (1|population:block), data = ff_data)

stems <- stems %>% inner_join(dplyr::select(env_data, population, Lat))
fit_ns_Lat <- lmer(num_of_stems ~ Lat + (1|population) + (1|block) + (1|population:block), data=stems)

stem_data <- stem_data %>% inner_join(dplyr::select(env_data, population, Lat))
fit_fruits_Lat <- lmer(fruits ~ Lat + (1|population) + (1|population:block:plant), data=stem_data)
fit_bf_Lat <- lmer(bds_flow ~ Lat + (1|population) + (1|population:block:plant), data=stem_data) #leave out (1|population:block), (1|block), variance is essentially zero for these effects
fit_forks_Lat <- lmer(forks ~ Lat + (1|population) + (1|population:block) + (1|population:block:plant), data=stem_data)
fit_stemd_Lat <- fit_stem_diam <- lmer(log(diam_stem) ~ Lat + (1|population) + (1|block) + (1|population:block) + (1|population:block:plant), data=stem_data)
fit_capsd_Lat <- lmer(diam_caps ~ Lat + (1|population) + (1|population:block) + (1|population:block:plant), data=stem_data)

yield_df <- dplyr::select(env_data, population, site, source, Lat) %>% inner_join(yield_df)
fit_yield_Lat <- lmer(log(EST_YIELD) ~ Lat + (1|population) + (1|block) + (1|population:block), data = yield_df)

fit_num_seeds_Lat <- lmer(EST_seeds_per_plant ~ Lat + (1|population) + (1|block) + (1|population:block), data=yield_df)

fit_hf_Lat <- lmer(HYP_fecundity ~ Lat + (1|population) + (1|block), data=yield_df)
summary(fit_hf_Lat)

# Make lists and storage for the for() loop
fit_list <- c(fit_sw_Lat, fit_ff_Lat, fit_ns_Lat, fit_fruits_Lat, fit_bf_Lat, fit_forks_Lat, fit_stemd_Lat, fit_capsd_Lat, fit_yield_Lat, fit_num_seeds_Lat)
trait_list <- c("Seed weight", "Fruit_fill", "Stems", "Fruits", "Buds_flowers", "Forks", "log_Stem_dia", "Capsule_dia", "log_Est_yield", "Est_num_seeds")
plot_list = list()
# Loop through each model, calculating population means and confidence intervals of regression lines
for (i in 1:length(fit_list)) {
  fit <- fit_list[[i]]
  pred_df <- make_pred_df(fit) #get population means
  
  # Obtain confidence interval for regression line
  newd <- data.frame(Lat = seq(min(geo_data$Lat, na.rm=T), max(geo_data$Lat, na.rm=T), length.out=100))
  lmm_boots <- bootMer(fit, predict_fun, nsim = 100)
  pred_ci <- cbind(newd, confint(lmm_boots))
  
  # Plot population means vs latitude
  plot <- ggplot(data=pred_df) +
    geom_abline(intercept=fixef(fit)[1], slope=fixef(fit)[2], lty=2) +
    geom_point(mapping=aes(x=Lat, y=pop_b0), color="royalblue2", alpha=0.5) +
    geom_linerange(mapping=aes(x=Lat, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se), color="royalblue2", alpha=0.5) +
    geom_ribbon(data=pred_ci, aes(x=Lat, ymin=`2.5 %`, ymax=`97.5 %`), alpha=0.25) +
    labs(x="Lat", y=paste(trait_list[[i]])) +
    theme_minimal()
  
  plot_list[[i]] <- plot
}
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggdraw(add_sub(p, "Latitude", vpadding=grid::unit(0,"lines"),y=6, x=0.53, vjust=4.5))

png("plots/traits_vs_Lat.png", width=11, height=9, res=300, units="in")
p
dev.off()



#' #### Model selection of traits vs climate predictors. bio1, bio12, bio10, bio18, bio3, bio15, bio4 ####

# Merge trait and climate data
sw_clim_data <- inner_join(sd_wt_data, clim_df)
sw_clim_data[,10:28] <- scale(sw_clim_data[,10:28]) #scale predictors

fit_swc1 <- lmer(sd_wt_50_ct ~ bio01*bio12+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc2 <- lmer(sd_wt_50_ct ~ bio01+bio12+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc3 <- lmer(sd_wt_50_ct ~ bio01+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc4 <- lmer(sd_wt_50_ct ~ bio01+bio12 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc5 <- lmer(sd_wt_50_ct ~ bio12+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc6 <- lmer(sd_wt_50_ct ~ bio01 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc7 <- lmer(sd_wt_50_ct ~ bio12 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc8 <- lmer(sd_wt_50_ct ~ bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)

fit_swc9 <- lmer(sd_wt_50_ct ~ bio10*bio18+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc10 <- lmer(sd_wt_50_ct ~ bio10+bio18+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc11 <- lmer(sd_wt_50_ct ~ bio10+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc12 <- lmer(sd_wt_50_ct ~ bio10+bio18 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc13 <- lmer(sd_wt_50_ct ~ bio18+bio03 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc14 <- lmer(sd_wt_50_ct ~ bio10 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)
fit_swc15 <- lmer(sd_wt_50_ct ~ bio18 + (1|population) + (1|block) + (1|population:block), data=sw_clim_data, REML = FALSE)

models <- list(fit_swc1, fit_swc2, fit_swc3, fit_swc4, fit_swc5, fit_swc6, fit_swc7, fit_swc8, fit_swc9, fit_swc10, fit_swc11, fit_swc12, fit_swc13, fit_swc14, fit_swc15)
model_names <- c("fit_swc1", "fit_swc2", "fit_swc3", "fit_swc4", "fit_swc5", "fit_swc6", "fit_swc7", "fit_swc8", "fit_swc9","fit_swc10", "fit_swc11", "fit_swc12", "fit_swc13", "fit_swc14", "fit_swc15")
aictab(cand.set = models, modnames = model_names)

# mod selection using BIC for geo variables
fit_LaxLoxE <- lmer(sd_wt_50_ct ~ Lat_s*Long_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LaxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + 
                               Lat_s:Long + Lat_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                               Lat_s:Long_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                              Lat_s:Elev_m_s + Long:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LaxE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                                    Lat_s:Long_s + Lat_s:Elev_m_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                          Lat_s:Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LaLoE_LaxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                         Lat_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LaLoE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                         Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# two factor models
fit_LaxLo <- lmer(sd_wt_50_ct ~ Lat_s*Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LaxE <- lmer(sd_wt_50_ct ~ Lat_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LoxE <-lmer(sd_wt_50_ct ~ Long_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LaLo <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LaE <-lmer(sd_wt_50_ct ~ Lat_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_LoE <- lmer(sd_wt_50_ct ~ Long_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# single factor models
fit_Lat <- lmer(sd_wt_50_ct ~ Lat_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_Long <- lmer(sd_wt_50_ct ~ Long_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
fit_Elev <- lmer(sd_wt_50_ct ~ Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

models <- list(fit_LaxLoxE, fit_LaLoE, fit_LaLoE_LaxLo_LaxE, fit_LaLoE_LaxLo_LoxE, fit_LaLoE_LaxE_LoxE, fit_LaLoE_LaxLo_LaxE_LoxE, fit_LaLoE_LaxLo, fit_LaLoE_LaxE, fit_LaLoE_LoxE, fit_LaxLo, fit_LaxE, fit_LoxE, fit_LaLo, fit_LaE, fit_LoE, fit_Lat, fit_Long, fit_Elev)

model_names <- c('LaxLoxE', 'LaLoE', 'LaLoE_LaxLo_LaxE', 'LaLoE_LaxLo_LoxE', 'LaLoE_LaxE_LoxE', 'LaLoE_LaxLo_LaxE_LoxE', 'LaLoE_LaxLo', 'LaLoE_LaxE', 'LaLoE_LoxE', 'LaxLo', 'LaxE', 'LoxE', 'LaLo', 'LaE', 'LoE', 'Lat', 'Long', 'Elev')

bic_sw <- bictab(cand.set = models, modnames = model_names) #Lat is best model
step(fit_LaxLoxE, k=log(nobs(fit_LaxLoxE)))


#### PCA/RDA for a multivariate approach? helpful: http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
install.packages("vegan")
devtools::install_github("gavinsimpson/ggvegan")
library(vegan)
library(ggvegan)

my_clim_rda <- rda(clim_df[4:22], scale = T)
biplot(my_clim_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
summary(my_clim_rda) #PC1 explains 49%, PC2 explains 18.6%, PC3 explains 11.1%

# using base R
my_clim_pca <- prcomp(clim_df[4:22], scale = T) #same as vegan RDA above
summary(my_pca)

# regress traits vs clim PCA
ptm <- pop_trait_means %>% filter(!population %in% c('APPAR', '37'))
sw_clim_fit <- lm(ptm$seed_weight ~ my_clim_rda$CA$u[,1])
sw_clim_fit <- lm(ptm$seed_weight ~ my_pca$x[,1])
#prcomp() and rda() give same result. no association b/w seed weight and first PC of the clim data
summary(sw_clim_fit)
plot(ptm$seed_weight ~ my_clim_rda$CA$u[,1])
plot(ptm$seed_weight ~ my_pca$x[,1]) #only dif compared to RDA is PCA1 isnt scaled?

# try combining with climate predictors?
vegan_df <- inner_join(pop_trait_means, geo_clim_df)
vegan_traits <- vegan_df[,2:9]
vegan_geo <- vegan_df[,11:13]
vegan_clim <- vegan_df[,14:32] # geo variables 11:13, climate vars are 14:32
my_rda <- rda(vegan_traits, vegan_geo, scale = T)
my_rda
plot(my_rda, type='n', scaling=1)
orditorp(my_rda, display='sp', cex=0.5, scaling=1, col='blue')
text(my_rda, display='cn', col='red')
