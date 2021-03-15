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
library(lmerTest)
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

# PCA of environmental variables
rownames(geo_clim_df) <- geo_clim_df[,1]
my_env_rda <- rda(geo_clim_df[3:24], scale = T) #geo and clim vars
summary(my_env_rda)
my_env_pca <- prcomp(geo_clim_df[3:24], scale = T) #using base R, same as vegan RDA above
summary(my_env_pca)

biplot(my_env_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_env_rda, dis="sites", cex=0.5)

png("plots/env_pca.png", width=9, height=9, res=300, units="in")
my_env_rda
dev.off()

autoplot(my_env_rda, arrows = TRUE, geom = "text", legend = "none") #alternate plotting option

# Plot of proportion variance explained
data.frame(summary(eigenvals(my_env_rda)))[2,1:12] %>%
  pivot_longer(1:12, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  ggplot(aes(x=PC, y=Proportion_Explained)) +
    geom_col() +
    labs(title="PCA of 22 geographic and climate predictors")

# plot the eigenvalues
screeplot(my_env_rda)

# Get loadings of the env vars on each PC
BioClim_codes <- read.csv("BioClim_codes.csv")
env_PC1_loadings <- data.frame(scores(my_env_rda, choices=1, display = "species", scaling = 0)) %>%
  tibble::rownames_to_column("var") %>%
  arrange(desc(abs(PC1))) %>%
  full_join(dplyr::select(BioClim_codes, var, description)) %>%
  relocate(description, .after = var)
write.csv(env_PC1_loadings, file = "results_summaries/env_PC1_loadings.csv")

env_PC2_loadings <- data.frame(scores(my_env_rda, choices=2, display = "species", scaling = 0)) %>%
  tibble::rownames_to_column("var") %>%
  arrange(desc(abs(PC2))) %>%
  full_join(dplyr::select(BioClim_codes, var, description)) %>%
  relocate(description, .after = var)
write.csv(env_PC2_loadings, file = "results_summaries/env_PC2_loadings.csv")

env_PC3_loadings <- data.frame(scores(my_env_rda, choices=3, display = "species", scaling=0)) %>% #default is scaling=2 (scale by species) is what the summary() reports
  tibble::rownames_to_column("var") %>%
  arrange(desc(abs(PC3))) %>%
  full_join(dplyr::select(BioClim_codes, var, description)) %>%
  relocate(description, .after = var)
write.csv(env_PC3_loadings, file = "results_summaries/env_PC3_loadings.csv")

# Get site (source/population) scores to use in trait-env model selection
env_PC_scores <- data.frame(scores(my_env_rda, choices=1:3, display = "sites", scaling=0)) %>%
  tibble::rownames_to_column("source")
                              

# PCA of just climate vars
my_clim_rda <- rda(geo_clim_df[6:24], scale = T)
summary(my_clim_rda)
biplot(my_clim_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_clim_rda, dis="sites", cex=0.5)
summary(eigenvals(my_clim_rda))

data.frame(summary(eigenvals(my_clim_rda)))[2,1:12] %>% 
  pivot_longer(1:12, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()

summary(my_clim_rda)

my_clim_rda$CA$v #access PC scores (loadings) 
scores(my_clim_rda, choices=1:2, display = "species", scaling = 1) #alternate way to access scores

my_clim_rda$CA$u[,1:2] #loadings of sources i.e. 'sites'
clim_PCA_scores <- data.frame(scores(my_clim_rda, choices=1:2, display = "sites", scaling=2)) %>% #scaling=2, i.e. scale by species, is the default, is what the summary() reports
  tibble::rownames_to_column("source")





#' #### Model selection of traits vs env PCs ####

traits <- c("sd_wt_50_ct", "good_fill", "num_of_stems", "fruits", "bds_flow", "forks", "log_diam_stem", "diam_caps", "log_EST_YIELD")

datasets <- list(sd_wt_data, ff_data, stems, stem_data, stem_data, stem_data, stem_data, stem_data, yield_df)

fixefs <- c('PC1*PC2*PC3', 'PC1 + PC2 + PC3', 'PC1*PC2', 'PC1*PC3', 'PC2*PC3', 'PC1 + PC2', 'PC1 + PC3', 'PC2 + PC3', 'PC1', 'PC2', 'PC3') #List of each predictor combination to include in our model selection

results <- list()
for ( i in 1:length(traits) ){ #loop through each trait
  trait <- traits[i]
  data <- datasets[[i]] %>% inner_join(env_PC_scores)
  
  # Different random effects structure for 'stem_data' traits that have repeated measures from the same plants
  if( traits[i] %in% c("fruits", "bds_flow", "forks", "log_diam_stem", "diam_caps")){
    ranefs <- list('(1|population)', '(1|population:block:plant)') #removed two random effects terms bc of singular fits
  } else {
    ranefs <- list('(1|population)', '(1|block)', '(1|population:block)')
  }
  randstr <- paste(ranefs, collapse=" + ") #Set up the random effects structure for right side of model
  
  fits <- list() #Storage for each model fit of trait i
  for ( j in 1:length(fixefs) ){ #Loop through the different predictor sets we want to compare
    form <- reformulate(c(fixefs[j], randstr), response=trait) #Set the model formula
    print(form)
    fit <- lmer(form, REML = FALSE, data = data) # Finally, fit the actual model!!
    
    fits[[j]] <- fit 
  }
  bic <- bictab(cand.set = fits, modnames = fixefs)
  results[[i]] <- bic
}
names(results) <- traits
results

# Brent suggested including the env PCs (climate + geog) AND the geographic predictors. 

#' #### Exploring latitudinal clines ####

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
sd_wt_data <- sd_wt_data %>% inner_join(dplyr::select(env_data,source,population, Lat, Lat_s, Long, Long_s, Elev_m, Elev_m_s)) %>%
  inner_join(clim_PCA_scores)

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

fit_hf_Lat <- lmer(HYP_fecundity ~ Lat + (1|population) + (1|block), data=yield_df) #'hypothetical fecundity' = (num fruits + num buds/flowers) * fruit fill
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
