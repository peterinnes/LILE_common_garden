#' ---
#' output: github_document
#' ---

#' environment/climate PCA, misc trait environment models
#' Peter Innes
#' 10.5.20 / last updated 6.23.21

#+ results=FALSE, message=FALSE, warning=FALSE


devtools::install_github("jimhester/archive") #for CHELSA climate data
remotes::install_github("MirzaCengic/climatedata") #for CHELSA climate data

library(dplyr)
#library(magrittr)
library(ggplot2)
#library(interactions) #for interaction plots
library(raster)
library(sp)
library(rgdal)
#library(remotes)
library(climatedata)
#library(AICcmodavg)
library(lme4)
library(lmerTest)
#library(arm) #for se.ranef()
library(vegan)
library(ggrepel)
library(ggvegan)

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
  filter(!is.na(Lat) | !is.na(Long)) %>% #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)
  mutate(Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) # scale predictors

BioClim_codes <- read.csv("BioClim_codes.csv") #this file matches the vague bioclim codes (e.g. bio01) with actual descriptions (e.g. Mean Annual Temp)

#### DOWNLOAD climate data from the CHELSA database ####
# CHELSA has the bioclim data at high resolution (30 arc sec, ~1km()
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

# Make new df with everything scaled.
geo_clim_scaled_df <- geo_clim_df
geo_clim_scaled_df[3:24] <- apply(geo_clim_scaled_df[3:24], 2, function(x){ scale(x) } )

# Check correlations b/w geographic predictors. No significant correlations.
plot(geo_data$Elev_m_s ~ geo_data$Long_s)
plot(geo_data$Elev_m_s ~ geo_data$Lat_s)
plot(geo_data$Lat ~ geo_data$Long)
cor.test(geo_data$Long, geo_data$Elev_m)
cor.test(geo_data$Lat, geo_data$Elev_m)
cor.test(geo_data$Long, geo_data$Lat)

#### PCA of environmental variables (i.e. geography AND climate) ####

my_env_rda <- rda(geo_clim_df[3:24], scale = T) #PCA of scaled geo and clim vars (skip column 1 and 2 which have source/population ID)
summary(my_env_rda)
summary(eigenvals(my_env_rda))[2,1:12] #percent variance explained
# Get site (source/population) PC scores for use in trait-env model selection
env_PC_scores <- data.frame(scores(my_env_rda, choices=1:3, display = "sites", scaling=0)) %>%
  tibble::rownames_to_column("source")

# Get Loadings of the env vars on each PC
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

# Using base R, same as vegan RDA above
my_env_pca <- prcomp(geo_clim_df[3:24], scale = T) 
my_env_pca$rotation
scores(my_env_pca)
biplot(my_env_pca)

# Plot
pops <- fortify(my_env_rda, display='sites')
envs <- fortify(my_env_rda, display='species')
env_pca_plot <- ggplot() +
  geom_point(data=pops, aes(x = PC1, y = PC2), size=2, alpha=0.5) +
  #geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_text_repel(data=pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5) +
  geom_segment(data=envs, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text_repel(data=envs, 
            aes(x=PC1,y=PC2,label=Label), 
            color="red", size=4) +
  labs(x='PC1 (44.5%)', y='PC2 (19.1%)', size=4) +
  theme_minimal()

png("plots/env_pca.png", width=8, height=6, res=300, units="in")
env_pca_plot
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

# Brent wants me to check bio vars 03, 04, 08, 11, 18?
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio03)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio04)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio08)
plot(env_PC_scores$PC2, geo_clim_scaled_df$bio11)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio18)

#### PCA of just climate vars ####
# Thinking we want use these if we also want to include geographic variables as separate predictors in model selection.
my_clim_rda <- rda(geo_clim_df[6:24], scale = T)

summary(my_clim_rda)
summary(eigenvals(my_clim_rda))

# Clim var loadings
clim_PC_loadings <- data.frame(scores(my_clim_rda, choices=1:3, display = "species", scaling = 0)) %>%
  tibble::rownames_to_column("var") %>%
  full_join(dplyr::select(BioClim_codes, var, description)) %>%
  relocate(description, .after = var)
# Criterion for what constitutes a large loading:
clim_PC_loading_cutoff <- sqrt(1/ncol(geo_clim_df[6:24])) #loading of a single variable if each variable contributed equally; sum of squares of all loadings for an individual principal components must sum to 1.

# Population scores
clim_PC_scores <- data.frame(scores(my_clim_rda, choices=1:3, display = "sites")) %>% #scaling=2, i.e. scale by species, is the default, is what the summary() reports. matches the plotting coordinates
  tibble::rownames_to_column("population")
#my_clim_rda$CA$u[,1:2] #alternate way to access loadings of sources i.e. 'sites'. same as setting scaling=0 above

# plot using ggrepel and ggplot
pops <- fortify(my_clim_rda, display='sites')
clims <- fortify(my_clim_rda, display='species')
clim_pca_plot <- ggplot() +
  geom_point(data=pops, aes(x = PC1, y = PC2), size=2, alpha=0.5) +
  #geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_text_repel(data=pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5) +
  geom_segment(data=clims, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text_repel(data=clims, 
                  aes(x=PC1,y=PC2,label=Label), 
                  color="red", size=4) +
  labs(x='PC1 (49%)', y='PC2 (18.6%)', size=4) +
  theme_minimal()
clim_pca_plot

autoplot(my_clim_rda, arrows = TRUE, geom = "text", legend = "none") #alternate plotting option

# Save plot
png("plots/clim_pca.png", width=8, height=6, res=300, units="in")
clim_pca_plot
dev.off()

data.frame(summary(eigenvals(my_clim_rda)))[2,1:12] %>% 
  pivot_longer(1:12, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()


#### Model selection of traits vs env/climate PCs?? ####
traits <- c("Seed_weight", "Capsule_fill", "Stems_per_plant", "Est_Capsules_per_stem", "Indeterminancy_index", "Forks_per_stem", "Stem_diam", "Capsule_diam", "Est_yield", "Est_fecundity")

datasets <- list(sd_wt_data, ff_data, stems, stem_data, stem_data_DI, stem_data, stem_data, stem_data, yield_df, yield_df)

#List of each predictor combination to include in our model selection. Add interactions of geog/clim?
fixefs <- c('Lat + Long + Elev_m + PC1 + PC2',
            'Lat + Long + Elev_m + PC1',
            'Lat + Long + Elev_m + PC2',
            
            'Lat + Long + PC1 + PC2',
            'Lat + Elev_m + PC1 + PC2',
            'Long + Elev_m + PC1 + PC2',
            
            'Lat + Long + Elev_m',
            
            'Lat + Long + PC1',
            'Lat + Elev_m + PC1',
            'Long + Elev_m + PC1',
            
            'Lat + Long + PC2',
            'Lat + Elev_m + PC2',
            'Long + Elev_m + PC2',
            
            'Lat + PC1 + PC2',
            'Long + PC1 + PC2',
            'Elev_m + PC1 + PC2',
            
            'Lat + PC1',
            'Long + PC1',
            'Elev_m + PC1',
            
            'Lat + PC2',
            'Long + PC2',
            'Elev_m + PC2',
            
            'Lat*Long',
            'Lat*Elev_m',
            'Long*Elev_m',
            
            'Lat + Long',
            'Lat + Elev_m',
            'Long + Elev_m',
            
            'PC1*PC2',
            'PC1 + PC2',
            
            'Lat',
            'Long',
            'Elev_m',
            
            'PC1',
            'PC2',
            
            'bio01',
            'bio02',
            'bio03',
            'bio04',
            'bio05',
            'bio06',
            'bio07',
            'bio08',
            'bio09',
            'bio10',
            'bio11',
            'bio12',
            'bio13',
            'bio14',
            'bio15',
            'bio16',
            'bio17',
            'bio18',
            'bio19'
            
            )

results <- list()
for ( i in 1:length(traits) ){ #loop through each trait
  trait <- traits[i]
  data <- datasets[[i]] %>%
    inner_join(clim_PC_scores) %>%
    inner_join(geo_clim_scaled_df)
  
  # Different random effects structure for 'stem_data' traits that have repeated measures from the same plants
  if( traits[i] %in% c("sqr_fruits", "sqr_bds_flow", "sqr_forks", "log_diam_stem", "diam_caps")){
    ranefs <- list('(1|population)', '(1|population:block:plant)') #removed two random effects terms bc of singular fits
  } else {
    ranefs <- list('(1|population)', '(1|block)', '(1|population:block)')
  }
  randstr <- paste(ranefs, collapse=" + ") #Set up the random effects structure for right side of model
  
  fits <- list() #Storage for each model of trait i
  for ( j in 1:length(fixefs) ){ #Loop through the different predictor sets we want to compare
    form <- reformulate(c(fixefs[j], randstr), response=trait) #Set the model formula
    print(form)
    fit <- lmer(form, REML = FALSE, data = data) # Finally, fit the actual model, using maximum likelihood.
    
    fits[[j]] <- fit 
  }
  bic <- bictab(cand.set = fits, modnames = fixefs)
  results[[i]] <- bic
}
names(results) <- traits
results #log_EST_YIELD ~ Elev_m_s doesn't converge. Don't think this is really an issue. It converges w/ REML=T though. weird?

# write bic results tables to files
for ( i in 1:length(results)){
  df <- data.frame(results[i])
  write.csv(format(df, digits=4), file=paste0("results_summaries/BIC_results/", names(results)[i], "_BIC_results.csv"))
}


#### Univariate tests for latitudinal clines? ####
# make_pred_df is a function that takes a LMM as its argument and returns a data frame with estimated group means (intercepts). I use it to find the population trait means of trait~latitude models, which have population as a random effect.
make_pred_df <- function(fit){
  pred_df <- data.frame(coef(fit)$population,
                        se.ranef(fit)$population[,1])
  names(pred_df) <- c("pop_b0", "b1", "pop_b0_se")
  pred_df <- pred_df %>%
    tibble::rownames_to_column("population") %>%
    inner_join(dplyr::select(env_data, population, Lat))
  # Calculate means (intercepts?) for each source 
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

# missing lat models for all other traits

# Make lists and storage for the for() loop
fit_list <- c(fit_sw_Lat)
trait_list <- c("Seed weight")
plot_list = list()
# Loop through each model, calculating population means and confidence intervals of regression lines
for (i in 1:length(fit_list)) {
  fit <- fit_list[[i]]
  pred_df <- make_pred_df(fit) #get population means
  
  # Obtain confidence interval for regression line
  newd <- data.frame(Lat = seq(min(geo_data$Lat, na.rm=T), max(geo_data$Lat, na.rm=T), length.out=100))
  lmm_boots <- bootMer(fit, predict_fun, nsim = 100)
  pred_ci <- cbind(newd, confint(lmm_boots))
  
  # Plot population means (transform scale) vs latitude 
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

# Plot population means (response scale) vs latitude
plot_response_scale <- ggplot(pred_df, aes(Lat, (pop_b0)^2)) + 
  geom_smooth() +
  geom_point()
