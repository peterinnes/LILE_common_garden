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
library(magrittr)
library(ggplot2)
#library(interactions) #for interaction plots
library(rgdal)
library(raster)
library(sp)
library(rgdal)
library(remotes)
library(climatedata)
library(AICcmodavg)
library(lme4)
library(lmerTest)
library(arm) #for se.ranef()
library(vegan)
library(ggrepel)
library(ggvegan)

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
  filter(!is.na(Lat) | !is.na(Long)) %>% #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)
  mutate(Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m)) # scale predictors

BioClim_codes <- read.csv("BioClim_codes.csv") #this file matches the vague bioclim codes with actual descriptions

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

geo_clim_df <- inner_join(geo_data[1:5], clim_df) 
rownames(geo_clim_df) <- geo_clim_df[,1] #set source number to the rownames, otherwise we lose these labels in the PCA below.
head(geo_clim_df)

geo_clim_scaled_df <- geo_clim_df
geo_clim_scaled_df[3:24] <- apply(geo_clim_scaled_df[3:24], 2, function(x){ scale(x) } )

#' #### Checking for association b/w environmental variables ####
# Check correlations b/w geographic predictors. No significant correlations, that's good.
plot(geo_data$Elev_m_s ~ geo_data$Long_s)
plot(geo_data$Elev_m_s ~ geo_data$Lat_s)
plot(geo_data$Lat ~ geo_data$Long)
cor.test(geo_data$Long, geo_data$Elev_m)
cor.test(geo_data$Lat, geo_data$Elev_m)
cor.test(geo_data$Long, geo_data$Lat)

# PCA of environmental variables
my_env_rda <- rda(geo_clim_df[3:24], scale = T) #PCA of scaled geo and clim vars (skip column 1 and 2 which have source/population ID)
summary(my_env_rda)
# Get site (source/population) PC scores for use in trait-env model selection
env_PC_scores <- data.frame(scores(my_env_rda, choices=1:3, display = "sites", scaling=0)) %>%
  tibble::rownames_to_column("source")

# Using base R, same as vegan RDA above
my_env_pca <- prcomp(geo_clim_df[3:24], scale = T) 
scores(my_env_pca)
biplot(my_env_pca)


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

# plot using ggrepel and ggplot
pops <- fortify(my_env_rda, display='sites')
clims <- fortify(my_env_rda, display='species')
env_pca_plot <- ggplot() +
  geom_point(data=pops, aes(x = PC1, y = PC2), shape=1, size=2, alpha=0.75) +
  #geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_text_repel(data=pops, aes(x = PC1, y = PC2, label=Label), size=3) +
  geom_segment(data=clims, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.5, arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=clims, 
            aes(x=PC1,y=PC2,label=Label,
                hjust="inward",vjust=0.75*(1-sign(PC2))), 
            color="red", size=3, alpha=0.5) +
  theme_minimal()

# Plot of proportion variance explained
data.frame(summary(eigenvals(my_env_rda)))[2,1:12] %>%
  pivot_longer(1:12, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  ggplot(aes(x=PC, y=Proportion_Explained)) +
    geom_col() +
    labs(title="PCA of 22 geographic and climate predictors")

# plot the eigenvalues
screeplot(my_env_rda)

# Get loadings (scores) of the env vars on each PC
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

# Brent wants me to check bio vars 03, 04, 08, 11, 18
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio03)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio04)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio08)
plot(env_PC_scores$PC2, geo_clim_scaled_df$bio11)
cor.test(env_PC_scores$PC2, geo_clim_scaled_df$bio18)




# PCA of just climate vars. Thinking we would use these if we also want to include geographic variables as separate predictors in model selection.
my_clim_rda <- rda(geo_clim_df[6:24], scale = T)
clim_PC_scores <- data.frame(scores(my_clim_rda, choices=1:3, display = "sites", scaling=0)) %>% #scaling=2, i.e. scale by species, is the default, is what the summary() reports
  tibble::rownames_to_column("source")
my_clim_rda$CA$u[,1:2] #alternate way to access loadings of sources i.e. 'sites'

summary(my_clim_rda)
summary(eigenvals(my_clim_rda))

biplot(my_clim_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_clim_rda, dis="sites", cex=0.5)
autoplot(my_clim_rda, arrows = TRUE, geom = "text", legend = "none") #alternate plotting option

data.frame(summary(eigenvals(my_clim_rda)))[2,1:12] %>% 
  pivot_longer(1:12, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()

# clim var loadings/scores
clim_PC_loadings <- data.frame(scores(my_clim_rda, choices=1:2, display = "species", scaling = 1)) %>%
  tibble::rownames_to_column("var") %>%
  full_join(dplyr::select(BioClim_codes, var, description)) %>%
  relocate(description, .after = var)

# PCA of just temperature vars (bio1-11)
my_temp_rda <- rda(geo_clim_df[6:16], scale = T)

autoplot(my_temp_rda, rows = TRUE, geom = "text", legend = "none")

data.frame(summary(eigenvals(my_temp_rda)))[2,1:11] %>% 
  pivot_longer(1:11, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()

temp_PC_scores <- data.frame(scores(my_temp_rda, choices=1:3, display = "sites", scaling=0)) %>% #scaling=2, i.e. scale by species, is the default, is what the summary() reports
  tibble::rownames_to_column("source") %>%
  rename(temp1=PC1, temp2=PC2, temp3=PC3)

# PCA of just precip variables (bio12-19)
my_precip_rda <- rda(geo_clim_df[17:24])

autoplot(my_precip_rda, rows = TRUE, geom = "text", legend = "none")

data.frame(summary(eigenvals(my_precip_rda)))[2,1:8] %>% 
  pivot_longer(1:8, names_to = "PC", values_to = "Proportion_Explained") %>%
  mutate(PC=factor(PC, levels = PC)) %>%
  # Plot proportion explained
  ggplot(aes(x=PC, y=Proportion_Explained)) + 
  geom_col()

precip_PC_scores <- data.frame(scores(my_precip_rda, choices=1, display = "sites", scaling=0)) %>%
  tibble::rownames_to_column("source") %>%
  rename(precip=PC1)


#' #### Model selection of traits vs env PCs ####

traits <- c("sd_wt_50_ct", "good_fill", "num_of_stems", "sqr_fruits", "sqr_bds_flow", "sqr_forks", "log_diam_stem", "diam_caps", "log_EST_YIELD")

datasets <- list(sd_wt_data, ff_data, stems, stem_data, stem_data, stem_data, stem_data, stem_data, yield_df)

#List of each predictor combination to include in our model selection. Add interactions of geog/clim?
fixefs <- c('Lat_s + Long_s + Elev_m_s + temp1 + temp2 + precip',
            'Lat_s + Long_s + Elev_m_s + temp1 + precip',
            'Lat_s + Long_s + Elev_m_s + temp2 + precip',
            'Lat_s + Long_s + Elev_m_s + temp1',
            'Lat_s + Long_s + Elev_m_s + temp2',
            'Lat_s + Long_s + Elev_m_s + precip',
            'Lat_s + Long_s + temp1',
            'Lat_s + Elev_m_s + temp1',
            'Long_s + Elev_m_s + temp1',
            
            'Lat_s + Long_s + precip',
            'Lat_s + Elev_m_s + precip',
            'Long_s + Elev_m_s + precip',
            
            'Lat_s + temp1 + precip',
            'Long_s + temp1 + precip',
            'Elev_m_s + temp1 + precip',
            
            'Lat_s + Long_s + temp2',
            'Lat_s + Elev_m_s + temp2',
            'Long_s + Elev_m_s + temp2',
            
            'Lat_s + temp2 + precip',
            'Long_s + temp2 + precip',
            'Elev_m_s + temp2 + precip',
            
            'Lat_s + Long_s + Elev_m_s',
            
            'Lat_s*Long_s',
            'Lat_s*Elev_m_s',
            'Long_s*Elev_m_s',
            
            'Lat_s + Long_s',
            'Lat_s + Elev_m_s',
            'Long_s + Elev_m_s',
            
            'Lat_s',
            'Long_s',
            'Elev_m_s',
            'temp1*precip',
            'temp1 + precip',
            'temp1',
            'precip',
            'temp2*precip',
            'temp2 + precip',
            'temp2'
            
            )

results <- list()
for ( i in 1:length(traits) ){ #loop through each trait
  trait <- traits[i]
  data <- datasets[[i]] %>%
    inner_join(temp_PC_scores) %>%
    inner_join(precip_PC_scores) %>%
    inner_join(dplyr::select(geo_data, source, Lat_s, Long_s, Elev_m_s))
  
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

# Brent suggested includingAND the geographic predictors. 

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
