# trait-environment correlations
# Peter Innes
# 10.5.20
devtools::install_github("jimhester/archive")
devtools::install_github("mirzacengic/climatedata")
install.packages("rgdal")
install.packages("AICcmodavg")
library(dplyr)
library(ggplot2)
library(raster)
library(sp)
library(remotes)
library(climatedata)
library(AICcmodavg)
library(magrittr)
library(lme4)
library(arm) #for se.ranef()

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) 
env_data$source %<>% as.factor
env_data$population %<>% as.factor

env_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>%
  filter(!is.na(Lat) | !is.na(Long)) #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)
  
# scale the predictors
env_data$Lat_s <- scale(env_data$Lat)
env_data$Long_s <- scale(env_data$Long)
env_data$Elev_m_s <- scale(env_data$Elev_m)

# check correlations b/w geographic predictors. no correlations.
plot(env_data$Elev_m_s ~ env_data$Long_s)
plot(env_data$Elev_m_s ~ env_data$Lat_s)
plot(env_data$Lat ~ env_data$Long)
cor.test(env_data$Long, env_data$Elev_m)
cor.test(env_data$Lat, env_data$Elev_m)
cor.test(env_data$Long, env_data$Lat)

#### get climate data from WorldClim
# combine temp and precip to give a drought 'potential evaporation and transpiration' there is a formula for this?
# maybe try 'growing season length' from MODIS. Josh Gray made the most recent?


?raster::getData
# pull data from WorldClim
r <- raster::getData('worldclim', var='bio', res=5)
# get data from the CHELSA database, which has the bioclim data at higher resolution (30 arc sec, ~1km, compared to raster::getData command above, which )
chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

#r <- r[[c(1,9,10,12)]] #select columns (climate variables) that we want to keep
#names(r) <- c("ann_mean_temp","mean_temp_driestQ","mean_temp_warmestQ","annual_prec")

coords <- data.frame(Long=env_data$Long, Lat=env_data$Lat,
                     row.names = env_data$population) %>% na.omit()
points <- SpatialPoints(coords, proj4string = chelsa@crs)

values <- raster::extract(chelsa,points) #previously raster::extract(r,points)

clim_df <- cbind.data.frame(coordinates(points),values) %>%
  tibble::rownames_to_column("population")
env_df <- inner_join(env_data, clim_df) 

#### model selectionn--LMM version 
# 3-way factor models
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

# no predictor
fit_no_pred <- lmer(sd_wt_50_ct ~ (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

# LMM model comparison results: fit_LaxLo is best
LaxLo_ranefs <- coef(fit_LaxLo)$population[1] %>%
  rename("fit_LaxLo"="(Intercept)") %>%
  tibble::rownames_to_column("population")
ranefs <- full_join(ranefs, LaxLo_ranefs)
ranefs <- inner_join(ranefs, dplyr::select(env_data, population, Lat, Long))

#### sw vs temp and precip, fork vs temp and precip

sw_clim_data <- sd_wt_data %>% left_join(dplyr::select(clim_df, population, CHELSA_bio10_01,CHELSA_bio10_09, CHELSA_bio10_12, CHELSA_bio10_17, CHELSA_bio10_10, CHELSA_bio10_18 )) %>%
  rename("MAT"="CHELSA_bio10_01", "AP"="CHELSA_bio10_12", "TDQ"="CHELSA_bio10_09", "PDQ"="CHELSA_bio10_17", "TWQ"="CHELSA_bio10_10", "PWQ"="CHELSA_bio10_18") %>%
  mutate(MAT_s=scale(MAT), AP_s=scale(AP), TDQ_s=scale(TDQ), PDQ_s=scale(PDQ), TWQ_s=scale(TWQ), PWQ_s=scale(PWQ))

fit_sw_clim <- lmer(sd_wt_50_ct ~ MAT_s*AP_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #only MAT is significant
fit_sw_clim2 <- lmer(sd_wt_50_ct ~ TDQ_s*PDQ_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #TDQ and interaction are significant.
fit_sw_clim3 <- lmer(sd_wt_50_ct ~ TWQ_s*PWQ_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #TWQ and interaction are significant
summary(fit_sw_clim) 

forks_clim_data <- stem_data %>% left_join(dplyr::select(clim_df, population, CHELSA_bio10_01,CHELSA_bio10_09, CHELSA_bio10_12, CHELSA_bio10_17, CHELSA_bio10_10, CHELSA_bio10_18 )) %>%
  rename("MAT"="CHELSA_bio10_01", "AP"="CHELSA_bio10_12", "TDQ"="CHELSA_bio10_09", "PDQ"="CHELSA_bio10_17", "TWQ"="CHELSA_bio10_10", "PWQ"="CHELSA_bio10_18") %>%
  mutate(MAT_s=scale(MAT), AP_s=scale(AP), TDQ_s=scale(TDQ), PDQ_s=scale(PDQ), TWQ_s=scale(TWQ), PWQ_s=scale(PWQ))

fit_forks_clim <- lmer(forks ~ TDQ_s*PDQ_s + (1|population) + (1|population:block) + (1|population:block:plant), data = forks_clim_data)
summary(fit_forks_clim)


# extract population-level means from REML fit using coef()
# First need a simple data frame of population coordinates to be used later on in model equation.
clims <- clim_df %>%
  na.omit() %>%
  dplyr::select(population, CHELSA_bio10_09, CHELSA_bio10_17, CHELSA_bio10_10, CHELSA_bio10_18) %>%
  rename("TDQ"="CHELSA_bio10_09", "PDQ"="CHELSA_bio10_17", "TWQ"="CHELSA_bio10_10", "PWQ"="CHELSA_bio10_18") %>%
  mutate(TDQ_s=scale(TDQ), PDQ_s=scale(PDQ), TWQ_s=scale(TWQ), PWQ_s=scale(PWQ)) %>%
  arrange(population)

ml_pred_df <- data.frame(coef(fit_sw_clim3)$population,
                         se.ranef(fit_sw_clim3)$population[,1])
names(ml_pred_df) <- c("pop_b0", "b1", "b2", "b3", "pop_b0_se")
ml_pred_df <- ml_pred_df %>%
  tibble::rownames_to_column("population") %>%
  inner_join(clims)

# Calculate intercepts for each population 
ml_pred_df$pop_b0 <- ml_pred_df$pop_b0 + ml_pred_df$b1*ml_pred_df$TWQ_s + ml_pred_df$b2*ml_pred_df$PWQ_s + ml_pred_df$b3*ml_pred_df$TWQ_s*ml_pred_df$PWQ_s

# Plot
plot_sw_clim3 <- ggplot(data=ml_pred_df) +
  geom_abline(intercept=fixef(fit_sw_clim3)[1], slope=fixef(fit_sw_clim3)[3], col="blue", lty=2) +
  geom_point(mapping=aes(x=TWQ_s, y=pop_b0)) +
  geom_linerange(mapping=aes(x=TWQ_s, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se)) +
  labs(x="Mean Temp Warmest Quarter (scaled)", y="Estimated intercept in population l")

#### aggregated number of stems vs geo
plot(back_trans_mean_stem_num ~ Lat, data=trait_env_df)
fit_ns_geo <- lm(back_trans_mean_stem_num ~ Lat + Long + Elev_m, data=trait_env_df)
plot(fit_ns_geo)

# keep all the data per Brett's recommendation to not aggregate.
ns <- stem_data %>% dplyr::select(population,trt,block,row,plant,num_of_stems) %>% unique()
ns_env_df <- inner_join(ns, env_df)
plot(num_of_stems~Lat, data=ns_env_df)
fit_sn_all <- lmer(num_of_stems ~ Lat + (1|population) + (1|block) + (1|population:block), data=ns_env_df)
summary(fit_sn_all)
plot(fit_sn_all) #does not look good

coef(fit_sn_all)$population #also wonky. Way too high. I think there are a few points that shift these population means to be much greater than the fixed-effects model.

fit_clim <- lm(mean_sd_wt_50_ct ~ bio11, data=trait_env_df)
summary(fit_clim)

# significant slopes with bio2,bio3,bio4,bio11. Still working on re-doing this to calcultate pearson correlation coefficients, then correct for multiple testing.
bio <- 0
for (var in names(trait_env_df[,7:25])) {
  bio <- bio + 1
  test <- cor.test
  fit <- lm(mean_sd_wt_50_ct ~ get(var), data=trait_env_df)
  s <- summary(fit)
  if (s$coefficients[2,4] < 0.05){
    print(bio)
    print(s$coefficients)  
  } 
}
cor.test(trait_env_df$mean_sd_wt_50_ct, trait_env_df[,7:25])


# plots

# define a boxplot panel function, then plot pairwise correlations to check for collinearity: https://www.flutterbys.com.au/stats/tut/tut7.3a.html
panel.bxp <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 2))
  boxplot(x, add = TRUE, horizontal = T)
}
pairs(~mean_sd_wt_50_ct + bio1 + bio12, data = trait_env_df, lower.panel = panel.smooth,
      diag.panel = panel.bxp, upper.panel = NULL, gap = 0)

cor.test(trait_env_df$bio1,trait_env_df$bio12) #mean annual temp and annual precip are strongly correlated

plot(mean_sd_wt_50_ct ~ Elev_m, data=trait_env_df)
quartz()
plot(mean_sd_wt_50_ct ~ Lat, data=trait_env_df)
plot(mean_sd_wt_50_ct ~ bio3, data=trait_env_df)

#### messing around with PCA/RDA as a multivariate approach.
install.packages("vegan")
library(vegan)
pop_trait_means <- read.csv("data/pop_trait_means.csv", header = T)

my_trait_rda <- rda(pop_trait_means[2:9], scale = T) #scale everything bc they are vastly different units. Also don't use EST_YIELD (14th col) since that was a composite index of the other values
biplot(my_trait_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(my_trait_rda, dis="sites", labels=pop_trait_means$population, cex=0.5 )

summary(my_trait_rda)
#ordihull(my_rda, group = pop_trait_means$population) 

my_clim_rda <- rda(clim_df[4:22], scale = T)
biplot(my_clim_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))

# using base R
my_pca <- prcomp(clim_df[4:22], scale = T)

# regress traits vs clim PCA
pop_trait_means <- pop_trait_means %>% filter(!population %in% c('APPAR', '37'))
sw_clim_fit <- lm(pop_trait_means$seed_weight ~ my_clim_rda$CA$u[,1])
sw_clim_fit <- lm(pop_trait_means$number_of_stems ~ my_pca$x[,1])
#prcomp() and rda() give same result. no association b/w seed weight and first PC of the clim data
summary(sw_clim_fit)
plot(pop_trait_means$forks_per_stem ~ my_pca$x[,1])

# try combining with climate predictors?
vegan_df <- inner_join(pop_trait_means, env_df)
vegan_traits <- vegan_df[,2:9]
vegan_geo <- vegan_df[,11:13]
vegan_clim <- vegan_df[,14:32] # geo variables 11:13, climate vars are 14:32
my_rda <- rda(vegan_traits, vegan_geo, scale = T)
my_rda
plot(my_rda, type='n', scaling=1)
orditorp(my_rda, display='sp', cex=0.5, scaling=1, col='blue')
text(my_rda, display='cn', col='red')

#### k means clustering. not working currently 1.13.21
library(cluster)
library(factoextra)
library(purrr)
set.seed(123)
# function to compute total within-cluster sum of square
df <- pop_trait_means
wss <- function(k) {
  kmeans(df, k, nstart = 1 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 42 (numer of pops)
k_values <- 1:length(unique(df$population))

# extract wss for  clusters
wss_values <- map_dbl(k_values, wss)

plot(k_values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# alternatively, this function wraps everything above into one
fviz_nbclust(pop_trait_means, kmeans, method = "wss")

# continue with k=4
clusters <- kmeans(df, 4, nstart = 25)
print(clusters)
fviz_cluster(clusters, data = df) # I don't really like this plot
 


#### plot for bev sears grant
newd <- data.frame(Lat = seq(min(trait_env_df$Lat), max(trait_env_df$Lat), length.out=100))
pred_w_ci <- cbind(newd,predict(fit_geo, newd, interval = "confidence"))
colnames(pred_w_ci)[2] <- c("mean_sd_wt_50_ct") #have to change column name to match the dataframe of actual values to plot CI bands in ggplot.

ggplot(data = trait_env_df, aes(x=Lat,y=mean_sd_wt_50_ct)) +
  geom_point() +
  geom_abline(slope=fit_geo$coefficients[2], intercept=fit_geo$coefficients[1], alpha=0.5) +
  geom_ribbon(data=pred_w_ci, aes(ymin=lwr, ymax=upr), stat = "identity", alpha=0.1) + # confidence interval
  theme_classic() +
  xlab("Latitude") +
  ylab("Mean 50-count seed weight (g)")


