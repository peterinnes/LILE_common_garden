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

env_data <- read.csv("LILE_seed_collection_spreadsheet.csv", header=T) 
env_data$source %<>% as.factor
env_data$population %<>% as.factor

env_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!is.na(Lat) | !is.na(Long)) #keep only rows that have coordinates. still need to get missing lat/long data from Stan.

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

#gather everything together with seed weight
trait_env_df <- inner_join(sd_wt_means, env_df) %>%
  filter(!population %in% c('MAPLE_GROVE')) #exclude Maple Grove because seed size may have been influenced by the selection/breeding process, and also UT_1 is the same source population as Maple Grove, also exlcude repeat collections for UT_4 and UT_5

#trait_env_df$ann_mean_temp <- trait_env_df$ann_mean_temp/10 #temp vars are scaled up by factor of 10 in WorldClim dataset
#trait_env_df$mean_temp_warmestQ <- trait_env_df$mean_temp_warmestQ/10
#trait_env_df$mean_temp_driestQ <- trait_env_df$mean_temp_driestQ/10

#### model selection
# AICc fxn from Brett Melbourne
AICc <- function(fitmod) {
  ll <- logLik(fitmod)
  k <-  attr(ll, "df")
  n <- attr(ll,"nobs")
  return( -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / (n - k - 1) )
}
# 3-way factor models
fit_LaxLoxE <- lm(mean_sd_wt_50_ct ~ Lat*Long*Elev_m, data=trait_env_df)

fit_LaLoE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m, data=trait_env_df)

fit_LaLoE_LaxLo_LaxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                             Lat:Long + Lat:Elev_m, data=trait_env_df)
fit_LaLoE_LaxLo_LoxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                             Lat:Long + Long:Elev_m, data=trait_env_df)
fit_LaLoE_LaxE_LoxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                            Lat:Elev_m + Long:Elev_m, data=trait_env_df)
fit_LaLoE_LaxLo_LaxE_LoxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                                  Lat:Long + Lat:Elev_m + Long:Elev_m, data=trait_env_df)
  
fit_LaLoE_LaxLo <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                        Lat:Long, data=trait_env_df)
fit_LaLoE_LaxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                       Lat:Elev_m, data=trait_env_df)
fit_LaLoE_LoxE <- lm(mean_sd_wt_50_ct ~ Lat + Long + Elev_m +
                       Long:Elev_m, data=trait_env_df)

# two factor models
fit_LaxLo <- lm(mean_sd_wt_50_ct ~ Lat*Long, data=trait_env_df)
fit_LaxE <- lm(mean_sd_wt_50_ct ~ Lat*Elev_m, data=trait_env_df)
fit_LoxE <-lm(mean_sd_wt_50_ct ~ Long*Elev_m, data=trait_env_df)
fit_LaLo <- lm(mean_sd_wt_50_ct ~ Lat + Long, data=trait_env_df)
fit_LaE <-lm(mean_sd_wt_50_ct ~ Lat + Elev_m, data=trait_env_df)
fit_LoE <- lm(mean_sd_wt_50_ct ~ Long + Elev_m, data=trait_env_df)

# single factor models
fit_Lat <- lm(mean_sd_wt_50_ct ~ Lat, data=trait_env_df)
fit_Long <- lm(mean_sd_wt_50_ct ~ Long, data=trait_env_df)
fit_Elev <- lm(mean_sd_wt_50_ct ~ Elev_m, data=trait_env_df)

models <- list(fit_LaxLoxE, fit_LaLoE, fit_LaLoE_LaxLo_LaxE, fit_LaLoE_LaxLo_LoxE, fit_LaLoE_LaxE_LoxE, fit_LaLoE_LaxLo_LaxE_LoxE, fit_LaLoE_LaxLo, fit_LaLoE_LaxE, fit_LaLoE_LoxE, fit_LaxLo, fit_LaxE, fit_LoxE, fit_LaLo, fit_LaE, fit_LoE, fit_Lat, fit_Long, fit_Elev, fit_no_pred)

model_names <- c('LaxLoxE', 'LaLoE', 'LaLoE_LaxLo_LaxE', 'LaLoE_LaxLo_LoxE', 'LaLoE_LaxE_LoxE', 'LaLoE_LaxLo_LaxE_LoxE', 'LaLoE_LaxLo', 'LaLoE_LaxE', 'LaLoE_LoxE', 'LaxLo', 'LaxE', 'LoxE', 'LaLo', 'LaE', 'LoE', 'Lat', 'Long', 'Elev', 'No_predictor')

aictab(cand.set = models, modnames = model_names) #LaxLo and Lat are essentially the same. LaxLoAICc is -156.16; Lat AICc is -156.12 when UT_4b and UT_5b excluded. When UT_4a and UT_5a excluded, Lat is -155.40, LaxLo is 155.39.
summary(fit_LaxLo)
summary(fit_Lat)

#### model comparison--LMM version 
# 3-way factor models
  fit_LaxLoxE <- lmer(sd_wt_50_ct ~ Lat_s*Long_s*Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
  
  fit_LaLoE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)
  
fit_LaLoE_LaxLo_LaxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s + 
                               Lat:Long + Lat:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxLo_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                             Lat_s:Long_s + Long_s:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

fit_LaLoE_LaxE_LoxE <- lmer(sd_wt_50_ct ~ Lat_s + Long_s + Elev_m_s +
                            Lat:Elev_m_s + Long:Elev_m_s + (1|population) + (1|block) + (1|population:block), data=sd_wt_data, REML = FALSE)

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

# LMM model comparison results: fit_LaxLo is best. agrees with simpler approach
LaxLo_ranefs <- coef(fit_LaxLo)$population[1] %>%
  rename("fit_LaxLo"="(Intercept)") %>%
  tibble::rownames_to_column("population")
ranefs <- full_join(ranefs, LaxLo_ranefs)
ranefs <- inner_join(ranefs, dplyr::select(env_data, population, Lat, Long))


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

#significant slopes with bio2,bio3,bio4,bio11. Still working on re-doing this to calcultate pearson correlation coefficients, then correct for multiple testing.
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
vegan_df <- na.omit(yield_df)

my_rda <- rda(vegan_df[,6:13], scale = T) #scale everything bc they are vastly different units. Also don't use EST_YIELD (14th col) since that was a composite index of the other values

biplot(my_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordihull(my_rda,
         group = vegan_df$population)

# next try including climate predictors?
vegan_df <- inner_join(vegan_df, env_df)
vegan_traits <- vegan_df[,6:13]
vegan_clim <- vegan_df[,16:18] # geo variables 16:18, climate vars are 19:37
my_rda <- rda(vegan_traits, vegan_clim, scale = T)
my_rda
quartz()
plot(my_rda, type='n', scaling=1)
orditorp(my_rda, display='sp', cex=0.5, scaling=1, col='blue')
text(my_rda, display='cn', col='red')

#### k means clustering
library(cluster)
library(factoextra)
library(purrr)
set.seed(123)
# function to compute total within-cluster sum of square
df <- vegan_df[,6:13]
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 42 (numer of pops)
k_values <- 1:length(unique(vegan_df$population))

# extract wss for  clusters
wss_values <- map_dbl(k_values, wss)

plot(k_values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# alternatively, this function wraps everything above into one
fviz_nbclust(vegan_df[,6:13], kmeans, method = "wss")

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


