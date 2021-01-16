#' ---
#' output: github_document
#' ---

#' trait-environment correlations
#' Peter Innes
#' 10.5.20 / last updated 1.15.21

#+ results=FALSE, message=FALSE, warning=FALSE
devtools::install_github("jimhester/archive")
devtools::install_github("mirzacengic/climatedata")
install.packages("rgdal")
install.packages("AICcmodavg")
library(dplyr)
library(ggplot2)
library(interactions)
library(raster)
library(sp)
library(remotes)
library(climatedata)
library(AICcmodavg)
library(magrittr)
library(lme4)
library(arm) #for se.ranef()

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population), Lat_s=scale(Lat), Long_s=scale(Long), Elev_m_s=scale(Elev_m))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Lat_s,Long,Long_s,Elev_m,Elev_m_s) %>%
  filter(!source %in% c(2,5,22,32,38)) %>%
  filter(!is.na(Lat) | !is.na(Long)) #keep only pops that have coordinates (missing coords for source 37, and Appar doesn't have coords)

# check correlations b/w geographic predictors. No significant correlations, that's good.
plot(geo_data$Elev_m_s ~ geo_data$Long_s)
plot(geo_data$Elev_m_s ~ geo_data$Lat_s)
plot(geo_data$Lat ~ geo_data$Long)
cor.test(geo_data$Long, geo_data$Elev_m)
cor.test(geo_data$Lat, geo_data$Elev_m)
cor.test(geo_data$Long, geo_data$Lat)

#### get climate data from the CHELSA database, which has the bioclim data at high resolution (30 arc sec, ~1km) 
chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

coords <- data.frame(Long=geo_data$Long, Lat=geo_data$Lat,
                     row.names = geo_data$population) %>% na.omit()
points <- SpatialPoints(coords, proj4string = chelsa@crs)

values <- raster::extract(chelsa,points) #previously raster::extract(r,points)

clim_df <- cbind.data.frame(coordinates(points),values) %>%
  tibble::rownames_to_column("population")
geo_clim_df <- inner_join(geo_data, clim_df) 

# define a boxplot panel function, then plot pairwise correlations to check for collinearity: https://www.flutterbys.com.au/stats/tut/tut7.3a.html
panel.bxp <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 2))
  boxplot(x, add = TRUE, horizontal = T)
}
#pairs(~mean_sd_wt_50_ct + , data = geo_clim_df, lower.panel = panel.smooth,
      #diag.panel = panel.bxp, upper.panel = NULL, gap = 0)

#### model selection with geographic predictors. 
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

models <- list(fit_LaxLoxE, fit_LaLoE, fit_LaLoE_LaxLo_LaxE, fit_LaLoE_LaxLo_LoxE, fit_LaLoE_LaxE_LoxE, fit_LaLoE_LaxLo_LaxE_LoxE, fit_LaLoE_LaxLo, fit_LaLoE_LaxE, fit_LaLoE_LoxE, fit_LaxLo, fit_LaxE, fit_LoxE, fit_LaLo, fit_LaE, fit_LoE, fit_Lat, fit_Long, fit_Elev)

model_names <- c('LaxLoxE', 'LaLoE', 'LaLoE_LaxLo_LaxE', 'LaLoE_LaxLo_LoxE', 'LaLoE_LaxE_LoxE', 'LaLoE_LaxLo_LaxE_LoxE', 'LaLoE_LaxLo', 'LaLoE_LaxE', 'LaLoE_LoxE', 'LaxLo', 'LaxE', 'LoxE', 'LaLo', 'LaE', 'LoE', 'Lat', 'Long', 'Elev')

# LMM model comparison results: fit_LaxLo is best
aictab(cand.set = models, modnames = model_names)

# hypothesis testing for seed weight vs temp and precip 
sw_clim_data <- sd_wt_data %>% left_join(dplyr::select(clim_df, population, CHELSA_bio10_01,CHELSA_bio10_09, CHELSA_bio10_12, CHELSA_bio10_17, CHELSA_bio10_10, CHELSA_bio10_18 )) %>%
  rename("MAT"="CHELSA_bio10_01", "AP"="CHELSA_bio10_12", "TDQ"="CHELSA_bio10_09", "PDQ"="CHELSA_bio10_17", "TWQ"="CHELSA_bio10_10", "PWQ"="CHELSA_bio10_18") %>%
  mutate(MAT_s=scale(MAT), AP_s=scale(AP), TDQ_s=scale(TDQ), PDQ_s=scale(PDQ), TWQ_s=scale(TWQ), PWQ_s=scale(PWQ))

fit_sw_clim <- lmer(sd_wt_50_ct ~ MAT_s*AP_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #only MAT is significant
fit_sw_clim2 <- lmer(sd_wt_50_ct ~ TDQ_s*PDQ_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #TDQ and interaction are significant.
fit_sw_clim3 <- lmer(sd_wt_50_ct ~ TWQ_s*PWQ_s + (1|population) + (1|block) + (1|population:block), data = sw_clim_data) #TWQ and interaction are significant
summary(fit_sw_clim) 

# forks per stem vs climate
forks_clim_data <- stem_data %>% left_join(dplyr::select(clim_df, population, CHELSA_bio10_01,CHELSA_bio10_09, CHELSA_bio10_12, CHELSA_bio10_17, CHELSA_bio10_10, CHELSA_bio10_18 )) %>%
  rename("MAT"="CHELSA_bio10_01", "AP"="CHELSA_bio10_12", "TDQ"="CHELSA_bio10_09", "PDQ"="CHELSA_bio10_17", "TWQ"="CHELSA_bio10_10", "PWQ"="CHELSA_bio10_18") %>%
  mutate(MAT_s=scale(MAT), AP_s=scale(AP), TDQ_s=scale(TDQ), PDQ_s=scale(PDQ), TWQ_s=scale(TWQ), PWQ_s=scale(PWQ))

fit_forks_clim <- lmer(forks ~ TDQ_s*PDQ_s + (1|population) + (1|population:block) + (1|population:block:plant), data = forks_clim_data)
summary(fit_forks_clim)


#### extract population-level means from REML fit using coef()
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

# Plot pop-level means vs temp
plot_sw_clim3 <- ggplot(data=ml_pred_df) +
  geom_abline(intercept=fixef(fit_sw_clim3)[1], slope=fixef(fit_sw_clim3)[3], col="blue", lty=2) +
  geom_point(mapping=aes(x=TWQ_s, y=pop_b0)) +
  geom_linerange(mapping=aes(x=TWQ_s, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se)) +
  labs(x="Mean Temp Warmest Quarter (scaled)", y="Estimated intercept in population l")

# Visualize the interaction effect by turning precip into categorical variable (split into 3 groups)
ml_pred_df$cat_PWQ_s <- as.factor(as.numeric(cut_number(ml_pred_df$PWQ_s, 3)))
plot_sw_clim3x <- ggplot(data=ml_pred_df) + 
  geom_point(mapping=aes(x=TWQ_s, y=pop_b0, color=cat_PWQ_s)) +
  geom_linerange(mapping=aes(x=TWQ_s, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se))

interact_plot(fit_sw_clim3, pred=TWQ_s, modx=PWQ_s, interval=TRUE)
sim_slopes(fit_sw_clim3, pred = TWQ_s, modx = PWQ_s, johnson_neyman = FALSE)

# Geo interaction plot
interact_plot(fit_LaxLo, pred = Lat_s, modx = Long_s, interval=T, int.width = 0.95)
  

#### PCA/RDA as a multivariate approach. helpful: http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
install.packages("vegan")
devtools::install_github("gavinsimpson/ggvegan")
library(vegan)
library(ggvegan)
pop_trait_means <- read.csv("data/pop_trait_means.csv", header = T)
temp <- pop_trait_means[,-1]
rownames(temp) <- pop_trait_means[,1]
pop_trait_means <- temp

my_trait_rda <- rda(pop_trait_means, scale = T) #scale everything bc they are vastly different units. Also don't use EST_YIELD (14th col) since that was a composite index of the other values
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

# extract df of PC scores to manually build a plot
pops_rda <- fortify(my_trait_rda, display='sites')
traits_rda <- fortify(my_trait_rda, display='species')

ggplot() +
  geom_point(data=pops_rda, aes(x = PC1, y = PC2)) +
  geom_text(data=pops_rda, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=3) +
  geom_segment(data=traits_rda, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=traits_rda, 
            aes(x=PC1,y=PC2,label=Label,
                hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
            color="red", size=4)

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


#### k means clustering. not working currently 1.13.21
library(cluster)
library(factoextra)
library(purrr)

df <- pop_trait_means
# Function to compute total within-cluster sum of square
wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

wssplot(pop_trait_means)

# continue with k=4
clusters <- kmeans(df, 4, nstart = 25)
print(clusters)
fviz_cluster(clusters, data = df) # I don't really like this plot
