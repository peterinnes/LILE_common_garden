#' ---
#' output: github_document
#' ---

#' trait-environment correlations
#' Peter Innes
#' 10.5.20 / last updated 1.15.21

#+ results=FALSE, message=FALSE, warning=FALSE
#devtools::install_github("jimhester/archive")
#devtools::install_github("mirzacengic/climatedata")
#install.packages("rgdal")
#install.packages("AICcmodavg")
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
colnames(clim_df)[4:22] <- lapply(colnames(clim_df)[4:22], gsub, pattern = "CHELSA_bio10_", replacement = "bio") #simplify column names

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

#' ########
#' #### Hypothesis testing for latitudinal clines
#' ########

# This is a function that takes a LMM as its argument and returns a data frame with estimated group means (intercepts). I use it to find the population trait means of trait~latitude models, which have population as a random effect.
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
sd_wt_data <- sd_wt_data %>% inner_join(dplyr::select(env_data,source,population, Lat)) 
fit_sw_Lat <- lmer(sd_wt_50_ct ~ Lat + (1|population) + (1|block) + (1|population:block),
                   data=sd_wt_data)

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

yield_df <- yield_df %>% inner_join(dplyr::select(env_data, population, source, Lat))
fit_yield_Lat <- lmer(log(EST_YIELD) ~ Lat + (1|population) + (1|block) + (1|population:block), data = yield_df)

# Make lists and storage for the for() loop
fit_list <- c(fit_sw_Lat, fit_ff_Lat, fit_ns_Lat, fit_fruits_Lat, fit_bf_Lat, fit_forks_Lat, fit_stemd_Lat, fit_capsd_Lat, fit_yield_Lat)
trait_list <- c("Seed weight", "Fruit fill", "Stems", "Fruits", "Buds/flowers", "Forks", "log Stem dia", "Capsule dia", "log Est. yield")
plot_list = list()
# Loop through each model, calculating population means and confidence intervals of regression lines
for (i in 1:length(fit_list)) {
  fit <- fit_list[[i]]
  pred_df <- make_pred_df(fit) #get population means
  
  # Obtain confidence interval for regression line
  newd <- data.frame(Lat = seq(min(geo_data$Lat, na.rm=T), max(geo_data$Lat, na.rm=T), length.out=100))
  lmm_boots <- bootMer(fit, predict_fun, nsim = 100)
  pred_ci <- cbind(newd, confint(lmm_boots))
  
  plot <- ggplot(data=pred_df) +
    geom_abline(intercept=fixef(fit)[1], slope=fixef(fit)[2], lty=2) +
    geom_point(mapping=aes(x=Lat, y=pop_b0), color="royalblue2", alpha=0.5) +
    geom_linerange(mapping=aes(x=Lat, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se), color="royalblue2", alpha=0.5) +
    geom_ribbon(data=pred_ci, aes(x=Lat, ymin=`2.5 %`, ymax=`97.5 %`), alpha=0.25) +
    labs(x="Lat", y=paste(trait_list[[i]]))
  
  plot_list[[i]] <- plot
}
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggdraw(add_sub(p, "Latitude", vpadding=grid::unit(0,"lines"),y=6, x=0.53, vjust=4.5))

png("traits_vs_Lat.png", width=8, height=7, res=300, units="in")
p
dev.off()


#' ########
#' #### Model selection of traits vs climate predictors. bio1, bio12, bio10, bio18, bio3, bio15, bio4
#' ########

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

#### PCA/RDA as a multivariate approach. helpful: http://dmcglinn.github.io/quant_methods/lessons/multivariate_models.html
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
