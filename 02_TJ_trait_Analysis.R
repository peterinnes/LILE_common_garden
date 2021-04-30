# Analyzing Tom Jones' data
# 1.30.21
# Peter Innes

library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(emmeans)
library(vegan)
library(ggvegan)

# Key to convert Tom 'Entry' number to Stan 'source' number/population ID.

TvS_key <- read.csv("info_from_TomJ/tomJ_data/StanVsTomAccessions.csv", header=T) %>%
  mutate(source=factor(source), Entry=factor(Entry))

#### Stem and capsule data (Capsules per stem, capsules per plant, stems per plant) ####
tj_stems_caps <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) %>% #Spreadsheet has two headers. so we gotta skip the first one
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))

# summarise. Entry 5 has 10 reps when each entry should only have 8. Data mis-entered
sc_summary <- tj_stems_caps %>%
  group_by(Entry) %>%
  summarise(n=n())

# Fix data entry error, change plots 134 and 265 from Entry 5 to Entry 6
tj_stems_caps$Entry[134] <- 6
tj_stems_caps$Entry[265] <- 6
Plot_Entry_key <- dplyr::select(tj_stems_caps, Plot, Entry) #Save the corrected Plot-Entry combinations in order to fix other trait data with same issue

tj_stems_caps <- tj_stems_caps %>% 
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38)) #filter out mistaken Appar sources

names(tj_stems_caps)[12] <- "surv_4_27_13" #Rename last column since we lost the column name when initially skipping first line. This column has number of plants surviving on 4.27.13 (out of 10 plants initially planted I believe).

dim(tj_stems_caps)
head(tj_stems_caps)

# Linear models
# Capsules per stem. log transform.
fit_CPS <- lmer(log(Capsules/sub_stems) ~ -1 + Entry + (1|Rep), data = tj_stems_caps) 
summary(fit_CPS)
plot(fit_CPS)
qqnorm(resid(fit_CPS))

# Caps per plant: Capsules/sub_stems * ttl_stems = est_ttl_capsules. Then divide est_ttl_capules by the number of surviving plants in the plot. Using the survival data from the tj_stems_caps spreadsheet--this closely matches survival counts taken later in August, closer to Harvest (but this Aug data is incomplete so we don't want to use it directly)
# square root transform, this is what Tom used.
tj_stems_caps$CPP_tr <- sqrt(tj_stems_caps$est_ttl_capsules/tj_stems_caps$surv_4_27_13)


fit_CPP <- lmer(CPP_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_CPP)
qqnorm(resid(fit_CPP))
summary(fit_CPP)

# total capsules (i.e. capsules per plot. incorporates survival)
tj_stems_caps$ttl_caps_tr <- sqrt(tj_stems_caps$est_ttl_capsules)
fit_ttl_caps <- lmer(ttl_caps_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_ttl_caps)

# Stems per plant. Use same sqrt transform as for capsules per plant
tj_stems_caps$spp_tr <- sqrt(tj_stems_caps$ttl_stems/tj_stems_caps$surv_4_27_13) 
fit_SPP <- lmer(spp_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)

plot(fit_SPP)
qqnorm(resid(fit_SPP))
summary(fit_SPP) # How to estimate capsules per plant? est_ttl_capules / surviv_8_30_13?

# Stems per plot (total stems. incorporates survival)
tj_stems_caps$ttl_stems_tr <- sqrt(tj_stems_caps$ttl_stems)
fit_ttl_stems <- lmer(ttl_stems_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_ttl_stems)

#### Biomass data ####
tj_biomass <- read.csv("data/TomJ_BiomassFlax2013_14_15_ANALYZE_THIS.csv", skip=1, header = T) %>%
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))
names(tj_biomass) <- c("Plot", "Rep", "Entry", "num_planted", "ttl_weight_2013", "ttl_stem_2013", "sub_weight_2013", "ttl_weight_2014", "DW_2015", "survivorship_4_27_13")
tj_biomass[4:10] <- as.integer(unlist(tj_biomass[4:10]))

# Summary. Entry 12 has 10 reps; should only have 8; Entries 1, 13, 15, 17, 25 have 9 reps. Need to replace the Plot vs Entry assignments from biomass data with that from the stems/caps data.
bm_summary <- tj_biomass %>%
  group_by(Entry) %>%
  summarise(n=n())

tj_biomass$Entry <- Plot_Entry_key$Entry #replace Plot-Entry matches with corrected version from stem/caps data. 

# Join with Stan's 'source' IDs, filter out mistaken Appar accessions.
tj_biomass <- tj_biomass %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38))

dim(tj_biomass)

bm_summary2 <- full_join(tj_biomass, dplyr::select(tj_stems_caps, Plot, entry = Entry)) %>%
  group_by(entry) %>%
  summarise(n=n())

# Linear models
# Total weight (biomass per plot, which incorporates survival)
tj_biomass$ttl_weight_2013_tr <- sqrt(tj_biomass$ttl_weight_2013) #square root transform

fit_biomass <- lmer(ttl_weight_2013_tr ~ -1 + Entry + (1|Rep), data=tj_biomass)
plot(fit_biomass)
qqnorm(resid(fit_biomass))

# Biomass per plant
tj_biomass$BPP <- tj_biomass$ttl_weight_2013 / tj_biomass$survivorship_4_27_13 
tj_biomass$BPP[tj_biomass$BPP == 0] <- NA #convert Zeros in BPP to NA
tj_biomass$BPP_tr <- sqrt(tj_biomass$BPP) #square root transform

fit_bpp <- lmer(BPP_tr ~ -1 + Entry + (1|Rep), data=tj_biomass)
plot(fit_bpp)
qqnorm(resid(fit_bpp))

biomass_lsmeans <- as.data.frame(lsmeans(fit_bpp, "Entry")) %>%
  inner_join(dplyr::select(TvS_key, source, Entry)) %>%
  mutate(btr_lsmean=(lsmean^2))

biomass_lsmeans <- as.data.frame(lsmeans(fit_bpp, "Entry"))
biomass_lsmeans[c(2,3,5,6)] <- lapply(biomass_lsmeans[c(2,3,5,6)], function(x) x^2)
                                   
#### Canopy height and plant diameter data ####
ht_dia_data <- read.csv("data/TomJ_flax_avg_ht_dia_2013_RAW_DATA.csv", header = T, na.strings = '.') %>%
  mutate(Plot=as.factor(Plot), Rep=as.factor(Rep), Entry=as.factor(Entry)) 
ht_dia_data <- head(ht_dia_data, -3) #delete 3 empty rows at end

#View(ht_dia_data %>%
#       group_by(Entry) %>%
#       summarise(n=n()))

ht_dia_data$Entry <- Plot_Entry_key$Entry #replace Plot-Entry matches with corrected version from stem/caps data.

ht_dia_data <- ht_dia_data %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38)) #filter out mistaken Appar sources

# Separate height data from diameter and convert each to long format
ht_data <- ht_dia_data %>% dplyr::select(Plot, Rep, Entry, ht_pl1, ht_pl2, ht_pl3, ht_pl4, ht_pl5)
ht_data <- pivot_longer(ht_data, names_to = "Plant", values_to = "height", ht_pl1:ht_pl5) 

dia_data <- ht_dia_data %>% dplyr::select(Plot, Rep, Entry, dia.1:dia.10)
dia_data <- pivot_longer(dia_data, names_to = "Plant", values_to = "dia", dia.1:dia.10) 

# How many plants in each plot? We will use this number to corroborate survivorship, in order to estimate per plant capsule and stem numbers. Diam and Height were reportedly measured a week before harvest in October. So they should be the most accurate? # Should also compare to April and Aug surv data. Could corroborate with caps data as well. 
dia_plants_per_plot <- dia_data %>%
  na.omit() %>%
  group_by(Plot) %>%
  summarise(num_plants_diam = n())

# Linear models
fit_ht <- lmer(height ~ -1 + Entry + (1|Plot), data = ht_data)
summary(fit_ht) #adding random effect of plot because we have measurements of multiple individual plants per plot. Rep (block) random effect is zero so we will remove.
plot(fit_ht)

fit_dia <- lmer(dia ~ -1 + Entry + (1|Rep) + (1|Plot), data = dia_data)
summary(fit_dia)
plot(fit_dia)
qqnorm(resid(fit_dia))

#### Survival data ####
tj_surviv <- read.csv("data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T) %>%
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))
# rename columns and keep only those we need
tj_surviv <- tj_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13, X5.1.14, X8.12.14)
names(tj_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13", "surviv_5_1_14", "surviv_8_12_14") 
tj_surviv <- head(tj_surviv, -1) #delete last row which had totals

# Entry 5 also has 10 reps here in this dataset, need to correct.
#View(tj_surviv %>%
#       group_by(Entry) %>%
#       summarise(n=n()))

tj_surviv$Entry <- Plot_Entry_key$Entry #replace Plot-Entry matches with corrected version from stem/caps data.

# join w/ source/site info, remove mistaken Appar sources
tj_surviv <- tj_surviv %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38))
# Convert characters to integers
tj_surviv[5:10] <- as.integer(unlist(tj_surviv[5:10])) 

dim(tj_surviv)

# Gather all evidence for survival numbers at first harvest to compare and decide what values to use.
surv <- full_join(dplyr::select(tj_surviv, Plot, surviv_4_27_13, surviv_8_30_13), dia_plants_per_plot) %>% #survival as noted in survival spreadsheet joined w/ num plants measured for diam
  full_join(dplyr::select(tj_stems_caps, Plot, surv_4_27_13)) %>% #survival from stem/caps spreadsheet
  rename(caps_surv_4_27_13 = surv_4_27_13) %>% 
  full_join(dplyr::select(tj_biomass, Plot, survivorship_4_27_13)) %>% #survival from biomass spreadsheet (same as that in stem/caps spreadsheet)
  rename(biomass_surv_4_27_13 = survivorship_4_27_13) %>%
relocate(surviv_4_27_13, .before = num_plants_diam) %>%
  relocate(num_plants_diam, .after = surviv_8_30_13) %>%
  relocate(caps_surv_4_27_13, .after = surviv_4_27_13)


# Convert "." to 0 for column 5 and onwards. Not sure we actually want to do this, in some cases the '.' seems to denote no survival, in other cases it seems to denote missing data. May have to manually corroborate survival.
#na_code <- "."
#tj_surviv [5:9] <- apply(tj_surviv[5:9], 2, function(x){ as.numeric(ifelse(x %in% na_code, 0, x)) } )
#
#tj_surviv[5:10] <- as.integer(unlist(tj_surviv[4:9]))
#tj_surviv[1:3] <- as.factor(unlist(tj_surviv[1:3]))
#

# Linear models
# First year survival.
fit_surviv <- lmer(surviv_8_30_13/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv)
plot(fit_surviv)

# Second year survival. 
fit_surviv2 <- lmer(surviv_8_12_14/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv2)
plot(fit_surviv2)

#### Summarize ls-means of all traits ####

#tj_fit_list <- c(fit_CPS, fit_CPP, fit_SPP, fit_bpp, fit_ht, fit_dia, fit_surviv)
#tj_trait_list <- c("Capsules_per_stem", "Capsules_per_plant", "Stems_per_plant", "Biomass_per_plant",  "Plant_height", "Plant_diameter", "Survival")

tj_fit_list <- c(fit_CPS, fit_ttl_caps, fit_ttl_stems, fit_biomass, fit_ht, fit_dia)
tj_trait_list <- c("Capsules_per_stem", "Capsules_per_plot", "Stems_per_plot", "Biomass_per_plot",  "Plant_height", "Plant_diameter")

tj_results <- list() #list to store means and confidence intervals
tj_esp_list <- list() #list to store effect size plots

# confidence intervals taken from the lsmeans() fxn here bc it is fast. will pronbably want to use confint() for final pub instead.
for (i in 1:length(tj_fit_list) ){
  
  fit <- tj_fit_list[[i]]
  lsmeans <- as.data.frame(lsmeans(fit, "Entry"))
  
  # Conditional statements to back-transform lsmeans, SEs, and CIs of appropriate traits for plotting
  # Square root transformed traits
  if( tj_trait_list[i] %in% c("Capsules_per_plot", "Stems_per_plot", "Biomass_per_plot") ){
    plot_lsmeans <- as.data.frame(lsmeans(fit, "Entry"))
    plot_lsmeans[c(2,3,5,6)] <- lapply(plot_lsmeans[c(2,3,5,6)], function(x) x^2)
  
  # Log transformed, so esponentiate 
  } else if( tj_trait_list[i]=="Capsules_per_stem" ) {
    plot_lsmeans <- as.data.frame(lsmeans(fit, "Entry"))
    plot_lsmeans[c(2,3,5,6)] <- lapply(plot_lsmeans[c(2,3,5,6)], function(x) exp(x))
  
  # Untransformed traits
  } else {
    plot_lsmeans <- as.data.frame(lsmeans(fit, "Entry"))
  }
  plot_lsmeans <- plot_lsmeans %>% #join with 'source' ID from Stan
    inner_join(dplyr::select(TvS_key, Entry, source)) %>%
    relocate(source, .after = Entry)
  # Create effect size plot
  esp <- ggplot(data=plot_lsmeans, aes(x=lsmean, y=reorder(source, lsmean), xmin=lower.CL, xmax=upper.CL)) +
    geom_point() +
    geom_errorbar() +
    ylab("Population") +
    xlab(tj_trait_list[[i]]) +
    theme(axis.text.y = element_text(size = 6))
  
  tj_esp_list[[i]] <- esp #Store plot
  
  names(lsmeans)[2] <- tj_trait_list[[i]] #Change 'lsmean' column name to trait name before storing in results
  tj_results[[i]] <- lsmeans #store means and confidence intervals
  
  # sort dataframe and write to csv for summary
  lsmeans <- lsmeans %>% arrange(-lsmeans[2]) #sort descending trait value to make more readable
  #write.csv(means_ci, file=paste0("results_summaries/", names(means_ci)[4], "_summary.csv"))
}
names(tj_results) <- tj_trait_list
head(tj_results, 7L)

tj_esp_grid <- cowplot::plot_grid(plotlist = tj_esp_list, ncol = 3)
png("plots/tj_traits_esp.png", width=12, height=9, res=300, units="in")
x11()
tj_esp_grid
dev.off()

#### Trait correlations & PCA ####
# Get population-level trait means
tj_pop_trait_means <- data.frame(Entry=tj_results[[1]]$Entry)
for ( i in 1:length(tj_results) ){
  tj_pop_trait_means <- cbind(tj_pop_trait_means, tj_results[[i]][2]) #add trait means to growing df
}
scaled_tj_ptm <- scale(tj_pop_trait_means[3:8], center = T, scale = T)
#write.csv(tj_pop_trait_means, file="data/tj_pop_trait_means.csv", row.names = FALSE)

# PCA with source as rownames/labels
tj_pop_trait_means <- tj_pop_trait_means %>% #join with 'source' ID from Stan
  inner_join(dplyr::select(TvS_key, Entry, source)) %>%
  relocate(source, .after = Entry)

temp <- tj_pop_trait_means[,-1]
rownames(temp) <- temp[,1]
temp <- temp[,-1]

tj_trait_rda <- rda(temp, scale = T) #scale and center everything bc they are vastly different units.
biplot(tj_trait_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(tj_trait_rda, dis="sites", cex=0.5)

# ggvegan version
autoplot(tj_trait_rda, arrows = TRUE, geom = "text", legend = "none") #basic


#### Trait correlations. Use scaled/centered population means? ####

#### Latitudinal clines ####
fit_list <- c()
tj_trait_list <- c("Capsules_per_stem", "Capsules_per_plot", "Stems_per_plot", "Biomass_per_plot",  "Plant_height", "Plant_diameter")
datasets <- list(tj_stems_caps, tj_stems_caps, tj_stems_caps, tj_biomass, ht_data, dia_data)
plot_list = list()

# Loop through each model, calculating population means and confidence intervals of regression lines
for (i in 1:length(fit_list)) {
  fit <- fit_list[[i]]
  data <- datasets[[i]] %>%
    inner_join(dplyr::select(TvS_key, Entry, source)) %>%
    inner_join(dplyr::select(env_data, source, Lat))
  
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
    labs(x="Lat", y=paste(tj_trait_list[[i]])) +
    theme_minimal()
  
  plot_list[[i]] <- plot
}
p <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggdraw(add_sub(p, "Latitude", vpadding=grid::unit(0,"lines"),y=6, x=0.53, vjust=4.5))


#### PC transfer distance test for local adaptation ####
library(climatedata)
library(sp)

milville_coords <- data.frame(Long=-111.816, Lat=41.656)
chelsa <- get_chelsa(type = "bioclim", layer = 1:19, period = c("current"))

milville_point <- SpatialPoints(milville_coords, proj4string = chelsa@crs)

milville_value <- data.frame(raster::extract(chelsa, milville_point)) #previously raster::extract(r,points)
colnames(milville_value) <- lapply(colnames(milville_value), gsub, pattern = "CHELSA_bio10_", replacement = "bio") #simplify column names

milville_clim <- cbind.data.frame(milville_coords, milville_value) %>%
  mutate(source="MILVILLE_GARDEN", population=NA, Elev_m=1407)

loc_adapt_df <- full_join(geo_clim_df, milville_clim)
rownames(loc_adapt_df) <- loc_adapt_df[,1] #set source number to the rownames, otherwise we lose these labels in the PCA below.

loc_adapt_rda <- rda(loc_adapt_df[3:24], scale = T)
summary(loc_adapt_rda)
# Get site (source/population) PC scores for use in trait-env model selection
loc_adapt_PC_scores <- data.frame(scores(loc_adapt_rda, choices=1:3, display = "sites", scaling=0)) %>%
  tibble::rownames_to_column("source")

milville_garden_pc1 <- filter(loc_adapt_PC_scores, source=="MILVILLE_GARDEN")$PC1
milville_garden_pc2 <- filter(loc_adapt_PC_scores, source=="MILVILLE_GARDEN")$PC2

dist_from_garden <- data.frame(source=loc_adapt_PC_scores$source, pc_trd=(milville_garden_pc1 - loc_adapt_PC_scores$PC1))

biomass_vs_dist_df <- inner_join(dist_from_garden, dplyr::select(biomass_lsmeans, source, biomass=btr_lsmean)) 

ggplot(data=biomass_vs_dist_df, aes(x=pc_trd, y= biomass)) +
  geom_point() +
  geom_text(aes(label = source))

biomass_vs_dist_df$pc_trd2 <- biomass_vs_dist_df$pc_trd^2
biomass_trd_fit <- lm(biomass ~ pc_trd + pc_trd2, data=biomass_vs_dist_df)
summary(biomass_trd_fit)
