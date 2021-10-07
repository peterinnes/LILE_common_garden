# Analyzing Tom Jones' data
# 1.30.21
# Peter Innes

library(raster)
library(sp)
library(rgdal)
#library(remotes)
library(climatedata)
library(raster)
library(sp)
library(rgdal)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(emmeans)
library(pbkrtest)
library(multcomp)
library(multcompView)
library(sjstats)
library(vegan)
#library(remotes)
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)

options(contrasts = c("contr.sum","contr.poly"))

#### Env/Climate data ####
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

#### Stem and capsule data (Capsules per stem, capsules per plant, stems per plant) ####
# Key to convert Tom 'Entry' number to Stan 'source' number/population ID.

TvS_key <- read.csv("info_from_TomJ/tomJ_data/StanVsTomAccessions.csv", header=T) %>%
  mutate(source=factor(source), Entry=factor(Entry))

mv_stems_caps <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) %>% #Spreadsheet has two headers. so we gotta skip the first one
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))

# summarise. Entry 5 has 10 reps when each entry should only have 8. Data mis-entered
sc_summary <- mv_stems_caps %>%
  group_by(Entry) %>%
  summarise(n=n())

# Fix data entry error, change plots 134 and 265 from Entry 5 to Entry 6
mv_stems_caps$Entry[134] <- 6
mv_stems_caps$Entry[265] <- 6
Plot_Entry_key <- dplyr::select(mv_stems_caps, Plot, Entry) #Save the corrected Plot-Entry combinations in order to fix other trait data with same issue

mv_stems_caps <- mv_stems_caps %>% 
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38)) #filter out mistaken Appar sources

names(mv_stems_caps)[12] <- "surv_4_27_13" #Rename last column since we lost the column name when initially skipping first line. This column has number of plants surviving on 4.27.13 (out of 10 plants initially planted I believe).

dim(mv_stems_caps)
head(mv_stems_caps)
mv_stems_caps <- inner_join(mv_stems_caps, dplyr::select(env_data, source, population)) %>%
  relocate(population, .after = source)

# Linear models

# 1. Capsules per STEM. log transform.
mv_stems_caps$CPS <- mv_stems_caps$Capsules/mv_stems_caps$sub_stems
fit_CPS <- lmer(log(CPS) ~ population + (1|Rep), data = mv_stems_caps) 
#plot(fit_CPS) THIS IS CAUSING R TO ABORT first encountered 6.11.21 ... whish was after system freeze, hard shut down, and subsequent updates to the OS.
qqnorm(resid(fit_CPS))

#fit_CPS2 <- lmer(log(Capsules/sub_stems) ~ (1|population) + (1|Rep), data = mv_stems_caps) #also fit model with entry (population) as random effect
coef(fit_CPS2)$population

# Caps per plant: Capsules/sub_stems * ttl_stems = est_ttl_capsules. Then divide est_ttl_capules by the number of surviving plants in the plot. Using the survival data from the mv_stems_caps spreadsheet--this closely matches survival counts taken later in August, closer to Harvest (but this Aug data is incomplete so we don't want to use it directly)
# square root transform, this is what Tom used.
#fit_CPP <- lmer(sqrt(CPP_tr) ~ -1 + Entry + (1|Rep), data = mv_stems_caps)


# 2. Total capsules per PLOT (incorporates survival). square root transform
fit_ttl_caps <- lmer(sqrt(est_ttl_capsules) ~ population + (1|Rep), data = mv_stems_caps)
#fit_ttl_caps2 <- lmer(sqrt(est_ttl_capsules) ~ (1|population) + (1|Rep), data = mv_stems_caps)


## Stems per plant. Use same sqrt transform as for capsules per plant
#mv_stems_caps$spp_tr <- sqrt(mv_stems_caps$ttl_stems/mv_stems_caps$surv_4_27_13) 
#fit_SPP <- lmer(spp_tr ~ -1 + Entry + (1|Rep), data = mv_stems_caps)

# 3. Stems per PLOT (incorporates survival). square root-transform
fit_ttl_stems <- lmer(sqrt(ttl_stems) ~ population + (1|Rep), data = mv_stems_caps)
#fit_ttl_stems2 <- lmer(sqrt(ttl_stems) ~ (1|population) + (1|Rep), data = mv_stems_caps)


#### Biomass data ####
mv_biomass <- read.csv("data/TomJ_BiomassFlax2013_14_15_ANALYZE_THIS.csv", skip=1, header = T) %>%
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))
names(mv_biomass) <- c("Plot", "Rep", "Entry", "num_planted", "ttl_weight_2013", "ttl_stem_2013", "sub_weight_2013", "ttl_weight_2014", "DW_2015", "survivorship_4_27_13")
mv_biomass[4:10] <- as.integer(unlist(mv_biomass[4:10]))

# Summarise to make sure plots are labeled w/ correct 'Entry' (accession aka population). Entry 12 has 10 reps; should only have 8; Entries 1, 13, 15, 17, 25 have 9 reps. Need to correct the Plot vs Entry assignments from biomass data with that from the stems/caps data.
bm_summary <- mv_biomass %>%
  group_by(Entry) %>%
  summarise(n=n())
bm_summary2 <- full_join(mv_biomass, dplyr::select(mv_stems_caps, Plot, entry = Entry)) %>%
  group_by(entry) %>%
  summarise(n=n())

mv_biomass$Entry <- Plot_Entry_key$Entry #replace Plot-Entry matches with corrected version from stem/caps data. 

# Join with Stan's 'source' IDs, filter out mistaken Appar accessions.
mv_biomass <- mv_biomass %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38)) %>%
  inner_join(dplyr::select(env_data, source, population))

dim(mv_biomass) #should be same as the stem and caps data, 254 rows.

# Linear models
# 4. 2013 Biomass per PLOT (incorporates survival). square root transform.
fit_2013biomass <- lmer(sqrt(ttl_weight_2013) ~ population + (1|Rep), data=mv_biomass)
#fit_2013biomass2 <- lmer(sqrt(ttl_weight_2013) ~ (1|population) + (1|Rep), data=mv_biomass)

# 5. 2014 Biomass per PLOT (incorporates survival). square root transform.
fit_2014biomass <- lmer(sqrt(ttl_weight_2014) ~ population + (1|Rep), data=mv_biomass)
#fit_2014biomass2 <- lmer(sqrt(ttl_weight_2014) ~ (1|population) + (1|Rep), data=mv_biomass)
 
## Biomass per plant (2013). Again, skipping per-plant models bc survival data is unreliable
#mv_biomass$BPP <- mv_biomass$ttl_weight_2013 / mv_biomass$survivorship_4_27_13 
#mv_biomass$BPP[mv_biomass$BPP == 0] <- NA #convert Zeros in BPP to NA
#fit_bpp <- lmer(sqrt(BPP) ~ -1 + Entry + (1|Rep), data=mv_biomass)

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
  filter(!source %in% c(2,5,22,32,38)) %>% #filter out mistaken Appar sources
  inner_join(dplyr::select(env_data, source, population))


# Separate height data from diameter and convert each to long format
ht_data <- ht_dia_data %>% dplyr::select(Plot, Rep, population, ht_pl1, ht_pl2, ht_pl3, ht_pl4, ht_pl5)
ht_data <- pivot_longer(ht_data, names_to = "Plant", values_to = "height", ht_pl1:ht_pl5) 

dia_data <- ht_dia_data %>% dplyr::select(Plot, Rep, population, dia.1:dia.10)
dia_data <- pivot_longer(dia_data, names_to = "Plant", values_to = "dia", dia.1:dia.10) 

# Summarize plot averages just for fun
ht_plot_means <- ht_data %>% group_by(Plot, population) %>% na.omit() %>% summarise(mean=mean(height))

dia_plot_means <- dia_data %>% group_by(Plot, population) %>% na.omit() %>% summarise(mean=mean(dia))

## How many plants in each plot? Want use this number to corroborate survivorship, in order to estimate per plant capsule and stem numbers. Diam and Height were reportedly measured a week before harvest in October. So they should be the most accurate? # Should also compare to April and Aug surv data. Could corroborate with caps data as well?
#dia_plants_per_plot <- dia_data %>%
#  na.omit() %>%
#  group_by(Plot) %>%
#  summarise(num_plants_diam = n())

# Linear models
# 6. PLANT height
fit_ht <- lmer(height ~ population + (1|Plot), data = ht_data) #adding random effect of plot because we have measurements of multiple individual plants per plot. variance due to Rep (block) random effect is essentially zero so we will remove.
#fit_ht2 <- lmer(height ~ population + (1|Rep), data = ht_data) #model with Rep instead of Plot RE, used to calculate coefficient of variation

# 7. PLANT diameter
fit_dia <- lmer(dia ~ population + (1|Rep) + (1|Plot), data = dia_data)
#fit_dia2 <- lmer(dia ~ population + (1|Rep), data = dia_data) #model without Plot RE, used to calculate coefficient of variation

#### Survival data ####
## Skip analysis of survival except for use in transfer distance test because we do not have survival measured at exact date of harvest. Also, survival is ultimately taken into account for plot-level measurements of number of stems, biomass, etc.
mv_surviv <- read.csv("data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T) %>% 
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))
# rename columns and keep only those we need
mv_surviv <- mv_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13, X5.1.14, X8.12.14)
names(mv_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13", "surviv_5_1_14", "surviv_8_12_14") 
mv_surviv <- head(mv_surviv, -1) #delete last row which had totals

# Entry 5 also has 10 reps here in this dataset, need to correct.
#View(mv_surviv %>%
#       group_by(Entry) %>%
#       summarise(n=n()))

mv_surviv$Entry <- Plot_Entry_key$Entry #replace Plot-Entry matches with corrected version from stem/caps data.

# join w/ source/population info, remove mistaken Appar sources
mv_surviv <- mv_surviv %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  inner_join(dplyr::select(env_data, source, population)) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38))
head(mv_surviv)
# Convert characters to integers
mv_surviv[5:10] <- as.integer(unlist(mv_surviv[5:10])) 

test <- dplyr::inner_join(mv_surviv, mv_biomass)
cor.test(test$surviv_8_12_14, test$ttl_weight_2014)
           
#dim(mv_surviv)

## Gather all evidence for survival numbers at first harvest to compare and decide what values to use.
#surv <- full_join(dplyr::select(mv_surviv, Plot, surviv_4_27_13, surviv_8_30_13), dia_plants_per_plot) %>% #survival as noted in survival spreadsheet joined w/ num plants measured for diam
#  full_join(dplyr::select(mv_stems_caps, Plot, surv_4_27_13)) %>% #survival from stem/caps spreadsheet
#  rename(caps_surv_4_27_13 = surv_4_27_13) %>% 
#  full_join(dplyr::select(mv_biomass, Plot, survivorship_4_27_13)) %>% #survival from biomass spreadsheet (same as that in stem/caps spreadsheet)
#  rename(biomass_surv_4_27_13 = survivorship_4_27_13) %>%
#  relocate(surviv_4_27_13, .before = num_plants_diam) %>%
#  relocate(num_plants_diam, .after = surviv_8_30_13) %>%
#  relocate(caps_surv_4_27_13, .after = surviv_4_27_13)

# Linear models
# first over-winter survival.
fit_surviv <- lmer(surviv_4_27_13/planted_6_12_12 ~ population + (1|Rep), data = mv_surviv)
#summary(fit_surviv)
#plot(fit_surviv)

## third year survival. 
#fit_surviv2 <- lmer(surviv_8_12_14/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = mv_surviv)
#summary(fit_surviv2)
#plot(fit_surviv2)

#### Gather and plot ls-means of all traits ####
# Confidence intervals, taken from the lsmeans() aka emmeans() function of package emmeans 
mv_fit_list <- list(fit_ttl_caps, fit_CPS, fit_ttl_stems, fit_2013biomass, fit_2014biomass, fit_ht, fit_dia)
mv_trait_list <- c("Capsules_per_plot","Capsules_per_stem", "Stems_per_plot",  "Biomass_per_plot_2013","Biomass_per_plot_2014", "Height", "Plant_diameter")

mv_results <- list() #list to store means and confidence intervals
mv_results_bt <- list() ##list to store back-trasnformed means and confidence intervals 
#mv_esp_list <- list() #list to store effect size plots
mv_emm_list <- list() #list to store emm plot
#cld_list <- list() #list to store compact letter displays of each model

for (i in 1:length(mv_fit_list) ){
  
  fit <- mv_fit_list[[i]]
  lsmeans <- as.data.frame(lsmeans(fit, "population")) #as data frame
  emm1 <- lsmeans(fit, "population") #same as lsmeans above but dont convert to df
  emm2 <- as.data.frame(lsmeans(fit, "population", type="response")) #separate object for backtransformed lsmeans
  
  # Plotting means and CIs with emmeans:::plot.emmGrid. in order to reorder the populations on y-axis, we need to edit the .plot.srg function in emmeans package:
  # trace(emmeans:::.plot.srg, edit=T). Edit lines 239-240, change aes_() to aes() and delete tilde from in front of x and y variables. Then use reorder() on the y variable.
  emm_plot <- plot(emm1, type = "response", comparisons = T, colors = c("salmon", "blue", "black")) +
    theme_minimal() +
    xlab(mv_trait_list[i]) +
    ylab("")
  mv_emm_list[[i]] <- emm_plot #store plot
  
  # Compact letter display
  contrasts <- emmeans::emmeans(object=fit, type="response", pairwise ~ "population", adjust="tukey") #tests are on transformed scale but display on response scale
  cld <- emmeans:::cld.emmGrid(object=contrasts$emmeans, Letters=letters, sort=F)
  cld_df <- data.frame(cld)
  
  # Renaming columns and storing results
  names(lsmeans)[2] <- mv_trait_list[[i]] #Change 'lsmean' column name to trait name before storing in results
  mv_results[[i]] <- lsmeans #store means and confidence intervals
  lsmeans <- lsmeans %>% arrange(-lsmeans[2]) #sort descending trait value to make more readable
  # Same thing but with backtransformed clds/emms. This will be results table 2(?) in manuscript.
  names(cld_df)[2] <- mv_trait_list[[i]] 
  mv_results_bt[[i]] <- cld_df
  #emm2 <- emm2 %>% arrange(-emm2[2]) 
}
names(mv_results) <- mv_trait_list
names(mv_results_bt) <- mv_trait_list

# store emms in one dataframe with population as rowname. transformed variables are not back-transformed here.
mv_means_df <- data.frame(matrix(ncol = length(mv_trait_list), nrow = length(unique(mv_stems_caps$population))))
names(mv_means_df) <- mv_trait_list
rownames(mv_means_df) <- mv_results[[1]]$population
for (i in 1:length(mv_trait_list) ){
  mv_means_df[i]  <- mv_results[[i]][2]
}
mv_means_df <- mv_means_df %>% arrange(desc(Capsules_per_plot)) %>%
  #round(digits=2) %>%
  tibble::rownames_to_column("Accession") %>%
  relocate(Accession, .before = Capsules_per_plot)
#write.csv(mv_means_df, file="plots/millville_trait_means_table.csv", row.names = F)

# dataframe with back-transformed emmeans, without letters, for PCA/RDA.
mv_means_df2 <- data.frame(matrix(ncol = length(mv_trait_list), nrow = length(unique(mv_stems_caps$population))))
names(mv_means_df2) <- mv_trait_list
rownames(mv_means_df2) <- mv_results_bt[[1]]$population
for (i in 1:length(mv_trait_list) ){
  mv_means_df2[i] <- mv_results_bt[[i]][2]
}

# Store emms with clds in dataframe with column for population/accession. This will be for publication. 
mv_means_df3 <- data.frame(matrix(ncol=8, nrow=33))
names(mv_means_df3) <- names(mv_means_df)
mv_means_df3$Accession <- mv_results_bt[[1]]$population
for (i in 1:length(mv_trait_list) ){
  emm_sf <- data.frame(apply(mv_results_bt[[i]][c(2:6)], 1:2,
                             function(x) signif(x, 3))) %>%
    mutate(letter_group=mv_results_bt[[i]][7])# change sig figs
  mv_means_df3[i+1] <- apply(emm_sf[c(1,6)], 1, paste, collapse="") #combine emmeans and letters into single column
}
mv_means_df3 <- mv_means_df3 %>% arrange(desc(Capsules_per_plot))
#names(mv_means_df3)[2:8] <- c("Capsules per plot", "Capsules per stem", "Stems per plot", "2013 Biomass per plot (g)", "2014 Biomass per plot (g)", "Plant height (cm)", "Plant diameter (cm)")
write.csv(mv_means_df3, "plots/millville_trait_means_table.csv", row.names = F)

#### Cet coefficients of variation ####
# want to use cv() from sjstats. This function calculates the cv by first using performance::rmse() to find the root mean squared error of a mixed model. This is divided by the grand mean of the same model. problem is that it calculates the grand mean from untransformed response variables. So we will calculate cv manually, with the same performance::rmse() that the cv() function from sjstats uses.

mv_cvs <- c()
# Use the same models as above 
for (fit in mv_cv_fit_list) {
  #resid_sd <- data.frame(VarCorr(fit))$sdcor[2]
  rmse <- performance::rmse(fit)
  grand_mean <- fixef(fit)[1]
  cv <- rmse/grand_mean
  mv_cvs <- append(mv_cvs, cv)
}
names(mv_cvs) <- mv_trait_list


# Trying to figure out how to create a nice table
library(sjPlot)
tab_df(mv_means_df3, file="Milville_emms_cld.docx")

library(xtable)
print(xtable(mv_means_df3))
library(stargazer)
stargazer(mv_means_df3, type="text")
library(printr)
knitr::kable(mv_means_df3)
library(pander)
pander::pander(mv_means_df3)

# Join emm plots together 
mv_emm_grid <- cowplot::plot_grid(plotlist = mv_emm_list, ncol = 3) 

#layout <- "
#ABCD
#EFGH
#"
#mv_emm_patchwork <- mv_emm_list[[1]] + mv_emm_list[[2]] + mv_emm_list[[3]] + mv_emm_list[[4]] + mv_emm_list[[5]] + mv_emm_list[[6]] + mv_emm_list[[7]] + plot_layout(design = layout)

png("plots/millville_traits_comparisons.png", width=12, height=9, res=300, units="in")
#x11()
mv_emm_grid
dev.off()

#### Gather BLUPs/conditional modes (for use in PCA/RDA?) ####
#mv_fit_list2 <- c(fit_CPS2, fit_ttl_caps2, fit_ttl_stems2, fit_2013biomass2, #fit_2014biomass2, fit_ht2, fit_dia2)
#mv_trait_list <- c("Capsules_per_stem", "Capsules_per_plot", #"Stems_per_plot", "2013_Biomass_per_plot","2014_Biomass_per_plot", #"Height", "Plant_diameter")
#mv_blup_df <- data.frame(matrix(ncol = length(mv_trait_list), nrow = 33))
#names(mv_blup_df) <- mv_trait_list
#for (i in 1:length(mv_fit_list2) ){
#  fit <- mv_fit_list2[[i]]
#  blups <- coef(fit)$population
#  mv_blup_df[i] <- blups
#}
#rownames(mv_blup_df) <- rownames(blups)
#
##rownames(mv_blup_df) <- mv_blup_df$population
##mv_blup_df <- mv_blup_df[-1]
#
##write.csv(mv_pop_trait_means, file="data/mv_pop_trait_means.csv", row.names = FALSE)

#### Trait PCA of accession level EMMs ####
#mv_trait_pca <- rda(mv_blup_df, scale = T) #scale and center everything bc they are vastly different units.
#summary(mv_trait_pca)

#mv_trait_pca_noAppar <- rda(mv_blup_df[-2,], scale = T)
mv_trait_pca_noAppar2 <- rda(mv_means_df2[-2,], scale = T) #pca with emmeans instead of blups. exlcude Appar (row 2)
#scaled_mv_blups <- scale(mv_blup_df, center = T, scale = T)
#mv_trait_pca2 <- princomp(scaled_mv_blups) #just making sure the scaling function within the rda function is working as expected (centered and scaled)
#summary(mv_trait_pca2)
#biplot(mv_trait_pca)
#biplot(mv_trait_pca_noAppar)
#biplot(mv_trait_pca_noAppar2)

# Species (traits) loadings for first 3 PCs
mv_trait_PC_loadings <- round(data.frame(scores(mv_trait_pca_noAppar2, choices=1:3, display = "species", scaling = 0)), digits=3) %>%
  arrange(desc(abs(PC1)))
mv_trait_PC_loading_cutoff <- sqrt(1/ncol(mv_means_df2)) #.378 is loading of a single variable if each variable contributed equally; sum of squares of all loadings for an individual principal components must sum to 1. variables greater than this contribute substantially to a particular PC.
mv_PCA_eigenvals <- round(summary(eigenvals(mv_trait_pca_noAppar2))[,1:3], digits = 3)

write.csv(rbind(mv_trait_PC_loadings, mv_PCA_eigenvals), "plots/millville_trait_PC_loadings_eigenvals.csv")


mv_fort_pops <- fortify(mv_trait_pca_noAppar2, display='sites', scaling=0)
mv_fort_traits <- fortify(mv_trait_pca_noAppar2, display='species', scaling=0)

# PC1 v PC2
millville_trait_pca_plot <- ggplot() +
  geom_point(data=mv_fort_pops, aes(x = PC1, y = PC2), size=2, alpha=0.5) +
  geom_text(data=mv_fort_pops, aes(x = PC1, y = PC2, label=Label), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=mv_fort_pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=mv_fort_traits, aes(x=0, xend=PC1, y=0, yend=PC2), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.02,"npc"))) +
  geom_text_repel(data=mv_fort_traits, 
                  aes(x=PC1,y=PC2,label=Label), 
                  color="red", size=4) +
  labs(x=paste0("PC1 ","(",100*mv_PCA_eigenvals[2,1],"%)"), y=paste0("PC2 ", "(",100*mv_PCA_eigenvals[2,2],"%)"), title="Millville") +
  theme_bw() +
  theme(text = element_text(size = 14))

#png("plots/millville_traits_PCA.png", width=9, height=9, res=300, units="in")
#millville_trait_pca_plot
#dev.off()

# PC2 v PC3
millville_trait_pca_plot2 <- ggplot() +
  geom_point(data=mv_fort_pops, aes(x = PC2, y = PC3), size=2, alpha=0.5) +
  geom_text(data=mv_fort_pops, aes(x = PC2, y = PC3, label=Label), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=mv_fort_pops, aes(x = PC1, y = PC2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=mv_fort_traits, aes(x=0, xend=PC2, y=0, yend=PC3), 
               color="red", alpha=0.75, arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text_repel(data=mv_fort_traits, 
                  aes(x=PC2,y=PC3,label=Label), 
                  color="red", size=4) +
  labs(x=paste0("PC2 ","(",100*mv_PCA_eigenvals[2,2],"%)"), y=paste0("PC3 ", "(",100*mv_PCA_eigenvals[2,3],"%)"), title="Millville") +
  theme_bw() +
  theme(text = element_text(size = 14))

# ggvegan version
#autoplot(mv_trait_pca_noAppar2, arrows = TRUE, geom = "text", legend = "none") #basic

# Join Millville and Ephraim plots together
# library(cowplot)
fig2 <- plot_grid(millville_trait_pca_plot, ephraim_trait_pca_plot, labels=c("a)","b)"), ncol=1, nrow=2)
jpeg("plots/Fig2.jpg", width=17, height = 23, res=600, units = "cm")
fig2
dev.off()

#### Full RDA (all env predictors and all traits) of accession-level data ####
mv_means_df2$population <- rownames(mv_means_df2)
#fullRDA_df <- inner_join(mv_pop_trait_means, geo_clim_df) %>%
#  dplyr::select(-c(Entry,source))

mv_fullRDA_df <- inner_join(mv_means_df2, geo_clim_df)
rownames(mv_fullRDA_df) <- mv_fullRDA_df$population
mv_fullRDA_df <- dplyr::select(mv_fullRDA_df, -c(population, source))
mv_RDA_traits <- mv_fullRDA_df[1:7]
mv_RDA_preds <- mv_fullRDA_df[8:29]

mv_full_rda <- rda(mv_RDA_traits ~ ., data=mv_RDA_preds, scale = T)


summary(mv_full_rda)
RsquareAdj(mv_full_rda)

mv_rda_trait_loadings <- round(data.frame(scores(mv_full_rda, choices=1:4, display = "species", scaling = 0)), digits=3) %>%
  arrange(desc(abs(RDA1)))
mv_rda_trait_loading_cutoff <- sqrt(1/7) #7 traits. 0.3779645

# how to get env_loadings?? not sure 'bp' is correct here. 
mv_rda_env_loadings <- round(data.frame(scores(mv_full_rda, choices=1:4, display = "bp", scaling = 0)), digits=3) %>% 
  arrange(desc(abs(RDA1)))
#rda_env_loading_cutoff <- sqrt(1/22) #=.213 this is not how to assess importance of env loadings.
write.csv(mv_rda_env_loadings, "plots/Millville_rda_env_loadings.csv")

mv_rda_eigenvals <- round(summary(eigenvals(mv_full_rda, model = "constrained"))[,1:4], digits = 3) #constrained by climate
mv_rda_eigenvals_adj <- round(rbind(mv_rda_eigenvals["Eigenvalue",], data.frame(mv_rda_eigenvals[2:3,]) * RsquareAdj(mv_full_rda)[2]), digits = 3) #loadings adjusted by adjusted R squared.
rownames(mv_rda_eigenvals_adj)[1] <- "Eigenvalue"
write.csv(rbind(mv_rda_trait_loadings, mv_rda_eigenvals_adj), "plots/millville_rda_loadings_eigenvals.csv")

# Plotting
mv_rda.sp_sc0 <- scores(mv_full_rda, choices = 1:2, scaling=0, display="sp") #scaling 0
mv_rda.sp_sc1 <- scores(mv_full_rda, choices = 1:2, scaling=1, display="sp") #scaling 1
mv_rda.sp_sc <- data.frame(scores(mv_full_rda, choices = 1:2, scaling = 0, display="sp"))
mv_rda.env_sc <- data.frame(scores(mv_full_rda, choices = 1:2, scaling = 0, display = "bp"))
#mv_mul <- ordiArrowMul(scores(mv_full_rda, choices = 1:2, scaling = 0, display = "bp")) #multiplier for the coordinates of the head of the env vectors such that they fill set proportion of the plot region. This function is used in the default plot() function for rda objects in scaling=2.
#mv_rda.env_sc <- mv_rda.env_sc*mv_mul
mv_rda.site_sc <- data.frame(scores(mv_full_rda, choices = 1:2, scaling = 0, display = "wa"))


# fortify() from ggvegan converting rda scores into dataframes, same thing as what i've already done above.
#ggveg_rda_pops <- fortify(mv_full_rda, display="sites")
#ggveg_rda_traits <- fortify(mv_full_rda, display="species")
#ggveg_rda_env <- fortify(mv_full_rda, display="bp")

mv_rda_triplotgg <- ggplot() +
  geom_point(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2), size=2, alpha=0.5) +
  #geom_text(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2, label=rownames(mv_2da.site_sc)), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=mv_full_rda, aes(x = RDA1, y = RDA2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=mv_rda.sp_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(.02, "npc"))) +
  geom_text_repel(data=mv_rda.sp_sc, 
                  aes(x=RDA1,y=RDA2,label=rownames(mv_rda.sp_sc)), 
                  color="red", size=4) +
  geom_segment(data=mv_rda.env_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", alpha=0.5, arrow=arrow(length=unit(.02,"npc"))) +
  #geom_text_repel(data=mv_rda.env_sc, 
                  #aes(x=RDA1,y=RDA2,label=rownames(mv_rda.env_sc)), 
                  #color="blue", size=4) +
  annotate("text", x = -.5, y = -.19, label = "Lat", color='blue') +
  annotate("text", x = .3, y = .54, label = "Temperature", color='blue') +
  annotate("text", x = -.09, y = -.4, label = "Precipitation", color='blue') +
  annotate("text", x = -.1, y = -.56, label = "Long", color='blue') +
  annotate("text", x = .5, y = -.21, label = "MDR", color='blue') +
  labs(x=paste0("RDA1 ","(",100*mv_rda_eigenvals_adj[2,1],"%)"), y=paste0("RDA2 ", "(",100*mv_rda_eigenvals_adj[2,2],"%)"), title="Millville") +
  theme_bw() +
  theme(text = element_text(size = 14))
  
mv_rda_triplotgg

# Combine with Ephraim RDA plot (code for Ephraim plot is in separate file)
fig3 <- plot_grid(mv_rda_triplotgg, eph_rda_triplotgg, labels=c("a)","b)"), ncol=1, nrow=2)
#
jpeg("plots/fig3.jpg", width=17, height=23, res=600, units="cm")
fig3
dev.off()

# Plot for Supp Mat with all predictor arrows labeled.
mv_rda_triplotgg_SUPP <- ggplot() +
  geom_point(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2), size=2, alpha=0.5) +
  #geom_text(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2, label=rownames(mv_2da.site_sc)), hjust=0, vjust=0, size=4, alpha=.5) +
  #geom_text_repel(data=mv_full_rda, aes(x = RDA1, y = RDA2, label=Label), size=4, alpha=0.5, max.overlaps = 11) +
  geom_segment(data=mv_rda.sp_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(.02, "npc"))) +
  geom_text_repel(data=mv_rda.sp_sc, 
                  aes(x=RDA1,y=RDA2,label=rownames(mv_rda.sp_sc)), 
                  color="red", size=4) +
  geom_segment(data=mv_rda.env_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", alpha=0.5, arrow=arrow(length=unit(.02,"npc"))) +
  geom_text_repel(data=mv_rda.env_sc, 
  aes(x=RDA1,y=RDA2,label=rownames(mv_rda.env_sc)), 
  color="blue", alpha=0.5, size=4) +
  labs(x=paste0("RDA1 ","(",100*mv_rda_eigenvals_adj[2,1],"%)"), y=paste0("RDA2 ", "(",100*mv_rda_eigenvals_adj[2,2],"%)"), title="Millville") +
  theme_bw() +
  theme(text = element_text(size = 14))


figS3 <- plot_grid(mv_rda_triplotgg_SUPP, eph_rda_triplotgg_SUPP, labels=c("a)","b)"), ncol=1, nrow=2)
jpeg("plots/figS3.jpg", width=17, height=23, res=600, units="cm")
figS3
dev.off()

# Alt/base RDA plots for comparison. Custom plot above should have same arrow positions as the base plot() function. The autoplot() function from ggvegan is slightly different--different scaling I think?
plot(mv_full_rda, display = c("sp", "wa", "bp"))
arrows(0,0,mv_rda.sp_sc[,1], mv_rda.sp_sc[,2], length=0, lty=8, col="red")
autoplot(mv_full_rda, geom="text", layers = c("species", "sites", "biplot"), legend= "none", scaling=2) 

mv_rda_triplot <- autoplot(mv_full_rda, geom="text", layers = c("species", "sites", "biplot"), legend= "none", scaling=2) 

# Plot with fewer predictor arrows, for visualization purposes (still based on full model)
# *NEEDS WORK*. remove some collinear env predictor arrows simply for visualization for main body figure. 
#keep.mv_rda.env <- c("Lat", "Long", "bio01", "bio02", "bio08", "bio17", "bio18")
#mv_rda.env_sc.trimmed <- subset(mv_rda.env_sc, rownames(mv_rda.env_sc) %in% keep.mv_rda.env) 
#
#mv_rda_triplotgg_TRIMMED <- ggplot() +
#  geom_point(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2), size=2, alpha=0.5) +
#  #geom_text(data=mv_rda.site_sc, aes(x = RDA1, y = RDA2, label=rownames(mv_rda.site_sc)), #hjust=0, vjust=0, size=4, alpha=.5) +
#  #geom_text_repel(data=mv_full_rda, aes(x = RDA1, y = RDA2, label=Label), size=4, alpha=0.5, #max.overlaps = 11) +
#  geom_segment(data=mv_rda.sp_sc, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
#               color="red", arrow=arrow(length=unit(.02, "npc"))) +
#  geom_text_repel(data=mv_rda.sp_sc, 
#                  aes(x=RDA1,y=RDA2,label=rownames(mv_rda.sp_sc)), 
#                  color="red", size=4) +
#  geom_segment(data=mv_rda.env_sc.trimmed, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
#               color="blue", arrow=arrow(length=unit(.02,"npc"))) +
#  geom_text_repel(data=mv_rda.env_sc.trimmed, 
#                  aes(x=RDA1,y=RDA2,label=rownames(mv_rda.env_sc.trimmed)), 
#                  color="blue", size=4) +
#  labs(x="RDA1", y="RDA2", title="Millville") +
#  theme_bw() +
#  theme(text = element_text(size = 14))
#
## Save plot on its own
#png("plots/millville_RDA.png", width=9, height=9, res=300, units="in")
#mv_rda_triplotgg
#dev.off()
#
#fig3_trimmed <- plot_grid(mv_rda_triplotgg_TRIMMED, eph_rda_triplotgg_TRIMMED, labels=c("a)","b#)"), ncol=1, nrow=2)
#
#jpeg("plots/fig3_TRIMMED.jpg", width=17, height=23, res=600, units="cm")
#fig3_trimmed
#dev.off()



# tests of significance, variance partition, and variable selection--something we should do? Not doing variable selection.
set.seed(9)
# Global test of the RDA result
anova.cca(mv_full_rda, permutations=1000)
# Marginal effects of the terms (each marginal term analysed in a model with all other variables). 
anova.cca(mv_full_rda, by="margin", permutations=1000)
anova.cca(mv_full_rda, by="terms", permutations=1000) 
anova.cca(mv_full_rda, by="axis", permutations=1000)

# Variance partition
mv_RDA_geog <- mv_RDA_preds[1:3]
mv_RDA_clim <- mv_RDA_preds[4:22]
RsquareAdj(mv_full_rda)
varpart(mv_RDA_traits, mv_RDA_geog, mv_RDA_clim, scale = T, permutations = 1000)
showvarparts(2)

# Check for collinearity of explanatory variables with Variance Inflation Factors (VIF)?
#vif.cca(mv_full_rda)

# Variable selection with vegan's ordistep() and packfor's forward.sel()
#mv_mod0 <- rda(mv_RDA_traits ~ 1, data=mv_RDA_preds, scale = T)
#mv_mod1 <- rda(mv_RDA_traits ~ ., data=mv_RDA_preds, scale = T)
#set.seed(10)
#step.forward <- ordistep(mv_mod0, scope=formula(mv_mod1), direction="forward", pstep=1000, Pin = 0.1)

#R2step <- ordiR2step(mv_mod0, scope=formula(mv_mod1), direction="both", Pin = 0.05, #R2permutations = 1000, R2scope = T, permutations = how(nperm = 499)) 

#library(packfor)
#R2a.full <- RsquareAdj(mv_full_rda)$adj.r.squared
#forward.sel(mv_RDA_traits, mv_RDA_preds, adjR2thresh =  R2a.full) 

# Parsimonious model using variables selected above?
#mv_pars_rda <- rda(mv_RDA_traits ~ , data=mv_RDA_preds, scale = T)
#autoplot(mv_pars_rda, arrows=FALSE, geom="text", legend= "none")
#summary(mv_pars_rda)
#RsquareAdj(mv_pars_rda)
#vif.cca(mv_pars_rda)
#anova.cca(mv_pars_rda, by="axis", permutations=1000)



#### Trait Pearson correlations ####
# Use scaled/centered population means as used in PCA?
# 
library(ggcorrplot)
mv_trait_corr_df <- mv_means_df2[-2,] #exclude Appar bc different species

scaled_mv_trait_corr <- scale(mv_trait_corr_df, center = T, scale = T)

mv_corr_mat <- round(cor(mv_trait_corr_df, method=c("pearson"), use = "complete.obs"),4)

mv_p_mat <- cor_pmat(mv_trait_corr_df)
head(mv_p_mat)
quartz()
mv_corr_plot <- ggcorrplot(mv_corr_mat, hc.order = TRUE,type = "lower", lab = TRUE, p.mat=mv_p_mat, insig = "blank", lab_size = 4, tl.cex = 10, show.legend = FALSE, title = "Millville") + #Using default sig level of .05
  theme(text = element_text(size = 14))

jpeg("plots/millville_trait_corrplot.jpg", width=17, height=10, res=600, units="cm")
mv_corr_plot
dev.off()

both_corr_plots <- plot_grid(mv_corr_plot, eph_corr_plot, ncol=1, nrow=2)

jpeg("plots/corr_plots.jpg", width=17, height=23, res=600, units="cm")
both_corr_plots
dev.off()

#### Latitudinal clines? ####
#fit_list <- c()
#mv_trait_list <- c("Capsules_per_stem", "Capsules_per_plot", "Stems_per_plot", #"Biomass_per_plot",  "Height", "Plant_diameter")
#datasets <- list(mv_stems_caps, mv_stems_caps, mv_stems_caps, mv_biomass, ht_data, #dia_data)
#plot_list = list()
#
## Loop through each model, calculating population means and confidence intervals of #regression lines
#for (i in 1:length(fit_list)) {
#  fit <- fit_list[[i]]
#  data <- datasets[[i]] %>%
#    inner_join(dplyr::select(TvS_key, Entry, source)) %>%
#    inner_join(dplyr::select(env_data, source, Lat))
#  
#  pred_df <- make_pred_df(fit) #get population means
#  
#  # Obtain confidence interval for regression line
#  newd <- data.frame(Lat = seq(min(geo_data$Lat, na.rm=T), max(geo_data$Lat, na.rm=T), #length.out=100))
#  lmm_boots <- bootMer(fit, predict_fun, nsim = 100)
#  pred_ci <- cbind(newd, confint(lmm_boots))
#  
#  # Plot population means vs latitude
#  plot <- ggplot(data=pred_df) +
#    geom_abline(intercept=fixef(fit)[1], slope=fixef(fit)[2], lty=2) +
#    geom_point(mapping=aes(x=Lat, y=pop_b0), color="royalblue2", alpha=0.5) +
#    geom_linerange(mapping=aes(x=Lat, ymin=pop_b0-pop_b0_se,ymax=pop_b0+pop_b0_se), color#="royalblue2", alpha=0.5) +
#    geom_ribbon(data=pred_ci, aes(x=Lat, ymin=`2.5 %`, ymax=`97.5 %`), alpha=0.25) +
#    labs(x="Lat", y=paste(mv_trait_list[[i]])) +
#    theme_minimal()
#  
#  plot_list[[i]] <- plot
#}
#p <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#ggdraw(add_sub(p, "Latitude", vpadding=grid::unit(0,"lines"),y=6, x=0.53, vjust=4.5))
