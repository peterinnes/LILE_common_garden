# Analyzing Tom Jones' data
# 1.30.21
# Peter Innes

library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)

# To do: filter out Appar sources, convert Tom source number to Stan source number/population ID.

TvS_key <- read.csv("info_from_TomJ/tomJ_data/StanVsTomAccessions.csv", header=T) %>%
  mutate(source=factor(source), Entry=factor(Entry))

#### Stem and capsule data (Capsules per stem, capsules per plant, stems per plant) ####
tj_stems_caps <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) %>% #Spreadsheet has two headers. so we gotta skip the first one
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))

# summarise. Entry 5 has 10 reps when each entry should only have 8. Data mis-entered?
sc_summary <- tj_stems_caps %>%
  group_by(Entry) %>%
  summarise(n=n())

# Fix data entry error, change plots 134 and 265 from Entry 5 to Entry 6
tj_stems_caps$Entry[134] <- 6
tj_stems_caps$Entry[265] <- 6
Plot_Entry_key <- dplyr::select(tj_stems_caps, Plot, Entry)

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
exp(ls_means(fit_CPS)$Estimate)
summary(fit_CPS)
plot(fit_CPS)
qqnorm(resid(fit_CPS))

# Caps per plant: Capsules/sub_stems * ttl_stems = est_ttl_capsules. Then divide est_ttl_capules by the number of surviving plants in the plot. Using the survival data from the tj_stems_caps spreadsheet--this closely matches survival counts taken later in August, closer to Harvest (but this Aug data is incomplete so we don't want to use it directly)
# square root transform, this is what Tom originally used. It looks good.
tj_stems_caps$CPP_tr <- (sqrt(tj_stems_caps$est_ttl_capsules/tj_stems_caps$surv_4_27_13) - 1) / .5
fit_CPP <- lmer(CPP_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_CPP)
qqnorm(resid(fit_CPP))
summary(fit_CPP)

# Stems per plant. Use same sqrt transform as for capsules per plant
tj_stems_caps$spp_tr <- (sqrt(tj_stems_caps$ttl_stems/tj_stems_caps$surv_4_27_13) - 1) / .5
fit_SPP <- lmer(spp_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)

plot(fit_SPP)
qqnorm(resid(fit_SPP))
summary(fit_SPP) # How to estimate capsules per plant? est_ttl_capules / surviv_8_30_13?


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

# Join with source info, filter out mistaken Appar accessions.
tj_biomass <- tj_biomass %>%
  inner_join(dplyr::select(TvS_key, Entry, source) %>% na.omit()) %>%
  relocate(source, .after = Entry) %>%
  filter(!source %in% c(2,5,22,32,38))

dim(tj_biomass)

bm_summary2 <- full_join(tj_biomass, dplyr::select(tj_stems_caps, Plot, entry = Entry)) %>%
  group_by(entry) %>%
  summarise(n=n())

# Linear models
# Biomass per plant
tj_biomass$BPP <- tj_biomass$ttl_weight_2013 / tj_biomass$survivorship_4_27_13 
tj_biomass$BPP[tj_biomass$BPP == 0] <- NA #convert Zeros in BPP to NA
tj_biomass$BPP_tr <- (sqrt(tj_biomass$BPP) - 1) / 0.5 #square root transform

fit_biomass <- lmer(BPP_tr ~ -1 + Entry + (1|Rep), data=tj_biomass)
plot(fit_biomass)
qqnorm(resid(fit_biomass))

#### Canopy height and plant diameter data ####
ht_dia_data <- read.csv("data/TomJ_flax_avg_ht_dia_2013_RAW_DATA.csv", header = T, na.strings = '.') %>%
  mutate(Plot=as.factor(Plot), Rep=as.factor(Rep), Entry=as.factor(Entry)) 
ht_dia_data <- head(ht_dia_data, -3) #delete 3 empty rows at end

View(ht_dia_data %>%
       group_by(Entry) %>%
       summarise(n=n()))

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
summary(fit_ht) #adding random effect of plot because we have measurements of multiple individual plants per plot. Rep (block) random effect is zero so we will remove that re.
plot(fit_ht)

fit_dia <- lmer(dia ~ -1 + Entry + (1|Rep) + (1|Plot), data = dia_data)
summary(fit_dia)
plot(fit_dia)
qqnorm(resid(fit_dia))

#### Survival data ####
# not sure how reliable this data is. 
tj_surviv <- read.csv("data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T) %>%
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep))
# rename columns and keep only those we need
tj_surviv <- tj_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13, X5.1.14, X8.12.14)
names(tj_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13", "surviv_5_1_14", "surviv_8_12_14") 
tj_surviv <- head(tj_surviv, -1) #delete last row which had totals

# Entry 5 also has 10 reps here in this dataset, need to correct.
View(tj_surviv %>%
       group_by(Entry) %>%
       summarise(n=n()))

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
surv <- full_join(dplyr::select(tj_surviv, Plot, surviv_4_27_13, surviv_8_30_13), dia_plants_per_plot) %>%
  full_join(dplyr::select(tj_stems_caps, Plot, surv_4_27_13)) %>% #survival from stem/caps spreadsheet
  rename(caps_surv_4_27_13 = surv_4_27_13) %>% 
  full_join(dplyr::select(tj_biomass, Plot, survivorship_4_27_13)) %>%
  rename(biomass_surv_4_27_13 = survivorship_4_27_13)
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
# leaving out fit_surviv for now.
tj_fit_list <- c(fit_CPS, fit_CPP, fit_SPP, fit_biomass, fit_ht, fit_dia)
tj_trait_list <- c("Capsules_per_stem", "Capsules_per_plant", "Stems_per_plant", "Biomass_per_plant", "Plant_height", "Plant_diameter")
tj_results <- list() #list to store means and confidence intervals
tj_esp_list <- list() #list to store effect size plots

# confidence intervals taken from the lsmeans() fxn here bc it is fast. will pronbably want to use confint() for final pub instead.
for (i in 1:length(tj_fit_list) ){
  fit <- tj_fit_list[[i]]
  
  lsmeans <- as.data.frame(lsmeans(fit, "Entry"))
  #names(means)[2] <- "trait_value"
  
  # Create effect size plot
  esp <- ggplot(data=lsmeans, aes(x=lsmean, y=reorder(Entry, lsmean), xmin=lower.CL, xmax=upper.CL)) +
    geom_point() +
    geom_errorbar() +
    ylab("Population") +
    xlab(tj_trait_list[[i]]) +
    theme(axis.text.y = element_text(size = 6))
  
  tj_esp_list[[i]] <- esp #Store plot
  
  names(lsmeans)[2] <- tj_trait_list[[i]] #Change to actual trait name before storing in results
  tj_results[[i]] <- lsmeans #store means and confidence intervals
  
  # Tweak dataframe and write to csv for summary to send to Scott J et al
  lsmeans <- lsmeans %>% arrange(-lsmeans[2]) #sort descending trait value to make more readable
  #write.csv(means_ci, file=paste0("results_summaries/", names(means_ci)[4], "_summary.csv"))
}
names(tj_results) <- tj_trait_list
head(tj_results)

#### Trait correlations & PCA ####
library(vegan)
library(ggvegan)

tj_pop_trait_means <- data.frame(Entry=tj_results[[1]]$Entry)
for ( i in 1:length(tj_results) ){
  tj_pop_trait_means <- cbind(tj_pop_trait_means, tj_results[[i]][2]) #add trait means to growing df
}
#write.csv(tj_pop_trait_means, file="data/tj_pop_trait_means.csv", row.names = FALSE)

# PCA with source as rownames/labels
tj_pop_trait_means <- tj_pop_trait_means %>% 
  inner_join(dplyr::select(TvS_key, Entry, source)) %>%
  relocate(source, .after = Entry)

temp <- tj_pop_trait_means[,-1]
rownames(temp) <- temp[,1]
temp <- temp[,-1]

tj_trait_rda <- rda(temp, scale = T) #scale everything bc they are vastly different units.
biplot(tj_trait_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
ordilabel(tj_trait_rda, dis="sites", cex=0.5)

# ggvegan version
autoplot(tj_trait_rda, arrows = TRUE, geom = "text", legend = "none") #basic

# Trait correlations. Need to back transform the lsmeans in order to use pearson correlations. sqr back transform = ((tr / 2) + 1)^2.
