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

#### Stem and capsule data ####
tj_stems_caps <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) %>% #Spreadsheet has two headers. so we gotta skip the first one
  mutate(Plot=factor(Plot), Entry=factor(Entry), Rep=factor(Rep)) %>%
  inner_join(dplyr::select(TvS_key, source, Entry)) %>%
  relocate(source, .before=Capsules) %>%
  inner_join(dplyr::select(env_data, source, site)) %>%
  relocate(site, .before=Capsules)
head(tj_stems_caps)

names(tj_stems_caps)[11] <- "surv_4_27_13" #Rename last column since we lost the column name when initially skipping first line. This column has number of plants surviving on 4.27.13 (out of 10 plants initially planted I believe).

# Need to merge survival data from Aug 2013 from other spreadsheet.
tj_stems_caps <- tj_stems_caps %>%
  left_join(dplyr::select(tj_surviv, Plot, Rep, Entry, surviv_8_30_13))

head(arrange(tj_stems_caps, Entry, Rep), 20L)
tj_stems_caps$Entry <- as.factor(tj_stems_caps$Entry)


# Caps per stem
fit_CPS <- lmer(log(Capsules/sub_stems) ~ -1 + Entry + (1|Rep), data = tj_stems_caps) 
exp(ls_means(fit_CPS)$Estimate)
summary(fit_CPS)
plot(fit_caps_per_stem)
qqnorm(resid(fit_caps_per_stem))

# Caps per plant: Capsules/sub_stems * ttl_stems = est_ttl_capsules. Then divide est_ttl_capules by the number of surviving plants in the plot. There might be an issue here with missing data from the 8.20.13 survival counts. 
# square root transform, this is what Tom originally used. It looks good.
tj_stems_caps$CPP_tr <- (sqrt(tj_stems_caps$est_ttl_capsules/tj_stems_caps$surviv_8_30_13) - 1) / .5
fit_CPP <- lmer(CPP_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_CPP)
qqnorm(resid(fit_CPP))
summary(fit_CPP)

# Stems per plant. Use same sqrt transform as before
tj_stems_caps$spp_tr <- (sqrt(tj_stems_caps$ttl_stems/tj_stems_caps$surviv_8_30_13) - 1) / .5
fit_SPP <- lmer(spp_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)

plot(fit_SPP)
qqnorm(resid(fit_SPP))
summary(fit_SPP) # How to estimate capsules per plant? est_ttl_capules / surviv_8_30_13?


#### Biomass data ####
tj_biomass <- read.csv("data/TomJ_BiomassFlax2013_14_15_ANALYZE_THIS.csv", skip=1, header = T)
names(tj_biomass) <- c("Plot", "Rep", "Entry", "num_planted", "ttl_weight_2013", "ttl_stem_2013", "sub_weight_2013", "ttl_weight_2014", "DW_2015", "survivorship_4_27_13")
head(tj_biomass)

tj_biomass[4:10] <- as.integer(unlist(tj_biomass[4:10]))
tj_biomass[1:3] <- as.factor(unlist(tj_biomass[1:3]))

tj_biomass$BPP <- tj_biomass$ttl_weight_2013 / tj_biomass$survivorship_4_27_13 #biomass per plant
tj_biomass$BPP_tr <- (sqrt(tj_biomass$BPP) - 1) / 0.5 #square root transform
fit_biomass <- lmer(BPP_tr ~ -1 + Entry + (1|Rep), data=tj_biomass)
plot(fit_biomass)
qqnorm(resid(fit_biomass))

#### Survival data ####
tj_surviv <- read.csv("data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T)

# rename columns and keep only those we need
tj_surviv <- tj_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13, X5.1.14, X8.12.14)
names(tj_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13", "surviv_5_1_14", "surviv_8_12_14") 

tj_surviv <- head(tj_surviv, -4) #delete last 4 empty rows

# Convert "." to 0 for column 5 and onwards.
na_code <- "."
tj_surviv [5:9] <- apply(tj_surviv[5:9], 2, function(x){ as.numeric(ifelse(x %in% na_code, 0, x)) } )

tj_surviv[4:9] <- as.integer(unlist(tj_surviv[4:9]))
tj_surviv[1:3] <- as.factor(unlist(tj_surviv[1:3]))

# First year survival
fit_surviv <- lmer(surviv_8_30_13/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv)
plot(fit_surviv)

# Second year survival
fit_surviv2 <- lmer(surviv_8_12_14/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv2)
plot(fit_surviv2)


#### Canopy height and plant diameter data ####
ht_dia_data <- read.csv("data/TomJ_flax_avg_ht_dia_2013_RAW_DATA.csv", header = T, na.strings = '.') %>% mutate(Plot=as.factor(Plot), Rep=as.factor(Rep), Entry=as.factor(Entry))
ht_data <- ht_dia_data %>% dplyr::select(Plot, Rep, Entry, ht_pl1, ht_pl2, ht_pl3, ht_pl4, ht_pl5)
ht_data <- pivot_longer(ht_data, names_to = "Plant", values_to = "height", ht_pl1:ht_pl5)


dia_data <- ht_dia_data %>% dplyr::select(Plot, Rep, Entry, dia.1:dia.10)
dia_data <- pivot_longer(dia_data, names_to = "Plant", values_to = "dia", dia.1:dia.10)

# How many plants in each plot? We will use this number to corroborate survivorship, in order to estimate per plant capsule and stem numbers. Diam and Height were reportedly measured a week before harvest in October. So they should be the most accurate?
plants_per_plot <- dia_data %>%
  na.omit() %>%
  group_by(Plot) %>%
  summarise(num_plants = n())

# Should also compare to April and Aug surv data. Could corroborate with caps data as well. 



fit_ht <- lmer(height ~ -1 + Entry + (1|Rep) + (1|Plot), data = ht_data)
summary(fit_ht)
plot(fit_ht)

fit_dia <- lmer(dia ~ -1 + Entry + (1|Rep) + (1|Plot), data = dia_data)
summary(fit_dia)
plot(fit_dia)



