# Analyzing Tom Jones' data
# 1.30.21
# Peter Innes

library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)

#### Stem and capsule data ####
tj_stems_caps <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) #Spreadsheet has two headers. so we gotta skip the first one
head(tj_stems_caps) 

names(tj_stems_caps)[11] <- "num_surv" #Rename last column since we lost the column name when initially skipping first line. This column has number of plants surviving on 4.27.13 (out of 10 plants initially planted I believe).
head(arrange(tj_stems_caps, Entry, Rep), 20L)
tj_stems_caps$Entry <- as.factor(tj_stems_caps$Entry)

# Caps per stem
fit_CPS <- lmer(log(Capsules/sub_stems) ~ -1 + Entry + (1|Rep), data = tj_stems_caps) 
exp(ls_means(fit_caps_per_stem)$Estimate)
plot(fit_caps_per_stem)
qqnorm(resid(fit_caps_per_stem))

# Caps per plant
# square root transform. THis is the transform Tom used
tj_stems_caps$cpp_tr <- (sqrt(tj_stems_caps$est_ttl_capsules/tj_stems_caps$num_surv) - 1) / .5
fit_CPP <- lmer(cpp_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)
plot(fit_CPP)
qqnorm(resid(fit_CPP))
summary(fit_CPP)

# Stems per plant. Use same sqrt transform as before
tj_stems_caps$spp_tr <- (sqrt(tj_stems_caps$ttl_stems/tj_stems_caps$num_surv) - 1) / .5
fit_SPP <- lmer(spp_tr ~ -1 + Entry + (1|Rep), data = tj_stems_caps)

plot(fit_SPP)
qqnorm(resid(fit_SPP))
summary(fit_SPP) # How to estimate capsules per plant? est_ttl_capules / surviv_8_30_13?

#### Biomass data ####
tj_biomass <- read.csv("data/TomJ_BiomassFlax2013_14_15_ANALYZE_THIS.csv", skip=1, header = T)
names(tj_biomass) <- c("Plot", "Rep", "Entry", "num_planted", "ttl_weight_2013", "ttl_stem_2013", "sub_weight_2013", "ttl_weight_2014", "DW_2015", "survivorship_4_27_13")
head(tj_biomass)

#### Survival data ####
tj_surviv <- read.csv("raw_data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T)
tj_surviv <- tj_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13, X5.1.14, X8.12.14)
names(tj_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13", "surviv_5_1_14", "surviv_8_12_14") #clean up and keep only the columns we need
tj_surviv <- head(tj_surviv, -4) #delete last 4 empty rows

# Convert "." to 0
na_code <- "."
tj_surviv [4:9] <- apply(tj_surviv[4:9], 2, function(x){ as.numeric(ifelse(x %in% na_code, 0, x)) } )

# Change variable types. For some reason survival columns were read-in as characters
tj_surviv[4:9] <- as.integer(unlist(tj_surviv[4:9]))
tj_surviv[1:3] <- as.factor(unlist(tj_surviv[1:3]))

# First year survival
fit_surviv <- lmer(surviv_4_27_13/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv)

# Second year survival
fit_surviv2 <- lmer(surviv_5_1_14/planted_6_12_12 ~ -1 + Entry + (1|Rep), data = tj_surviv)
summary(fit_surviv2)


#### Height data ####
ht_data <- read.csv("data/TomJ_flax_avg_ht_dia_2013_RAW_DATA.csv", header = T, na.strings = '.') %>% mutate(Plot=as.factor(Plot), Rep=as.factor(Rep), Entry=as.factor(Entry))
ht_data <- ht_data %>% dplyr::select(Plot, Rep, Entry, ht_pl1, ht_pl2, ht_pl3, ht_pl4, ht_pl5)
ht_data <- gather(ht_data, key = "Plant", value = "height", ht_pl1:ht_pl5)

head(ht_data)
str(ht_data)
fit_ht <- lmer(height ~ -1 + Entry + (1|Rep) + (1|Plot), data = ht_data)
summary(fit_ht)

#### 
