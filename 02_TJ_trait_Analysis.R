# Analyzing Tom Jones' data
# 1.30.21
# Peter Innes

# Stem and capsule data
tj_stems <- read.csv("data/TomJ_2013_stems_caps_ht_diam_surv_ANALYZE_THIS.csv", skip=1, na.strings = ".", header = T) #Spreadsheet has two headers >:( so we gotta skip the first one
head(tj_stems) 

# Biomass data
tj_biomass <- read.csv("data/TomJ_BiomassFlax2013_14_15_ANALYZE_THIS.csv", skip=1, header = T)
names(tj_biomass) <- c("Plot", "Rep", "Entry", "num_planted", "ttl_weight_2013", "ttl_stem_2013", "sub_weight_2013", "ttl_weight_2014", "DW_2015", "survivorship_4_27_13")

# Survival data
tj_surviv <- read.csv("raw_data/TomJ_linum_lewisii_survivorship_flwr_full (ANALYZED IN FILE).csv", skip=1, header = T)
tj_surviv <- tj_surviv %>% dplyr::select(Plot, Rep, Entry, X6.12.12, X10.27.12, X4.27.13, X8.30.13)
names(tj_surviv) <- c("Plot", "Rep", "Entry", "planted_6_12_12", "surviv_10_27_12", "surviv_4_27_13", "surviv_8_30_13") #clean up and keep only the columns we need 
# Convert "." to 0
na_code <- "."
tj_surviv [4:7] <- apply(tj_surviv[4:7], 2, function(x){ as.numeric(ifelse(x %in% na_code, 0, x)) } )
tj_surviv <- head(tj_surviv, -4) #delete last 4 rows

# Stems/caps model
names(tj_stems)[11] <- "no_survival" #Rename last column since we lost the column name when initially skipping first line. This column has number of plants surviving on 4.27.13 (out of 10 plants initially planted I believe).
head(arrange(tj_stems, Entry, Rep), 20L)
tj_stems$Entry <- as.factor(tj_stems$Entry)

fit_caps_per_stem <- lmer(log(Capsules/sub_stems) ~ -1 + Entry + (1|Rep), data = tj_stems) 
exp(ls_means(fit_caps_per_stem)$Estimate)
plot(fit_caps_per_stem)
qqnorm(resid(fit_caps_per_stem))

# How to estimate capsules per plant? est_ttl_capules / surviv_8_30_13?

# Surviv model
tj_surviv$perc_surviv_10_27_12 <-  tj_surviv$surviv_10_27_12 / tj_surviv$planted_6_12_12 
tj_surviv$perc_surviv_4_27_13 <-  tj_surviv$surviv_4_27_13 / tj_surviv$planted_6_12_12 

fit_surviv_10_27_12 <- lmer(perc_surviv_10_27_12 ~ -1 + Entry + (1|Rep), data=tj_surviv)
ls_means(fit_surviv)
summary(fit_surviv)

fit_surviv_4_27_13 <- lmer(tj_surviv$perc_surviv_4_27_13 ~ -1 + Entry + (1|Rep), data=tj_surviv)
summary(fit_surviv_4_27_13)

# Height data
ht_data <- read.csv("data/TomJ_flax_avg_ht_dia_2013_RAW_DATA.csv", header = T, na.strings = '.') %>% mutate(Plot=as.factor(Plot), Rep=as.factor(Rep), Entry=as.factor(Entry))
ht_data <- ht_data %>% dplyr::select(Plot, Rep, Entry, ht_pl1, ht_pl2, ht_pl3, ht_pl4, ht_pl5)
ht_data <- gather(ht_data, key = "Plant", value = "height", ht_pl1:ht_pl5)

head(ht_data)
str(ht_data)
fit_ht <- lmer(height ~ -1 + Entry + (1|Rep) + (1|Plot), data = ht_data)
summary(fit_ht)
