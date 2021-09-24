# Oil traits analysis
# Peter Innes
# 5.15.21

library(dplyr)
library(ggplot2)
library(dichromat)
library(lme4)
library(lmerTest)
library(cowplot)
library(sjstats)

#### Oil composition ####

gc_data <- read.csv("data/GC_data/BrantGC/Brant_Flax_GC_Complete.csv", header = T, na.strings = c(".","n.a.")) %>%
  dplyr::select(Sample.Number, source=Accession, block=Plot, row=Range, BL, Palmitic, Palmitoleic, Stearic, Oleic, Linoleic, Alphalinolenic, Arachidic, Gondoic) 
gc_data <- gc_data %>% mutate(source=as.factor(source), block=as.factor(block), row=as.factor(row), BL=as.factor(BL)) %>% 
  full_join(dplyr::select(env_data, source, population, Lat, Long)) %>%
  filter(!source %in% c(2,5,22,32,38)) #exclude mistaken Appar collections
# multiply fatty acid percentages by 10 to convert to grams per kilogram fatty acid
gc_data[6:13] <- apply(gc_data[6:13], 2, function(x) x*10)

# quick EDA
ala_hist <- ggplot(data=gc_data, aes(x=Alphalinolenic), stat = "identity") +
  geom_histogram()

gc_data %>% 
  mutate(yjit=jitter(0*Alphalinolenic)) %>%
  ggplot() +
  geom_point(aes(x=Alphalinolenic, col=block, y=yjit)) +
  facet_wrap(facets = ~ source)

# linear models. 
fit_ala <- lmer(Alphalinolenic ~ population + (1|BL), data=gc_data)
#fit_ala2 <- lmer(Alphalinolenic ~ (1|population) + (1|block) + (1|population:block), data=gc_data)
#cbind(coef(fit_ala2)$population, emmeans(fit_ala, specs = "population")) 

fit_linoleic <- lmer(Linoleic ~ population + (1|BL), data=gc_data)
#fit_linoleic2 <- lmer(Linoleic ~ (1|population) + (1|block) + (1|population:block), data=gc_data)

fit_palmitic <- lmer(Palmitic ~ population + (1|BL), data=gc_data)
fit_stearic <- lmer(Stearic ~ population + (1|BL), data=gc_data)
fit_oleic <- lmer(Oleic ~ population + (1|BL), data=gc_data)

# 
#fit_oilComp_Geog <- lmer(Linoleic ~ Lat + Long + (1|population) + (1|block), data=gc_data)
#summary(fit_oilComp_Geog)
#anova(fit_oilComp_Geog)
#plot(fit_oilComp_Geog)
#coef(fit_oilComp_geog)$population

# get emms
ala_emm <- data.frame(emmeans(fit_ala, specs = "population")) %>% dplyr::select(population, Alphalinolenic=emmean)
linoleic_emm <- data.frame(emmeans(fit_linoleic, specs = "population")) %>% dplyr::select(population, Linoleic=emmean)
palmitic_emm <- data.frame(emmeans(fit_palmitic, specs="population")) %>% dplyr::select(population, Palmitic=emmean)
stearic_emm <- data.frame(emmeans(fit_stearic, specs="population")) %>% dplyr::select(population, Stearic=emmean)
oleic_emm <- data.frame(emmeans(fit_oleic, specs="population")) %>% dplyr::select(population, Oleic=emmean)

# plot emms
#ala_emm_plot <- plot(emmeans(fit_ala, specs = "population"), comparisons = T) +
  #theme_minimal() +
  #labs(x="% Alphalinolenic acid", y="Population")

# plot composition as stacked bar graph
oil_comp_df <- merge(ala_emm, linoleic_emm) %>%
  merge(palmitic_emm) %>%
  merge(stearic_emm) %>%
  merge(oleic_emm)

oil_comp_df_long <- pivot_longer(oil_comp_df, names_to = "Fatty_acid", values_to = "Composition", cols = 2:6 ) %>% inner_join(dplyr::select(env_data, population, Lat, Long))

# FIG4. stacked bar graph with accessions ordered by latitude. Composition is percentage by weight of total fatty acids
oil_comp_plot <- oil_comp_df_long %>% ggplot(aes(x=reorder(population, Lat), y=Composition, fill=Fatty_acid)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  labs(x="Accession", y="Fatty acid composition (g/kg)",fill="Fatty acid") +
  theme(text = element_text(size=14), legend.title = element_blank()) +
  scale_fill_manual(values = colorschemes$DarkRedtoBlue.12[c(1,3,5,9,11)])
oil_comp_plot

jpeg("plots/oil_composition_plot.png", width = 17, height = 10, res = 600, units = "cm")
oil_comp_plot
dev.off()


#### Oil Content ####

nmr_data <- read.csv("data/perennial_flax_paperNMR.csv", header=T) %>%
  dplyr::select(source=ID, Rep, Entry_number, Adj_Oil)
nmr_data <- nmr_data %>%
  mutate(source=as.factor(source), Rep=as.factor(Rep)) %>%
  full_join(dplyr::select(env_data, source, population)) %>%
  filter(!source %in% c(2,5,22,32,38)) %>%
  full_join(geo_data)
nmr_data$Adj_Oil <- nmr_data$Adj_Oil*10 #convert to grams per kilogram

nmr_hist <- ggplot(data=nmr_data, aes(x=Adj_Oil), stat = "identity") +
  geom_histogram()

fit_oil_cont <- lm(Adj_Oil ~ population, data=nmr_data)
#fit_oil_cont2 <- lmer(Adj_Oil ~ (1|population), data=nmr_data)
#fit_oilCont_Geog <- lmer(Adj_Oil ~ Long + (1|population), data=nmr_data)
#summary(fit_oilCont_Geog)

#plot(fit_oil_cont)

#emm_oil_cont <- lsmeans(fit_oil_cont, "population")
#oc_contrasts <- emmeans::emmeans(object=fit_oil_cont, pairwise ~ "population", adjust="tukey")
#oc_cld <- emmeans:::cld.emmGrid(object=contrasts$emmeans, Letters=letters)

#oc_cld_df <- data.frame(cld)
#oc_cld_df[c(2:6)] <- apply(cld_df[c(2:6)], 1:2, function(x) round(x, digits = 2))
#oc_cld_df$emm_letter <- apply(cld_df[c(2,7)], 1, paste, collapse="")

oil_cont_emm_plot <- plot(emm_oil_cont, type = "response", comparisons = T, colors = c("salmon", "blue", "black")) +
  theme_minimal() +
  xlab(tj_trait_list[i]) +
  ylab("")

# Get emms, clds, for table
#### Gather and plot ls-means of all traits. Confidence intervals taken from the lsmeans() aka emmeans() function of package emmeans ####
oil_fit_list <- list(fit_oil_cont, fit_ala, fit_linoleic, fit_oleic, fit_palmitic, fit_stearic)
oil_trait_list <- c("Oil_content", "Alphalinolenic", "Linoleic", "Oleic", "Palmitic", "Stearic")

oil_results <- list() #list to store means and confidence intervals
oil_results_cld <- list() ##list to store back-trasnformed means and confidence intervals 
#oil_esp_list <- list() #list to store effect size plots
oil_emm_list <- list() #list to store emm plot
#cld_list <- list() #list to store compact letter displays of each model

for (i in 1:length(oil_fit_list) ){
  
  fit <- oil_fit_list[[i]]
  lsmeans <- as.data.frame(lsmeans(fit, "population")) #as data frame
  emm1 <- lsmeans(fit, "population") #same as lsmeans above but dont convert to df
  emm2 <- as.data.frame(lsmeans(fit, "population", type="response")) #separate object for backtransformed lsmeans
  
  # Plotting means and CIs with emmeans:::plot.emmGrid. in order to reorder the populations on y-axis, we need to edit the .plot.srg function in emmeans package:
  # trace(emmeans:::.plot.srg, edit=T). Edit lines 239-240, change aes_() to aes() and delete tilde from in front of x and y variables. Then use reorder() on the y variable.
  emm_plot <- plot(emm1, type = "response", comparisons = T, colors = c("salmon", "blue", "black")) +
    theme_minimal() +
    xlab(oil_trait_list[i]) +
    ylab("")
  oil_emm_list[[i]] <- emm_plot #store plot
  
  # Compact letter display
  contrasts <- emmeans::emmeans(object=fit, type="response", pairwise ~ "population", adjust="tukey") #tests are on transformed scale but display on response scale
  cld <- emmeans:::cld.emmGrid(object=contrasts$emmeans, Letters=letters, sort=F)
  cld_df <- data.frame(cld)
  cld_df[c(2:6)] <- apply(cld_df[c(2:6)], 1:2, function(x) signif(x,3))
  cld_df$emm_letter <- apply(cld_df[c(2,7)], 1, paste, collapse="") #combine emmeans and letters into single column
  
  # Renaming columns and storing results
  names(lsmeans)[2] <- oil_trait_list[[i]] #Change 'lsmean' column name to trait name before storing in results
  oil_results[[i]] <- lsmeans #store means and confidence intervals
  lsmeans <- lsmeans %>% arrange(-lsmeans[2]) #sort descending trait value to make more readable
  # Same thing but with emms and clds together. This will be results table in manuscript.
  names(cld_df)[8] <- oil_trait_list[[i]] 
  oil_results_cld[[i]] <- cld_df
  #emm2 <- emm2 %>% arrange(-emm2[2]) 
}
names(oil_results) <- oil_trait_list
names(oil_results_cld) <- oil_trait_list

# store emms in one dataframe with population as rowname
oil_means_df <- data.frame(matrix(ncol = length(oil_trait_list), nrow = length(oil_results[[1]]$population)))
names(oil_means_df) <- oil_trait_list
rownames(oil_means_df) <- oil_results[[1]]$population
for (i in 1:length(oil_trait_list) ){
  oil_means_df[i]  <- oil_results[[i]][2]
}

# Store emms with clds in dataframe with column for population. this will be table in manuscript.
oil_means_df2 <- data.frame(matrix(ncol = length(oil_trait_list), nrow = length(oil_results[[1]]$population)))
names(oil_means_df2) <- oil_trait_list
rownames(oil_means_df2) <- oil_results_cld[[1]]$population
for (i in 1:length(oil_trait_list) ){
  oil_means_df2[i] <- oil_results_cld[[i]][8]
}
oil_means_df2 <- oil_means_df2 %>% arrange(desc(Oil_content)) %>%
  tibble::rownames_to_column("Accession") %>%
  relocate(Accession, .before = Oil_content)
write.csv(oil_means_df2, "plots/Oil_means_table.csv", row.names = F)

#### Cet coefficients of variation ####
oil_cvs <- c()
# Use the same models as above 
for (fit in oil_fit_list) {
  #resid_sd <- data.frame(VarCorr(fit))$sdcor[2]
  rmse <- performance::rmse(fit)
  grand_mean <- mean(insight::get_response(fit))
  cv <- rmse/grand_mean
  oil_cvs <- append(oil_cvs, cv)
}
names(oil_cvs) <- oil_trait_list

#oil_cv_fit_list <-  list(fit_ala, fit_linoleic, fit_oleic, fit_palmitic, fit_stearic)
#for (fit in oil_cv_fit_list) {
#  resid_sd <- data.frame(VarCorr(fit))$sdcor[2]
#  grand_mean <- fixef(fit)[1]
#  oil_cvs <- append(oil_cvs, coef.var(resid_sd, grand_mean))
#}
#names(oil_cvs) <- oil_trait_list
