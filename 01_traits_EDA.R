#' ---
#' output: github_document
#' ---

#' Exploratory Data analysis of LILE traits
#' author: Peter Innes
#' date: 1.13

#+ results=FALSE, message=FALSE, warning=FALSE
library(dplyr)
library(ggplot2)
library(lme4)
library(arm)
library(magrittr)
library(tidyr)

#' Read-in collection data
env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>%
  dplyr::select(source,population,Lat,Long,Elev_m)
env_data$source %<>% as.factor
env_data$population %<>% as.factor

#' Read-in seed weight data
sd_wt_data <- read.csv("data/cleaned_LILE_yield_data_2013_seed_wt.csv", header = T)
sd_wt_data$source %<>% as.factor
sd_wt_data$block %<>% as.factor
sd_wt_data <- filter(sd_wt_data, is.na(notes)) %>% #filter out observations with fewer than 50 seeds (described in notes column). All other obs are weights of 50 seeds. 
  dplyr::select(!notes) %>% #don't need notes column anymore
  filter(!source %in% c(2,5,22,32,38)) %>% # exclude these sources bc they were found to be mostly 'Appar', which is already represented (source 41)
  left_join(dplyr::select(env_data,source,population,Lat))

# Read-in fruit fill data.
ff_data <- read.csv("data/cleaned_LILE_yield_data_2013_fruit_fill.csv", header=T)
ff_data$source %<>% as.factor
ff_data$block %<>% as.factor

Appar <- ff_data %>% #list of individual Appar plants to exclude. 6 plants from source 22 are listed as Appar. Brent says go ahead and exclude this entire source.
  filter(notes==c("Appar")) %>%
  filter(trt=="B") %>%
  dplyr::select(source,trt,block,row,plot,plant)
Appar

ff_data <- ff_data %>% 
  filter(trt=="B") %>%  #filter out 'trt A' (non-study/non-harvested plants). 
  filter(!source %in% c(2,5,22,32,38)) %>% #exclude mistaken Appar sources
  left_join(dplyr::select(env_data,source,population,Lat))  #add population ID, latitude

# Read-in stem data
stem_data <- read.csv("data/cleaned_LILE_yield_data_2013_stem_and_fruit.csv", header=T)
stem_data$source %<>% as.factor
stem_data$block %<>% as.factor
stem_data <- stem_data %>% 
  filter(trt=="B") %>%
  filter(!source %in% c(2,5,22,32,38)) %>%
  left_join(dplyr::select(env_data,source,population,Lat)) 

# gather traits to estimate Yield. The method here is to multiply the trait values within each accession at the lowest level possible, since we lack individual plant data for seed weight (the seed weight values are pooled at the 'plot' level—population within block) Also, we have to take averages, at the plant level, of the fruit per stem and buds/flowers per stem traits, since we have those counts for multiple stems (up to 20) per plant. Also of note is several cases there are two plants per block, due to sampling methods top 8 vigourous plants across all blocks selected as the 'trt B' study plants.
a <- stem_data %>%
  dplyr::select(source,population,block,row,plot,plant,num_of_stems) %>% 
  unique() #%>% 
#group_by(source,block) %>%
#summarise(num_stems=mean(num_of_stems))
b <- stem_data %>%
  group_by(source,population,block,row,plot,plant) %>%
  summarise(fruit_per_stem=mean(na.omit(fruits)), #mean number fruits, number buds and flowers, number forks, capsule diam, and stem diam, for individual plants.
            bds_flws_per_stem=mean(na.omit(bds_flow)),
            forks=mean(na.omit(forks)),
            caps_diam=mean(na.omit(diam_caps)),
            stem_diam=mean(na.omit(diam_stem))) 
c <- ff_data %>%
  group_by(source,population,block,row,plot,plant) %>%
  dplyr::select(source,population,block,row,plot,plant,good_fill) #%>%
#na.omit()
d <- sd_wt_data %>%
  group_by(source,population,block, row, plot) %>%
  summarise(sd_wt_50_ct=mean(sd_wt_50_ct)) #take average at the pooled population:block level
yield_df <- full_join(a,b)
yield_df <- full_join(yield_df,c)
yield_df <- full_join(yield_df,d)

# next, impute missing seed weight values so that we can estimate yield: For plants from plots without seed weight data, use the mean of seed weight data from all other plots of its population. 
for ( i in 1:281 ){
  pop <- yield_df$population[i]
  pop_mn <- yield_df %>% filter(population==pop) %>% 
    dplyr::select(source,row,plot,sd_wt_50_ct) %>%
    unique() %>%
    summarise(mean(na.omit(sd_wt_50_ct)))
  pop_mn %<>% as.numeric()
  if( is.na(yield_df$sd_wt_50_ct[i]) ){
    print(yield_df[i,])
    yield_df$sd_wt_50_ct[i] <- pop_mn
  }
}

yield_df <- yield_df %>%
  na.omit() %>%
  mutate(EST_YIELD = num_of_stems 
         * (fruit_per_stem + bds_flws_per_stem) 
         * good_fill 
         * (sd_wt_50_ct/50)) %>%
  group_by(population,block) %>%
  arrange(as.character(population))
write.csv(yield_df,file="data/yield_data.csv", row.names = FALSE)


#' #### Exploratory data analysis
#' Histogram of seed weight
sd_wt_data %>%
  ggplot() +
  geom_histogram(mapping=aes(x=sd_wt_50_ct,y=stat(density)),bins=30)

#' Violin plots of seed weight by population, and by block, separately
ggplot(data=sd_wt_data, aes(x=reorder(population, sd_wt_50_ct), y=sd_wt_50_ct)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))
ggplot(data=sd_wt_data, aes(x=reorder(block, sd_wt_50_ct), y=sd_wt_50_ct)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) #does not appear to be any block effect. good. 

#' seed weight by population and block, simultaneously. A couple things jump out from this plot: NV 3 has way more variation than most of the other populations, UT_7 seems to have a couple outliers, but looks like it is just due to a single block
sd_wt_data %>%
  mutate(yjit=jitter(0*sd_wt_50_ct)) %>%
  ggplot() +
  geom_point(mapping=aes(x=sd_wt_50_ct, col=block, y=yjit),shape=1,alpha=0.5) +
  facet_wrap(facets = ~ population) +
  ylim(-0.1,0.1)

#' double check sample sizes to verify balanced design. Most all pops have n=32, a few populations have n=28. Pop 38 has n=12—this may be of concern. But, looking at the previous plots, pop 38 does not appear to have elevated uncertainty.
sd_wt_data %>%
  group_by(population) %>%
  summarise(sample_size=n()) %>%
  arrange(-sample_size) %>%
  print(n=Inf)

#' Based on initial histogram I suspect a normal distribution will be appropriate for modeling seed weight. Let's check the normal fit for just a few populations.

#' summary stats for 5 counties
#+ results=FALSE, message=FALSE, warning=FALSE
set.seed(4)
sw5pops <- sd_wt_data %>% 
  dplyr::select(population,block,sd_wt_50_ct) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(sd_wt_50_ct), sd=sd(sd_wt_50_ct), n=n(), min=min(sd_wt_50_ct), max=max(sd_wt_50_ct)) %>%
  sample_n(5)
sw5pops

#' Normal fitted for the 5 pops
norm_df <- NULL
for ( i in 1:5 ) {
  x <- seq(sw5pops$min[i],sw5pops$max[i],length.out = 100)
  y <- dnorm(x, sw5pops$mean[i], sw5pops$sd[i])
  norm_df <- rbind(norm_df,data.frame(x,y,population=sw5pops$population[i]))
}
rm(x,y) #clean up
head(norm_df)

#' plot observed data with expected distribution overlaid, by population. The normal appears to fit some populations better than others.
#+ results=FALSE, message=FALSE, warning=FALSE
sd_wt_data %>%
  group_by(population) %>%
  filter(population%in%sw5pops$population) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=sd_wt_50_ct, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=sd_wt_50_ct), col="blue") +
  geom_line(data=norm_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ population)

#' Now for the whole data set, by block
#+ results=FALSE, message=FALSE, warning=FALSE
sw_byblock <- sd_wt_data %>% 
  dplyr::select(population,block,sd_wt_50_ct) %>%
  na.omit() %>%
  group_by(block) %>%
  summarise(mean=mean(sd_wt_50_ct), sd=sd(sd_wt_50_ct), n=n(), min=min(sd_wt_50_ct), max=max(sd_wt_50_ct))

#' Normal fitted for 8 blocks
norm_df <- NULL
for (i in 1:length(sw_byblock$block)) {
  x <- seq(sw_byblock$min[i],sw_byblock$max[i],length.out = 100)
  y <- dnorm(x, sw_byblock$mean[i], sw_byblock$sd[i])
  norm_df <- rbind(norm_df, data.frame(x,y, block=sw_byblock$block[i]))
}
rm(x,y)
head(norm_df)

#' Plot observed data and fitted normal. Looks good.
#+ results=FALSE, message=FALSE, warning=FALSE
sd_wt_data %>%
  ggplot() +
  geom_histogram(mapping=aes(x=sd_wt_50_ct, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=sd_wt_50_ct), col="blue") +
  geom_line(data=norm_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ block)

#' EDA of fruit fill
# summary stats for 5 counties
#+ results=FALSE, message=FALSE, warning=FALSE
set.seed(7)
ff5pops <- ff_data %>% 
  dplyr::select(population,block,good_fill) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(good_fill), sd=sd(good_fill), n=n(), min=min(good_fill), max=max(good_fill)) %>%
  sample_n(5)
ff5pops

#' Normal fitted for the 5 pops
norm_df <- NULL
for ( i in 1:5 ) {
  x <- seq(ff5pops$min[i],ff5pops$max[i],length.out = 100)
  y <- dnorm(x, ff5pops$mean[i], ff5pops$sd[i])
  norm_df <- rbind(norm_df,data.frame(x,y,population=ff5pops$population[i]))
}
rm(x,y) #clean up
head(norm_df)

#' plot observed vs fitted. normal distro appears to fit sufficiently well.
ff_data %>% group_by(population) %>%
  filter(population%in%ff5pops$population) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=good_fill, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=good_fill), col="blue") +
  geom_line(data=norm_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ population)

#' #### EDA of number of stems per plant

#' Histograms of each trait. For this project, just interested in num_of_stems. The distribution of num_of_stems is skewed right, as expected for count data. Poisson or negative binomial could be a good fit, depending on how the mean and variance compare.
#+ results=FALSE, message=FALSE, warning=FALSE
stem_data %>%
  dplyr::select(fruits,bds_flow,forks,diam_caps,diam_stem, num_of_stems) %>%
  gather(key="trait", value="trait_value") %>%
  ggplot() +
  geom_histogram(mapping=aes(x=trait_value,y=stat(density)), bins=75) +
  facet_wrap(facets = ~ trait, scales="free")

#' Will just look at num_of_stems moving forward 
stem_data <- stem_data %>% dplyr::select(population,block,row,plot,plant,num_of_stems) %>%
  unique() #unique values only because original spreadsheet had both plant-level and stem-level data. We just want number of stems per plant, which is plant-level.

#' Summary statistics, by population. Variance appears much larger than the mean, so it's over-dispersed. Poisson might not be appropriate. negative binomial instead? 
ns_summary <- stem_data %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(num_of_stems),sd=sd(num_of_stems),var=sd(num_of_stems)^2,n=n()) 
head(ns_summary) 
table(ns_summary$n) #Mostly balanced design, with n=8 for 34 populations, n=7 for 1 pop, n=5 for 2 pops, n=4 for 2 pops.

#' Histograms by population. Small sample size at the population level could make this difficult to model—the distributions at this level are mostly uniform.
#+  results=FALSE, message=FALSE, warning=FALSE
stem_data %>%
  ggplot() +
  geom_histogram(mapping=aes(x=num_of_stems,y=stat(density))) +
  facet_wrap(facets = ~ population, scales="free")

#' Violin plots of stem_number by population, and by block, separately. Populations with the largest num_of_stems value also appear to have the largest variance. Blocks 1,2,3 appear to have more outliers in large stem number values compared to the other blocks. For stem data, there is just one plant per block.
#+ results=FALSE, message=FALSE, warning=FALSE
ggplot(data=stem_data, aes(x=reorder(population, num_of_stems), y=num_of_stems)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))
#+ results=FALSE, message=FALSE, warning=FALSE
ggplot(data=stem_data, aes(x=block, y=num_of_stems)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) 

#' Number of stems by population and block, simultaneously. Population is faceted, block is color-coded. This plot reveals some peculilarities in the data: not every population is represented in every block. Some populations have multiple plants from same block. But, it is still close to a balanced design.
#+ results=FALSE, message=FALSE, warning=FALSE
stem_data %>%
  mutate(yjit=jitter(0*num_of_stems)) %>%
  ggplot() +
  geom_point(mapping=aes(x=num_of_stems, col=block, y=yjit),shape=1,alpha=0.5) +
  facet_wrap(facets = ~ population) +
  ylim(-0.1,0.1)

#' Now to decide what distribution to use in modeling number of stems. It appears to be a highly variable trait—variance is much larger than the mean, nearly across the board. This coincides with the distributions being uniform at the population level, for most of the populations. I'm not really sure how to proceed at this point. Negative binomial seems reasonable at the scale of the entire experiment, but when estimating means at the population level, would it still be appropriate? \n


#' EDA for estimated yield.
yield_summ <- yield_df %>%
  group_by(source,population) %>%
  na.omit() %>%
  summarise(mean=mean(EST_YIELD), se=sd(EST_YIELD)/sqrt(n()), n=n())

yield_df %>%
  ggplot() +
  geom_histogram(mapping=aes(x=EST_YIELD,y=stat(density)),bins=30)

ggplot(data=yield_df, aes(x=reorder(population, EST_YIELD), y=EST_YIELD)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))

ggplot(data=yield_df, aes(x=reorder(block, EST_YIELD), y=EST_YIELD)) +
  geom_violin(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) #does not appear to be any block effect. good. 

yield_df %>% 
  na.omit() %>%
  mutate(yjit=jitter(0*EST_YIELD)) %>%
  ggplot() +
  geom_point(mapping=aes(x=EST_YIELD, col=block, y=yjit),shape=1,alpha=0.5) +
  facet_wrap(facets = ~ population) +
  ylim(-0.1,0.1)

yield_df %>% na.omit() %>%
  group_by(population) %>%
  summarise(sample_size=n()) %>%
  arrange(-sample_size) %>%
  print(n=Inf)

# check fit of different distros
set.seed(39)
yield_df <- mutate(yield_df, log_EST_YIELD=log(ifelse(EST_YIELD==0,0.1,EST_YIELD)))
yield5pops <- yield_df %>% 
  dplyr::select(population,block,EST_YIELD,log_EST_YIELD) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(EST_YIELD), sd=sd(EST_YIELD), min=min(EST_YIELD), max=max(EST_YIELD), n=n(), logmean=mean(log_EST_YIELD), logsd=sd(log_EST_YIELD), logmin=min(log_EST_YIELD), logmax=max(log_EST_YIELD)) %>%
  sample_n(25)
yield5pops

# normal fitted for the 5 pops
norm_df <- NULL
for ( i in 1:25 ) {
  x <- seq(yield5pops$min[i],yield5pops$max[i],length.out = 100)
  y <- dnorm(x, yield5pops$mean[i], yield5pops$sd[i])
  norm_df <- rbind(norm_df,data.frame(x,y,population=yield5pops$population[i]))
}
rm(x,y) #clean up
head(norm_df)

#' gamma fitted for the 5 pops
gamma_df <- NULL
for ( i in 1:25 ) {
  x <- seq(yield5pops$min[i],yield5pops$max[i],length.out = 100)
  y <- dgamma(x, (yield5pops$mean[i]/yield5pops$sd[i])^2, yield5pops$mean[i]/yield5pops$sd[i]^2)
  gamma_df <- rbind(gamma_df,data.frame(x,y,population=yield5pops$population[i]))
}
rm(x,y) #clean up
head(gamma_df)

# lognormal fitted. this distro appears to fit best.
lognorm_df <- NULL
for ( i in 1:25 ) {
  x <- seq(yield5pops$logmin[i],yield5pops$logmax[i],length.out = 100)
  y <- dnorm(x, yield5pops$logmean[i], yield5pops$logsd[i])
  lognorm_df <- rbind(lognorm_df,data.frame(x,y,population=yield5pops$population[i]))
}
rm(x,y) #clean up
head(lognorm_df)

yield_df %>%
  group_by(population) %>%
  filter(population%in%yield5pops$population) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=EST_YIELD, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=EST_YIELD), col="blue") +
  geom_line(data=norm_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ population)

yield_df %>%
  group_by(population) %>%
  filter(population%in%yield5pops$population) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=EST_YIELD, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=EST_YIELD), col="blue") +
  geom_line(data=gamma_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ population)

yield_df %>%
  group_by(population) %>%
  filter(population%in%yield5pops$population) %>%
  ggplot() +
  geom_histogram(mapping=aes(x=log_EST_YIELD, y=stat(density)), bins=30) +
  geom_density(mapping=aes(x=log_EST_YIELD), col="blue") +
  geom_line(data=lognorm_df, mapping=aes(x=x,y=y), col="red") +
  facet_wrap(facets = ~ population)

#' EDA misc scraps
# histograms
stem_data %>%
  select(fruits,bds_flow,forks,diam_caps,diam_stem, num_of_stems) %>%
  gather(key="trait", value="trait_value") %>%
  ggplot() +
  geom_histogram(mapping=aes(x=trait_value,y=stat(density)), bins=75) +
  facet_wrap(facets = ~ trait, scales="free")

sd_wt_data %>%
  ggplot() +
  geom_histogram(mapping=aes(x=sd_wt_50_ct,y=stat(density)),bins=30)

quartz()
ff_data %>%
  select(good_sds,bad_sds,tot_sds,good_fill,tot_fill) %>%
  gather(key="trait", value="trait_value") %>%
  ggplot() +
  geom_histogram(mapping=aes(x=trait_value,y=stat(density)), bins=30) +
  facet_wrap(facets = ~ trait, scales="free")

# data summaries
sw_summ <-
  sd_wt_data %>% 
  dplyr::select(population,block,sd_wt_50_ct) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(sd_wt_50_ct), sd=sd(sd_wt_50_ct), n=n(), min=min(sd_wt_50_ct), max=max(sd_wt_50_ct))

sn_summ <- stem_data %>%
  dplyr::select(source,population,block,row,plot,plant,num_of_stems) %>%
  unique() %>%
  na.omit() %>%
  group_by(population, source) %>%
  summarise(mean=mean(num_of_stems),se=sd(num_of_stems)/sqrt(n()), cv=100*(sd(num_of_stems)/mean(num_of_stems)), n=n()) #variance is much larger than the mean, so its overdispersed. poisson might not be appropriate. negative binomial instead?
View(sn_summ)

fruit_summ <- stem_data %>%
  select(population,block,row,plant,stem_no,fruits) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(fruits), sd=sd(fruits), var=sd(fruits)^2, n=n()) #also overdispersed. poisson might not be appropriate. negative binomial instead?

# compare 'trt A (non-study)' with 'trt B (study)'
ff_data %>%
  group_by(trt) %>%
  summarise(mean=mean(na.omit(good_fill)), sd=sd(na.omit(good_fill)), n=n())
#ggplot(aes(x=trt, y=good_fill)) +
#geom_violin()


forks_summ <- stem_data %>%
  select(population,block,row,plant,stem_no,forks) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(forks), sd=sd(forks), var=sd(forks)^2, n=n()) #mostly underdispersed, but some overdispersion (at population level). Also, not all pops have fork numbers from full 160 stems (8 blocks x 1 plant x 20 stems)

# check for block effect in forks data
stem_data %>%
  group_by(block) %>%
  ggplot(aes(x=block, y=forks)) +
  geom_violin()

ff_summ <- ff_data %>%
  select(population,block,good_fill) %>%
  na.omit() %>%
  group_by(population) %>%
  summarise(mean=mean(good_fill), sd=sd(good_fill), var=sd(good_fill)^2, n=n()) #fruit fill is underdispersed—variance is less than the mean.

#### Box plots of all the traits
# estimated yield
ggplot(data=yield_df, aes(x=reorder(population,EST_YIELD), y=EST_YIELD)) +
  geom_boxplot()

# seed weight
ggplot(data=sd_wt_data, aes(x=population, y=sd_wt_50_ct)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))

# fruit fill
ggplot(data=filter(ff_data, trt=="B"), aes(x=population, y=good.fill)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))

# fruit per stem
ggplot(data=filter(stem_data, trt=="B"), aes(x=population, y=fruits)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) 

# number of stems
ggplot(data=filter(stem_data, trt=="B"), aes(x=population, y=num_of_stems)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))

# forks per stem
ggplot(data=filter(stem_data, trt=="B"), aes(x=population, y=forks)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) 

# stem diameter
ggplot(data=filter(stem_data, trt=="B"), aes(x=population, y=diam.stem)) + # big outlier in source 6. 8mm stem diameter?? 
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6)) 

# capsule diameter
ggplot(data=filter(stem_data, trt=="B"), aes(x=population, y=diam.caps)) +
  geom_boxplot(alpha=0.5) + 
  theme(axis.text.x = element_text(angle=45, size=6))


