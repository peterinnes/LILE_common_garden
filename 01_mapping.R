# mapping LILE project sample locations
# Peter Innes
# 10.5.20

library(dplyr)
library(ggplot2)
library(raster)
library(sp)
library(sf)
library(maps)
library(tmap)

env_data <- read.csv("LILE_seed_collection_spreadsheet.csv", header=T) 
env_data$source %<>% as.factor
env_data$population %<>% as.factor

env_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!is.na(Lat) | !is.na(Long)) #keep only rows that have coordinates

env_data$State

#### mapping 
# basemap
# using maps package
states <- sf::st_as_sf(map('state', c('colorado','idaho','utah', 'nevada'), plot = FALSE, fill = TRUE))

LILE_map <- ggplot() +
  geom_sf(data = states) +
  geom_point(data = coords, aes(x=Long, y=Lat, color=Elev_m), alpha=0.75, size=3) +
  scale_color_gradient(low="forestgreen", high="royalblue1") +
  theme_classic()
LILE_map

tm_shape(states) +
  tm_raster("elevation", palette = terrain.colors(10))

png("LILE_map.png", width=8, height=7, res=300, units="in")
LILE_map
dev.off()