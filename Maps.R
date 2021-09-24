# mapping LILE project sample locations
# Peter Innes
# 10.5.20

library(dplyr)
library(ggplot2)
library(raster)
library(sp)
library(rgdal)
library(sf)
library(tmap)
library(maps)
library(elevatr)
library(maptools)
library(rgeos)

env_data <- read.csv("data/LILE_seed_collection_spreadsheet.csv", header=T) %>% 
  mutate(source=as.factor(source), population=as.factor(population))

geo_data <- env_data %>% dplyr::select(source,population,Lat,Long,Elev_m) %>%
  filter(!source %in% c(2,5,22,32,38)) %>% #remove mistaken/duplicate Appar
  filter(!is.na(Lat) | !is.na(Long)) #only keep sources with lat/long data (exludes Appar and source 37)

coords <- data.frame(Long=geo_data$Long, Lat=geo_data$Lat, Elev_m=geo_data$Elev_m, population=geo_data$population, row.names = geo_data$population) %>% na.omit()

# coordinates of the two gardens
ephraim_coords <- data.frame(Long=-111.5782, Lat=39.3706) #elevation 5532'
milville_coords <- data.frame(Long=-111.816, Lat=41.656)
garden_coords <- rbind(ephraim_coords, milville_coords)
garden_coords$Garden <- c("Ephraim", "Milville")

# convert coords df to a SpatialPointsDataframe so we can map the sites
coords_sp <- SpatialPointsDataFrame(coords = coords[,c(1,2)],
                                    data = coords,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

garden_coords_sp <-  SpatialPointsDataFrame(coords = garden_coords[,c(1,2)],
                                            data = garden_coords,
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#### Mapping ####
us <- raster::getData("GADM", country="USA", level=1) #get US States shapefile
iwstates <- c("Colorado", "Utah", "Idaho", "Nevada") #list of states I want to plot

elevation.sub <- raster("data/lile_elevation.gri")

# Subset US shapefile by desired states
state.sub <- us[as.character(us@data$NAME_1) %in% iwstates, ]

# Use elevatr package to get elevation raster data
elev_model <- get_elev_raster(state.sub, z = 9) #zoom=9

# Crop and mask the elevation raster to our subset of sates
elevation.sub <- raster::crop(elev_model, extent(state.sub))
elevation.sub <- mask(elevation.sub, state.sub)

writeRaster(elevation.sub, "data/lile_elevation.raster")


# basic map
plot(elevation.sub)
plot(state.sub, add = TRUE)

# plot with tmap for more control over details 
LILE_map_2.0 <- tm_shape(elevation.sub) +
  tm_raster(title = "Elevation (m)") +
  tm_legend(legend.position = c("right","top")) +
  tm_shape(state.sub) +
  tm_borders() +
  tm_shape(coords_sp) +
  tm_text("population", size = .75, auto.placement = T) +
  tm_symbols(shape=1, size = 0.4, col = "black") +
  tm_layout(frame = F, legend.title.size = 1.5, legend.text.size = 1) +
  tm_shape(garden_coords_sp) +
  tm_text("Garden", col="slateblue1", auto.placement = T) +
  tm_symbols(shape=2, size=0.6, col="slateblue1") +
  tm_grid(lines = F, labels.size = 1) +
  tm_xlab("Longitude", size=1.5) +
  tm_ylab("Latitude", size=1.5)

LILE_map_2.0

# basic map without elevation raster
states <- sf::st_as_sf(map('state', c('colorado','idaho','utah', 'nevada'), plot = FALSE, fill = TRUE))
LILE_map <- ggplot() +
  geom_sf(data = states) +
  geom_point(data = coords, aes(x=Long, y=Lat, color=Elev_m), alpha=0.75, size=3) +
  scale_color_gradient(low="forestgreen", high="royalblue1") +
  theme_classic()
LILE_map

png("plots/LILE_map_2.0.png", width=17, height=17, res=600, units="cm")
LILE_map_2.0
dev.off()
