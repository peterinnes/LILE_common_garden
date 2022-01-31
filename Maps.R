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
library(tmaptools)
library(maps)
library(elevatr)
library(maptools)
library(rgeos)
library(grid)

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
non_contiguous <- c("Hawaii", "Alaska")

# Subset US shapefile by desired states
state.sub <- us[as.character(us@data$NAME_1) %in% iwstates, ]

# load in the saved data 
elevation.sub <- raster("data/lile_elevation.gri")

# Use elevatr package to get elevation raster data. Don't run if .gri file is already saved.
elev_model <- get_elev_raster(state.sub, z = 6) #zoom=9
elevation.sub <- raster::crop(elev_model, extent(state.sub)) #Crop and mask the elevation raster to our subset of sates
elevation.sub <- raster::mask(elevation.sub, state.sub) #mask

writeRaster(elevation.sub, "data/lile_elevation.raster")


# basic map
plot(elevation.sub)
plot(state.sub, add = TRUE)

# plot with tmap for more control over details. Careful, actually loading the map takes a TON of memory and may crash the R sessions.  
set.seed(15)
LILE_map_2.0 <- tm_shape(elevation.sub) +
  tm_raster(title = "Elevation (m)") +
  tm_shape(state.sub) +
  tm_borders() +
  tm_shape(coords_sp) +
  tm_text("population", size = .5, auto.placement = T) +
  tm_symbols(shape=1, size = 0.2, col = "black") +
  tm_shape(garden_coords_sp) +
  tm_text("Garden", col="slateblue1", auto.placement = T) +
  tm_symbols(shape=2, size=0.6, col="slateblue1") +
  tm_grid(lines = F, labels.inside.frame = F, labels.size = 1) +
  tm_ylab("Latitude", size=1) +
  tm_xlab("Longitude", size=1) +
  tm_scale_bar(text.size=1, position = c(.625,0.325)) +
  tm_layout(frame = F, legend.title.size = 1, legend.text.size = .75,
            legend.just=c("right", "top"), legend.position = c(1,1),
            legend.width = -.275, legend.height = -.275,
            legend.frame = T)

LILE_map_2.0

# bounding box of main map 
main_xy <- st_bbox(state.sub)
asp <- (main_xy$ymax - main_xy$ymin)/(main_xy$xmax - main_xy$xmin)

# Inset map
inset_state.sub <- us[!as.character(us@data$NAME_1) %in% non_contiguous, ]
inset_map <- tm_shape(inset_state.sub, projection = 2163) +
  tm_polygons() + 
  tm_layout(inner.margins = c(.01,.01,.01,.01), outer.margins=c(0,0,0,0)) +
  tm_shape(state.sub, projection=2163) +
  tm_borders(col = "red", lwd=1)

# apect ratio of the 
inset_xy <- st_bbox(inset_state.sub)
asp2 <- (inset_xy$xmax - inset_xy$xmin)/(inset_xy$ymax - inset_xy$ymin)
w <- .225
h <- asp2 * w


vp <- viewport(x=.878, y=0.878, width = w, height=h, just=c("right", "top"))

# FIG 1
tmap_save(LILE_map_2.0,filename="Figure_1_collection_map.jpg",
          dpi=600, insets_tm=inset_map, insets_vp=vp,
          height=asp*17, width=17, units="cm")

## basic map without elevation raster
#states <- sf::st_as_sf(map('state', c('colorado','idaho','utah', 'nevada'), plot = FALSE, #fill = TRUE))
#LILE_map <- ggplot() +
#  geom_sf(data = states) +
#  geom_point(data = coords, aes(x=Long, y=Lat, color=Elev_m), alpha=0.75, size=3) +
#  scale_color_gradient(low="forestgreen", high="royalblue1") +
#  theme_classic()
#LILE_map

png("plots/LILE_map_2.0.png", width=17, height=17, res=600, units="cm")
LILE_map_2.0
dev.off()
