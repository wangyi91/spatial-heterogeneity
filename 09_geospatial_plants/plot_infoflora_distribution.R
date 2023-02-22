
library(rgdal)
library(ggmap)
library(ggrepel)
library(ggplot2)


plot_infoflora_distribution <- function(df, tab){
  
  # df is a dataframe with columns "x_wgs84", "y_wgs84", "obs_nb" and "species", "genus"
  # tab is a dataframe combining reads and long/lat of a taxon
  
  
  # if to transfer coordination system, transfer data.frame to SpatialPointsDataFrame-class
  #df_spd <- SpatialPoints(df[,c("x_wgs84", "y_wgs84")], proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # create a data frame variable with coordinates and info we want to plot".
  #tmp = as.data.frame(df_spd)
  #colnames(tmp) = c("long", "lat")
  #tmp$obs_nb = df$obs_nb
  #tmp$species = df$species
  #tmp$genus = df$genus
  
  
  
  
  # Use ggmap to plot coordinates on Google Maps.
  # Set my API Key. Here you need paste your own API.
  ggmap::register_google(key = "MY_API_KEY")
  
  # get map of Lake Constance and catchment area
  #map <- get_googlemap(center = c(lon = 9.190567, lat = 46.8),
  #                     zoom = 8, scale = 2, size = c(1280, 1280),
  #                     maptype ='terrain',
  #                     color = 'color',
  #                     path = "&style=feature:all|element:labels|visibility:off") 
  
  map <- get_map(location = c(lon = 10.2, lat = 46.8),
                       zoom = 8, scale = 2, size = c(1280, 1280),
                       maptype ='terrain-background',
                       color = 'color', source="stamen")
                       
  
  # Plot coordinates and nb observations on map
  ggmap(map, darken = c(0.35, "white")) + 
    geom_point(data=df, aes(x = x_wgs84, y = y_wgs84, size=obs_nb), colour=alpha("#c54b8c", 0.7), shape=16) + 
    geom_point(data=tab, aes(x = long, y = lat, size=2*get(colnames(tab[1]))),
               colour=alpha("#5d76cb", 0.3), shape=15) +
    scale_size_identity(trans="sqrt",guide="legend") +
    labs(size="Number of observations") +
    theme(legend.position="bottom") +
    ggtitle(bquote(paste(~italic(.(df$genus))~italic(.(df$species)), sep=" "))) + 
    xlab("°E") +
    ylab("°N") +
    theme(
      plot.title = element_text(colour = "#87a96b", size=25)
    )
}
