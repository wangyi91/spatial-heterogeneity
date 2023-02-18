# Code to generate the file "grid_elevations_info_flora" 
# Needs objects "subif" from code "plante_elevation_plot.R"

# Tutorial for elevatr: https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html#Get_Raster_Elevation_Data

library(elevatr) # get elevation data
library(raster) # raster to points
library(plyr)
library(dplyr)



# get grid information with only columns "gridcell_id","x","y"
unique_grid <- subif[,c("gridcell_id","x_wgs84","y_wgs84")] %>% unique %>% mutate(x=x_wgs84,y=y_wgs84) %>% dplyr::select(-x_wgs84,-y_wgs84)


# function to get median elevations for a list of info flora grids

get_infoflora_grid_elevation <- function(df) {
  # df is dataframe with three columns: gridcell_id, x, y
  # return median elevation in each grid: gridcell_id, median_elevation
  
  prj_dd <- "EPSG:4326" # the EPSG identifier of WGS84
  
  out <- data.frame(gridcell_id=NA, median_elevation=NA)
  i=0; n=nrow(df)
  
  for (id in df$gridcell_id) {
    elevation_df <- get_elev_raster(df[df$gridcell_id==id, c("x","y")], prj = prj_dd, expand=0.01, z = 14)
    
    # get the median elevation of grid
    median_grid_elev <- median(data.frame(rasterToPoints(elevation_df))[,3])
    
    out <- out %>% add_row(gridcell_id=id, median_elevation=median_grid_elev)
    i=i+1
    print(i/n)
  }
  
  return(out[2:nrow(out),])
}


# generate grid elevations
grid_elevations <- get_infoflora_grid_elevation(unique_grid)

saveRDS(grid_elevations, file = "./00_R-repository/grid_elevations_info_flora")
