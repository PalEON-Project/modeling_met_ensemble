#------------------------------------
# Script Information
#------------------------------------
# Purpose: Calculate available water content for a given area 
# Creator: Drew Duckett, 1 February 2017
# Contact: duckettdj@g.cofc.edu
#------------------------------------

#------------------------------------
# Description
#------------------------------------
# Function (ssurgo_extraction_v5 to download soil data from the SSURGO database  
#   and calculate the average available water content for the area weighted by 
#   polygon size

# Input:
# 1. min and max latitude and longitude of target area (for creating spatial polygon)
# 2. label (for labeling files downloaded)

# Output:
# 1. a single number, the weighted average available water content in cm
# Note: the awc is from the top 25 cm of soil

#------------------------------------
# Workflow Overview
#------------------------------------
# 1) Create extent object with min/max lat/long input. Convert to spatial polygon
# 2) Use spatial polygon to download SSURGO data, which includes spatial and tabular
# 3) Project the downloaded shapefile (aea) to calculate the area of each polygon 
# 4) Match the polygon areas to awc data in tabular format using MUKEYs and calculate
#     weighted average
#------------------------------------



ssurgo_extraction_v5 <- function(lat_min, lat_max, lon_min, lon_max, label){
  
  #load required libraries
  library(FedData)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(sp)
  
  #Get extent and download data
  area_extent <- extent(lon_min, lon_max, lat_min, lat_max) #create extent object
  area_polygon <- polygon_from_extent(area_extent, proj4string = "+proj=longlat +datum=WGS84") #turn extent into polygon
  ssurgo_extent <- get_ssurgo(template = area_polygon, label = label, extraction.dir = paste0("./EXTRACTIONS/", label, "/SSURGO")) #download ssurgo data for polygon
  
  #import shapefile
  prewd <- getwd() #get wd
  prewd <- strsplit(prewd, "/") #split the wd string
  prewd <- prewd[[1]][length(prewd[[1]])] #get last element to be used to reset wd
  wd <- paste0("~/", prewd, "/EXTRACTIONS/", label, "/SSURGO") #create path for wd
  setwd(wd) #set wd to import shapefile
  shape_label <- paste0("", label, "_SSURGO_Mapunits.shp") #create label for reading in shapefile
  ssurgo_shape <- readOGR(dsn = shape_label) #import shape file
  
  #project data to albers equal area and calculate polygon areas
  ssurgo_shape <- spTransform(ssurgo_shape, CRS = CRS("+proj=aea +ellps=WGS84 +datum=WGS84 +units=m")) #project using albers equal area - need to match proj4string to original shapefile except for +proj
  area_list <- gArea(ssurgo_shape, byid = TRUE) #calculate area for each polygon
  names(area_list) <- ssurgo_extent[["spatial"]]$MUKEY #name list indices as mukeys
  
  #calculate total area
  area_sum <- sum(area_list) #sum entire area
  
  #match areas to mukeys
  test_list <- list() #create empty list
  m = 1 #set counter
  for (mukey in unique(ssurgo_extent[["spatial"]]$MUKEY)){ #for each unique mukey
    
    index_list <- which(names(area_list) == mukey) #create list of lists where each list contains indices for each unique mukey
    test_list[m] <- list(area_list[index_list]) #create list of lists where each list contains all area measurements for a single mukey
    
    m = m + 1
  }
  names(test_list) <- unique(ssurgo_extent[["spatial"]]$MUKEY) #name list indices as mukeys
  
  #calculate area of each mukey and proportion of area compared to total area
  list_sums <- sapply(test_list, sum) #sum each list (calculate total area for each unique mukey)
  list_prop <- list_sums / area_sum #get proportion of area for each mukey
  
  #extract awc values and calculate weighted average
  awc <- ssurgo_extent[["tabular"]]$muaggatt$aws025wta #extract all awc values
  names(awc) <- ssurgo_extent[["tabular"]]$muaggatt$mukey #name list indices as mukeys
  awc_final <- sum((awc * list_prop), na.rm = TRUE) #weight awc values by proportion of area and sum to get weighted average
  
  return(awc_final)
}
