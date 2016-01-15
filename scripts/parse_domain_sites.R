# ------------------------------------------------------------------------------------
# This file extracts and stores the following two products from CMIP5 Met drivers:
# 1) The PalEON domain from CMIP5 Met drivers 
# 2) Individual site-variables (not yet implemented)
#
# This script is set to work with 3 data sets that will be used for the paleon drivers
# 1) Experiment: past1000,   Resolution: Daily
# 2) Experiment: historical, Resolution: Daily
# 3) Experiment: historical, Resolution: 3-hourly
#
# Christy Rollinson, crollinson@gmail.com
# 10 April, 2015
# ------------------------------------------------------------------------------------

#---------------------------------------
#Load libraries & set directories
#---------------------------------------
library(raster)
library(ncdf4)

# directories (or parts) for the file structure
setwd("/projectnb/dietzelab/paleon/ED_runs/met_ensemble/data")
in.dir <- "full_raw"
out.domain <- "paleon_domain"
out.sites <- "paleon_sites"
test.dir <- "~/Desktop/CMIP5_climate_daily"

# The data structure; currently not separating driver types since it all needs to get done
models <- c("MRI-CGCM3")
experiments <- c("p1000", "historical")
time.res <- c("3hr", "day")
cmip5.spatref <- CRS("+proj=longlat")

#---------------------------------------

#---------------------------------------
# Read in & parse the different data sets for the Domain
#---------------------------------------
paleon.domain <- extent(c(360-98.6, 360-66.1, 36.5, 49.75))

files <- dir(file.path(in.dir, models[m], "p1000"), ".nc")
# A test on a single file on p1000
for(i in 1:length(files)){
	nc <- stack(file.path(in.dir, models[m], "p1000", files[i]))
	proj4string(nc) <- cmip5.spatref
#	nc

	domain <- crop(nc, paleon.domain, filename=file.path(out.domain, models[m], "p1000", paste0(substr(files[1], 1, nchar(files[1])-20), ".nc")), format="CDF", overwrite=T, bylayer=T, suffix=names(nc))
	}
#---------------------------------------

# #---------------------------------------
# # Read in & parse the different data sets for Individual Sites
# #---------------------------------------
# sites <- data.frame(Site=c( "PHA",  "PHO",  "PUN",  "PBL",  "PDL", "PMB" ), 
					# Lon =c(-72.18, -68.73, -89.53, -94.58, -95.17, -82.83) + 360,
					# Lat =c( 42.54,  45.25,  46.22,  46.28,  47.17,  43.61)
					# )
# coordinates(sites) <- sites[,c("Lon", "Lat")]

# paleon.domain <- extent(c(360-98.6, 360-66.1, 36.5, 49.75))
# #files <- dir(file.path(test.dir, models[m], "p1000"), ".nc")
# files <- dir(file.path(test.dir, "paleon_domain"), ".nc")
					
# for(i in 1:length(files)){
	# for(s in 1:nrow(sites)){
# #		nc <- stack(file.path(out.domain, models[m], "p1000", files[i]))
		# nc <- stack(file.path(test.dir, "paleon_domain", files[i]))
		# nc
		
		# crop()
	# }
# }
# #---------------------------------------
