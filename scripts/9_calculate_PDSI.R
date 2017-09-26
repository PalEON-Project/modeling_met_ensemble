##' @param path.in File path to where the CF-standard .nc files should be housed
##' @param years.pdsi - which years to calculate PDSI for; if NULL (default), all available years will be used
##' @param years.calib - years to calibrate the PDSI against;
##' @param watcap - vector of length 2 indicating the water holding capacity of the upper and lower layers of the soil
##' @return list with data frame layers for temperature, precipitaiton, daylength, and PDSI with dims=c(years, months)

path.in = "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/aggregated/month/CCSM4/CCSM4_001.01"
years.pdsi = NULL
years.calib = c(1931, 1990)
site.lat <- 42.54
site.lon <- -72.18

# ----------
# Extract & calculate our soil water values
# ----------
source("calc.awc.R")
source("pdsi1.R")
source("pdsix.R")
source("PE.thornthwaite.R")
source("soilmoi1.R")

# path.soil <- "~/PalEON_CR/ED_PalEON/MIP2_Region/phase2_env_drivers_v2/soil/"
path.soil <- "~/Dropbox/PalEON_CR/env_regional/phase2_env_drivers_v2/soil" 

sand.t <- ncdf4::nc_open(file.path(path.soil, "paleon_soil_t_sand.nc"))
sand.s <- ncdf4::nc_open(file.path(path.soil, "paleon_soil_s_sand.nc"))
clay.t <- ncdf4::nc_open(file.path(path.soil, "paleon_soil_t_clay.nc"))
clay.s <- ncdf4::nc_open(file.path(path.soil, "paleon_soil_s_clay.nc"))
depth  <- ncdf4::nc_open(file.path(path.soil, "paleon_soil_soil_depth.nc"))

lon <- ncdf4::ncvar_get(sand.t, "longitude")
lat <- ncdf4::ncvar_get(sand.t, "latitude")

x.ind <- which(lon-0.25<=site.lon & lon+0.25>=site.lon)
y.ind <- which(lat-0.25<=site.lat & lat+0.25>=site.lat)

sand1 <- ncdf4::ncvar_get(sand.t, "t_sand", c(x.ind, y.ind), c(1,1))
sand2 <- ncdf4::ncvar_get(sand.s, "s_sand", c(x.ind, y.ind), c(1,1))
clay1 <- ncdf4::ncvar_get(clay.t, "t_clay", c(x.ind, y.ind), c(1,1))
clay2 <- ncdf4::ncvar_get(clay.s, "s_clay", c(x.ind, y.ind), c(1,1))
depth2 <- ncdf4::ncvar_get(depth, "soil_depth", c(x.ind, y.ind), c(1,1))

awc1 <- calc.awc(sand1, clay1)
awc2 <- calc.awc(sand2, clay2)

wcap1 <- awc1*ifelse(depth2>30, 30, depth2-1) * 1/2.54 # 30 cm top depth * 1 in / 2.54 cm
wcap2 <- awc2*ifelse(depth2>30, depth2-30, 1) * 1/2.54 # remaining depth * 1 in / 2.54 cm

watcap <- c(wcap1, wcap2)
