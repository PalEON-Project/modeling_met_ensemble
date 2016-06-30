library(RCurl)
library(ncdf4)
library(XML)
test <- getURI("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb.nc?var=Pressure&latitude=42.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&accept=netcdf&point=true")
summary(test)
test

library(XML)
test <- xmlParse("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=Pressure&latitude=45.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=xml&point=true")
test <- xmlParse("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=LW_radiation_flux_downwards_surface&var=Precipitation_hourly_total&var=Pressure&var=SW_radiation_flux_downwards_surface&var=N2-m_above_ground_Specific_humidity&var=N2-m_above_ground_Temperature&var=N10-m_above_ground_Meridional_wind_speed&var=N10-m_above_ground_Zonal_wind_speed&latitude=45.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=xml&point=true")
summary(test)
test2 <- xmlToList(test)
summary()


for(i in 1:length(test2$point)){
  names(test2$point)[i] <- test2$point[[i]]$.attrs[1]
}

dat.out <- data.frame(var=names(test2$point))
for(i in 1:length(test2$point)){
  dat.out[i,"value"] <- as.numeric(test2$point[[i]]$text)
}

val <- as.numeric(test2$point$Pressure$text)
names(test2$point)
names(test2$point)
test2$point$data
test2$point[[4]]


test <- read_xml("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=Pressure&latitude=45.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=xml&point=true")
xmlToList
test2 <- xmlRoot(test)
test2
test3 <- xmlSApply(test2, xmlValue)

summary(test$doc)


http://thredds.daac.ornl.gov/thredds/dodsC/ornldaac/1220/mstmip_driver_global_hd_climate_

/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb
test <- getURL("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/?var=Pressure&var=N2-m_above_ground_Temperature&latitude=42.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=netcdf&point=true")
test2 <- nc_open(test)



test <- getURLContent("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=Pressure&spatial=all&north=52.9375&west=-124.9375&east=-67.0625&south=25.0625&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&horizStride=&addLatLon=true")
test2 <- nc_open(con <- textConnection(test))
test

test2 <- nc_open(test)

http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb/dataset.html
NLDAS_FORA0125_MC.ACLIM12.002.grb.nc

test <- getURLContent("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=Pressure&spatial=all&north=52.9375&west=-124.9375&east=-67.0625&south=25.0625&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&horizStride=&addLatLon=true")
test2 <- nc_open(test)

test <- "http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb.nc"
test

test2 <- getURL(test)
test2 <- nc_open(test)

test <- getURL("http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=Pressure&var=N2-m_above_ground_Temperature&latitude=42.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=netcdf&point=true")

http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb?var=LW_radiation_flux_downwards_surface&var=Precipitation_hourly_total&var=Pressure&var=SW_radiation_flux_downwards_surface&var=N2-m_above_ground_Specific_humidity&var=N2-m_above_ground_Temperature&var=N10-m_above_ground_Meridional_wind_speed&var=N10-m_above_ground_Zonal_wind_speed&latitude=42.54&longitude=-72.18&temporal=all&time_start=1980-12-01T00%3A00%3A00Z&time_end=1980-12-01T00%3A00%3A00Z&time=1980-12-01T00%3A00%3A00Z&vertCoord=&accept=netcdf&point=true, 
http://hydro1.sci.gsfc.nasa.gov/thredds/ncss/grid/NLDAS_FORA0125_MC.002/NLDAS_FORA0125_MC.ACLIM12.002.grb/pointDataset.html

NLDAS_FORA0125_MC.ACLIM12.002.grb.nc