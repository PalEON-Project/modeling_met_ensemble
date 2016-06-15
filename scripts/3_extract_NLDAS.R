# Script to extract NLDAS data for one site
source("download.NLDAS.R")
source("download.NLDAS_raster.R")

start.time=Sys.time()
dat.pha <- download.NLDAS(outfolder="~/Desktop/", start_date="1986-05-16", end_date="1986-05-17", site_id=2, lat.in=42.54, lon.in=-72.18)
end.time=Sys.time()


start.time
end.time
