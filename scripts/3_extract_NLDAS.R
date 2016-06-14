# Script to extract NLDAS data for one site
source("download.NLDAS.R")

start.time=Sys.time()
dat.pha <- download.NLDAS(outfolder="~/Desktop", start_date="1986-01-01", end_date="1986-12-31", site_id=1, lat.in=42.54, lon.in=-72.18)
end.time=Sys.time()


start.time
end.time
