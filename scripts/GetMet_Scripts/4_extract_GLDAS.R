# Script to extract GLDAS data for one site
source("download.GLDAS.R")

dir.PHA <- "~/Desktop/LDAS/"
start.time <- Sys.time()
download.GLDAS(outfolder=dir.PHA, start_date="2001-01-01", end_date="2001-12-31", site_id="Harvard", lat.in=42.54, lon.in=-72.18)
end.time <- Sys.time()

start.time; end.time