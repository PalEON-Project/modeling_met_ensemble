# Script to extract NLDAS data for one site
source("download.NLDAS.R")

dir.PHA <- "/projectnb/dietzelab/paleon/met_ensemble/data/paleon_sites/Harvard"
download.NLDAS(outfolder=dir.PHA, start_date="1980-01-01", end_date="2015-12-31", site_id="Harvard", lat.in=42.54, lon.in=-72.18)
