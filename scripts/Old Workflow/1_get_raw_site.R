# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Extract raw met for an individual site
# Creator: Christy Rollinson, 1 July 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Extract raw met for a given site for model at different spatial & temporal scales 
# using get.raw function.
#  -- Funciton Arguments:
#     1. wd.base   = base path to the github repository; there should be subfolders of data & scripts
#     2. site.name = what you want to call the site you're extracting
#     3. lat       = latitude of the site you're extracting (N = positive, 0-90 degrees)
#     4. lon       = longitude of hte site you're processing (W = negative, -180 - 0 degrees)
#     5. ldas.type = whether to extract 3-hourly GLDAS or 1-hourly NLDAS data for the training data
#                     ** If the site is outside the lower 48, you must use GLDAS
#     6. GCM       = Which GCM you want to extract
#  -- Extracts/Processes 3 types of met data:
#     1. LDAS (NLDAS or GLDAS)
#         - Native Temporal Extent: 1980-2015
#         - Driver Usage: 1980-2015, raw
#     2. CRUNCEP
#         - Native Temporal Extent: 1901-2010
#         - Driver Usage: 1901-1979, bias-corrected, temp. downscaled
#     3. GCM
#         A. p1000 simulation
#            - Native Temporal Extent: 850-1850
#            - Driver Usage: 850-1849, bias-corrected, temp. downscaled
#         B. historical simulation
#            - Native Temporal Extent: 1850-2005
#            - Driver Usage: 1850-1900, bias-corrected, temp. downscaled
#  ** Met from each dataset will be saved in {wd.base}/data/paleon_sites/{site.name}
#  -- Function will check to see if each type of data has been done yet before processing
#  -- See get_point_raw.R for internal workflow & more details
# -----------------------------------
wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
setwd(wd.base)

site.name = "HARVARD"
site.lat  = 42.5
site.lon  = -72.18
path.out = "~/Desktop/Research/met_ensembles/data/paleon_sites"


# Download NLDAS (note: this is not a pecan script & requires fill LDAS stored somewhere locally)
source("scripts/extract_local_NLDAS.R")
ldas.type = "NLDAS"
path.nldas = "/Volumes/Celtis/Meteorology/LDAS/NLDAS_FORA0125_H.002/netcdf/"
extract.local.NLDAS(outfolder=file.path(path.out, site.name, "NLDAS"), in.path=path.nldas, 
                    start_date="1980-01-01", end_date="2015-12-31", 
                    site_id=site.name, lat.in=site.lat, lon.in=site.lon)

# Note: This keeps breaking every 5-10 years; so I'm having to go real slow at it
path.pecan <- "~/Desktop/Research/pecan/modules/data.atmosphere/R/"
source(file.path(path.pecan, "download.CRUNCEP_Global.R"))
download.CRUNCEP(outfolder=file.path(path.out, site.name, "CRUNCEP"), 
                 start_date="2007-01-01", end_date=paste0("2010-12-31"), 
                 site_id=site.name, lat.in=site.lat, lon.in=site.lon)

# Extract from the GCMs:
source("scripts/extract_local_CMIP5.R")
path.cmip5 = "/Volumes/Celtis/Meteorology/CMIP5/"
GCM.scenarios = c("p1000", "historical")
GCM.list  = c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
# GCM.list="CCSM4"
for(GCM in GCM.list){
  for(scenario in GCM.scenarios){
    if(scenario=="p1000"){
      cmip5.start = "0850-01-01"
      cmip5.end   = "1849-12-31"
    } else if (scenario == "historical"){
      cmip5.start = "1850-01-01"
      cmip5.end   = "2005-12-31"
    } else {
      stop("Scenario not implemented yet")
    }
    
    print(paste(GCM, scenario, sep=" - "))
    # huss_day_MIROC-ESM_past1000_r1i1p1_08500101-10091231.nc
    extract.local.CMIP5(outfolder = file.path(path.out, site.name, GCM, scenario), in.path = file.path(path.cmip5, GCM, scenario), 
                        start_date = cmip5.start, end_date = cmip5.end, 
                        site_id = site.name, lat.in = site.lat, lon.in = site.lon, 
                        model = GCM, scenario = scenario, ensemble_member = "r1i1p1")   
  } # end GCM.scenarios
} # End GM lop

