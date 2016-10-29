# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Create statistical models to predict subdaily meteorology from daily means
# Creator: Christy Rollinson, 28 October 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Apply the statistical models from step 3 to convert the daily, bias-corrected met 
# files from step 2 (daily means) and predict subdaily values.  This gets done by 
# filtering backwards in time starting with the present (where the trianing data is).
#
# There are ways to improve this and speed it up, but hopefully this works for now.
# We whould also probably think about applying this filter approach to the bias-
# correction step to avoid abrupt and unreasonable jumps in climate.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Load libraries, set up file paths, etc
# ----- Loop through by ensemble member ----------
#    1. Load and format prediction data (1 file from step 2)
#       1.1 Load output file from bias correction (bias ensemble member)
#       1.2 select year we're working with (single file for output)
#    ----- Loop through by year ----------
#      2. Predict subdaily values for whole year, filtering backwards in time
#      3. Write annual output into .nc files 
#         - separate file for each year/ensemle member; 
#         - all met vars in one annual file (similar to pecan met structure)
#    ----- recycle steps 2 & 3 for all years in file ----------
# ----- recycle step 1 for all files for ensemble member ----------
# -----------------------------------


# -----------------------------------
# 0. Load libraries, set up file paths, etc
# -----------------------------------
# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(MASS)
# library(lubridate)
library(ggplot2)
# library(tictoc)
rm(list=ls())


dir.base <- "/projectnb/dietzelab/paleon/met_ensemble/data/met_ensembles/HARVARD/"
# mod.out <- "../data/met_ensembles/HARVARD/subday_models"
# mod.out <- "~/Desktop/met_ensembles/HARVARD/subday_models"
# fig.dir <- file.path(mod.out, "model_qaqc")

# if(!dir.exists(mod.out)) dir.create(mod.out, recursive = T)
# if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)

GCM.list = c("CCSM4", "MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1")
# -----------------------------------

for(GCM in GCM.list){
  path.gcm <- file.path(dir.base, GCM, "day")
  setwd(path.gcm)
  # All of the daily ensembles should be zipped to save space
  ens.day <- dir(".", ".tar.bz2")
  
  for(i in 1:length(ens.day)){
    # Extract the base name of this ensemble member
    ens.name <- substr(ens.day[i], 1, nchar(ens.day[i])-8)

    # uncompress the ensemble member
    system(paste0("tar -jxvf ", ens.day[i])) 
    
    # Get list of files
    files.day <- dir(ens.name, ".nc")
    
    # -----------------------------------
    # 1. Load bias-correction file
    # NOTE: Need to start with the last and work to the first
    # -----------------------------------
    for(j in length(files.day):1){
      ncT <- nc_open(file.path(ens.name, files.day[j]))
      
      dat.nc <- data.frame(time        = ncvar_get(ncT, "time"   ),
                           tmax.day    = ncvar_get(ncT, "tmax"   ),
                           tmin.day    = ncvar_get(ncT, "tmin"   ),
                           precipf.day = ncvar_get(ncT, "precipf"),
                           swdown.day  = ncvar_get(ncT, "swdown" ),
                           lwdown.day  = ncvar_get(ncT, "lwdown" ),
                           qair.day    = ncvar_get(ncT, "qair"   ),
                           press.day   = ncvar_get(ncT, "press"  ),
                           wind.day    = ncvar_get(ncT, "wind"   )
                           )
      nc_close(ncT)
      
      dat.nc$date <- date(dat.nc$time, )
      
    }
    # -----------------------------------
      
  }
}