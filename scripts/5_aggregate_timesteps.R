# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Aggregate ensemble members to monthly and/or daily timesteps
# Creator: Christy Rollinson, 12 June 2018
# Contact: crollinson@mortonarb.org
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Script to take hourly ensemble members and aggregate them to monthly time step
# for PDSI calculation.  This script can/should be adapted to save daily time step
# at a later point in time if any models or data users want daily data.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Set up file paths, constants, etc
# 1. Extract Hourly data & aggregate to monthly time step
# 2. Save output (.csv?)
# -----------------------------------
rm(list=ls())

# -----------------------------------
# 0. Set up file paths, constants, etc
# -----------------------------------
# Setting up the path to the meteorology
path.met <- "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/"
path.out <- "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/monthly_all"
if(!dir.exists(path.out)) dir.create(path.out)

# Getting a list of what GCMs we have ensembles for
gcms <- dir(path.met)
gcms <- gcms[!gcms %in% c("subday_models", "monthly_all")]
# -----------------------------------

# -----------------------------------
# 1. Extract Hourly data & aggregate to monthly time step
#    - uses script built by drew to do the aggregation.
# -----------------------------------
library(ncdf4)
library(lubridate)
source("met_summaries.R")

# getting an estimate of how many files assuming we have even effort across GCMs
ens.gcm <- dir(file.path(path.met, gcms[1], "1hr"))
ens.gcm <- ens.gcm[!ens.gcm %in% c("subdaily_qaqc")] # Make sure to omit this folder
files.ens <- dir(file.path(path.met, gcms[1], "1hr", ens.gcm[1]))

pb <- txtProgressBar(min = 0, max = length(gcms)*length(ens.gcm)*length(files.ens), style = 3)
stat = 0
for(GCM in gcms){
  ens.gcm <- dir(file.path(path.met, GCM, "1hr"))
  ens.gcm <- ens.gcm[!ens.gcm %in% c("subdaily_qaqc")] # Make sure to omit this folder
  for(ENS in ens.gcm){
    files.ens <- dir(file.path(path.met, GCM, "1hr", ENS))
    

    for(fnow in unique(files.ens)){
      stat=stat+1
      ncT <- nc_open(file.path(path.met, GCM, "1hr", ENS, fnow))
      
      # Get daily & monthly summaries
      out.day <- daily_sums(ncT)
      out.mon <- monthly_sums(out.day)
      # out.mon$dayfact <- test.mon$hrs_sun/12*dpm/30
      
      nc_close(ncT)
      
      if(fnow==files.ens[1]){
        dat.ens <- out.mon
      } else {
        dat.ens <- rbind(dat.ens, out.mon)
      }
      
      setTxtProgressBar(pb, stat)
    } # End file loop; have 1 ensemble member done
    
    
    if(GCM == gcms[1] & ENS==ens.gcm[1]){
      tair   <- data.frame(dat.ens$tair_mean)
      precip <- data.frame(dat.ens$precip_tot)
      daylen <- data.frame(dat.ens$hrs_sun)
      names(tair) <- names(precip) <- names(daylen) <- ENS
      
      yrs <- str_pad(850:2015, width=4, pad="0")
      mos <- str_pad(1:12, width=2, pad="0")
      row.ind <- paste(rep(yrs, each=12), rep(mos, length.out=length(yrs)*12), sep="-")
      
      row.names(tair) <- row.names(precip) <- row.names(daylen) <- row.ind
      
    } else {
      tair  [,ENS] <- dat.ens$tair_mean
      precip[,ENS] <- dat.ens$precip_tot
      daylen[,ENS] <- dat.ens$daylen
    }
      
    
  } # End ensemble member loop; have 1 GCM done
} # End extraction for all GCMs/ensemble members


# Save data in this "wide" format
write.csv(tair  , file.path(path.out, "tair.csv"     ), row.names=T)
write.csv(precip, file.path(path.out, "precip.csv"   ), row.names=T)
write.csv(daylen, file.path(path.out, "daylength.csv"), row.names=T)

# out.mon$dayfact <- test.mon$hrs_sun/12*dpm/30

# -----------------------------------

