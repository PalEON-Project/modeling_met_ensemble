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
#       1.2 select year we're working with
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
library(lubridate)
library(ggplot2)
library(stringr)
library(tictoc)
library(parallel)
# library(tictoc)
rm(list=ls())
set.seed(0017)

# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble/"
wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble/"
setwd(wd.base)

# Load the scripts that do all the heavy lifting
source("scripts/temporal_downscale.R")
source("scripts/temporal_downscale_functions.R")


# dat.base <- "/projectnb/dietzelab/paleon/met_ensemble/data/met_ensembles/HARVARD/"
# dat.train <- read.csv("/projectnb/dietzelab/paleon/met_ensemble/data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")

# dat.base <- "~/Desktop/met_ensembles/HARVARD/"
dat.base <- "~/Desktop/met_bias_day/data/met_ensembles/HARVARD/"

dat.train <- read.csv(file.path(wd.base, "data/paleon_sites/HARVARD/NLDAS_1980-2015.csv"))

# if(!dir.exists(mod.out)) dir.create(mod.out, recursive = T)
# if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)

# Hard-coding numbers for Harvard
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18

# GCM.list = c("CCSM4", "MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1")
GCM.list = "MIROC-ESM"
ens.hr  <- 3 # Number of hourly ensemble members to create
n.day <- 5 # Number of daily ensemble members to process
yrs.plot <- c(2015, 1985, 1920, 1875, 1800, 1000, 850)
years.sim=2015:1900
cores.max = 12

# Defining variable names, longname & units
vars.info <- data.frame(name    =c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind"),
                        name.cf = c("air_temperature", 
                                    "precipitation_flux",
                                    "surface_downwelling_shortwave_flux_in_air",
                                    "surface_downwelling_longwave_flux_in_air",
                                    "air_pressure",
                                    "specific_humidity",
                                    "wind"
                                    ),
                        longname=c("2 meter mean air temperature", 
                                   "cumulative precipitation (water equivalent)",
                                   "incident (downwelling) showtwave radiation",
                                   "incident (downwelling) longwave radiation",
                                   'Pressure at the surface',
                                   'Specific humidity measured at the lowest level of the atmosphere',
                                   'Wind speed' 
                                   ),
                        units= c("K", "kg m-2 s-1", "W m-2", "W m-2", "Pa", "kg kg-1", "m s-1")
                        )
# Make a few dimensions we can use
dimY <- ncdim_def( "lon", units="degrees", longname="latitude", vals=site.lat )
dimX <- ncdim_def( "lat", units="degrees", longname="longitude", vals=site.lon )
# -----------------------------------

# NOTE: all precip needs to be converted precip back to kg/m2/s from kg/m2/day
# This gets done when formatting things for downscaling

for(GCM in GCM.list){
  # tic()
  # Set the directory where the output is & load the file
  path.gcm <- file.path(dat.base, GCM, "day")
  dat.day <- dir(path.gcm, ".Rdata")
  
  load(file.path(path.gcm, dat.day)) # Loads dat.out.full

  # Set & create the output directory
  path.out <- file.path(dat.base, "test_ensembles", GCM, "1hr")
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # -----------------------------------
  # 1. Format output so all ensemble members can be run at once
  # NOTE: Need to start with the last and work to the first
  # -----------------------------------
  if(is.null(years.sim)){
    years.sim <- max(dat.out.full$tmax$sims$year):min(dat.out.full$tmax$sims$year)
  }
  
  tot.ens <- which(substr(names(dat.out.full$tmax$sims),1,1)=="X")
  ens.day <- sample(1:length(tot.ens), n.day, replace=F) # For now, randomly choose which ensemble members to downscale
  # Initialize the lags
  lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally
  for(e in ens.day){
    tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
    tmax.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tmax"]
    tmin.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tmin"]
    precipf.init <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"precipf"]/(60*60*24)
    swdown.init  <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"swdown"]
    lwdown.init  <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"lwdown"]
    press.init   <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"press"]
    qair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"qair"]
    wind.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"wind"]
    
    lags.init[[paste0("X", e)]][["tair"   ]] <- data.frame(array(tair.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["tmax"   ]] <- data.frame(array(tmax.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["tmin"   ]] <- data.frame(array(tmin.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["precipf"]] <- data.frame(array(precipf.init, dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["swdown" ]] <- data.frame(array(swdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["lwdown" ]] <- data.frame(array(lwdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["press"  ]] <- data.frame(array(press.init  , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["qair"   ]] <- data.frame(array(qair.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["wind"   ]] <- data.frame(array(wind.init   , dim=c(1, ens.hr)))
  }
  
  for(y in years.sim){
    # Create a list with just the data from that year
    dat.yr <- list()
    dat.yr$tmax    <- dat.out.full$tmax   $sims[dat.out.full$tmax   $sims$year==y,]
    dat.yr$tmin    <- dat.out.full$tmin   $sims[dat.out.full$tmin   $sims$year==y,]
    dat.yr$precipf <- dat.out.full$precipf$sims[dat.out.full$precipf$sims$year==y,]
    dat.yr$swdown  <- dat.out.full$swdown $sims[dat.out.full$swdown $sims$year==y,]
    dat.yr$lwdown  <- dat.out.full$lwdown $sims[dat.out.full$lwdown $sims$year==y,]
    dat.yr$qair    <- dat.out.full$qair   $sims[dat.out.full$qair   $sims$year==y,]
    dat.yr$press   <- dat.out.full$press  $sims[dat.out.full$press  $sims$year==y,]
    dat.yr$wind    <- dat.out.full$wind   $sims[dat.out.full$wind   $sims$year==y,]
    
    # Setting up the 'next' dataset
    dat.nxt <- list()
    rows.next <-  which(dat.out.full$tmax$sims$time>=min(dat.yr$tmax$time)-1 & dat.out.full$tmax$sims$time<=max(dat.yr$tmax$time)-1)
    dat.nxt$tmax    <- dat.out.full$tmax   $sims[rows.next,]
    dat.nxt$tmin    <- dat.out.full$tmin   $sims[rows.next,]
    dat.nxt$precipf <- dat.out.full$precipf$sims[rows.next,]
    dat.nxt$swdown  <- dat.out.full$swdown $sims[rows.next,]
    dat.nxt$lwdown  <- dat.out.full$lwdown $sims[rows.next,]
    dat.nxt$press   <- dat.out.full$press  $sims[rows.next,]
    dat.nxt$qair    <- dat.out.full$qair   $sims[rows.next,]
    dat.nxt$wind    <- dat.out.full$wind   $sims[rows.next,]
    
    # If this is the first year in the dataset (850), there is no "next" value, 
    # so we need to add it in to have the right dimensions
    if(y == min(dat.out.full$tmax$sims$year)){
      dat.nxt$tmax    <- rbind(dat.nxt$tmax   [1,], dat.nxt$tmax   )
      dat.nxt$tmin    <- rbind(dat.nxt$tmin   [1,], dat.nxt$tmin   )
      dat.nxt$precipf <- rbind(dat.nxt$precipf[1,], dat.nxt$precipf)
      dat.nxt$swdown  <- rbind(dat.nxt$swdown [1,], dat.nxt$swdown )
      dat.nxt$lwdown  <- rbind(dat.nxt$lwdown [1,], dat.nxt$lwdown )
      dat.nxt$press   <- rbind(dat.nxt$press  [1,], dat.nxt$press  )
      dat.nxt$qair    <- rbind(dat.nxt$qair   [1,], dat.nxt$qair   )
      dat.nxt$wind    <- rbind(dat.nxt$wind   [1,], dat.nxt$wind   )
    }
    
    dat.ens <- list() # a new list for each ensemble member as a new layer
    df.hour <- data.frame(hour=0:23)
    
    # Create a list layer for each ensemble member
    for(e in ens.day){
      dat.ens[[paste0("X", e)]] <- data.frame(dataset      =dat.yr $tmax   $dataset,
                                              year         =dat.yr $tmax   $year,
                                              doy          =dat.yr $tmax   $doy,
                                              date         =dat.yr $tmax   $time,
                                              tmax.day     =dat.yr $tmax   [,paste0("X", e)],
                                              tmin.day     =dat.yr $tmin   [,paste0("X", e)],
                                              precipf.day  =dat.yr $precipf[,paste0("X", e)]/(60*60*24),
                                              swdown.day   =dat.yr $swdown [,paste0("X", e)],
                                              lwdown.day   =dat.yr $lwdown [,paste0("X", e)],
                                              press.day    =dat.yr $press  [,paste0("X", e)],
                                              qair.day     =dat.yr $qair   [,paste0("X", e)],
                                              wind.day     =dat.yr $wind   [,paste0("X", e)],
                                              next.tmax    =dat.nxt$tmax   [,paste0("X", e)],
                                              next.tmin    =dat.nxt$tmin   [,paste0("X", e)],
                                              next.precipf =dat.nxt$precipf[,paste0("X", e)]/(60*60*24),
                                              next.swdown  =dat.nxt$swdown [,paste0("X", e)],
                                              next.lwdown  =dat.nxt$lwdown [,paste0("X", e)],
                                              next.press   =dat.nxt$press  [,paste0("X", e)],
                                              next.qair    =dat.nxt$qair   [,paste0("X", e)],
                                              next.wind    =dat.nxt$wind   [,paste0("X", e)]
                                             )
      dat.ens[[paste0("X", e)]]$time.day <- as.numeric(difftime(dat.ens[[paste0("X", e)]]$date, "2016-01-01", tz="GMT", units="day"))
      dat.ens[[paste0("X", e)]] <- merge(dat.ens[[paste0("X", e)]], df.hour, all=T)
      
      dat.ens[[paste0("X", e)]]$date <- strptime(paste(dat.ens[[paste0("X", e)]]$year, dat.ens[[paste0("X", e)]]$doy+1, dat.ens[[paste0("X", e)]]$hour, sep="-"), "%Y-%j-%H", tz="GMT")
      dat.ens[[paste0("X", e)]]$time.hr <- as.numeric(difftime(dat.ens[[paste0("X", e)]]$date, "2016-01-01", tz="GMT", units="hour"))
      dat.ens[[paste0("X", e)]] <- dat.ens[[paste0("X", e)]][order(dat.ens[[paste0("X", e)]]$time.hr, decreasing=F),]
    }
    
    # Set up the time dimension for this year
    hrs.now <- as.numeric(difftime(dat.ens[[paste0("X", ens.day[1])]]$date, "0850-01-01", tz="GMT", units="hour"))
    dim.t <- ncdim_def(name = "time",
                       units = paste0("hours since 0850-01-01 00:00:00:"),
                       vals = hrs.now, # calculating the number of months in this run
                       calendar = "standard", unlim = TRUE)
    
    
    # -----------------------------------
    # 2. Predict met vars for each ensemble member
    # Note: Using a loop for each ensemble member for now, but this will get 
    #       parallelized to speed it up soon, but we'll prototype in parallel
    # -----------------------------------
    # cores.use <- min(cores.max, length(dat.ens))
    # ens.sims  <- mclapply(dat.ens, predict.subdaily, mc.cores=cores.use, n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.init=lags.init[[paste0("X", e)]], dat.train=dat.train)
    ens.sims <- list()
    for(e in unique(ens.day)){
      # # Do the prediction
      ens.sims[[paste0("X", e)]] <- predict.subdaily(dat.mod=dat.ens[[paste0("X", e)]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.init=lags.init[[paste0("X", e)]], dat.train=dat.train)
      # qair.max <- quantile(as.matrix(ens.sims[[paste0("X", e)]]$qair[,c(1:ens.hr)]), 0.99)
      # ens.sims[[paste0("X", e)]]$qair[ens.sims[[paste0("X", e)]]$qair>qair.max] <- qair.max

      # If this is one of our designated QAQC years, makes some graphs
      if(y %in% yrs.plot){
        day.name <- paste0(site.name, "_", GCM, "_1hr_", str_pad(e, 3, pad=0))
        fig.ens <- file.path(path.out, "subdaily_qaqc", day.name)
        if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
        
        for(v in names(ens.sims[[paste0("X", e)]])){
          graph.predict(dat.mod=dat.ens[[paste0("X", e)]], dat.ens=ens.sims[[paste0("X", e)]], var=v, fig.dir=fig.ens)
        }
      }
      
      
      # Update the initial lags for next year
      for(v in names(ens.sims[[paste0("X", e)]])){
        lags.init[[paste0("X",e)]][[v]] <- data.frame(ens.sims[[paste0("X", e)]][[v]][length(ens.sims[[paste0("X", e)]][[v]]),])
      }
      
      # -----------------------------------
      # Write each year for each ensemble member into its own .nc file
      # -----------------------------------
      for(i in 1:ens.hr){
        ens.name <- paste0(site.name, "_", GCM, "_1hr_", str_pad(e, 3, pad=0), "-", str_pad(i, 3, pad=0))
        
        if(!dir.exists(file.path(path.out, ens.name))) dir.create(file.path(path.out, ens.name), recursive=T)
        path.out <- file.path(dat.base, GCM, "1hr")
        
        var.list <- list()
        dat.list <- list()
        for(v in names(ens.sims[[paste0("X", e)]])){
          var.cf = vars.info[vars.info$name==v, "name.cf"]
          var.list[[v]] <- ncvar_def(v, units=paste(vars.info[vars.info$name==v, "units"]), dim=list(dimX, dimY, dim.t), longname=paste(vars.info[vars.info$name==v, "longname"]))
          dat.list[[v]] <- array(ens.sims[[paste0("X", e)]][[v]][,i], dim=c(1,1,length(hrs.now)))
        }
        
        # Naming convention: [SITE]_[GCM]_1hr_[bias_ens_member]-[subday_ens_member]_[YEAR].nc
        nc <- nc_create(file.path(path.out, ens.name, paste0(ens.name,"_", str_pad(y, 4, pad=0), ".nc")), var.list)
        for(v in 1:length(var.list)) {
          ncvar_put(nc, var.list[[v]], dat.list[[v]])
        }
        nc_close(nc)    
      }
      # -----------------------------------
      
    } # End ensemble member prediction for 1 year
    # -----------------------------------
  } # End Year Loop
  # -----------------------------------
  
  # Do some clean-up to save space
  # dir.compress <- dir(path.out, GCM)
  
  # setwd(path.out)
  # for(ens in dir.compress){
  #   system(paste0("tar -jcvf ", ens, ".tar.bz2 ", ens)) # Compress the folder
  #   system(paste0("rm -rf ", ens)) # remove the uncompressed folder
  # }
  # setwd(wd.base)
  toc()
} # End GCM loop
