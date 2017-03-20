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

wd.base <- "~/Dropbox/PalEON_CR/met_ensemble/"
# wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble/"
setwd(wd.base)

# Load the scripts that do all the heavy lifting
source("scripts/temporal_downscale.R")
source("scripts/temporal_downscale_functions.R")


# dat.base <- "/projectnb/dietzelab/paleon/met_ensemble/data/met_ensembles/HARVARD/"
dat.base <- "~/Desktop/Research/met_ensembles/data/met_ensembles/VCM/"
# dat.base <- "~/Desktop/met_bias_day/data/met_ensembles/HARVARD/"

# dat.train <- read.csv(file.path(wd.base, "data/paleon_sites/HARVARD/NLDAS_1980-2015.csv"))
dat.train <- read.csv(file.path(wd.base, "data/paleon_sites/VCM/Ameriflux_2007-2014.csv"))
# dat.train$hour <- dat.train$hour + minute(dat.train$date)/60

# aggregate to hour
dat.train <- aggregate(dat.train[,c("tair", "precipf", "swdown", "lwdown", "press", "qair", "uas", "vas", "wind")],
                       by=dat.train[,c("year", "doy", "hour")],
                       FUN=mean)
dat.train$date <- as.POSIXct(paste(dat.train$year, dat.train$doy, dat.train$hour, sep="-"), "%Y-%j-%H", tz="GMT")
summary(dat.train)


# Calculate the timestep in hours to make things easier
timestep <- unique(minute(dat.train$date)/60)+1 # add 1 so hourly data = 1 hr
if(length(timestep)>1) timestep <- mean(diff(timestep)) # if we have more than one minute mark, do the difference


df.hour <- data.frame(hour=unique(dat.train$hour)) # match this to whatever your "hourly" timestep is



# Hard-coding numbers for Harvard
site.name="VCM"
site.lat=35.89
site.lon=-106.53

# GCM.list = c("CCSM4", "MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1")
# GCM.list = "MIROC-ESM"
ens.hr  <- 3 # Number of hourly ensemble members to create
n.day <- 40 # Number of daily ensemble members to process
yrs.plot <- c(2015, 1985, 1920, 1901)
# years.sim=2015:1900
years.sim=NULL
cores.max = 8

# Set up the appropriate seed
set.seed(0017)
seed.vec <- sample.int(1e6, size=500, replace=F)

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

# for(GCM in GCM.list){
  GCM="Ameriflux"
  # tic()
  # Set the directory where the output is & load the file
  path.gcm <- file.path(dat.base, "day")
  # dat.day <- dir(path.gcm, ".Rdata")
  # load(file.path(path.gcm, dat.day)) # Loads dat.out.full

  # Set & create the output directory
  path.out <- file.path(dat.base, paste0(timestep, "hr2"))
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # -----------------------------------
  # 1. Format output so all ensemble members can be run at once
  # NOTE: Need to start with the last and work to the first
  # -----------------------------------
  # Figure out which daily ensembles we should pull from
  day.dirs <- dir(path.gcm,  paste0(site.name, "_", GCM, "_day_"))
  ens.day.all <- substr(day.dirs, nchar(day.dirs)-2, nchar(day.dirs))
  
  hrs.dir <- dir(file.path(path.out),  paste0(site.name, "_", GCM, "_1hr_"))
  ens.day.done <- unique(substr(hrs.dir, nchar(hrs.dir)-6, nchar(hrs.dir)-4))
  
  ens.list <- ens.day.all[!(ens.day.all %in% ens.day.done)]
  
  seed <- seed.vec[length(hrs.dir)+1] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed
  set.seed(seed)
  ens.day <- sample(ens.list, n.day, replace=F) # For now, randomly choose which ensemble members to downscale

  # If we haven't specified a years subset, extract it from the available files
  if(is.null(years.sim)){
    files.avail <- dir(file.path(path.gcm, day.dirs[1]))
    years.files <- as.numeric(substr(files.avail,nchar(files.avail)-6, nchar(files.avail)-3))
    years.sim <- max(years.files):min(years.files)
  }
  
  # Initialize the lags
  lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally
  for(e in ens.day){
    path.day <- file.path(path.gcm,  paste0(site.name, "_", GCM, "_day_",e))
    file.ens <- dir(path.day, paste(str_pad(years.sim[1], 4, pad=0)))
    
    nc.now <- nc_open(file.path(path.day, file.ens))
    # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
    nc.time <- ncvar_get(nc.now, "time")
    tmax.init    <- ncvar_get(nc.now, "tmax")[length(nc.time)]
    tmin.init    <- ncvar_get(nc.now, "tmin")[length(nc.time)]
    precipf.init <- ncvar_get(nc.now, "precipf")[length(nc.time)]/(60*60*24)
    swdown.init  <- ncvar_get(nc.now, "swdown")[length(nc.time)]
    lwdown.init  <- ncvar_get(nc.now, "lwdown")[length(nc.time)]
    press.init   <- ncvar_get(nc.now, "press")[length(nc.time)]
    qair.init    <- ncvar_get(nc.now, "qair")[length(nc.time)]
    wind.init    <- ncvar_get(nc.now, "wind")[length(nc.time)]
    nc_close(nc.now)
    
    lags.init[[paste0("X", e)]][["tair"   ]] <- data.frame(array(mean(c(tmax.init, tmin.init)), dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["tmax"   ]] <- data.frame(array(tmax.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["tmin"   ]] <- data.frame(array(tmin.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["precipf"]] <- data.frame(array(precipf.init, dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["swdown" ]] <- data.frame(array(swdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["lwdown" ]] <- data.frame(array(lwdown.init , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["press"  ]] <- data.frame(array(press.init  , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["qair"   ]] <- data.frame(array(qair.init   , dim=c(1, ens.hr)))
    lags.init[[paste0("X", e)]][["wind"   ]] <- data.frame(array(wind.init   , dim=c(1, ens.hr)))
  }

  # Initialize a progress bar
  pb.index=1
  pb <- txtProgressBar(min=1, max=length(years.sim)*length(ens.day), style=3)
  
  for(y in years.sim){
    dat.ens <- list() # a new list for each ensemble member as a new layer
    ens.sims <- list()
    # Create a list layer for each ensemble member
    # NOTE: NEED TO CHECK THE TIME STAMPS HERE TO MAKE SURE WE'RE PULLING AND ADDING IN THE RIGHT PLACE
    for(e in ens.day){
      setTxtProgressBar(pb, pb.index) #update progress bar
      
      path.day <- file.path(path.gcm,  paste0(site.name, "_", GCM, "_day_",e))
      file.ens <- dir(path.day, str_pad(y, 4, pad=0))
      
      nc.now <- nc_open(file.path(path.day, file.ens))
      dat.yr <- data.frame(time    = ncvar_get(nc.now, "time"   ),
                           tmax    = ncvar_get(nc.now, "tmax"   ),
                           tmin    = ncvar_get(nc.now, "tmin"   ),
                           precipf = ncvar_get(nc.now, "precipf")/(60*60*24),
                           swdown  = ncvar_get(nc.now, "swdown" ),
                           lwdown  = ncvar_get(nc.now, "lwdown" ),
                           press   = ncvar_get(nc.now, "press"  ),
                           qair    = ncvar_get(nc.now, "qair"   ),
                           wind    = ncvar_get(nc.now, "wind"   )
                           )
      nc_close(nc.now)
      
      # Do some stuff to get the right time variables for dat.yr
      dat.yr$year <- y
      dat.yr$date <- as.Date(dat.yr$time, origin="850-01-01")
      dat.yr$doy <- yday(dat.yr$date)
      
      # Create the data frame for the "next" values
      dat.nxt <- dat.yr
      # Shift everyting up by a day to get the preview of the next day to get processed
      # Note: Because we work backwards through time, Jan 1 is the "next" day to get processed with Jan 2
      dat.nxt[2:(nrow(dat.nxt)), c("tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind")] <- dat.nxt[1:(nrow(dat.nxt)-1), c("tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
      
      # Need to add in the "next" value for january (dec 31 of next year) 
      # Note: if we're past the end of our daily data, the best we can do is leave things as is (copy the last day's value)
      if(y>min(years.sim)){
        file.nxt <- dir(path.day, str_pad(y-1, 4, pad=0))
        
        nc.nxt <- nc_open(file.path(path.day, file.nxt))
        # tair.init    <- dat.out.full$met.bias[dat.out.full$met.bias$year==max(years.sim) & dat.out.full$met.bias$doy==364,"tair"]
        nxt.time <- ncvar_get(nc.nxt, "time")
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmax"   ] <- ncvar_get(nc.nxt, "tmax")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"tmin"   ] <- ncvar_get(nc.nxt, "tmin")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"precipf"] <- ncvar_get(nc.nxt, "precipf")[length(nxt.time)]/(60*60*24)
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"swdown" ] <- ncvar_get(nc.nxt, "swdown")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"lwdown" ] <- ncvar_get(nc.nxt, "lwdown")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"press"  ] <- ncvar_get(nc.nxt, "press")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"qair"   ] <- ncvar_get(nc.nxt, "qair")[length(nxt.time)]
        dat.nxt[dat.nxt$time==min(dat.nxt$time),"wind"   ] <- ncvar_get(nc.nxt, "wind")[length(nxt.time)]
        
        nc_close(nc.nxt)
      } 
      
      dat.ens[[paste0("X", e)]] <- data.frame(ens.day      =as.factor(paste0("X", e)),
                                              year         =dat.yr $year   ,
                                              doy          =dat.yr $doy    ,
                                              date         =dat.yr $date   ,
                                              tmax.day     =dat.yr $tmax   ,
                                              tmin.day     =dat.yr $tmin   ,
                                              precipf.day  =dat.yr $precipf,
                                              swdown.day   =dat.yr $swdown ,
                                              lwdown.day   =dat.yr $lwdown ,
                                              press.day    =dat.yr $press  ,
                                              qair.day     =dat.yr $qair   ,
                                              wind.day     =dat.yr $wind   ,
                                              next.tmax    =dat.nxt$tmax   ,
                                              next.tmin    =dat.nxt$tmin   ,
                                              next.precipf =dat.nxt$precipf,
                                              next.swdown  =dat.nxt$swdown ,
                                              next.lwdown  =dat.nxt$lwdown ,
                                              next.press   =dat.nxt$press  ,
                                              next.qair    =dat.nxt$qair   ,
                                              next.wind    =dat.nxt$wind
                                             )
      
      dat.ens[[paste0("X", e)]]$time.day <- as.numeric(difftime(dat.ens[[paste0("X", e)]]$date, "2016-01-01", tz="GMT", units="day"))
      dat.ens[[paste0("X", e)]] <- merge(dat.ens[[paste0("X", e)]], df.hour, all=T)
      
      # A very ugly hack way to get the minutes in there
      dat.ens[[paste0("X", e)]]$minute <- abs(dat.ens[[paste0("X", e)]]$hour-(round(dat.ens[[paste0("X", e)]]$hour, 0)))*60
      
      # timestep <- max(dat.ens[[paste0("X", e)]]$minute)/60 # Calculate the timestep to make the date function work
      # dat.ens[[paste0("X", e)]]$date <- strptime(paste(dat.ens[[paste0("X", e)]]$year, dat.ens[[paste0("X", e)]]$doy, round(dat.ens[[paste0("X", e)]]$hour-timestep/2),dat.ens[[paste0("X", e)]]$minute, sep="-"), "%Y-%j-%H-%M", tz="GMT")
      dat.ens[[paste0("X", e)]]$date <- strptime(paste(dat.ens[[paste0("X", e)]]$year, dat.ens[[paste0("X", e)]]$doy, dat.ens[[paste0("X", e)]]$hour,dat.ens[[paste0("X", e)]]$minute, sep="-"), "%Y-%j-%H-%M", tz="GMT")
      dat.ens[[paste0("X", e)]]$time.hr <- as.numeric(difftime(dat.ens[[paste0("X", e)]]$date, "2016-01-01", tz="GMT", units="hour")) #+ minute(dat.train$date)/60
      dat.ens[[paste0("X", e)]] <- dat.ens[[paste0("X", e)]][order(dat.ens[[paste0("X", e)]]$time.hr, decreasing=F),]
      
      # Do the modeling here rather than in parallel
      # ens.sims[[paste0("X", e)]] <- predict.subdaily(dat.ens[[paste0("X", e)]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
      
      pb.index <- pb.index + 1 # Advance our progress bar  # Doing it here rather than at end of writing files because that should be fast
    } # End ensembles setup
    
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
    cores.use <- min(cores.max, length(dat.ens))
    ens.sims  <- mclapply(dat.ens, predict.subdaily, mc.cores=cores.use, n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
    
    # ens.sims <- predict.subdaily(dat.ens[[1]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.init=lags.init[[1]], dat.train=dat.train)
    
    # If this is one of our designated QAQC years, makes some graphs
    # Now doing this for the whole GCM
    if(y %in% yrs.plot){
      dat.plot <- data.frame()
      ens.plot <- list()
      for(i in names(ens.sims[[1]])){
        ens.plot[[i]] <- data.frame(matrix(nrow=nrow(ens.sims[[1]][[i]]), ncol=0))
      }
      for(e in names(ens.sims)){
        dat.plot <- rbind(dat.plot, dat.ens[[e]])
        for(i in names(ens.sims[[e]])){
          ens.plot[[i]] <- cbind(ens.plot[[i]], ens.sims[[e]][[i]])
        }
      }
      day.name <- paste0(site.name, "_", GCM, "_", timestep, "hr")
      fig.ens <- file.path(path.out, "subdaily_qaqc", day.name)
      if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
      # dat.ens.full <- dat.ens
      for(v in names(ens.plot)){
        graph.predict(dat.mod=dat.plot, dat.ens=ens.plot, var=v, yr=y, fig.dir=fig.ens)
      }
    }
    
    
    # ens.sims <- list()
    for(e in unique(ens.day)){
      # # Do the prediction
      # ens.sims[[paste0("X", e)]] <- predict.subdaily(dat.mod=dat.ens[[paste0("X", e)]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.list=lags.init, lags.init=NULL, dat.train=dat.train)
      # qair.max <- quantile(as.matrix(ens.sims[[paste0("X", e)]]$qair[,c(1:ens.hr)]), 0.99)
      # ens.sims[[paste0("X", e)]]$qair[ens.sims[[paste0("X", e)]]$qair>qair.max] <- qair.max

      # Update the initial lags for next year
      for(v in names(ens.sims[[paste0("X", e)]])){
        lags.init[[paste0("X",e)]][[v]] <- data.frame(ens.sims[[paste0("X", e)]][[v]][length(ens.sims[[paste0("X", e)]][[v]]),])
      }
      
      # -----------------------------------
      # Write each year for each ensemble member into its own .nc file
      # -----------------------------------
      for(i in 1:ens.hr){
        ens.name <- paste0(site.name, "_", GCM, "_", timestep, "hr", "_", str_pad(e, 3, pad=0), "-", str_pad(i, 3, pad=0))
        
        if(!dir.exists(file.path(path.out, ens.name))) dir.create(file.path(path.out, ens.name), recursive=T)
        # path.out <- file.path(dat.base, GCM, "1hr")
        
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
# } # End GCM loop
