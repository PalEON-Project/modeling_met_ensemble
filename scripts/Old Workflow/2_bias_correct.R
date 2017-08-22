# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Bias-correction to create a smooth daily met product from multiple sources 
#          of varying temporal frequencies and extents
# Creator: Christy Rollinson, 1 July 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Bias-correct raw met data & take monthly variables to daily time step
# The end state of this script is continuous, smoothly daily output from 850-2010+ 
# that can be used as is or fed into the day to subday script to get hourly drivers
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Set up file structure, etc.
# 1. Read in & format the different datasets
#    - aggregate all variable to daily for the bias-correction
#    - do some exploratory graphing along the way
# 2. Define training window for all datasets
#    - defines met.train (e.g. nldas data) and met.raw (e.g. CRUNCEP)
#    - make sure either years line up (best) or to calculate the met means
# 3. Doing the DOY bias-corrections : 
#     - Met.train ~ s(doy) + Met.Raw
#     - I *think* this can be used to get monthly to daily data as well
#     ** Don't forget to use the residuals to add in the "error" in the model
#        -- assuming all goes well this will just be rnorm(mean=mean(resid), sd=sd(resid))
#     - These residuals will be what's added back into the bias-correction step
# 4. Generate an ensemble of the predicted mean day of year
#     - use covariance matrix and make array of n ensemble members & run raw data
# 5. Add in anomaly
#    - use anomaly from original climate mean rather than residual from predicted model
# 6. Save an .Rdata file with an ensemble of bias-corrected daily data
#    - this ensemble file can either be read into the temporal downscaling or fed into
#      a script that writes to netcdf format for models that take daily drivers
#
# NOTE: In the GCMs, short- and longwave radiation are generally only available at the 
#       monthly timestep.  Do the month -> Day downscaling of these LAST leveraging the
#       covariances with other met variables
# -----------------------------------

rm(list=ls())

# -----------------------------------
# 0. Set up file structure, load packages, etc
# -----------------------------------
# Load libraries
library(ncdf4)
library(mgcv); library(ggplot2)
library(stringr)
library(lubridate)

# Set the working directory
wd.base <- "~/Dropbox/PalEON_CR/met_ensemble/"
out.base <- "~/Desktop/met_ensembles/"
setwd(wd.base)

# Setting some important file paths
path.pecan <- "~/Desktop/Research/pecan"

# Defining a site name -- this can go into a function later
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18
GCM.list=c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
# GCM.list=NULL
# GCM.list=c("bcc-csm1-1")
LDAS="NLDAS"
# ens=1:50
ens=1:10
n=length(ens)
# n=10 # Number of ensemble members
# Make some vectors with met sets to make things easier
# met.subday <- c(LDAS, "CRUNCEP")
met.subday <- c("CRUNCEP", LDAS, "Ameriflux")



# Set up the appropriate seed
set.seed(1159)
seed.vec <- sample.int(1e6, size=500, replace=F)
seed <- seed.vec[min(ens)] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed

# -----------------------------------


# -----------------------------------
# 1. generate a daily training dataset to get us started
# 
# this will end up being a 1-member "ensemble" of the NLDAS dataset
# -----------------------------------

df.var <- data.frame(CF.name = c("air_temperature", "air_temperature_maximum", "air_temperature_minimum", 
                                 "surface_downwelling_longwave_flux_in_air",
                                 "air_pressure", "surface_downwelling_shortwave_flux_in_air", 
                                 "eastward_wind", "northward_wind", "wind_speed", "specific_humidity", "precipitation_flux"), 
                     units = c("Kelvin", "Kelvin", "Kelvin", "W/m2", "Pascal", "W/m2", "m/s", "m/s", "m/s", "g/g", "kg/m2/s"))

nc.info <- data.frame(CF.name = c("air_temperature_minimum", "air_temperature_maximum", "precipitation_flux", 
        "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air", 
        "air_pressure", "specific_humidity", "wind_speed"), 
        longname = c("2 meter minimum air temperature", "2 meter maximum air temperature", 
        "cumulative precipitation (water equivalent)", "incident (downwelling) showtwave radiation", 
        "incident (downwelling) longwave radiation", "air_pressureure at the surface", 
        "Specific humidity measured at the lowest level of the atmosphere", 
        "Wind speed"), 
        units = c("K", "K", "kg m-2 s-1", "W m-2", "W m-2", "Pa", 
        "kg kg-1", "m s-1"))

path.ldas <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/NLDAS/"
files.train <- dir(path.ldas)

outfolder <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/NLDAS_day/"
dir.create(outfolder, recursive=T)
for(i in 1:length(files.train)){
	
	# Figure out what year we're working with
	yr.now <- as.numeric(strsplit(files.train[i], "[.]")[[1]][2])
	nday <- ifelse(leap_year(yr.now), 366, 365)

	dat.day <- list()
	
	# Open the file so we can query from it
	ncT <- nc_open(file.path(path.ldas, files.train[i]))
	
	# Extract som plot dimensions
	lat.nc <- ncvar_get(ncT, "latitude")
	lon.nc <- ncvar_get(ncT, "longitude")
	
	time.nc <- ncvar_get(ncT, "time") 
	time.day <- apply(matrix(time.nc, ncol=nday), 2, mean) # get the daily time stamps
	
	# Extract plot info & aggregate to daily resolution
	for(v in names(ncT$var)){
		dat.hr <- matrix(ncvar_get(ncT, v), ncol=nday)
		if(v == "air_temperature"){
			dat.day[["air_temperature_minimum"]] <- apply(dat.hr, 2, min)
			dat.day[["air_temperature_maximum"]] <- apply(dat.hr, 2, max)
		} else if(v %in% c("eastward_wind", "northward_wind")) {
			wind.e <- matrix(ncvar_get(ncT, "eastward_wind"), ncol=nday)
			wind.n <- matrix(ncvar_get(ncT, "northward_wind"), ncol=nday)
			wind <- sqrt(wind.e^2 + wind.n^2)
			dat.day[["wind_speed"]] <- apply(wind, 2, mean)
		} else {
			dat.day[[v]] <- apply(dat.hr, 2, mean)
		}
	}
	
	# Create a daily .nc file for each year
	dim.lat <- ncdim_def(name='latitude', units='degree_north', vals=lat.nc, create_dimvar=TRUE)
    dim.lon <- ncdim_def(name='longitude', units='degree_east', vals=lon.nc, create_dimvar=TRUE)
    dim.time <- ncdim_def(name='time', units="sec", vals=time.day, create_dimvar=TRUE, unlim=TRUE)
    nc.dim=list(dim.lat,dim.lon,dim.time)

    var.list = list()
    for(v in names(dat.day)){
      var.list[[v]] = ncvar_def(name=v, units=as.character(nc.info[nc.info$CF.name==v, "units"]), dim=nc.dim, missval=-999, verbose=verbose)
    }

	loc.file <- file.path(outfolder, paste("NLDAS_day", str_pad(yr.now, width=4, side="left",  pad="0"), "nc", sep = "."))
    loc <- nc_create(filename = loc.file, vars = var.list, verbose = verbose)

    for (v in names(dat.day)) {
    	ncvar_put(nc = loc, varid = as.character(v), vals = dat.day[[v]])
    }
    nc_close(loc)	
}

# -----------------------------------

# -----------------------------------
# Working through the datasets to be bias-corrected
# General Workflow
# 1. Align Data:
#    1.1. Common temporal resolution 
#         - match training dataset (daily)
#    1.2. Pair Ensemble members
#         A. n train = n out --> pair by member ID
#         B. n train < n out -->
#            --> if n out is multiple of n train --> save multiple out for each in
#            --> if n out not multiple of n train --> randomly pick which members to save multiples for
#         C. n train > n out --> randomly subset members
# 2. Debias Met 
#    2.1. Climatology Debias: Annual means + seasonal cycle
#    2.2. Anomaly Debias: adjust anomaly variance/distribution to try to maintain variance
# 3. Save Met
#    - write years of debiased (or training) met to keep with new ensemble & member ID
# -----------------------------------
GCM.list
n.ens = 10
ens.ID = "TEST"

# 1. Align CRU 6-hourly with LDAS daily
source(file.path(path.pecan, "modules/data.atmosphere/R", "align_met.R"))
train.path <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/NLDAS_day"
source.path <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/CRUNCEP"

# For first round, we only want a single in & out
met.out <- align.met(train.path, source.path, n.ens=1, seed=201708, pair.mems = FALSE)

# Calculate wind speed if it's not already there
if(!"wind_speed" %in% names(met.out$dat.source)){
	met.out$dat.source$wind_speed <- sqrt(met.out$dat.source$eastward_wind^2 + met.out$dat.source$northward_wind^2)
}

# 2. Pass the training & source met data into the bias-correction functions; this will get written to the ensemble
source(file.path(path.pecan, "modules/data.atmosphere/R", "debias_met_regression.R"))
debias.met.regression(train.data=met.out$dat.train, source.data=met.out$dat.source, n.ens=10, vars.debias=NULL, CRUNCEP=TRUE,
                      pair.anoms = TRUE, pair.ens = FALSE, uncert.prop="mean", resids = FALSE, seed=Sys.Date(),
                      outfolder="~/Desktop/Research/met_ensembles/data/met_ensembles/TEST", 
                      yrs.save=1901:1979, ens.name="TEST2", ens.mems=NULL, lat.in=site.lat, lon.in=site.lon,
                      save.diagnostics=TRUE, path.diagnostics="~/Desktop/Research/met_ensembles/data/met_ensembles/TEST/bias_correct_qaqc",
                      parallel = FALSE, n.cores = NULL, overwrite = TRUE, verbose = FALSE) 



# -----------------------------------















# source("scripts/debias.gcm.R")
# debias.gcm(GCM=GCM.list[1], LDAS=LDAS, wd.base=wd.base, out.base=out.base, site.name=site.name, site.lat=site.lat, site.lon=site.lon, n=n)


for(GCM in GCM.list){
  
  print(GCM)
  path.dat <- file.path(wd.base, "data/paleon_sites", site.name)
  path.out <- file.path(out.base, "data/met_ensembles", site.name, GCM, "day")
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)  
  
  met.done <- dir(path.out, ".csv")
  # met.all <- c(paste0(GCM, ".hist"), paste0(GCM, ".p1000"), "CRUNCEP", LDAS, "Ameriflux")
  # met.sets <- c("CRUNCEP", LDAS, "Ameriflux")
  met.sets <- c(paste0(GCM, ".p1000"), paste0(GCM, ".hist"), "CRUNCEP", LDAS)
  
  # -----------------------------------
  # 1. Read in & format the different datasets
  # -----------------------------------
  
  # Find the appropriate file name for each
  # file.flux <- dir(path.dat, "Ameriflux")
  file.ldas <- dir(path.dat, LDAS)
  file.cru <- dir(path.dat, "CRUNCEP")
  file.hist <- dir(path.dat, paste0(GCM, "_historical"))
  file.p1k <- dir(path.dat, paste0(GCM, "_p1000"))
  
  # flux     <- read.csv(file.path(path.dat, file.flux))
  ldas     <- read.csv(file.path(path.dat, file.ldas))
  cruncep  <- read.csv(file.path(path.dat, file.cru))
  gcm.hist <- read.csv(file.path(path.dat, file.hist))
  gcm.p1k  <- read.csv(file.path(path.dat, file.p1k))
  
  # Adding an hour field to the gcm; setting as noon (middle of window) for simplicity
  gcm.p1k$hour  <- 12.00
  gcm.hist$hour <- 12.00
  
  # Adding minute stamps to everything
  ldas$minute <- 30
  cruncep$minute <- 30
  gcm.hist$minute <- 30
  gcm.p1k$minute <- 30
  
  # Adding date stamps to datasets
  ldas$date <- paste0(ldas$year, "-", ldas$doy+1, " ", ldas$hour, ":", ldas$minute)
  ldas$date <- strptime(ldas$date, format="%Y-%j %H:%M")
  summary(ldas)
  
  cruncep$date <- paste0(cruncep$year, "-", cruncep$doy+1, " ", cruncep$hour, ":", cruncep$minute)
  cruncep$date <- strptime(cruncep$date, format="%Y-%j %H:%M")
  summary(cruncep)
  
  gcm.hist$date <- paste0(gcm.hist$year, "-", gcm.hist$doy+1, " ", gcm.hist$hour, ":", gcm.hist$minute)
  gcm.hist$date <- strptime(gcm.hist$date, format="%Y-%j %H:%M")
  summary(gcm.hist)
  
  gcm.p1k$date <- paste0(gcm.p1k$year, "-", gcm.p1k$doy+1, " ", gcm.p1k$hour, ":", gcm.p1k$minute)
  gcm.p1k$date <- strptime(gcm.p1k$date, format="%Y-%j %H:%M")
  summary(gcm.p1k)
  
  # # right now the ameriflux .csv doesn't match the others (script updated, but not re-run)
  # flux$dataset <- as.factor("Ameriflux")
  # flux$minute <- minute(flux$date)
  # summary(flux)
  
  # flux <- flux[,names(ldas)]
  
  
  # making sure all datasets have tair, tmax, and tmin
  # flux$tmax    <- flux$tmin    <- flux$tair
  ldas$tmax    <- ldas$tmin    <- ldas$tair
  cruncep$tmax <- cruncep$tmin <- cruncep$tair
  gcm.p1k$tair  <- apply(gcm.p1k [,c("tmax", "tmin")], 1, FUN=mean)
  gcm.hist$tair <- apply(gcm.hist[,c("tmax", "tmin")], 1, FUN=mean)
  
  # ******* TYPO CORRECTIONS!! *******
  # **** NOTE: WILL NEED TO CANCEL THIS OUT UPON NEW EXTRACTION
  # make first day of year start on 0
  # ldas$doy    <- ldas$doy - 1
  # cruncep$doy <- cruncep$doy -1
  # ldas$precipf <- ldas$precipf*.1
  # flux$doy <- flux$doy-1
  
  # summary(flux)
  summary(ldas)
  summary(cruncep)
  summary(gcm.p1k)
  summary(gcm.hist)
  
  
  
  # vars.met <- c("tair", "tmax", "tmin", "precipf", "press", "qair", "wind", "swdown", "lwdown")
  vars.met <- c("tair", "tmax", "tmin", "qair", "swdown", "press", "lwdown", "wind", "precipf")
  cols.bind <- c("dataset", "year", "doy", "hour", "minute", vars.met)
  # met.all <- rbind(flux[,cols.bind], ldas[,cols.bind], cruncep[,cols.bind])
  met.all <- rbind(ldas[,cols.bind], cruncep[,cols.bind], gcm.hist[,cols.bind], gcm.p1k[,cols.bind])
  met.all$Date <- as.Date(met.all$doy, origin=as.Date(paste(met.all$year, "01", "01", sep="-")))
  # met.all$Date <- as.Date(paste0(met.all$doy, " ", met.all$hour, ":", met.all$minute), format="%j %H:%M", origin=as.Date(paste(met.all$year, "01", "01", sep="-")))
  # met.all$Date2 <- as.POSIXct(paste0(met.all$Date, " ", met.all$hour, ":", met.all$minute), format="%Y-%m-%d %H:%M")
  summary(met.all)
  # ----------------
  # Creating a daily met record
  # ----------------
  met.day <- aggregate(met.all[,vars.met], by=met.all[,c("dataset", "year", "doy")], FUN=mean)
  
  # getting tmax & tmin for ldas & cru
  met.day2 <- aggregate(met.all[met.all$dataset %in% met.subday,"tair"], 
                        by=met.all[met.all$dataset %in% met.subday,c("dataset", "year", "doy")], 
                        FUN=max)
  names(met.day2)[4] <- "tmax"
  met.day2$tmin <- aggregate(met.all[met.all$dataset %in% met.subday,"tair"], 
                             by=met.all[met.all$dataset %in% met.subday,c("dataset", "year", "doy")], 
                             FUN=min)[,4]
  summary(met.day2)
  
  # merging tmax & tmin back into met.day 
  met.day[met.day$dataset %in% met.subday, "tmax"] <- aggregate(met.all[met.all$dataset %in%met.subday,"tair"], 
                                                                by=met.all[met.all$dataset %in% met.subday,c("dataset", "year", "doy")], 
                                                                FUN=max)[,4]
  met.day[met.day$dataset %in% met.subday, "tmin"] <- aggregate(met.all[met.all$dataset %in% met.subday,"tair"], 
                                                                by=met.all[met.all$dataset %in% met.subday,c("dataset", "year", "doy")], 
                                                                FUN=min)[,4]
  summary(met.day)
  
  
  # Make precipf in total kg/m2/day rather than per s to make life easier
  met.day$precipf <- met.day$precipf*60*60*24
  summary(met.day)
  # ----------------
  
  # ----------------
  # Aggregating to year and doy in long format for exploratory
  # ----------------
  # Getting a scaffold to build off of
  met.year0 <- aggregate(met.day[,vars.met], by=met.day[,c("dataset", "year")], FUN=mean)
  met.doy0  <- aggregate(met.day[,vars.met], by=met.day[,c("dataset", "doy")], FUN=mean)
  
  # Stacking to long format
  met.year <- stack(met.year0[,vars.met])
  names(met.year) <- c("value", "met")
  met.year[,c("dataset", "year")] <- met.year0[,c("dataset", "year")]
  met.year$lwr <- stack(aggregate(met.day[,vars.met], by=met.day[,c("dataset", "year")], FUN=quantile, 0.025)[,vars.met])[,1]
  met.year$upr <- stack(aggregate(met.day[,vars.met], by=met.day[,c("dataset", "year")], FUN=quantile, 0.975)[,vars.met])[,1]
  summary(met.year)
  
  met.doy <- stack(met.doy0[,vars.met])
  names(met.doy) <- c("value", "met")
  met.doy[,c("dataset", "doy")] <- met.doy0[,c("dataset", "doy")]
  met.doy$lwr <- stack(aggregate(met.day[,vars.met], by=met.day[,c("dataset", "doy")], FUN=quantile, 0.025)[,vars.met])[,1]
  met.doy$upr <- stack(aggregate(met.day[,vars.met], by=met.day[,c("dataset", "doy")], FUN=quantile, 0.975)[,vars.met])[,1]
  summary(met.doy)
  
  met.year$dataset <- factor(met.year$dataset, levels=met.sets)
  met.doy $dataset <- factor(met.doy $dataset, levels=met.sets)
  met.year$met     <- factor(met.year$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
  met.doy $met     <- factor(met.doy $met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
  
  # Putting in a dummy upper bound for precip to deal with some of the real oddballs
  precip.cutoff <- quantile(met.doy[met.doy$met=="precipf","upr"], 0.95)
  met.doy$upr2 <- ifelse(met.doy$met=="precipf" & met.doy$upr>=precip.cutoff, precip.cutoff, met.doy$upr)
  
  # GCM="Ameriflux"
  
  png(file.path(path.out, paste0(GCM, "_Raw_Year_0850-2015.png")), height=11, width=8.5, "in", res=180)
  print(
    ggplot(data=met.year) +
      facet_grid(met~., scales="free_y") +
      geom_line(aes(x=year, y=value, color=dataset)) +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(name="Annual Mean Value") +
      theme_bw() +
      theme(legend.position="top",
            legend.direction="horizontal")
  )
  dev.off()
  
  png(file.path(path.out, paste0(GCM, "_Raw_DOY_All.png")), height=11, width=8.5, "in", res=180)
  print(
    ggplot(data=met.doy) +
      facet_grid(met~., scales="free_y") +
      geom_ribbon(aes(x=doy, ymin=lwr, ymax=upr2, fill=dataset), alpha=0.3) +
      geom_line(aes(x=doy, y=value, color=dataset)) +
      scale_x_continuous(expand=c(0,0), name="Day of Year (Julian)") +
      scale_y_continuous(name="Daily Mean Value") +
      theme_bw() +
      theme(legend.position="top",
            legend.direction="horizontal")
  )
  dev.off()
  # ----------------
  # -----------------------------------
  
  # -----------------------------------
  # 3. Doing the DOY bias-corrections : 
  # -----------------------------------
  # bias.correct <- function(met.day, met.var, dat.train, yrs.cal, n){
  
  # The met vars we need (in the order we want to do them)
  # vars.met <- c("tmax", "tmin", "swdown", "lwdown", "precipf", "qair", "press", "wind")
  vars.met <- c("tmax", "tmin", "qair", "swdown", "press", "lwdown", "wind", "precipf")
  dat.cal = LDAS
  
  # Note This dataframe is only for the datasets to be bias-corrected!
  # yrs.cal = data.frame(dataset = c(LDAS, "CRUNCEP"),
  #                      cal.min = c(2007, 1980),
  #                      cal.max = c(2014, 2010)
  #                      )
  yrs.cal = data.frame(dataset = c("CRUNCEP", paste0(GCM, ".hist"), paste0(GCM, ".p1000")),
                       cal.min = c(1980,             1901,              1829),
                       cal.max = c(2010,             1921,              1849)
  )
  summary(met.all)
  # yrs.cal = yrs.cal[1,]
  
  # --------------------------
  # Looking at the covariance among different variables
  # Very very clunky way of plotting the linear model of several pairs
  # --------------------------
  {
    pdf(file.path(path.out, paste0(GCM, "_MetCovariance.pdf")))
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=tmin, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs Tmin") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=swdown, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs swdown") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=lwdown, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs lwdown") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=precipf, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs precipf") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=press, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs press") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=qair, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs Qair") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=tmax, y=wind, color=dataset, fill=dataset), method="lm") +
        ggtitle("Tmax vs Wind") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=swdown, y=precipf, color=dataset, fill=dataset), method="lm") +
        ggtitle("Swdown vs Precipf") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=precipf, y=qair, color=dataset, fill=dataset), method="lm") +
        ggtitle("Precipf vs Qair") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=precipf, y=press, color=dataset, fill=dataset), method="lm") +
        ggtitle("Precipf vs Press") +
        theme_bw()
    )
    print(
      ggplot(data=met.day) +
        stat_smooth(aes(x=precipf, y=wind, color=dataset, fill=dataset), method="lm") +
        ggtitle("Precipf vs Wind") +
        theme_bw()
    )
    dev.off()
  }
  
  source("scripts/bias_correct_day_ensloop.R")
  met.bias <- met.day
  dat.out.full <- bias.correct(met.bias=met.bias, vars.met=vars.met, dat.train=LDAS, GCM=GCM, yrs.cal=yrs.cal, n=n, path.out=path.out, seed=seed)
  # --------------------
  
  
  # --------------------
  # Looking at the new bias-corrected time series! 
  # --------------------
  library(grid)
  for(met.var in vars.met){
    print(met.var)
    dat.out <- dat.out.full[[met.var]]
    dat.yr.bias <- aggregate(dat.out$ci[,c("X", "anom.raw", "mean", "lwr", "upr")],
                             by=dat.out$ci[,c("dataset", "met", "year")],
                             FUN=mean)
    dat.yr.bias$dataset <- factor(dat.yr.bias$dataset, levels=met.sets)
    # dat.yr.bias$met     <- factor(dat.yr.bias$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
    summary(dat.yr.bias)
    summary(met.year[met.year$met==met.var,]) #raw annual data
    
    LDAS.use.yr <- dat.yr.bias[dat.yr.bias$dataset==LDAS,]
    CRU.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="CRUNCEP" & dat.yr.bias$year<min(LDAS.use.yr$year),]
    GCM.hist.use.yr <- dat.yr.bias[dat.yr.bias$dataset==paste0(GCM, ".hist") & dat.yr.bias$year<min(CRU.use.yr$year),]
    GCM.p1k.use.yr <- dat.yr.bias[dat.yr.bias$dataset==paste0(GCM, ".p1000") & dat.yr.bias$year<min(GCM.hist.use.yr$year),]
    met.final.yr <- rbind(LDAS.use.yr, CRU.use.yr, GCM.hist.use.yr, GCM.p1k.use.yr)
    
    # Original 
    summary(met.year)
    plot.orig <-   ggplot(data=met.final.yr[met.final.yr$met==met.var,]) +
      facet_grid(met~., scales="free_y") +
      geom_line(aes(x=year, y=X, color=dataset)) +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Annual Mean Value") +
      ggtitle("Raw (Annual)") +
      theme_bw() +
      theme(legend.position="top",
            legend.direction="horizontal")  
    
    plot.bias <-   ggplot(data=met.final.yr[met.final.yr$met==met.var,]) +
      facet_grid(met~., scales="free_y") +
      geom_ribbon(aes(x=year, ymin=lwr, ymax=upr, fill=dataset), alpha=0.3) +
      geom_line(aes(x=year, y=mean, color=dataset)) +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Annual Mean Value") +
      ggtitle("Bias-Corrected (Annual)") +
      theme_bw() +
      theme(legend.position="top",
            legend.direction="horizontal")  
    
    LDAS.use <- dat.out$ci[dat.out$ci$dataset==LDAS,]
    CRU.use <- dat.out$ci[dat.out$ci$dataset=="CRUNCEP" & dat.out$ci$year<min(LDAS.use$year),]
    GCM.hist.use <- dat.out$ci[dat.out$ci$dataset==paste0(GCM, ".hist") & dat.out$ci$year<min(CRU.use$year),]
    GCM.p1k.use <- dat.out$ci[dat.out$ci$dataset==paste0(GCM, ".p1000") & dat.out$ci$year<min(GCM.hist.use$year),]
    
    met.final.day <- rbind(LDAS.use, CRU.use, GCM.hist.use, GCM.p1k.use)
    met.final.day$splice <- NA
    met.final.day[met.final.day$year>=1847 & met.final.day$year<=1852,"splice"] <- "GCM.p1000-GCM.hist"
    met.final.day[met.final.day$year>=1898 & met.final.day$year<=1903,"splice"] <- "GCM.hist-CRUNCEP"
    met.final.day[met.final.day$year>=1977 & met.final.day$year<=1982,"splice"] <- "CRUNCEP-LDAS"
    met.final.day$splice <- as.factor(met.final.day$splice)
    met.final.day$year.frac <- met.final.day$year + met.final.day$doy/366
    summary(met.final.day)
    met.final.day$dataset <- factor(met.final.day$dataset, levels=c(paste0(GCM, ".p1000"), paste0(GCM, ".hist"), "CRUNCEP", LDAS))
    
    plot.splice <- ggplot(data=met.final.day[!is.na(met.final.day$splice),]) +
      facet_wrap(~splice, scales="free_x", ncol=1) +
      geom_ribbon(aes(x=year.frac, ymin=lwr, ymax=upr, fill=dataset), alpha=0.3) +
      geom_line(aes(x=year.frac, y=mean, color=dataset)) +
      scale_x_continuous(expand=c(0,0), name="Year") +
      scale_y_continuous(expand=c(0,0), name="Daily Mean Value") +
      ggtitle("Bias-Corrected (Splices)") +
      theme_bw() +
      theme(legend.position="top",
            legend.direction="horizontal")
    
    png(file.path(path.out, paste0(GCM, "_", met.var, "_bias-correction.png")), height=8.5, width=14, "in", res=180)
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1,3)))
    print(plot.orig  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(plot.bias  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(plot.splice, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
    dev.off()
    
  } # End QA/QC graphs
  # -----------------------------------
  
  # -----------------------------------
  # Get rid of years we don't want to use
  # -----------------------------------
  CRU.max <- min(dat.out.full$met.bias[dat.out.full$met.bias$dataset==LDAS,"year"])-1
  GCM.hist.max <- min(dat.out.full$met.bias[dat.out.full$met.bias$dataset=="CRUNCEP","year"])-1
  GCM.p1000.max <- min(dat.out.full$met.bias[dat.out.full$met.bias$dataset==paste0(GCM, ".hist"),"year"])-1
  
  for(v in names(dat.out.full)){
    if(v == "met.bias"){
      dat.out.full[[v]] <- dat.out.full[[v]][dat.out.full[[v]]$dataset==LDAS | 
                                               (dat.out.full[[v]]$dataset=="CRUNCEP" & dat.out.full[[v]]$year<=CRU.max) |
                                               (dat.out.full[[v]]$dataset==paste0(GCM, ".hist") & dat.out.full[[v]]$year<=GCM.hist.max) |
                                               (dat.out.full[[v]]$dataset==paste0(GCM, ".p1000") & dat.out.full[[v]]$year<=GCM.p1000.max),]
    } else {
      dat.out.full[[v]]$ci <- dat.out.full[[v]]$ci[dat.out.full[[v]]$ci$dataset==LDAS | 
                                                     (dat.out.full[[v]]$ci$dataset=="CRUNCEP" & dat.out.full[[v]]$ci$year<=CRU.max) |
                                                     (dat.out.full[[v]]$ci$dataset==paste0(GCM, ".hist") & dat.out.full[[v]]$ci$year<=GCM.hist.max) |
                                                     (dat.out.full[[v]]$ci$dataset==paste0(GCM, ".p1000") & dat.out.full[[v]]$ci$year<=GCM.p1000.max),]
      dat.out.full[[v]]$sims <- dat.out.full[[v]]$sims[dat.out.full[[v]]$sims$dataset==LDAS | 
                                                         (dat.out.full[[v]]$sims$dataset=="CRUNCEP" & dat.out.full[[v]]$sims$year<=CRU.max) |
                                                         (dat.out.full[[v]]$sims$dataset==paste0(GCM, ".hist") & dat.out.full[[v]]$sims$year<=GCM.hist.max) |
                                                         (dat.out.full[[v]]$sims$dataset==paste0(GCM, ".p1000") & dat.out.full[[v]]$sims$year<=GCM.p1000.max),]
      
      # Make sure the data is in chronological order so that it gets saved into the netcdf file properly
      dat.out.full[[v]]$ci   <- dat.out.full[[v]]$ci[order(dat.out.full[[v]]$ci$time),]
      dat.out.full[[v]]$sims <- dat.out.full[[v]]$sims[order(dat.out.full[[v]]$sims$time),]
    }
    
  } # End variable selection
  
  # --------------
  # Looking at the low-frequency trends in a continuous time series for 3 ensemble members
  # --------------
  library(zoo)
  sims.check <- paste0("X", sample(1:n, 3)) # randomly pick 3 ensemble members
  for(v in names(dat.out.full)[!names(dat.out.full)=="met.bias"]){
    dat.temp <- dat.out.full[[v]]$sims[,c("year", "doy", sims.check)]
    
    # aggregating up to the annual scale
    dat.temp <- aggregate(dat.temp[,sims.check], by=list(dat.temp$year), FUN=mean)
    names(dat.temp)[1] <- "year"
    # Making sure that data are ordered 850-2010
    dat.temp <- dat.temp[order(dat.temp$year),]
    
    dat.smooth <- data.frame(rollapply(dat.temp[,sims.check], width=10, fill=NA, FUN=mean))
    
    dat.check <- stack(dat.temp[,c(sims.check)])
    names(dat.check) <- c("annual", "EnsMem")
    dat.check$decadal <- stack(dat.smooth[,c(sims.check)])[,1]
    dat.check$year <- dat.temp$year
    dat.check$Met <- v
    summary(dat.check)
    
    if(v == names(dat.out.full)[1]){
      dat.final <- dat.check
    } else {
      dat.final <- rbind(dat.final, dat.check)
    }
  }
  
  
  dat.final$Met <- factor(dat.final$Met, levels=vars.met)
  png(file.path(path.out, paste0(GCM, "_Ensembles_Smoothed.png")), height=8.5, width=14, "in", res=180)
  print(
    ggplot(data=dat.final) +
      # facet_wrap(~Met, scales="free_y") +
      facet_grid(Met~., scales="free_y") +
      geom_line(aes(x=year, y=annual, color=EnsMem), size=0.5, alpha=0.4) +
      geom_line(aes(x=year, y=decadal, color=EnsMem)) +
      scale_x_continuous(name="Year (A.D.)", expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      guides(color=F) +
      geom_vline(xintercept=c(1850, 1901, 1980), linetype="dashed") +
      theme_bw()
  )
  dev.off()
  # --------------
  
  save(dat.out.full, file=file.path(path.out, paste0(GCM, "_", str_pad(min(ens), 3, pad=0), "-", str_pad(max(ens), 3, pad=0), "_day_alldata.Rdata")))
  # -----------------------------------
  
  # -----------------------------------
  # Saving each variable as a netcdf 
  # -----------------------------------
  # Breaking from previous versions and putting all variables into the same file in 100 year
  # chunks
  library(ncdf4)
  library(lubridate)
  library(stringr)
  
  # Defining variable names, longname & units
  vars.info <- data.frame(name=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"),
                          longname=c("2 meter mean air temperature", 
                                     "2 meter daily maximum air temperature", 
                                     "2 meter daily minimum air temperature",
                                     paste('The per unit area and time ',
                                           'precipitation representing the sum of convective rainfall, ',
                                           'stratiform rainfall, and snowfall', sep=''),
                                     paste('Incident (downwelling) radiation in ',
                                           'the shortwave part of the spectrum averaged over the time ',
                                           'step of the forcing data', sep=''),
                                     paste('Incident (downwelling) longwave ',
                                           'radiation averaged over the time step of the forcing data', sep=''),
                                     'Pressure at the surface',
                                     'Specific humidity measured at the lowest level of the atmosphere',
                                     'Wind speed' 
                          ),
                          units= c("K", "K", "K", "kg m-2 s-1", "W m-2", "W m-2", "Pa", "kg kg-1", "m s-1")
  )
  # Make a few dimensions we can use
  dimY <- ncdim_def( "lon", units="degrees", longname="latitude", vals=site.lat )
  dimX <- ncdim_def( "lat", units="degrees", longname="longitude", vals=site.lon )
  
  
  # yr.bins <- c(min(dat.out.full$met.bias$year), seq(min(dat.out.full$met.bias$year)+50, round(max(dat.out.full$met.bias$year),-2), by=100))
  yrs <- min(dat.out.full$met.bias$year):max(dat.out.full$met.bias$year)
  # GCM="Ameriflux"
  for(i in 1:n){
    # Make a directory for each ensemble member
    out.name <- paste0(site.name, "_", GCM, "_day_", str_pad(ens[i], 3, pad=0))
    new.dir <- file.path(path.out, out.name)
    if(!dir.exists(new.dir)) dir.create(new.dir, recursive=T)  
    
    for(j in 1:length(yrs)){
      date.start <- as.Date(paste0(yrs[j], "-01-01"))
      date.end   <- as.Date(paste0(yrs[j], "-12-31"))
      day.vec <- seq(date.start, date.end, by="day")
      day.vec <- julian(day.vec, origin=as.Date("850-01-01"))
      
      dim.t <- ncdim_def(name = "time",
                         units = paste0("days since 850"),
                         vals = day.vec, # calculating the number of months in this run
                         calendar = "standard", unlim = TRUE)
      
      var.list <- list()
      dat.list <- list()
      
      for(v in names(dat.out.full)[1:(length(dat.out.full)-1)]){
        var.list[[v]] <- ncvar_def(v, units=paste(vars.info[vars.info$name==v, "units"]), dim=list(dimX, dimY, dim.t), longname=paste(vars.info[vars.info$name==v, "longname"]))
        dat.list[[v]] <- array(dat.out.full[[v]]$sims[dat.out.full[[v]]$sims$year>=year(date.start) & dat.out.full[[v]]$sims$year<=year(date.end), paste0("X", i)], dim=c(1,1,length(day.vec)))
      }
      
      # put precip back in (kg m-2 s-1) to keep units consistent
      dat.list$precipf <- dat.list$precipf/(60*60*24) 
      
      # Naming convention: [SITE]_[GCM]_day_[member]_[YEAR].nc
      nc <- nc_create(file.path(new.dir, paste0(out.name, "_", str_pad(year(date.start), 4, pad=0), ".nc")), var.list)
      for(v in 1:length(var.list)) {
        ncvar_put(nc, var.list[[v]], dat.list[[v]])
      }
      nc_close(nc)    
    }
    # Compress the ensemble file  setwd(dir.out) # Go to the output directory so we don't get annoying file paths
    # setwd(path.out)
    # system(paste0("tar -jcvf ", out.name, ".tar.bz2 ", out.name)) # Compress the folder
    # system(paste0("rm -rf ", out.name)) # remove the uncompressed folder
    setwd(wd.base) # Go back to our base directory
    
    
  } # end splitting ensemble members
  # -----------------------------------
  
} # End GCM loop
