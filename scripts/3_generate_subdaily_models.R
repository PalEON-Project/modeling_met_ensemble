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
# Make statistical models that take the daily, bias-corrected met files that come out 
# of step 2 (daily means) and predict subdaily values (e.g. hourly or 3-hourly) using 
# the a training dataset (e.g. NLDAS, GLDAS)
#
# This script just generates and stores the models so that they can be applied and 
# filtered through the bias-corrected met.  There are many ways in which both the 
# models and approach can be sped, up (saving models & betas separately, etc.), but 
# this should hopefully just get it working for now.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Load Libraries, set up file directories
# 1. Load and format training data
#    1.0 Read data & Make time stamps
#    1.1 Coming up with the daily means that are what we can use as predictors
#    1.2 Setting up a 1-hour lag -- smooth transitions at midnight
#    1.3 Setting up a variable to 'preview' the next day's mean to help get smoother transitions
#    1.4 calculate tmin & tmax as departure from mean; order data
# 2. Train the models for each variable and save them to be read in as needed
#    2.1 Generating all the daily models, save the output as .Rdata files, then clear memory
# -----------------------------------

# ------------------------------------------
# 0. Load Libraries, set up file directories
# ------------------------------------------
# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(MASS)
# library(lubridate)
library(ggplot2)
# library(tictoc)
rm(list=ls())

mod.out <- "/projectnb/dietzelab/paleon/met_ensemble/data/met_ensembles/HARVARD/subday_models"
# mod.out <- "~/Desktop/met_ensembles/HARVARD/subday_models"
fig.dir <- file.path(mod.out, "model_qaqc")

if(!dir.exists(mod.out)) dir.create(mod.out, recursive = T)
if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)
# ------------------------------------------


# ------------------------------------------
# 1. Load and format training data
# ------------------------------------------
# ----------
# 1.0 Read data & Make time stamps
# ----------
{
  # Load the data
  dat.train <- read.csv("../data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")
  # dat.train$doy <- as.ordered(dat.train$doy)
  
  # order the data just o make life easier
  dat.train <- dat.train[order(dat.train$year, dat.train$doy, dat.train$hour, decreasing=T),]
  dat.train[1:25,]
  head(dat.train)
  summary(dat.train)
  
  # Add various types of time stamps to make life easier
  dat.train$date <- strptime(paste(dat.train$year, dat.train$doy+1, dat.train$hour, sep="-"), "%Y-%j-%H", tz="GMT")
  dat.train$time.hr <- as.numeric(difftime(dat.train$date, "2016-01-01", tz="GMT", units="hour"))
  dat.train$time.day <- as.numeric(difftime(dat.train$date, "2016-01-01", tz="GMT", units="day"))+1/24
  dat.train$time.day2 <- as.integer(dat.train$time.day)-1
  dat.train <- dat.train[order(dat.train$time.hr, decreasing=T),]
  # dat.train[1:25,]
  # head(dat.train)
  
  # For some reason certain days are getting an extra hour, so make sure it lines up right
  for(i in max(dat.train$time.day2):min(dat.train$time.day2)){
    rows.now <- which(dat.train$time.day2==i)
    if(length(rows.now)<=24) next
    
    dat.train[rows.now[25],"time.day2"] <- i-1
  }
}
# ----------

# ----------
# 1.1 Coming up with the daily means that are what we can use as predictors
# ----------
{
  train.day <- aggregate(dat.train[,c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")], 
                         by=dat.train[,c("year", "doy")],
                         FUN=mean)
  names(train.day)[3:9] <- c("tmean.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")
  train.day$tmax.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=max)$x
  train.day$tmin.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=min)$x
  summary(train.day)
  
  
  dat.train <- merge(dat.train[,], train.day, all.x=T, all.y=T)
  summary(dat.train)
}
# ----------

# ----------
# 1.2 Setting up a 1-hour lag -- smooth transitions at midnight
# NOTE: because we're filtering from the present back through the past, -1 will associate the closest 
#       hour that we've already done (midnight) with the day we're currently working on
# ----------
{
  vars.hour <- c("tair","precipf", "swdown", "lwdown", "press", "qair", "wind")
  vars.lag <- c("lag.tair", "lag.precipf", "lag.swdown", "lag.lwdown", "lag.press", "lag.qair", "lag.wind") 
  lag.day <- dat.train[dat.train$hour==0,c("year", "doy", "time.day2", vars.hour)]
  names(lag.day)[4:10] <- vars.lag
  # lag.day$lag.diff <- dat.train[dat.train$hour==6,"tair"] - lag.day$lag.day # Lag is the change in temp in the proximate 3 hours
  lag.day <- aggregate(lag.day[,vars.lag],
                       by=lag.day[,c("year", "doy", "time.day2")],
                       FUN=mean)
  lag.day$lag.tmin <- aggregate(dat.train[,c("tair")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=min)[,"x"] # Add in a lag for the next day's min temp
  lag.day$lag.tmax <- aggregate(dat.train[,c("tair")],  
                                by=dat.train[,c("year", "doy", "time.day2")],
                                FUN=max)[,"x"] # Add in a lag for the next day's min temp
  lag.day$time.day2 <- lag.day$time.day2-1
  head(lag.day)
  summary(lag.day)
  
  dat.train <- merge(dat.train, lag.day[,c("time.day2", vars.lag, "lag.tmin", "lag.tmax")], all.x=T)
}
# ----------


# ----------
# 1.3 Setting up a variable to 'preview' the next day's mean to help get smoother transitions
# NOTE: because we're filtering from the present back through the past, +1 will associate 
#       the mean for the next day we're going to model with the one we're currently working on
# ----------
{
  vars.day <- c("tmean.day", "tmax.day", "tmean.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")
  vars.next <- c("next.tmean", "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind") 
  
  next.day <- dat.train[c("year", "doy", "time.day2", vars.day)]
  names(next.day)[4:12] <- vars.next
  # next.day$next.diff <- dat.train[dat.train$hour==6,"tair"] - next.day$next.day # Lag is the change in temp in the proximate 3 hours
  next.day <- aggregate(next.day[,vars.next],
                        by=next.day[,c("year", "doy", "time.day2")],
                        FUN=mean)
  next.day$time.day2 <- next.day$time.day2+1  
  
  dat.train <- merge(dat.train, next.day[,c("time.day2", vars.next)], all.x=T)
}
# ----------


# ----------
# 1.4 calculate tmin & tmax as departure from mean; order data
# ----------
{
  # Lookign at max & min as departure from mean
  dat.train$max.dep <- dat.train$tmax.day - dat.train$tmean.day
  dat.train$min.dep <- dat.train$tmin.day - dat.train$tmean.day
  summary(dat.train)
  
  # Order the data just to help with my sanity when visualizing
  dat.train <- dat.train[order(dat.train$time.hr, decreasing=T),]
  summary(dat.train)
}
# ----------
# ------------------------------------------



# ------------------------------------------
# 2 Train the models for each variable and save them to be read in as needed
# ------------------------------------------
source("temporal_downscale_functions.R")

# ---------
# 2.1 Generating all the daily models, save the output as .Rdata files, then clear memory
# Note: Could save Betas as .nc files that we pull from as needed to save memory; but for now just leaving it in the .Rdata file for eas
# Note: To avoid propogating too much wonkiness in hourly data, any co-variates are at the daily level
# ---------
# Settings for the calculations
n.beta=1000
resids=F
paralell=F
n.cores=4

mod.tair.doy    <- model.tair   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.tair.doy, betas="betas", outfile=file.path(mod.out, "betas_tair.nc"))
save.model(model.out=mod.tair.doy, model="model", outfile=file.path(mod.out, "model_tair.Rdata"))
rm(mod.tair.doy)

mod.precipf.doy <- model.precipf(dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.precipf.doy, betas="betas", outfile=file.path(mod.out, "betas_precipf.nc"))
save.model(model.out=mod.precipf.doy, model="model", outfile=file.path(mod.out, "model_precipf.Rdata"))
rm(mod.precipf.doy)

mod.swdown.doy  <- model.swdown (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.swdown.doy, betas="betas", outfile=file.path(mod.out, "betas_swdown.nc"))
save.model(model.out=mod.swdown.doy, model="model", outfile=file.path(mod.out, "model_swdown.Rdata"))
rm(mod.swdown.doy)

mod.lwdown.doy  <- model.lwdown (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.lwdown.doy, betas="betas", outfile=file.path(mod.out, "betas_lwdown.nc"))
save.model(model.out=mod.lwdown.doy, model="model", outfile=file.path(mod.out, "model_lwdown.Rdata"))
rm(mod.lwdown.doy)

mod.press.doy   <- model.press  (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.press.doy, betas="betas", outfile=file.path(mod.out, "betas_press.nc"))
save.model(model.out=mod.press.doy, model="model", outfile=file.path(mod.out, "model_press.Rdata"))
rm(mod.press.doy)

mod.qair.doy    <- model.qair   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.qair.doy, betas="betas", outfile=file.path(mod.out, "betas_qair.nc"))
save.model(model.out=mod.qair.doy, model="model", outfile=file.path(mod.out, "model_qair.Rdata"))
rm(mod.qair.doy)

mod.wind.doy    <- model.wind   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
save.betas(model.out=mod.wind.doy, betas="betas", outfile=file.path(mod.out, "betas_wind.nc"))
save.model(model.out=mod.wind.doy, model="model", outfile=file.path(mod.out, "model_wind.Rdata"))
rm(mod.wind.doy)
# ---------

# ------------------------------------------
