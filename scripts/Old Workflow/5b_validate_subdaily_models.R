# NOTE: THIS SCRIPT DOES NOT CURRENTLY WORK --> NEED TO CHANGE THE INPUT FORMATTING!
#
#
#
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
library(tictoc)
rm(list=ls())

# mod.out <- "/projectnb/dietzelab/paleon/met_ensemble/data/met_ensembles/HARVARD/subday_models"
site.name = "GOOSE"
vers=".v1"
site.lat  = 43.068496
site.lon  = -73.287425

mod.out <- file.path("~/met_ensemble/data/met_ensembles", paste0(site.name, vers), "1hr/mods.tdm")
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
  dat.train <- read.csv(file.path("../data/paleon_sites", paste0(site.name, vers), "NLDAS_1980-2015.csv"))
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
parallel=T
n.cores=4
n.ens=ens.hr=10


dat.mod <- dat.train
dat.mod$time.day <- dat.train$time.day2
dat.mod <- dat.mod[!is.na(dat.mod$next.tmax),]

dat.sim <- list()


lags.init <- list() # Need to initialize lags first outside of the loop and then they'll get updated internally
lags.init[["tair"   ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"tair"]   , dim=c(1, ens.hr)))
lags.init[["tmax"   ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"tmax.day"]   , dim=c(1, ens.hr)))
lags.init[["tmin"   ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"tmin.day"]   , dim=c(1, ens.hr)))
lags.init[["precipf"]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"precipf"], dim=c(1, ens.hr)))/(60*60*24)
lags.init[["swdown" ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"swdown"] , dim=c(1, ens.hr)))
lags.init[["lwdown" ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"lwdown"] , dim=c(1, ens.hr)))
lags.init[["press"  ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"press"]  , dim=c(1, ens.hr)))
lags.init[["qair"   ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"qair"]   , dim=c(1, ens.hr)))
lags.init[["wind"   ]] <- data.frame(array(dat.mod[dat.mod$time.hr==max(dat.mod$time.hr),"wind"]   , dim=c(1, ens.hr)))



# ---------
# Tair
# ---------
{
  mod.tair.doy    <- model.tair   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["tair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$tair = -99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$tair)
      names(sim.lag) <- c("lag.tair", "ens")
      
      sim.lag$lag.tmin <- stack(lags.init$tmin)[,1]
      sim.lag$lag.tmax <- stack(lags.init$tmax)[,1]
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["tair"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$tair)))))
      names(sim.lag) <- c("lag.tair", "ens")
      sim.lag$lag.tmin <- stack(apply(dat.sim[["tair"]][dat.mod$time.day==(i+1),], 2, min))[,1]
      sim.lag$lag.tmax <- stack(apply(dat.sim[["tair"]][dat.mod$time.day==(i+1),], 2, max))[,1]
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.tair.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.tair.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.tair.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$tair), replace=T)
    
    for(j in 1:ncol(dat.sim$tair)){
      dat.sim[["tair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "tair")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="tair"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=tair), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=tair), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=tair), color="black") +
        geom_point(aes(x=date, y=tair), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
      
  }
  
  
  graph.resids(var="tair", dat.train=dat.train, model.var=mod.tair.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.tair.doy, betas="betas", outfile=file.path(mod.out, "betas_tair.nc"))
  save.model(model.out=mod.tair.doy, model="model", outfile=file.path(mod.out, "model_tair.Rdata"))
  rm(mod.tair.doy)
}
# ---------


# ---------
# precipf
# ---------
{
  mod.precipf.doy    <- model.precipf   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["precipf"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$precipf = -99999 # Dummy value so there's a column
    dat.temp$rain.prop = 99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$precipf)
      names(sim.lag) <- c("lag.precipf", "ens")
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["precipf"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$precipf)))))
      names(sim.lag) <- c("lag.precipf", "ens")
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    

    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.precipf.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.precipf.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.precipf.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    # Re-distribute negative probabilities -- add randomly to make more peaky
    if(max(dat.pred)>0){ # If there's no rain on this day, skip the re-proportioning
      tmp <- 1:nrow(dat.pred) # A dummy vector of the 
      for(j in 1:ncol(dat.pred)){
        if(min(dat.pred[,j])>=0) next
        rows.neg <- which(dat.pred[,j]<0)
        rows.add <- sample(tmp[!tmp %in% rows.neg],length(rows.neg), replace=T)
        
        for(z in 1:length(rows.neg)){
          dat.pred[rows.add[z],j] <- dat.pred[rows.add[z],j] - dat.pred[rows.neg[z],j]
          dat.pred[rows.neg[z],j] <- 0  
        }
      }
      dat.pred <- dat.pred/rowSums(dat.pred)
      dat.pred[is.na(dat.pred)] <- 0
    }
    # Convert precip into real units
    dat.pred <- dat.pred*(dat.temp$precipf.day*24)
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$precipf), replace=T)
    
    for(j in 1:ncol(dat.sim$precipf)){
      dat.sim[["precipf"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  dat.sim$precipf <- dat.sim$precipf*24
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "precipf")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="precipf"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=precipf), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=precipf), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=precipf), color="black") +
        geom_point(aes(x=date, y=precipf), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="precipf", dat.train=dat.train, model.var=mod.precipf.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.precipf.doy, betas="betas", outfile=file.path(mod.out, "betas_precipf.nc"))
  save.model(model.out=mod.precipf.doy, model="model", outfile=file.path(mod.out, "model_precipf.Rdata"))
  rm(mod.precipf.doy)
}
# ---------


# ---------
# swdown
# ---------
{
  mod.swdown.doy    <- model.swdown   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["swdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    
    # For SWDOWN, we only want to model daylight hours -- make sure this matches what's in the swdown function
    day.now = unique(dat.mod[dat.mod$time.day==i, "doy"])
    
    # Use the training data to figure out night/day
    hrs.day = unique(dat.train[dat.train$doy==day.now & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
    
    rows.now = which(dat.mod$time.day==i)
    rows.mod = which(dat.mod$time.day==i & dat.mod$hour %in% hrs.day)
    
    dat.temp <- dat.mod[rows.mod,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    
    dat.temp$swdown = -99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.swdown.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.swdown.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.swdown.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    
    dat.pred[dat.pred<0] <- 0
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$swdown), replace=T)
    for(j in 1:ncol(dat.sim$swdown)){
      dat.sim[["swdown"]][rows.mod,j] <- dat.pred[,cols.prop[j]]
    }
    
    # For night time hours, value shoudl be 0
    dat.sim[["swdown"]][rows.now[!rows.now %in% rows.mod],] <- 0
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "swdown")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="swdown"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=swdown), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=swdown), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=swdown), color="black") +
        geom_point(aes(x=date, y=swdown), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="swdown", dat.train=dat.train, model.var=mod.swdown.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.swdown.doy, betas="betas", outfile=file.path(mod.out, "betas_swdown.nc"))
  save.model(model.out=mod.swdown.doy, model="model", outfile=file.path(mod.out, "model_swdown.Rdata"))
  rm(mod.swdown.doy)
}
# ---------

# ---------
# lwdown
# ---------
{
  mod.lwdown.doy    <- model.lwdown   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["lwdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$lwdown = 99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$lwdown)
      names(sim.lag) <- c("lag.lwdown", "ens")
      
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["lwdown"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$lwdown)))))
      names(sim.lag) <- c("lag.lwdown", "ens")
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.lwdown.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.lwdown.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.lwdown.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    dat.pred <- dat.pred^2 # because squared to prevent negative numbers
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$lwdown), replace=T)
    
    for(j in 1:ncol(dat.sim$lwdown)){
      dat.sim[["lwdown"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "lwdown")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="lwdown"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=lwdown), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=lwdown), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=lwdown), color="black") +
        geom_point(aes(x=date, y=lwdown), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="lwdown", dat.train=dat.train, model.var=mod.lwdown.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.lwdown.doy, betas="betas", outfile=file.path(mod.out, "betas_lwdown.nc"))
  save.model(model.out=mod.lwdown.doy, model="model", outfile=file.path(mod.out, "model_lwdown.Rdata"))
  rm(mod.lwdown.doy)
}
# ---------


# ---------
# press
# ---------
{
  mod.press.doy    <- model.press   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["press"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$press = 99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$press)
      names(sim.lag) <- c("lag.press", "ens")
      
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["press"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$press)))))
      names(sim.lag) <- c("lag.press", "ens")
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.press.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.press.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.press.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$press), replace=T)
    
    for(j in 1:ncol(dat.sim$press)){
      dat.sim[["press"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "press")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="press"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=press), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=press), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=press), color="black") +
        geom_point(aes(x=date, y=press), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="press", dat.train=dat.train, model.var=mod.press.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.press.doy, betas="betas", outfile=file.path(mod.out, "betas_press.nc"))
  save.model(model.out=mod.press.doy, model="model", outfile=file.path(mod.out, "model_press.Rdata"))
  rm(mod.press.doy)
}
# ---------

# ---------
# qair
# ---------
{
  mod.qair.doy    <- model.qair   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["qair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$qair = 99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$qair)
      names(sim.lag) <- c("lag.qair", "ens")
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["qair"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$qair)))))
      names(sim.lag) <- c("lag.qair", "ens")
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.qair.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))

    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.qair.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    dat.pred <- exp(dat.pred) # because log-transformed
    
    # having a *very* hard time keeping humidity reasonable, even with transformations so truncating the max
    dat.pred[dat.pred>max(dat.train$qair)] <- max(dat.train$qair) 
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$qair), replace=T)
    
    for(j in 1:ncol(dat.sim$qair)){
      dat.sim[["qair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "qair")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="qair"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=qair), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=qair), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=qair), color="black") +
        geom_point(aes(x=date, y=qair), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="qair", dat.train=dat.train, model.var=mod.qair.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.qair.doy, betas="betas", outfile=file.path(mod.out, "betas_qair.nc"))
  save.model(model.out=mod.qair.doy, model="model", outfile=file.path(mod.out, "model_qair.Rdata"))
  rm(mod.qair.doy)
}
# ---------



# ---------
# wind
# ---------
{
  mod.wind.doy    <- model.wind   (dat.train=dat.train[,], resids=resids, parallel=parallel, n.cores=n.cores, n.beta=n.beta)
  
  # Doing the prediction 
  dat.sim[["wind"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
  pb <- txtProgressBar(min=min(abs(dat.mod$time.day)), max=max(abs(dat.mod$time.day)), style=3)
  tic()
  for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
    setTxtProgressBar(pb, abs(i))
    rows.now = which(dat.mod$time.day==i)
    dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                   "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                   "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
    dat.temp$wind = 99999 # Dummy value so there's a column
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(dat.mod$time.day)){
      sim.lag <- stack(lags.init$wind)
      names(sim.lag) <- c("lag.wind", "ens")
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[["wind"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$wind)))))
      names(sim.lag) <- c("lag.wind", "ens")
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    rows.beta <- sample(1:n.beta, n.ens, replace=T)
    Rbeta <- as.matrix(mod.wind.doy[[paste(day.now)]]$betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    # dimnames(Rbeta)[[2]] <- names(coef(mod.wind.doy[[paste(day.now)]]))
    
    dat.pred <- predict.met(newdata=dat.temp, 
                            model.predict=mod.wind.doy[[paste(day.now)]]$model, 
                            Rbeta=Rbeta, 
                            resid.err=F,
                            model.resid=NULL, 
                            Rbeta.resid=NULL, 
                            n.ens=n.ens)
    dat.pred <- dat.pred^2 # because squared to prevent negative numbers
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim$wind), replace=T)
    
    for(j in 1:ncol(dat.sim$wind)){
      dat.sim[["wind"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
  }
  toc()
  
  # Graphing
  fig.ens <- file.path(fig.dir, "model_validation", "wind")
  if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
  
  var="wind"
  dat.mod$var.pred <- apply(dat.sim[[var]],1,mean)
  dat.mod$var.025 <- apply(dat.sim[[var]],1,quantile, 0.025)
  dat.mod$var.975 <- apply(dat.sim[[var]],1,quantile, 0.975)
  
  for(y in max(dat.mod$year):min(dat.mod$year)){
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14) & dat.mod$year==y,]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14) & dat.mod$year==y,]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14) & dat.mod$year==y,]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14) & dat.mod$year==y,]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    png(file.path(fig.ens, paste0(var, y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        geom_point(aes(x=date, y=wind), color="black", size=0.1, alpha=0.5) +
        geom_line(aes(x=date, y=wind), color="black", size=0.2, alpha=0.5) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
    
    png(file.path(fig.ens, paste0(var, "_examples_", y,".png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.graph[dat.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_line(aes(x=date, y=wind), color="black") +
        geom_point(aes(x=date, y=wind), color="black", size=0.5) +
        geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=var.pred), color="blue") +
        geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
        scale_y_continuous(name=var) +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(paste0(var, ": ", y)) +
        theme_bw()
    )
    dev.off()
    
  }
  
  
  graph.resids(var="wind", dat.train=dat.train, model.var=mod.wind.doy, fig.dir=fig.dir)
  save.betas(model.out=mod.wind.doy, betas="betas", outfile=file.path(mod.out, "betas_wind.nc"))
  save.model(model.out=mod.wind.doy, model="model", outfile=file.path(mod.out, "model_wind.Rdata"))
  rm(mod.wind.doy)
}
# ---------

# ------------------------------------------
