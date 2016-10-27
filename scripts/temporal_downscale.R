# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(MASS)
library(lubridate)
library(ggplot2)
library(tictoc)
rm(list=ls())

fig.dir <- "../data/met_ensembles/HARVARD/subday_qaqc"
if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)

dat.train <- read.csv("../data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")
# dat.train$doy <- as.ordered(dat.train$doy)

# Now trying things in a predictive framework
dat.train <- dat.train[order(dat.train$year, dat.train$doy, dat.train$hour, decreasing=T),]
dat.train[1:25,]
head(dat.train)
summary(dat.train)

# dat.train <- dat.train[,c("dataset", "year", "doy", "hour", "tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
dat.train$date <- strptime(paste(dat.train$year, dat.train$doy+1, dat.train$hour, sep="-"), "%Y-%j-%H", tz="GMT")
dat.train$time.hr <- as.numeric(difftime(dat.train$date, "2016-01-01", tz="GMT", units="hour"))
dat.train$time.day <- as.numeric(difftime(dat.train$date, "2016-01-01", tz="GMT", units="day"))+1/24
dat.train$time.day2 <- as.integer(dat.train$time.day)-1
dat.train <- dat.train[order(dat.train$time.hr, decreasing=T),]
# dat.train[1:25,]
# head(dat.train)

# For some reason certain days are getting an extra hour
for(i in max(dat.train$time.day2):min(dat.train$time.day2)){
  rows.now <- which(dat.train$time.day2==i)
  if(length(rows.now)<=24) next
  
  dat.train[rows.now[25],"time.day2"] <- i-1
}

train.day <- aggregate(dat.train[,c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")], 
                       by=dat.train[,c("year", "doy")],
                       FUN=mean)
names(train.day)[3:9] <- c("tmean.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")
train.day$tmax.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=max)$x
train.day$tmin.day <- aggregate(dat.train[,c("tair")], by=dat.train[,c("year", "doy")], FUN=min)$x
summary(train.day)


dat.train <- merge(dat.train[,], train.day, all.x=T, all.y=T)
summary(dat.train)

# ----------
# Setting up a 1-hour lag -- smooth transitions at midnight
# NOTE: because we're filtering from the present back through the past, -1 will associate the closest 
#       hour that we've already done (midnight) with the day we're currently working on
# ----------
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
# ----------


# ----------
# Setting up a variable to 'preview' the next day's mean to help get smoother transitions
# NOTE: because we're filtering from the present back through the past, +1 will associate 
#       the mean for the next day we're going to model with the one we're currently working on
# ----------
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
# ----------


# Order the data just to help with my sanity when visualizing
dat.train <- dat.train[order(dat.train$time.hr, decreasing=T),]
summary(dat.train)

# Lookign at max & min as departure from mean
dat.train$max.dep <- dat.train$tmax.day - dat.train$tmean.day
dat.train$min.dep <- dat.train$tmin.day - dat.train$tmean.day
summary(dat.train)
# ------------------------------------------

# ------------------------------------------
# Set up data for modeling with uncertainty propogation
# ------------------------------------------
# Set up data
dat.mod <- dat.train[dat.train$year>=2014,]
dat.mod[, "lag.tair"] <- NA
dat.mod[, "mod.tair"] <- NA
# head(dat.mod)

n.ens=30
dat.sim <- list() # Each variable needs to be a layer in a list so we can propogate that uncertainty
# ------------------------------------------

# ------------------------------------------
# Modeling SWDOWN 
# ------------------------------------------
{
  # ---------
  # Generating all the daily models and storing it into an easy-to-find list
  # ---------
  source("temporal_downscale_functions.R")
  tic()
  # dat.train[dat.train$swdown==0, "swdown"] <- 1e-6
  # dat.train[dat.train$swdown==0, "lag.swdown"] <- 1e-6
  mod.swdown.doy <- model.swdown(dat.train=dat.train[,], resids=T, parallel=F, n.cores=4, n.beta=100)
  toc()
  length(mod.swdown.doy)
  # ---------
  
  # ---------
  # Looking at the residuals from the model
  # ---------
  {
    for(i in names(mod.swdown.doy)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      hrs.day = unique(dat.train[dat.train$doy==as.numeric(i) & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
      dat.train[dat.train$doy==as.numeric(i) & dat.train$hour %in% hrs.day, "resid"] <- resid(mod.swdown.doy[[i]]$model)
      dat.train[dat.train$doy==as.numeric(i) & dat.train$hour %in% hrs.day, "predict"] <- predict(mod.swdown.doy[[i]]$model)
      
      dat.train[dat.train$doy==as.numeric(i) & !(dat.train$hour %in% hrs.day), "predict"] <- 0
    }
    summary(dat.train)
  
    png(file.path(fig.dir, "SWDOWN_Resid_vs_Hour.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ hour, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
  
    png(file.path(fig.dir, "SWDOWN_Resid_vs_DOY.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ doy, data=dat.train, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
    dev.off()
  
    png(file.path(fig.dir, "SWDOWN_Resid_vs_Predict.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ predict, data=dat.train); abline(h=0, col="red")
    dev.off()
  
    png(file.path(fig.dir, "SWDOWN_Resid_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ swdown, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
  
    png(file.path(fig.dir, "SWDOWN_Predict_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(predict ~ swdown, data=dat.train, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
  
    # Looking at the daily maxes & mins
    day.stats <- aggregate(dat.train[,c("swdown", "resid", "predict")],
                           by=dat.train[,c("time.day2", "year", "doy")],
                           FUN=mean, na.rm=T)
    day.stats$mod.max <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=max)[,"x"]
    day.stats$mod.min <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=min)[,"x"]
  
    png(file.path(fig.dir, "SWDOWN_Predict_vs_SWmean.png"), height=8, width=8, units="in", res=180)
    plot(predict~ swdown, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
  
    png(file.path(fig.dir, "SWDOWN_ResidualHistograms.png"), height=6, width=8, units="in", res=180)
    # par(mfrow=c(3,1))
    par(mfrow=c(1,1))
    hist(day.stats$resid, main="Daily Mean Residuals")
    dev.off()
  }
  # ---------
  
  # ---------
  # Predict via Day filter
  # Note: this can be generalized to just run by DOY for all years at once since there's no memory in the system
  # ---------
  {
    # library(MASS)
    
    dat.sim[["swdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    dat.sim[["swdown"]][1,] <- dat.mod[1,"swdown"]
    
    pb <- txtProgressBar(min=min(dat.mod$doy), max=max(dat.mod$doy), style=3)
    set.seed(138)
    tic()
    for(i in min(dat.mod$doy):max(dat.mod$doy)){
      setTxtProgressBar(pb, abs(i))
      
      # For SWDOWN, we only want to model daylight hours -- make sure this matches what's in the swdown function
      day.now = i
      
      # Use the training data to figure out night/day
      hrs.day = unique(dat.train[dat.train$doy==day.now & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
      
      rows.now = which(dat.mod$doy==i)
      rows.mod = which(dat.mod$doy==i & dat.mod$hour %in% hrs.day)
      
      dat.temp <- dat.mod[rows.mod,c("time.day2", "year", "doy", "hour", "tmax.day", "tmin.day", "tmean.day", "max.dep", "min.dep", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")]
      dat.temp$swdown = 99999 # Dummy value so there's a column
      dat.temp <- dat.temp[complete.cases(dat.temp),]
      
      # dat.pred <- predict.met(newdata=dat.temp, model.predict=mod.swdown.doy[[paste(day.now)]]$model, betas=mod.swdown.doy[[paste(day.now)]]$betas, n.ens=n.ens)
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.swdown.doy[[paste(day.now)]]$model, 
                              betas=mod.swdown.doy[[paste(day.now)]]$betas, 
                              resid.err=T,
                              model.resid=mod.swdown.doy[[paste(day.now)]]$model.resid, 
                              betas.resid=mod.swdown.doy[[paste(day.now)]]$betas.resid, 
                              n.ens=n.ens)
      dat.pred[dat.pred<0] <- 0
      # dat.pred <- exp(dat.pred)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$swdown), replace=T)
      for(j in 1:ncol(dat.sim$swdown)){
        dat.sim[["swdown"]][rows.mod,j] <- dat.pred[,cols.prop[j]]
      }
      
      # For night time hours, value shoudl be 0
      dat.sim[["swdown"]][rows.now[!rows.now %in% rows.mod],] <- 0
      
      dat.mod[rows.now,"mod.swdown"] <- apply(dat.sim[["swdown"]][rows.now,],1, mean)
    }
    # dat.mod$mod.swdown.mean <- apply(dat.sim[["swdown"]], 1, mean)
    dat.mod$mod.swdown.025   <- apply(dat.sim[["swdown"]], 1, quantile, 0.025)
    dat.mod$mod.swdown.975   <- apply(dat.sim[["swdown"]], 1, quantile, 0.975)
    summary(dat.mod)
  }
  toc()
  # ---------
  
  # ---------
  # Graph the output
  # ---------
  {
    for(y in unique(dat.mod$year)){
      png(file.path(fig.dir, paste0("swdown_", y, "_year.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_ribbon(aes(x=date, ymin=mod.swdown.025, ymax=mod.swdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.swdown), color="blue") +
          geom_point(aes(x=date, y=mod.swdown), color="blue", size=0.5) +
          geom_line(aes(x=date, y=swdown), color="black", alpha=0.5) +
          geom_point(aes(x=date, y=swdown), color="black", size=0.3, alpha=0.5) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("swdown_", y, "_year_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_point(aes(x=swdown, y=mod.swdown), color="black", size=0.5) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
    }  
    
    
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    for(y in unique(dat.graph$year)){
      png(file.path(fig.dir, paste0("swdown_",y,"_examples.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_line(aes(x=date, y=swdown), color="black") +
          geom_point(aes(x=date, y=swdown), color="black", size=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.swdown.025, ymax=mod.swdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.swdown), color="blue") +
          geom_point(aes(x=date, y=mod.swdown), color="blue", size=0.5) +
          scale_y_continuous(name="shortwave radiation") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("swdown_", y, "_examp.es_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_point(aes(x=swdown, y=mod.swdown), color="black", size=0.5) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
    }
  }
  # ---------
  rm(mod.swdown.doy) # Clear out the model to save memory
}
# ------------------------------------------

# ------------------------------------------
# Modeling Temperature 
# ------------------------------------------
{
  # ---------
  # Generating all the daily models and storing it into an easy-to-find list
  # ---------
  source("temporal_downscale_functions.R")
  tic()
  mod.tair.doy <- model.tair(dat.train=dat.train[,], resids=T, parallel=F, n.cores=4, n.beta=100)
  toc()
  length(mod.tair.doy)
  # ---------
  
  # ---------
  # Looking at the residuals from the model
  # ---------
  {
    for(i in names(mod.tair.doy)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.tair), "resid"] <- resid(mod.tair.doy[[i]]$model)
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.tair), "predict"] <- predict(mod.tair.doy[[i]]$model)
    }
    summary(dat.train)
    
    png(file.path(fig.dir, "TAIR_Resid_vs_Hour.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ hour, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "TAIR_Resid_vs_DOY.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ doy, data=dat.train, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
    dev.off()
    
    png(file.path(fig.dir, "TAIR_Resid_vs_Predict.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ predict, data=dat.train); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "TAIR_Resid_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ tair, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "TAIR_Predict_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(predict ~ tair, data=dat.train, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    # Looking at the daily maxes & mins
    day.stats <- aggregate(dat.train[,c("tair", "resid", "predict")],
                           by=dat.train[,c("time.day2", "year", "doy")],
                           FUN=mean, na.rm=T)
    day.stats$mod.max <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=max)[,"x"]
    day.stats$mod.min <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=min)[,"x"]
    
    png(file.path(fig.dir, "TAIR_Predict_vs_SWmean.png"), height=8, width=8, units="in", res=180)
    plot(predict~ tair, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    png(file.path(fig.dir, "TAIR_ResidualHistograms.png"), height=6, width=8, units="in", res=180)
    # par(mfrow=c(3,1))
    par(mfrow=c(1,1))
    hist(day.stats$resid, main="Daily Mean Residuals")
    dev.off()
  }
  # ---------
  
  # ---------
  # Predict via Day filter
  # ---------
  {
    library(MASS)
    
    dat.sim[["tair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    dat.sim[["tair"]][1,] <- dat.mod[1,"tair"]
  
    pb <- txtProgressBar(min=abs(max(dat.mod$time.day2)), max=abs(min(dat.mod$time.day2)), style=3)
    set.seed(138)
    tic()
    for(i in max(dat.mod$time.day2):min(dat.mod$time.day2)){
      setTxtProgressBar(pb, abs(i))
      
      rows.now = which(dat.mod$time.day2==i)
      dat.temp <- dat.mod[rows.now,c("time.day2", "year", "doy", "hour", "tmax.day", "tmin.day", "tmean.day", "max.dep", "min.dep", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")]
      dat.temp$tair = -99999 # Dummy value so there's a column
      dat.temp <- dat.temp[complete.cases(dat.temp),]
      # day.now = dat.temp$doy
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day2)){
        sim.lag <- stack(data.frame(array(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour==0),"tair"], dim=c(1, ncol(dat.sim$tair)))))
        names(sim.lag) <- c("lag.tair", "ens")
        
        sim.lag$lag.tmin <- stack(data.frame(array(min(dat.mod[which(dat.mod$time.day2==i ),"tair"]), dim=c(1, ncol(dat.sim$tair)))))[,1]
        sim.lag$lag.tmax <- stack(data.frame(array(max(dat.mod[which(dat.mod$time.day2==i ),"tair"]), dim=c(1, ncol(dat.sim$tair)))))[,1]
  
        # sim.lag <- stack(data.frame(array(mean(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour<=2),"tair"]), dim=c(1, ncol(dat.sim)))))
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["tair"]][dat.mod$time.day2==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$tair)))))
        names(sim.lag) <- c("lag.tair", "ens")
        sim.lag$lag.tmin <- stack(apply(dat.sim[["tair"]][dat.mod$time.day2==(i+1),], 2, min))[,1]
        sim.lag$lag.tmax <- stack(apply(dat.sim[["tair"]][dat.mod$time.day2==(i+1),], 2, max))[,1]
        
        # sim.lag <- stack(apply(dat.sim[dat.mod$time.day2==(i+1)  & dat.mod$hour<=2,],2, mean))
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.tair.doy[[paste(day.now)]]$model, 
                              betas=mod.tair.doy[[paste(day.now)]]$betas, 
                              resid.err=T,
                              model.resid=mod.tair.doy[[paste(day.now)]]$model.resid, 
                              betas.resid=mod.tair.doy[[paste(day.now)]]$betas.resid, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$tair), replace=T)
      # pred.prop <- array(dim=c(1,ncol(dat.sim)))
      for(j in 1:ncol(dat.sim$tair)){
        dat.sim[["tair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
      dat.mod[rows.now,"mod.tair"] <- apply(dat.sim[["tair"]][rows.now,],1, mean)
      # dat.sim[i,] <- dat.pred
      # dat.mod[i,"mod.hr"] <- predict(mod.fill, dat.mod[i,])
      if(i>min(dat.mod$time.day2)){ 
        dat.mod[dat.mod$time.day2==i-1,"lag.tair" ] <- dat.mod[dat.mod$time.day2==i & dat.mod$hour==0,"mod.tair"] 
        dat.mod[dat.mod$time.day2==i-1,"lag.tmin" ] <- min(dat.mod[dat.mod$time.day2==i,"mod.tair"] )
        dat.mod[dat.mod$time.day2==i-1,"lag.tmax" ] <- max(dat.mod[dat.mod$time.day2==i,"mod.tair"] )
      }
      # if(i>min(dat.mod$time.day2)) dat.mod[dat.mod$time.day2==i-1,"lag.day"] <- mean(dat.mod[dat.mod$time.day2==i & dat.mod$hour<=2,"mod.tair"])
    }
    # dat.mod$mod.tair.mean <- apply(dat.sim[["tair"]], 1, mean)
    dat.mod$mod.tair.025   <- apply(dat.sim[["tair"]], 1, quantile, 0.025)
    dat.mod$mod.tair.975   <- apply(dat.sim[["tair"]], 1, quantile, 0.975)
    summary(dat.mod)
  }
  toc()
  # ---------
  
  # ---------
  # Graph the output
  # ---------
  {
  for(y in unique(dat.mod$year)){
    png(file.path(fig.dir, paste0("Tair_", y, "_year.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        # geom_point(aes(x=date, y=tmax, color=as.factor(doy)), size=0.25, alpha=0.5) +
        # geom_point(aes(x=date, y=tmin, color=as.factor(doy)), size=0.25, alpha=0.5) +
        geom_ribbon(aes(x=date, ymin=mod.tair.025, ymax=mod.tair.975), alpha=0.5, fill="blue") +
        geom_line(aes(x=date, y=mod.tair), color="blue") +
        geom_point(aes(x=date, y=mod.tair), color="blue", size=0.5) +
        geom_line(aes(x=date, y=tair), color="black") +
        geom_point(aes(x=date, y=tair), color="black", size=0.5) +
        # geom_vline(xintercept=seq(min(dat.mod$time.hr), max(dat.mod$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
        scale_x_datetime(expand=c(0,0)) +
        ggtitle(y) +
        theme_bw()
    )
    dev.off()
    
    png(file.path(fig.dir, paste0("Tair_", y, "_year_scatter.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=dat.mod[dat.mod$year==y,]) +
        geom_point(aes(x=tair, y=mod.tair), color="black", size=0.5) +
        geom_abline(slope=1, intercept=0, color="red") +
        ggtitle(y) +
        theme_bw()
    )
    dev.off()
    
  }  
    
    
  dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
  dat.graph1$season <- as.factor("winter")
  dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
  dat.graph2$season <- as.factor("spring")
  dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
  dat.graph3$season <- as.factor("summer")
  dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
  dat.graph4$season <- as.factor("fall")
  
  tair.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
  
  for(y in unique(tair.graph$year)){
    png(file.path(fig.dir, paste0("Tair_",y,"_examples.png")), height=8, width=10, units="in", res=220)
    print(
    ggplot(data=tair.graph[tair.graph$year==y,]) +
      facet_wrap(~season, scales="free") +
      geom_line(aes(x=date, y=tair), color="black") +
      geom_point(aes(x=date, y=tair), color="black", size=0.5) +
      geom_ribbon(aes(x=date, ymin=mod.tair.025, ymax=mod.tair.975), alpha=0.5, fill="blue") +
      geom_line(aes(x=date, y=mod.tair), color="blue") +
      geom_point(aes(x=date, y=mod.tair), color="blue", size=0.5) +
      scale_y_continuous(name="Hourly Air Temperature") +
      scale_x_datetime(expand=c(0,0)) +
      ggtitle(y) +
      theme_bw()
    )
    dev.off()
    png(file.path(fig.dir, paste0("Tair_",y,"_examples_scatter.png")), height=8, width=10, units="in", res=220)
    print(
      ggplot(data=tair.graph[tair.graph$year==y,]) +
        facet_wrap(~season, scales="free") +
        geom_point(aes(x=tair, y=mod.tair), color="black", size=0.5) +
        geom_abline(intercept=0, slope=1, color="red") +
        ggtitle(y) +
        theme_bw()
    )
    dev.off()
    
  }
  }
  # ---------
  
  rm(mod.tair.doy) # Clear out the model to save memory
}
# ------------------------------------------



# ------------------------------------------
# Modeling LWDOWN 
# ------------------------------------------
{
  # ---------
  # Generating all the daily models and storing it into an easy-to-find list
  # ---------
  source("temporal_downscale_functions.R")
  tic()
  mod.lwdown.doy <- model.lwdown(dat.train=dat.train[,], resids=T, parallel=F, n.cores=4, n.beta=100)
  toc()
  length(mod.lwdown.doy)
  # ---------
  
  # ---------
  # Looking at the residuals from the model
  # ---------
  {
    for(i in names(mod.lwdown.doy)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.lwdown), "resid"] <- resid(mod.lwdown.doy[[i]]$model)
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.lwdown), "predict"] <- predict(mod.lwdown.doy[[i]]$model)
    }
    summary(dat.train)
    
    png(file.path(fig.dir, "LWDOWN_Resid_vs_Hour.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ hour, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "LWDOWN_Resid_vs_DOY.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ doy, data=dat.train, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
    dev.off()
    
    png(file.path(fig.dir, "LWDOWN_Resid_vs_Predict.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ predict, data=dat.train); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "LWDOWN_Resid_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ lwdown, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "LWDOWN_Predict_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(predict ~ lwdown, data=dat.train, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    # Looking at the daily maxes & mins
    day.stats <- aggregate(dat.train[,c("lwdown", "resid", "predict")],
                           by=dat.train[,c("time.day2", "year", "doy")],
                           FUN=mean, na.rm=T)
    day.stats$mod.max <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=max)[,"x"]
    day.stats$mod.min <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=min)[,"x"]
    
    png(file.path(fig.dir, "LWDOWN_Predict_vs_SWmean.png"), height=8, width=8, units="in", res=180)
    plot(predict~ lwdown, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    png(file.path(fig.dir, "LWDOWN_ResidualHistograms.png"), height=6, width=8, units="in", res=180)
    # par(mfrow=c(3,1))
    par(mfrow=c(1,1))
    hist(day.stats$resid, main="Daily Mean Residuals")
    dev.off()
  }
  # ---------
  
  # ---------
  # Predict via Day filter
  # ---------
  {
    library(MASS)
    
    dat.sim[["lwdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    dat.sim[["lwdown"]][1,] <- dat.mod[1,"lwdown"]
    
    pb <- txtProgressBar(min=abs(max(dat.mod$time.day2)), max=abs(min(dat.mod$time.day2)), style=3)
    set.seed(138)
    tic()
    for(i in max(dat.mod$time.day2):min(dat.mod$time.day2)){
      setTxtProgressBar(pb, abs(i))
      
      rows.now = which(dat.mod$time.day2==i)
      dat.temp <- dat.mod[rows.now,c("time.day2", "year", "doy", "hour", "tmax.day", "tmin.day", "tmean.day", "max.dep", "min.dep", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")]
      dat.temp$lwdown = -99999 # Dummy value so there's a column
      dat.temp <- dat.temp[complete.cases(dat.temp),]
      # day.now = dat.temp$doy
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day2)){
        sim.lag <- stack(data.frame(array(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour==0),"lwdown"], dim=c(1, ncol(dat.sim$lwdown)))))
        names(sim.lag) <- c("lag.lwdown", "ens")
        
        # sim.lag$lag.tmin <- stack(data.frame(array(min(dat.mod[which(dat.mod$time.day2==i ),"lwdown"]), dim=c(1, ncol(dat.sim$lwdown)))))[,1]
        # sim.lag$lag.tmax <- stack(data.frame(array(max(dat.mod[which(dat.mod$time.day2==i ),"lwdown"]), dim=c(1, ncol(dat.sim$lwdown)))))[,1]
        
        # sim.lag <- stack(data.frame(array(mean(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour<=2),"lwdown"]), dim=c(1, ncol(dat.sim)))))
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["lwdown"]][dat.mod$time.day2==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$lwdown)))))
        names(sim.lag) <- c("lag.lwdown", "ens")
        # sim.lag$lag.tmin <- stack(apply(dat.sim[["lwdown"]][dat.mod$time.day2==(i+1),], 2, min))[,1]
        # sim.lag$lag.tmax <- stack(apply(dat.sim[["lwdown"]][dat.mod$time.day2==(i+1),], 2, max))[,1]
        
        # sim.lag <- stack(apply(dat.sim[dat.mod$time.day2==(i+1)  & dat.mod$hour<=2,],2, mean))
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      dat.temp$tair <- stack(dat.sim$tair[rows.now,])[,1]
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.lwdown.doy[[paste(day.now)]]$model, 
                              betas=mod.lwdown.doy[[paste(day.now)]]$betas, 
                              resid.err=F,
                              model.resid=mod.lwdown.doy[[paste(day.now)]]$model.resid, 
                              betas.resid=mod.lwdown.doy[[paste(day.now)]]$betas.resid, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$lwdown), replace=T)
      # pred.prop <- array(dim=c(1,ncol(dat.sim)))
      for(j in 1:ncol(dat.sim$lwdown)){
        dat.sim[["lwdown"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
      dat.mod[rows.now,"mod.lwdown"] <- apply(dat.sim[["lwdown"]][rows.now,],1, mean)
      # dat.sim[i,] <- dat.pred
      # dat.mod[i,"mod.hr"] <- predict(mod.fill, dat.mod[i,])
      if(i>min(dat.mod$time.day2)){ 
        dat.mod[dat.mod$time.day2==i-1,"lag.lwdown" ] <- dat.mod[dat.mod$time.day2==i & dat.mod$hour==0,"mod.lwdown"] 
        dat.mod[dat.mod$time.day2==i-1,"lag.tmin" ] <- min(dat.mod[dat.mod$time.day2==i,"mod.lwdown"] )
        dat.mod[dat.mod$time.day2==i-1,"lag.tmax" ] <- max(dat.mod[dat.mod$time.day2==i,"mod.lwdown"] )
      }
      # if(i>min(dat.mod$time.day2)) dat.mod[dat.mod$time.day2==i-1,"lag.day"] <- mean(dat.mod[dat.mod$time.day2==i & dat.mod$hour<=2,"mod.lwdown"])
    }
    # dat.mod$mod.lwdown.mean <- apply(dat.sim[["lwdown"]], 1, mean)
    dat.mod$mod.lwdown.025   <- apply(dat.sim[["lwdown"]], 1, quantile, 0.025)
    dat.mod$mod.lwdown.975   <- apply(dat.sim[["lwdown"]], 1, quantile, 0.975)
    summary(dat.mod)
  }
  toc()
  # ---------
  
  # ---------
  # Graph the output
  # ---------
  {
    for(y in unique(dat.mod$year)){
      png(file.path(fig.dir, paste0("lwdown_", y, "_year.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          # geom_point(aes(x=date, y=tmax, color=as.factor(doy)), size=0.25, alpha=0.5) +
          # geom_point(aes(x=date, y=tmin, color=as.factor(doy)), size=0.25, alpha=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.lwdown.025, ymax=mod.lwdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.lwdown), color="blue") +
          geom_point(aes(x=date, y=mod.lwdown), color="blue", size=0.5) +
          geom_line(aes(x=date, y=lwdown), color="black") +
          geom_point(aes(x=date, y=lwdown), color="black", size=0.5) +
          # geom_vline(xintercept=seq(min(dat.mod$time.hr), max(dat.mod$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("lwdown_", y, "_year_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_point(aes(x=lwdown, y=mod.lwdown), color="black", size=0.5) +
          geom_abline(slope=1, intercept=0, color="red") +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
    }  
    
    
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
    dat.graph4$season <- as.factor("fall")
    
    lwdown.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    for(y in unique(lwdown.graph$year)){
      png(file.path(fig.dir, paste0("lwdown_",y,"_examples.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=lwdown.graph[lwdown.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_line(aes(x=date, y=lwdown), color="black") +
          geom_point(aes(x=date, y=lwdown), color="black", size=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.lwdown.025, ymax=mod.lwdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.lwdown), color="blue") +
          geom_point(aes(x=date, y=mod.lwdown), color="blue", size=0.5) +
          scale_y_continuous(name="Hourly Air Temperature") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      png(file.path(fig.dir, paste0("lwdown_",y,"_examples_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=lwdown.graph[lwdown.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_point(aes(x=lwdown, y=mod.lwdown), color="black", size=0.5) +
          geom_abline(intercept=0, slope=1, color="red") +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
    }
  }
  # ---------
  
  rm(mod.lwdown.doy) # Clear out the model to save memory
}
# ------------------------------------------

# ------------------------------------------
# Modeling PRESS 
# ------------------------------------------
{
  # ---------
  # Generating all the daily models and storing it into an easy-to-find list
  # ---------
  source("temporal_downscale_functions.R")
  tic()
  mod.press.doy <- model.press(dat.train=dat.train[,], resids=T, parallel=F, n.cores=4, n.beta=100)
  toc()
  length(mod.press.doy)
  # ---------
  
  # ---------
  # Looking at the residuals from the model
  # ---------
  {
    for(i in names(mod.press.doy)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.press) & !is.na(dat.train$next.press), "resid"] <- resid(mod.press.doy[[i]]$model)
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.press) & !is.na(dat.train$next.press), "predict"] <- predict(mod.press.doy[[i]]$model)
    }
    summary(dat.train)
    
    png(file.path(fig.dir, "PRESS_Resid_vs_Hour.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ hour, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "PRESS_Resid_vs_DOY.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ doy, data=dat.train, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
    dev.off()
    
    png(file.path(fig.dir, "PRESS_Resid_vs_Predict.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ predict, data=dat.train); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "PRESS_Resid_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(resid ~ press, data=dat.train, cex=0.5); abline(h=0, col="red")
    dev.off()
    
    png(file.path(fig.dir, "PRESS_Predict_vs_Obs.png"), height=8, width=8, units="in", res=180)
    plot(predict ~ press, data=dat.train, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    # Looking at the daily maxes & mins
    day.stats <- aggregate(dat.train[,c("press", "resid", "predict")],
                           by=dat.train[,c("time.day2", "year", "doy")],
                           FUN=mean, na.rm=T)
    day.stats$mod.max <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=max)[,"x"]
    day.stats$mod.min <- aggregate(dat.train[,c("predict")],
                                   by=dat.train[,c("time.day2", "year", "doy")],
                                   FUN=min)[,"x"]
    
    png(file.path(fig.dir, "PRESS_Predict_vs_Mean.png"), height=8, width=8, units="in", res=180)
    plot(predict~ press, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
    dev.off()
    
    png(file.path(fig.dir, "PRESS_ResidualHistograms.png"), height=6, width=8, units="in", res=180)
    # par(mfrow=c(3,1))
    par(mfrow=c(1,1))
    hist(day.stats$resid, main="Daily Mean Residuals")
    dev.off()
  }
  # ---------
  
  # ---------
  # Predict via Day filter
  # ---------
  {
    library(MASS)
    
    dat.sim[["press"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    dat.sim[["press"]][1,] <- dat.mod[1,"press"]
    
    pb <- txtProgressBar(min=abs(max(dat.mod$time.day2)), max=abs(min(dat.mod$time.day2)), style=3)
    set.seed(138)
    tic()
    for(i in max(dat.mod$time.day2):min(dat.mod$time.day2)){
      setTxtProgressBar(pb, abs(i))
      
      rows.now = which(dat.mod$time.day2==i)
      dat.temp <- dat.mod[rows.now,c("time.day2", "year", "doy", "hour", "tmax.day", "tmin.day", "tmean.day", "max.dep", "min.dep", "precipf.day", "swdown.day", "press.day", "press.day", "qair.day", "wind.day", "next.press")]
      dat.temp$press = -99999 # Dummy value so there's a column
      dat.temp <- dat.temp[complete.cases(dat.temp),]
      # day.now = dat.temp$doy
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day2)){
        sim.lag <- stack(data.frame(array(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour==0),"press"], dim=c(1, ncol(dat.sim$press)))))
        names(sim.lag) <- c("lag.press", "ens")
        
        # sim.lag$lag.tmin <- stack(data.frame(array(min(dat.mod[which(dat.mod$time.day2==i ),"press"]), dim=c(1, ncol(dat.sim$press)))))[,1]
        # sim.lag$lag.tmax <- stack(data.frame(array(max(dat.mod[which(dat.mod$time.day2==i ),"press"]), dim=c(1, ncol(dat.sim$press)))))[,1]
        
        # sim.lag <- stack(data.frame(array(mean(dat.mod[which(dat.mod$time.day2==i & dat.mod$hour<=2),"press"]), dim=c(1, ncol(dat.sim)))))
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["press"]][dat.mod$time.day2==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$press)))))
        names(sim.lag) <- c("lag.press", "ens")
        # sim.lag$lag.tmin <- stack(apply(dat.sim[["press"]][dat.mod$time.day2==(i+1),], 2, min))[,1]
        # sim.lag$lag.tmax <- stack(apply(dat.sim[["press"]][dat.mod$time.day2==(i+1),], 2, max))[,1]
        
        # sim.lag <- stack(apply(dat.sim[dat.mod$time.day2==(i+1)  & dat.mod$hour<=2,],2, mean))
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      # dat.temp$tair <- stack(dat.sim$tair[rows.now,])[,1]
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.press.doy[[paste(day.now)]]$model, 
                              betas=mod.press.doy[[paste(day.now)]]$betas, 
                              resid.err=T,
                              model.resid=mod.press.doy[[paste(day.now)]]$model.resid, 
                              betas.resid=mod.press.doy[[paste(day.now)]]$betas.resid, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$press), replace=T)
      # pred.prop <- array(dim=c(1,ncol(dat.sim)))
      for(j in 1:ncol(dat.sim$press)){
        dat.sim[["press"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
      dat.mod[rows.now,"mod.press"] <- apply(dat.sim[["press"]][rows.now,],1, mean)
      # dat.sim[i,] <- dat.pred
      # dat.mod[i,"mod.hr"] <- predict(mod.fill, dat.mod[i,])
      if(i>min(dat.mod$time.day2)){ 
        dat.mod[dat.mod$time.day2==i-1,"lag.press" ] <- dat.mod[dat.mod$time.day2==i & dat.mod$hour==0,"mod.press"] 
        dat.mod[dat.mod$time.day2==i-1,"lag.tmin" ] <- min(dat.mod[dat.mod$time.day2==i,"mod.press"] )
        dat.mod[dat.mod$time.day2==i-1,"lag.tmax" ] <- max(dat.mod[dat.mod$time.day2==i,"mod.press"] )
      }
      # if(i>min(dat.mod$time.day2)) dat.mod[dat.mod$time.day2==i-1,"lag.day"] <- mean(dat.mod[dat.mod$time.day2==i & dat.mod$hour<=2,"mod.press"])
    }
    # dat.mod$mod.press.mean <- apply(dat.sim[["press"]], 1, mean)
    dat.mod$mod.press.025   <- apply(dat.sim[["press"]], 1, quantile, 0.025)
    dat.mod$mod.press.975   <- apply(dat.sim[["press"]], 1, quantile, 0.975)
    summary(dat.mod)
  }
  toc()
  # ---------
  
  # ---------
  # Graph the output
  # ---------
  {
    for(y in unique(dat.mod$year)){
      png(file.path(fig.dir, paste0("press_", y, "_year.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          # geom_point(aes(x=date, y=tmax, color=as.factor(doy)), size=0.25, alpha=0.5) +
          # geom_point(aes(x=date, y=tmin, color=as.factor(doy)), size=0.25, alpha=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.press.025, ymax=mod.press.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.press), color="blue") +
          geom_point(aes(x=date, y=mod.press), color="blue", size=0.5) +
          geom_line(aes(x=date, y=press), color="black") +
          geom_point(aes(x=date, y=press), color="black", size=0.5) +
          # geom_vline(xintercept=seq(min(dat.mod$time.hr), max(dat.mod$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("press_", y, "_year_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_point(aes(x=press, y=mod.press), color="black", size=0.5) +
          geom_abline(slope=1, intercept=0, color="red") +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
    }  
    
    
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
    dat.graph4$season <- as.factor("fall")
    
    press.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    for(y in unique(press.graph$year)){
      png(file.path(fig.dir, paste0("press_",y,"_examples.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=press.graph[press.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_line(aes(x=date, y=press), color="black") +
          geom_point(aes(x=date, y=press), color="black", size=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.press.025, ymax=mod.press.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.press), color="blue") +
          geom_point(aes(x=date, y=mod.press), color="blue", size=0.5) +
          scale_y_continuous(name="Hourly Air Temperature") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      png(file.path(fig.dir, paste0("press_",y,"_examples_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=press.graph[press.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_point(aes(x=press, y=mod.press), color="black", size=0.5) +
          geom_abline(intercept=0, slope=1, color="red") +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
    }
  }
  # ---------
  
  rm(mod.press.doy) # Clear out the model to save memory
}
# ------------------------------------------


# ------------------------------------------
# Remaining vars
# ------------------------------------------
# # These variables have a strong diurnal trend
# test.swdown <- gam(swdown ~ s(hour) -1, data=dat.train)
# plot(test.swdown)
# 
# test.tair <- gam(tair ~ s(hour) -1, data=dat.train)
# plot(test.tair)
# 
# test.lwdown <- gam(lwdown ~ s(hour) -1, data=dat.train)
# plot(test.lwdown)

test.qair <- gam(qair ~ s(hour) -1, data=dat.train)
plot(test.qair)

test.wind <- gam(wind ~ s(hour) -1, data=dat.train)
plot(test.wind)

# Those below don't have a strong diurnal trend
test.precipf <- gam(precipf ~ s(hour) -1, data=dat.train)
plot(test.precipf)

test.press <- gam(press ~ s(hour) -1, data=dat.train)
plot(test.press)
# ------------------------------------------


# ------------------------------------------
# ------------------------------------------
# ------------------------------------------


# ------------------------------------------
# ------------------------------------------
# ------------------------------------------


# ------------------------------------------
# ------------------------------------------
# ------------------------------------------


