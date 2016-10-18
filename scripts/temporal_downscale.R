# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(lubridate)
library(ggplot2)
library(tictoc)
rm(list=ls())

dat.train <- read.csv("../data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")
# dat.train$doy <- as.ordered(dat.train$doy)

# Now trying things in a predictive framework
dat.train <- dat.train[order(dat.train$year, dat.train$doy, dat.train$hour, decreasing=T),]
dat.train[1:25,]
head(dat.train)
summary(dat.train)

dat.tair <- dat.train[,c("dataset", "year", "doy", "hour", "tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
dat.tair$date <- strptime(paste(dat.tair$year, dat.tair$doy+1, dat.tair$hour, sep="-"), "%Y-%j-%H", tz="GMT")
dat.tair$time.hr <- as.numeric(difftime(dat.tair$date, "2016-01-01", tz="GMT", units="hour"))
dat.tair$time.day <- as.numeric(difftime(dat.tair$date, "2016-01-01", tz="GMT", units="day"))+1/24
dat.tair$time.day2 <- as.integer(dat.tair$time.day)-1
dat.tair <- dat.tair[order(dat.tair$time.hr, decreasing=T),]
# dat.tair[1:25,]
head(dat.tair)

# For some reason certain days are getting an extra hour
for(i in max(dat.tair$time.day2):min(dat.tair$time.day2)){
  rows.now <- which(dat.tair$time.day2==i)
  if(length(rows.now)<=24) next
  
  dat.tair[rows.now[25],"time.day2"] <- i-1
}

tair.day <- aggregate(dat.tair[,c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")], by=dat.tair[,c("year", "doy")], FUN=mean)
names(tair.day)[which(names(tair.day)=="tair")] <- "tmean"
tair.day$tmax <- aggregate(dat.tair[,c("tair")], by=dat.tair[,c("year", "doy")], FUN=max)$x
tair.day$tmin <- aggregate(dat.tair[,c("tair")], by=dat.tair[,c("year", "doy")], FUN=min)$x
summary(tair.day)

dat.tair <- merge(dat.tair[,c("dataset", "year", "doy", "hour", "date", "time.hr", "time.day", "time.day2", "tair")], tair.day, all.x=T, all.y=T)
summary(dat.tair)

# adding in the 1-hour lag
lag.hr <- dat.tair
names(lag.hr)[which(names(lag.hr)=="tair")] <- "lag.hr"
lag.hr$time.hr <- lag.hr$time.hr-1
head(lag.hr)

dat.tair <- merge(dat.tair, lag.hr[,c("time.hr", "lag.hr")], all.x=T)
dat.tair <- dat.tair[order(dat.tair$time.hr, decreasing=T),]
dat.tair[is.na(dat.tair$lag.hr),"lag.hr"] <- dat.tair[1,"tair"]
head(dat.tair)
summary(dat.tair)

# creating a variable for 11 pm the day before (to work at daily timestep)
lag.day <- dat.tair[dat.tair$hour==0,]
names(lag.day)[which(names(lag.day)=="tair")] <- "lag.day" # The previousnext morning's temperature
lag.day$lag.diff <- dat.tair[dat.tair$hour==6,"tair"] - lag.day$lag.day # Lag is the change in temp in the proximate 3 hours
lag.day <- aggregate(lag.day[,c("lag.day", "lag.diff", "tmean", "tmax", "tmin")],
                     by=lag.day[,c("year", "doy", "time.day2")],
                     FUN=mean)
lag.day$lag.min <- aggregate(dat.tair[,c("tair")],  
                             by=dat.tair[,c("year", "doy", "time.day2")],
                             FUN=min)[,"x"] # Add in a lag for the next day's min temp
lag.day$lag.max <- aggregate(dat.tair[,c("tair")],  
                             by=dat.tair[,c("year", "doy", "time.day2")],
                             FUN=max)[,"x"] # Add in a lag for the next day's min temp
# lag.day <- aggregate(lag.day[,c("lag.day", "tmean", "tmax", "tmin")],
#                      by=lag.day[,c("year", "doy", "time.day2")],
#                      FUN=function(x){max(x)-min(x)})
lag.day$time.day2 <- lag.day$time.day2-1
# lag.day[is.na(lag.day$lag.day), "lag.day"] <- dat.tair[1,"tair"]
head(lag.day)
summary(lag.day)

dat.tair <- merge(dat.tair, lag.day[,c("time.day2", "lag.day", "lag.diff", "lag.min", "lag.max")], all.x=T)
dat.tair <- dat.tair[order(dat.tair$time.hr, decreasing=T),]
# dat.tair[is.na(dat.tair$lag.day), "lag.day"] <- dat.tair[1,"tair"]
# dat.tair[is.na(dat.tair$lag.diff), "lag.diff"] <- mean(dat.tair[,"lag.diff"],na.rm=T)
# dat.tair <- dat.tair[complete.cases(dat.tair),]
head(dat.tair)
summary(dat.tair)

# Lookign at max & min as departure from mean
dat.tair$max.dep <- dat.tair$tmax - dat.tair$tmean
dat.tair$min.dep <- dat.tair$tmin - dat.tair$tmean
summary(dat.tair)
# ------------------------------------------

# dat.test <- dat.tair[dat.tair$hour<=24 & dat.tair$doy<=5,]
# 
# gam.test <- gam(tair ~ s(hour, by=as.factor(doy), k=4) + as.factor(doy), data=dat.test)
# plot(gam.test, pages=1)

# ------------------------------------------
# Generating all the daily models and storing it into an easy-to-find list
# ------------------------------------------
source("temporal_downscale_functions.R")
tic()
mod.tair.doy <- model.tair(dat.tair=dat.tair[,], parallel=T, n.cores=4, n.beta=250)
toc()
length(mod.tair.doy)
# summary(mod.tair.doy[[1]]$model)

# mod.tair.doy[[1]]$model$call
# ------------------------------------------

# ------------------------------------------
# Looking at the residuals from the model
# ------------------------------------------
{
  for(i in names(mod.tair.doy)){
    if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
    dat.tair[dat.tair$doy==as.numeric(i) & !is.na(dat.tair$lag.day), "resid"] <- resid(mod.tair.doy[[i]]$model)
    dat.tair[dat.tair$doy==as.numeric(i) & !is.na(dat.tair$lag.day), "predict"] <- predict(mod.tair.doy[[i]]$model)
  }
  summary(dat.tair)

  png("Resid_vs_Hour.png")
  plot(resid ~ hour, data=dat.tair, cex=0.5); abline(h=0, col="red")
  dev.off()
  png("Resid_vs_DOY.png")
  plot(resid ~ doy, data=dat.tair, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
  dev.off()
  # plot(resid ~ predict, data=dat.tair); abline(h=0, col="red")
  png("Resid_vs_Obs.png")
  plot(resid ~ tair, data=dat.tair, cex=0.5); abline(h=0, col="red")
  dev.off()
  png("Predict_vs_Obs.png")
  plot(predict ~ tair, data=dat.tair, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()

  # Looking at the daily maxes & mins
  day.stats <- aggregate(dat.tair[,c("tair", "tmax", "tmin", "resid", "predict")],
                         by=dat.tair[,c("time.day2", "year", "doy")],
                         FUN=mean)
  day.stats$mod.max <- aggregate(dat.tair[,c("predict")],
                                 by=dat.tair[,c("time.day2", "year", "doy")],
                                 FUN=max)[,"x"]
  day.stats$mod.min <- aggregate(dat.tair[,c("predict")],
                                 by=dat.tair[,c("time.day2", "year", "doy")],
                                 FUN=min)[,"x"]
  day.stats$max.resid <- day.stats$mod.max - day.stats$tmax
  day.stats$min.resid <- day.stats$mod.min - day.stats$tmin

  png("Predict_vs_Tmax.png")
  plot(mod.max ~ tmax, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()
  png("Predict_vs_Tmin.png")
  plot(mod.min ~ tmin, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()
  png("Predict_vs_Tmean.png")
  plot(predict ~ tair, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()

  png("ResidualHistograms.png", height=6, width=8, units="in", res=220)
  par(mfrow=c(3,1))
  hist(day.stats$resid, main="Daily Mean Residuals")
  hist(day.stats$min.resid, main="Daily Min Residuals")
  hist(day.stats$max.resid, main="Daily Max Residuals")
  par(mfrow=c(1,1))
  dev.off()
}
# ------------------------------------------


# ------------------------------------------
# Testing the run with propogating uncertainty
# ------------------------------------------
source("temporal_downscale_functions.R")

# Set up data
tair.mod <- dat.tair[dat.tair$year>=2010,]
tair.mod[, "lag.hr"] <- NA
tair.mod[, "lag.day"] <- NA
tair.mod[, "mod.hr"] <- NA
tair.mod[, "mod.day"] <- NA
# head(tair.mod)

# Initializing the lags
# tair.mod[tair.mod$time.hr==max(tair.mod$time.hr),"lag.hr"] <- tair.mod[tair.mod$time.hr==max(tair.mod$time.hr),"tair"]
# tair.mod[tair.mod$time.day2==max(tair.mod$time.day2),"lag.day"] <- mean(tair.mod[tair.mod$time.day2==max(tair.mod$time.day2) & tair.mod$hour==0,"tair"], na.rm=T)
# tair.mod[tair.mod$time.day2==max(tair.mod$time.day2),"lag.diff"] <- mean(dat.tair[,"lag.diff"], na.rm=T)
# head(tair.mod)
# tair.mod[1:27,]
# tair.mod[,c("year", "doy", "hour", "tair", "tair2", "mod3", "time")]

n.ens=50
dat.sim <- data.frame(array(dim=c(nrow(tair.mod), n.ens)))
dat.sim[1,] <- tair.mod[1,"tair"]
# dat.simstack <- 


# ---------
# Hour filter
# ---------
{
  # library(MASS)
  # # Do the calculation by Hour
  # pb <- txtProgressBar(min=1, max=nrow(tair.mod), style=3)
  # set.seed(138)
  # tic()
  # for(i in 1:nrow(tair.mod)){
  #   setTxtProgressBar(pb, i)
  # 
  #   dat.temp <- tair.mod[i,c("doy", "hour", "tmax", "tmin","tmean", "max.dep", "min.dep", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
  #   dat.temp$tair = -99999 # Dummy value so there's a column
  #   day.now = dat.temp$doy
  #   
  #   # Set up the lags
  #   if(i==1){
  #     sim.lag <- stack(data.frame(array(dat.sim[i,], dim=c(1, ncol(dat.sim)))))
  #   } else {
  #     sim.lag <- stack(data.frame(array(dat.sim[i-1,], dim=c(1, ncol(dat.sim)))))
  #   }
  #   names(sim.lag) <- c("lag.hr", "ens")
  #   dat.temp <- merge(dat.temp, sim.lag, all.x=T)
  # 
  # 
  #   dat.pred <- predict.met(newdata=dat.temp, mod.predict=mod.tair.doy[[paste(day.now)]]$model, betas=mod.tair.doy[[paste(day.now)]]$betas, n.ens=n.ens)
  #   
  #   # Randomly pick which values to save & propogate
  #   cols.prop <- sample(1:n.ens, ncol(dat.sim), replace=T)
  #   # pred.prop <- array(dim=c(1,ncol(dat.sim)))
  #   for(j in 1:ncol(dat.sim)){
  #     dat.sim[i,j] <- dat.pred[j,cols.prop[j]]
  #   }
  # 
  #   tair.mod[i,"mod.hr"] <- mean(as.numeric(dat.sim[i,]))
  #   # dat.sim[i,] <- dat.pred
  #   # tair.mod[i,"mod.hr"] <- predict(mod.fill, tair.mod[i,])
  #   if(i<nrow(tair.mod)) tair.mod[i+1,"lag.hr"] <- tair.mod[i,"mod.hr"]
  # }
  # tair.mod$mod.hr2 <- apply(dat.sim, 1, mean)
  # tair.mod$mod.hr.low <- apply(dat.sim, 1, quantile, 0.025)
  # tair.mod$mod.hr.hi <- apply(dat.sim, 1, quantile, 0.975)
  # summary(tair.mod)
  # toc()
}
# ---------

# ---------
# Day filter
# ---------
# {
  library(MASS)
  # Do the calculation by Hour
  pb <- txtProgressBar(min=abs(max(tair.mod$time.day2)), max=abs(min(tair.mod$time.day2)), style=3)
  set.seed(138)
  tic()
  for(i in max(tair.mod$time.day2):min(tair.mod$time.day2)){
    setTxtProgressBar(pb, abs(i))
    
    rows.now = which(tair.mod$time.day2==i)
    dat.temp <- tair.mod[rows.now,c("time.day2", "year", "doy", "hour", "tmax", "tmin","tmean", "max.dep", "min.dep", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
    dat.temp$tair = -99999 # Dummy value so there's a column
    dat.temp <- dat.temp[complete.cases(dat.temp),]
    # day.now = dat.temp$doy
    day.now = unique(dat.temp$doy)
    
    # Set up the lags
    if(i==max(tair.mod$time.day2)){
      sim.lag <- stack(data.frame(array(tair.mod[which(tair.mod$time.day2==i & tair.mod$hour==0),"tair"], dim=c(1, ncol(dat.sim)))))
      names(sim.lag) <- c("lag.day", "ens")
      
      sim.lag$lag.diff <- stack(data.frame(array(mean(tair.mod[which(tair.mod$time.day2==i & tair.mod$hour==0),"lag.diff"]), dim=c(1, ncol(dat.sim)))))[,1]
      sim.lag$lag.min <- stack(data.frame(array(min(tair.mod[which(tair.mod$time.day2==i ),"tair"]), dim=c(1, ncol(dat.sim)))))[,1]
      sim.lag$lag.max <- stack(data.frame(array(max(tair.mod[which(tair.mod$time.day2==i ),"tair"]), dim=c(1, ncol(dat.sim)))))[,1]

      # sim.lag <- stack(data.frame(array(mean(tair.mod[which(tair.mod$time.day2==i & tair.mod$hour<=2),"tair"]), dim=c(1, ncol(dat.sim)))))
    } else {
      sim.lag <- stack(data.frame(array(dat.sim[tair.mod$time.day2==(i+1)  & tair.mod$hour==0,], dim=c(1, ncol(dat.sim)))))
      names(sim.lag) <- c("lag.day", "ens")
      sim.lag$lag.diff <- stack(data.frame(array(dat.sim[tair.mod$time.day2==(i+1)  & tair.mod$hour==6,] - dat.sim[tair.mod$time.day2==(i+1)  & tair.mod$hour==0,], dim=c(1, ncol(dat.sim)))))[,1]
      sim.lag$lag.min <- stack(apply(dat.sim[tair.mod$time.day2==(i+1),], 2, min))[,1]
      sim.lag$lag.max <- stack(apply(dat.sim[tair.mod$time.day2==(i+1),], 2, max))[,1]
      
      # sim.lag <- stack(apply(dat.sim[tair.mod$time.day2==(i+1)  & tair.mod$hour<=2,],2, mean))
    }
    dat.temp <- merge(dat.temp, sim.lag, all.x=T)
    
    
    dat.pred <- predict.met(newdata=dat.temp, mod.predict=mod.tair.doy[[paste(day.now)]]$model, betas=mod.tair.doy[[paste(day.now)]]$betas, n.ens=n.ens)
    
    # Randomly pick which values to save & propogate
    cols.prop <- sample(1:n.ens, ncol(dat.sim), replace=T)
    # pred.prop <- array(dim=c(1,ncol(dat.sim)))
    for(j in 1:ncol(dat.sim)){
      dat.sim[rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
    }
    
    tair.mod[rows.now,"mod.day"] <- apply(dat.sim[rows.now,],1, mean)
    # dat.sim[i,] <- dat.pred
    # tair.mod[i,"mod.hr"] <- predict(mod.fill, tair.mod[i,])
    if(i>min(tair.mod$time.day2)){ 
      tair.mod[tair.mod$time.day2==i-1,"lag.day" ] <- tair.mod[tair.mod$time.day2==i & tair.mod$hour==0,"mod.day"] 
      tair.mod[tair.mod$time.day2==i-1,"lag.diff"] <- tair.mod[tair.mod$time.day2==i & tair.mod$hour==6,"mod.day"] - tair.mod[tair.mod$time.day2==i & tair.mod$hour==0,"mod.day"] 
      tair.mod[tair.mod$time.day2==i-1,"lag.min" ] <- min(tair.mod[tair.mod$time.day2==i,"mod.day"] )
      tair.mod[tair.mod$time.day2==i-1,"lag.max" ] <- max(tair.mod[tair.mod$time.day2==i,"mod.day"] )
    }
    # if(i>min(tair.mod$time.day2)) tair.mod[tair.mod$time.day2==i-1,"lag.day"] <- mean(tair.mod[tair.mod$time.day2==i & tair.mod$hour<=2,"mod.day"])
  }
  tair.mod$mod.day2 <- apply(dat.sim, 1, mean)
  tair.mod$mod.day.low <- apply(dat.sim, 1, quantile, 0.025)
  tair.mod$mod.day.hi <- apply(dat.sim, 1, quantile, 0.975)
  summary(tair.mod)
  toc()
# }
# ---------
for(y in unique(tair.mod$year)){
  png(paste0("Tair_", y, "_year.png"), height=8, width=10, units="in", res=220)
  print(
    ggplot(data=tair.mod[tair.mod$year==y,]) +
      # geom_point(aes(x=date, y=tmax, color=as.factor(doy)), size=0.25, alpha=0.5) +
      # geom_point(aes(x=date, y=tmin, color=as.factor(doy)), size=0.25, alpha=0.5) +
      geom_ribbon(aes(x=date, ymin=mod.day.low, ymax=mod.day.hi), alpha=0.5, fill="blue") +
      geom_line(aes(x=date, y=mod.day), color="blue") +
      geom_point(aes(x=date, y=mod.day), color="blue", size=0.5) +
      geom_line(aes(x=date, y=tair), color="black") +
      geom_point(aes(x=date, y=tair), color="black", size=0.5) +
      # geom_vline(xintercept=seq(min(tair.mod$time.hr), max(tair.mod$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
      scale_x_datetime(expand=c(0,0)) +
      theme_bw()
  )
  dev.off() 
}  
  
  
dat.graph1 <- tair.mod[tair.mod$doy>=32 & tair.mod$doy<=(32+14),]
dat.graph1$season <- as.factor("winter")
dat.graph2 <- tair.mod[tair.mod$doy>=123 & tair.mod$doy<=(123+14),]
dat.graph2$season <- as.factor("spring")
dat.graph3 <- tair.mod[tair.mod$doy>=214 & tair.mod$doy<=(213+14),]
dat.graph3$season <- as.factor("summer")
dat.graph4 <- tair.mod[tair.mod$doy>=305 & tair.mod$doy<=(305+14),]
dat.graph4$season <- as.factor("fall")

tair.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)

for(y in unique(tair.graph$year)){
  png(paste0("Tair_",y,"_examples.png"), height=8, width=10, units="in", res=220)
  print(
  ggplot(data=tair.graph[tair.graph$year==y,]) +
    facet_wrap(~season, scales="free") +
    # geom_point(aes(x=date, y=tmax), size=0.1, alpha=0.5, color="gray50") +
    # geom_point(aes(x=date, y=tmin), size=0.1, alpha=0.5, color="gray50") +
    geom_line(aes(x=date, y=tair), color="black") +
    geom_point(aes(x=date, y=tair), color="black", size=0.5) +
    geom_ribbon(aes(x=date, ymin=mod.day.low, ymax=mod.day.hi), alpha=0.5, fill="blue") +
    geom_line(aes(x=date, y=mod.day), color="blue") +
    geom_point(aes(x=date, y=mod.day), color="blue", size=0.5) +
    scale_y_continuous(name="Hourly Air Temperature") +
    scale_x_datetime(expand=c(0,0)) +
    theme_bw()
  )
  dev.off()
}
# ------------------------------------------

