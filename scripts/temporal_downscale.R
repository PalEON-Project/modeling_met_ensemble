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
lag.day <- dat.tair[dat.tair$hour>=21,]
names(lag.day)[which(names(lag.day)=="tair")] <- "lag.day"
lag.day <- aggregate(lag.day[,c("lag.day", "tmean", "tmax", "tmin")],
                     by=lag.day[,c("year", "doy", "time.day2")],
                     FUN=mean)
lag.day$time.day2 <- lag.day$time.day2-1
lag.day[is.na(lag.day$lag.day), "lag.day"] <- dat.tair[1,"tair"]
head(lag.day)
summary(lag.day)

dat.tair <- merge(dat.tair, lag.day[,c("time.day2", "lag.day")], all.x=T)
dat.tair <- dat.tair[order(dat.tair$time.hr, decreasing=T),]
dat.tair[is.na(dat.tair$lag.day), "lag.day"] <- dat.tair[1,"tair"]
head(dat.tair)
summary(dat.tair)

# Lookign at max & min as departure from mean
dat.tair$max.dep <- dat.tair$tmax - dat.tair$tmean
dat.tair$min.dep <- dat.tair$tmin - dat.tair$tmean
summary(dat.tair)

# ------------------------------------------
# Generating all the daily models and storing it into an easy-to-find list
# ------------------------------------------
source("temporal_downscale_functions.R")
tic()
mod.tair.doy <- model.tair(dat.tair=dat.tair[,], parallel=T, n.cores=8)
toc()
length(mod.tair.doy)
# summary(mod.tair.doy[[1]]$model)

# mod.tair.doy[[1]]$model$call
# ------------------------------------------



# ------------------------------------------
# Testing the run with propogating uncertainty
# ------------------------------------------
source("temporal_downscale_functions.R")

# Set up data
tair.mod <- dat.tair[dat.tair$year==2015,]
tair.mod[, "lag.hr"] <- NA
tair.mod[, "lag.day"] <- NA
tair.mod[, "mod.hr"] <- NA
tair.mod[, "mod.day"] <- NA
# head(tair.mod)

# Initializing the lags
tair.mod[tair.mod$time.hr==max(tair.mod$time.hr),"lag.hr"] <- tair.mod[tair.mod$time.hr==max(tair.mod$time.hr),"tair"]
# head(tair.mod)
# tair.mod[1:27,]
# tair.mod[,c("year", "doy", "hour", "tair", "tair2", "mod3", "time")]

n.ens=10
dat.sim <- data.frame(array(dim=c(nrow(tair.mod), n.ens)))
dat.sim[1,] <- tair.mod[1,"tair"]
# dat.simstack <- 

library(MASS)
# Do the calculation by Hour
pb <- txtProgressBar(min=1, max=nrow(tair.mod), style=3)
set.seed(138)
tic()
for(i in 1:nrow(tair.mod)){
  setTxtProgressBar(pb, i)

  dat.temp <- tair.mod[i,c("doy", "hour", "tmax", "tmin","tmean", "max.dep", "min.dep", "precipf", "swdown", "lwdown", "press", "qair", "wind")]
  dat.temp$tair = -99999 # Dummy value so there's a column
  day.now = dat.temp$doy
  
  # Set up the lags
  if(i==1){
    sim.lag <- stack(data.frame(array(dat.sim[i,], dim=c(1, ncol(dat.sim)))))
  } else {
    sim.lag <- stack(data.frame(array(dat.sim[i-1,], dim=c(1, ncol(dat.sim)))))
  }
  names(sim.lag) <- c("lag.hr", "ens")
  dat.temp <- merge(dat.temp, sim.lag, all.x=T)


  dat.pred <- predict.met(newdata=dat.temp, mod.predict=mod.tair.doy[[paste(day.now)]]$model, betas=mod.tair.doy[[paste(day.now)]]$betas, n.ens=n.ens)
  
  # Randomly pick which values to save & propogate
  cols.prop <- sample(1:n.ens, ncol(dat.sim), replace=T)
  # pred.prop <- array(dim=c(1,ncol(dat.sim)))
  for(j in 1:ncol(dat.sim)){
    dat.sim[i,j] <- dat.pred[j,cols.prop[j]]
  }

  tair.mod[i,"mod.hr"] <- mean(as.numeric(dat.sim[i,]))
  # dat.sim[i,] <- dat.pred
  # tair.mod[i,"mod.hr"] <- predict(mod.fill, tair.mod[i,])
  if(i<nrow(tair.mod)) tair.mod[i+1,"lag.hr"] <- tair.mod[i,"mod.hr"]
}
tair.mod$mod.hr2 <- apply(dat.sim, 1, mean)
tair.mod$mod.hr.low <- apply(dat.sim, 1, quantile, 0.025)
tair.mod$mod.hr.hi <- apply(dat.sim, 1, quantile, 0.975)
summary(tair.mod)
toc()

# tair.mod1 <- tair.mod

# ggplot(data=tair.mod[,]) +
#   # geom_point(aes(x=date, y=tmax, color=as.factor(doy)), size=0.25, alpha=0.5) +
#   # geom_point(aes(x=date, y=tmin, color=as.factor(doy)), size=0.25, alpha=0.5) +
#   geom_ribbon(aes(x=date, ymin=mod.hr.low, ymax=mod.hr.hi), alpha=0.5, fill="blue") +
#   geom_line(aes(x=date, y=mod.hr), color="blue") +
#   geom_point(aes(x=date, y=mod.hr), color="blue", size=0.5) +
#   geom_line(aes(x=date, y=tair), color="black") +
#   geom_point(aes(x=date, y=tair), color="black", size=0.5) +
#   # geom_vline(xintercept=seq(min(tair.mod$time.hr), max(tair.mod$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
#   # scale_x_continuous(name="Hourly Air Temperature", expand=c(0,0)) +
#   theme_bw()


dat.graph1 <- tair.mod[tair.mod$doy>=32 & tair.mod$doy<=(32+14),]
dat.graph1$season <- as.factor("winter")
dat.graph2 <- tair.mod[tair.mod$doy>=123 & tair.mod$doy<=(123+14),]
dat.graph2$season <- as.factor("spring")
dat.graph3 <- tair.mod[tair.mod$doy>=214 & tair.mod$doy<=(213+14),]
dat.graph3$season <- as.factor("summer")
dat.graph4 <- tair.mod[tair.mod$doy>=305 & tair.mod$doy<=(305+14),]
dat.graph4$season <- as.factor("fall")

tair.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)

png("Tair_2015_examples.png", height=8, width=10, units="in", res=220)
ggplot(data=tair.graph[,]) +
  facet_wrap(~season, scales="free") +
  # geom_point(aes(x=date, y=tmax), size=0.1, alpha=0.5, color="gray50") +
  # geom_point(aes(x=date, y=tmin), size=0.1, alpha=0.5, color="gray50") +
  geom_line(aes(x=date, y=tair), color="black") +
  geom_point(aes(x=date, y=tair), color="black", size=0.5) +

  geom_ribbon(aes(x=date, ymin=mod.hr.low, ymax=mod.hr.hi), alpha=0.5, fill="blue") +
  geom_line(aes(x=date, y=mod.hr), color="blue") +
  geom_point(aes(x=date, y=mod.hr), color="blue", size=0.5) +
  scale_y_continuous(name="Hourly Air Temperature")+
  theme_bw()
dev.off()

# ------------------------------------------

