# Script to prototype temporal downscaling
library(ncdf4)
library(mgcv)
library(lubridate)
library(ggplot2)
rm(list=ls())

dat.train <- read.csv("../data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")
# dat.train$doy <- as.ordered(dat.train$doy)

# Now trying things in a predictive framework
dat.train <- dat.train[order(dat.train$year, dat.train$doy, dat.train$hour, decreasing=T),]
dat.train[1:25,]
head(dat.train)
summary(dat.train)

dat.tair3 <- dat.train[,c("dataset", "year", "doy", "hour", "tair")]
dat.tair3$date <- strptime(paste(dat.tair3$year, dat.tair3$doy+1, dat.tair3$hour, sep="-"), "%Y-%j-%H", tz="GMT")
dat.tair3$time.hr <- as.numeric(difftime(dat.tair3$date, "2016-01-01", tz="GMT", units="hour"))
dat.tair3$time.day <- as.numeric(difftime(dat.tair3$date, "2016-01-01", tz="GMT", units="day"))+1/24
dat.tair3$time.day2 <- as.integer(dat.tair3$time.day)-1
dat.tair3 <- dat.tair3[order(dat.tair3$time.hr, decreasing=T),]
# dat.tair3[1:25,]
head(dat.tair3)

tair.day <- aggregate(dat.tair3[,c("tair")], by=dat.tair3[,c("year", "doy")], FUN=mean)
names(tair.day)[which(names(tair.day)=="x")] <- "tmean"
tair.day$tmax <- aggregate(dat.tair3[,c("tair")], by=dat.tair3[,c("year", "doy")], FUN=max)$x
tair.day$tmin <- aggregate(dat.tair3[,c("tair")], by=dat.tair3[,c("year", "doy")], FUN=min)$x
summary(tair.day)

dat.tair3 <- merge(dat.tair3, tair.day, all.x=T, all.y=T)
summary(dat.tair3)

# adding in the 1-hour lag
lag.hr <- dat.tair3
names(lag.hr)[which(names(lag.hr)=="tair")] <- "lag.hr"
lag.hr$time.hr <- lag.hr$time.hr-1
head(lag.hr)

dat.tair3 <- merge(dat.tair3, lag.hr[,c("time.hr", "lag.hr")], all.x=T)
dat.tair3 <- dat.tair3[order(dat.tair3$time.hr, decreasing=T),]
dat.tair3[is.na(dat.tair3$lag.hr),"lag.hr"] <- dat.tair3[1,"tair"]
head(dat.tair3)
summary(dat.tair3)

# creating a variable for 11 pm the day before (to work at daily timestep)
lag.day <- dat.tair3[dat.tair3$hour>=21,]
names(lag.day)[which(names(lag.day)=="tair")] <- "lag.day"
lag.day <- aggregate(lag.day[,c("lag.day", "tmean", "tmax", "tmin")],
                     by=lag.day[,c("year", "doy", "time.day2")],
                     FUN=mean)
lag.day$time.day2 <- lag.day$time.day2-1
lag.day[is.na(lag.day$lag.day), "lag.day"] <- dat.tair3[1,"tair"]
head(lag.day)
summary(lag.day)

dat.tair3 <- merge(dat.tair3, lag.day[,c("time.day2", "lag.day")], all.x=T)
dat.tair3 <- dat.tair3[order(dat.tair3$time.hr, decreasing=T),]
dat.tair3[is.na(dat.tair3$lag.day), "lag.day"] <- dat.tair3[1,"tair"]
head(dat.tair3)
summary(dat.tair3)


dat.tair4 <- dat.tair3[dat.tair3$year==2015 & dat.tair3$doy<31,c("dataset", "year", "doy", "hour", "date", "time.hr", "time.day", "time.day2", "tair", "tmax", "tmin", "tmean", "lag.hr", "lag.day")]
dat.tair4[, "lag.hr"] <- NA
dat.tair4[, "lag.day"] <- NA
dat.tair4[, "mod.hr"] <- NA
dat.tair4[, "mod.day"] <- NA
head(dat.tair4)

# Initializing the lags
dat.tair4[dat.tair4$time.hr==max(dat.tair4$time.hr),"lag.hr"] <- dat.tair4[dat.tair4$time.hr==max(dat.tair4$time.hr),"tair"]
# for(i in max(dat.tair4[dat.tair4$time.day2==max(dat.tair4$time.day2),"time.hr"]):min(dat.tair4[dat.tair4$time.day2==max(dat.tair4$time.day2),"time.hr"])){
#   dat.tair4[dat.tair4$time.hr==i-1,"lag.hr"] <- dat.tair4[dat.tair4$time.hr==i,"tair"]
# }
# dat.tair4[dat.tair4$time.hr!=max(dat.tair4$time.hr),"lag.hr"] <- dat.tair4[dat.tair4$time.hr==max(dat.tair4$time.hr),"tair"]
dat.tair4[dat.tair4$time.day2==max(dat.tair4$time.day2),"lag.day"] <- dat.tair4[dat.tair4$time.day2==max(dat.tair4$time.day2) & dat.tair4$hour==23,"tair"]
head(dat.tair4)
dat.tair4[1:27,]
# dat.tair4[,c("year", "doy", "hour", "tair", "tair2", "mod3", "time")]

n=30
dat.sim <- array(dim=c(nrow(dat.tair4), n))

library(MASS)
# Do the calculation by Hour
pb <- txtProgressBar(min=1, max=nrow(dat.tair4), style=3)
for(i in 1:nrow(dat.tair4)){
  setTxtProgressBar(pb, i)
  day.now = dat.tair4[i,"doy"]
  # mod.fill <- gam(tair ~ s(hour, by=as.factor(doy)) + tmax*tmin*tmean  -1, data=dat.tair3[dat.tair3$doy==day.now,])
  mod.fill <- lm(tair ~ as.factor(hour)*lag.hr*(tmax*tmin +tmean) - as.factor(hour)-1, data=dat.tair3[dat.tair3$doy==day.now,])
  # summary(mod.fill)
  # plot(mod.fill)
  mod.terms <- terms(mod.fill)
  # Terms <- delete.response(tt)
  # p <- mod.fill$rank
  # p1 <- seq_len(p)
  # piv <- if (p)  qr(mod.fill)$pivot[p1]
  # drop(X[, piv, drop = FALSE] %*% beta[piv])
  # 
  m <- model.frame(mod.terms, dat.tair4[i,], xlev = mod.fill$xlevels)
  Xp <- model.matrix(mod.terms, m, contrasts.arg = mod.fill$contrasts)
  # 
  mod.coef <- coef(mod.fill)
  mod.cov  <- vcov(mod.fill)
  piv <- as.numeric(which(!is.na(mod.coef)))
  # piv2 <- as.numeric(which(!is.na(mod.coef)))
  
  Rbeta <- mvrnorm(n=n, mod.coef[piv], mod.cov)
  # mod.resid <- resid(mod.fill)
  # 
  dat.sim[i,] <- Xp[,piv] %*% t(Rbeta) #+ rnorm(n, mean=mean(mod.resid), sd=sd(mod.resid))
  # test <- Xp[,piv] %*% t(Rbeta)
  # summary(out)
  
  # dat.tair4[i,"mod.hr"] <- mean(dat.sim[i,])
  dat.tair4[i,"mod.hr"] <- predict(mod.fill, dat.tair4[i,])
  if(i<nrow(dat.tair4)) dat.tair4[i+1,"lag.hr"] <- dat.tair4[i,"mod.hr"]
}
dat.tair4$mod.hr2 <- apply(dat.sim, 1, mean)
dat.tair4$mod.hr.low <- apply(dat.sim, 1, quantile, 0.025)
dat.tair4$mod.hr.hi <- apply(dat.sim, 1, quantile, 0.975)
summary(dat.tair4)

plot(tair ~ time.hr, data=dat.tair4[dat.tair4$doy<14,], type="l", lwd=2)
points(tair ~ time.hr, data=dat.tair4[dat.tair4$doy<14,], pch=19, cex=0.5)
abline(v=seq(min(dat.tair4$time.hr), max(dat.tair4$time.hr), by=24)-0.5, col="gray50", cex=0.5, lty="dashed")
# lines(mod.day ~ abs(time.hr), data=dat.tair4, pch=19, cex=0.5, col="red")
# points(mod.day ~ abs(time.hr), data=dat.tair4, pch=19, cex=0.5, col="red")
lines(mod.hr ~ time.hr, data=dat.tair4[dat.tair4$doy<14,], pch=19, cex=0.5, col="blue")
points(mod.hr ~ time.hr, data=dat.tair4[dat.tair4$doy<14,], pch=19, cex=0.5, col="blue")

png("Tair_test.png")
ggplot(data=dat.tair4[dat.tair4$doy<14,]) +
  geom_ribbon(aes(x=date, ymin=mod.hr.low, ymax=mod.hr.hi), alpha=0.5, fill="blue") +
  geom_line(aes(x=date, y=mod.hr), color="blue") +
  geom_point(aes(x=date, y=mod.hr), color="blue", size=0.5) +
  geom_line(aes(x=date, y=tair), color="black") +
  geom_point(aes(x=date, y=tair), color="black", size=0.5) +
  # geom_vline(xintercept=seq(min(dat.tair4$time.hr), max(dat.tair4$time.hr), by=24)-0.5, linetype="dashed", color="gray50") +
  # scale_x_continuous(limits=range(dat.tair4[dat.tair4$doy<14,"time.day"]), expand=c(0,0)) +
  theme_bw()
dev.off()

