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
library(mgcv); library(ggplot2)

# Set the working directory
wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble"
# wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

# Source to my gamm helper scripts
# Github: https://github.com/crollinson/R_Functions
path.gam <- "~/Desktop/Research/R_Functions"

# Defining a site name -- this can go into a function later
site.name="HARVARD"

path.dat <- file.path(wd.base, "data/paleon_sites", site.name)
path.out <- file.path(wd.base, "data/met_ensembles", site.name)
if(!dir.exists(path.out)) dir.create(path.out, recursive=T)  

met.done <- dir(path.out, ".csv")
# -----------------------------------

# -----------------------------------
# 1. Read in & format the different datasets
# -----------------------------------
#    - aggregate all variable to daily for the bias-correction
#    - do some exploratory graphing along the way
ldas     <- read.csv(file.path(path.dat, "NLDAS_1980-2015.csv"))
cruncep  <- read.csv(file.path(path.dat, "CRUNCEP_1901-2010.csv"))
gcm.p1k  <- read.csv(file.path(path.dat, "MIROC-ESM_historical_1850-2005.csv"))
gcm.hist <- read.csv(file.path(path.dat, "MIROC-ESM_p1000_850-1849.csv"))

# Adding an hour field to the gcm; setting as noon for simplicity
gcm.p1k$hour  <- 12.00
gcm.hist$hour <- 12.00

# making sure all datasets have tair, tmax, and tmin
ldas$tmax    <- ldas$tmin    <- ldas$tair
cruncep$tmax <- cruncep$tmin <- cruncep$tair
gcm.p1k$tair  <- apply(gcm.p1k [,c("tmax", "tmin")], 1, FUN=mean)
gcm.hist$tair <- apply(gcm.hist[,c("tmax", "tmin")], 1, FUN=mean)

# ******* TYPO CORRECTIONS!! *******
# **** NOTE: WILL NEED TO CANCEL THIS OUT UPON NEW EXTRACTION
# make first day of year start on 0
ldas$doy    <- ldas$doy - 1
cruncep$doy <- cruncep$doy -1
ldas$precipf <- ldas$precipf*.1

summary(ldas)
summary(cruncep)
summary(gcm.p1k)
summary(gcm.hist)

vars.met <- c("tair", "tmax", "tmin", "precipf", "press", "qair", "wind", "swdown", "lwdown")
cols.bind <- c("dataset", "year", "doy", "hour", vars.met)
met.all <- rbind(ldas[,cols.bind], cruncep[,cols.bind], gcm.p1k[,cols.bind], gcm.hist[,cols.bind])
met.all$Date <- as.Date(met.all$doy, origin=as.Date(paste(met.all$year, "01", "01", sep="-")))
summary(met.all)

# ----------------
# Creating a daily met record
# ----------------
met.day <- aggregate(met.all[,vars.met], by=met.all[,c("dataset", "year", "doy")], FUN=mean)

# getting tmax & tmin for ldas & cru
met.day2 <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                      by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
                      FUN=max)
names(met.day2)[4] <- "tmax"
met.day2$tmin <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                           by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
                           FUN=min)[,4]
summary(met.day2)

# merging tmax & tmin back into met.day 
met.day[met.day$dataset %in% c("NLDAS", "CRUNCEP"), "tmax"] <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                                                                         by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
                                                                         FUN=max)[,4]
met.day[met.day$dataset %in% c("NLDAS", "CRUNCEP"), "tmin"] <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                                                                         by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
                                                                         FUN=min)[,4]
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

met.year$dataset <- factor(met.year$dataset, levels=c("MIROC-ESM.p1000", "MIROC-ESM.hist", "CRUNCEP", "NLDAS"))
met.doy $dataset <- factor(met.doy $dataset, levels=c("MIROC-ESM.p1000", "MIROC-ESM.hist", "CRUNCEP", "NLDAS"))
met.year$met     <- factor(met.year$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
met.doy $met     <- factor(met.doy $met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))

# Putting in a dummy upper bound for precip to deal with some of the real oddballs
precip.cutoff <- quantile(met.doy[met.doy$met=="precipf","upr"], 0.95)
met.doy$upr2 <- ifelse(met.doy$met=="precipf" & met.doy$upr>=precip.cutoff, precip.cutoff, met.doy$upr)

png(file.path(path.out, "Met_Raw_Year_0850-2015.png"), height=11, width=8.5, "in", res=180)
ggplot(data=met.year) +
  facet_grid(met~., scales="free_y") +
  geom_line(aes(x=year, y=value, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
dev.off()

png(file.path(path.out, "Met_Raw_DOY_All.png"), height=11, width=8.5, "in", res=180)
ggplot(data=met.doy) +
  facet_grid(met~., scales="free_y") +
  geom_ribbon(aes(x=doy, ymin=lwr, ymax=upr2, fill=dataset), alpha=0.3) +
  geom_line(aes(x=doy, y=value, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Day of Year (Julian)") +
  scale_y_continuous(name="Daily Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
dev.off()
# ----------------
# -----------------------------------

# -----------------------------------
# 2. Define training window for all datasets
# -----------------------------------
vars.met <- c("tair", "tmax", "tmin", "precipf", "press", "qair", "wind", "swdown", "lwdown")

# Define the calibration window for each dataset
#  right now this is hard-coded, but should probably be soft-coded eventually
yrs.cal = data.frame(dataset=c("NLDAS", "CRUNCEP", "MIROC-ESM.hist", "MIROC-ESM.p1000"),
                     yr.min =c(   1980,      1980,             1901,              1829),
                     yr.max =c(   2010,      2010,             1921,              1849)
                     )
# -----------------------------------

# -----------------------------------
# 3. Doing the DOY bias-corrections : 
# -----------------------------------
#     - Met.train ~ s(doy) + Met.Raw
#     - I *think* this can be used to get monthly to daily data as well
#     ** Don't forget to use the residuals to add in the "error" in the model
#        -- assuming all goes well this will just be rnorm(mean=mean(resid), sd=sd(resid))
# Sourcing my gam processing scripts
source(file.path(path.gam, "Calculate_GAMM_Posteriors.R"))
library(MASS)
set.seed(321)
# Number of simulations  
n=10

# Creating a new data frame for the bias-corrected met means
met.bias <- met.day
summary(met.bias)

# Setting the variable we're working with for right now
var.now   = "tmax"

# -------------
# Creating the output list with LDAS as our base before we forget
# -------------
dat.train = "NLDAS"
yr.min    = yrs.cal[yrs.cal$dataset==dat.train, "yr.min"]
yr.max    = yrs.cal[yrs.cal$dataset==dat.train, "yr.max"]

# Calculating the anomalies from the current means
dat.temp <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min & met.bias$year<=yr.max, ]
dat.temp$Y <- dat.temp[,var.now]
lm.train <- lm(Y ~ as.factor(doy)-1, data=dat.temp)
train.anom <- met.bias[met.bias$dataset=="NLDAS", var.now] - predict(lm.train, newdata=met.bias[met.bias$dataset=="NLDAS", c("year", "doy")])
summary(train.anom)

train.ci <- data.frame(dataset="NLDAS", met=var.now, met.bias[met.bias$dataset=="NLDAS", c("year", "doy")], 
                       X=met.bias[met.bias$dataset=="NLDAS", var.now], 
                       anomaly=train.anom,
                       mean=met.bias[met.bias$dataset=="NLDAS", var.now], 
                       lwr=met.bias[met.bias$dataset=="NLDAS", var.now], 
                       upr=met.bias[met.bias$dataset=="NLDAS", var.now],
                       time=as.Date(met.bias[met.bias$dataset=="NLDAS", "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset=="NLDAS", "year"], "-01-01")))
                       )
train.sims <- data.frame(dataset="NLDAS", met=var.now, met.bias[met.bias$dataset=="NLDAS", c("year", "doy")], 
                         time=as.Date(met.bias[met.bias$dataset=="NLDAS", "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset=="NLDAS", "year"], "-01-01")))
                         )
train.sims[,paste0("X", 1:n)] <- met.bias[met.bias$dataset=="NLDAS", var.now]

dat.out <- list()
dat.out$ci <- train.ci
dat.out$sims <- train.sims
# -------------

# --------------------
# Working with CRUNCEP & tair first
# --------------------
# Define the training dataset
dat.train = "NLDAS"
dat.bias  = "CRUNCEP"
yr.min    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.min"]
yr.max    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.max"]

# establish a generlized data frames that can be used with any variables
# The calibration data
dat.temp <- data.frame(year = met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "year" ],
                       doy  = met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"  ],
                       Y    = met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min & met.bias$year<=yr.max, var.now],
                       X    = met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, var.now]
                       )
# The prediction Data
dat.pred <- data.frame(year = met.bias[met.bias$dataset==dat.bias, "year" ],
                       doy  = met.bias[met.bias$dataset==dat.bias, "doy"  ],
                       X    = met.bias[met.bias$dataset==dat.bias, var.now]
                       )

# Calculating the anomalies from the current means
lm.anom <- lm(X ~ as.factor(doy)-1, data=dat.temp)
dat.pred$anomaly <- dat.pred$X - predict(lm.anom, newdata=dat.pred)
summary(dat.pred)

# Plotting the correlation between dat.train & dat.bias
plot(Y ~ X, data=dat.temp) 
abline(a=0, b=1, col="red")
abline(lm(Y ~ X, data=dat.temp), col="red", lty="dashed")

ggplot(data=dat.temp) +
  geom_point(aes(x=X, y=Y, color=doy), size=0.5)

# Running a model
mod.bias <- gam(Y ~ s(doy) + X, data=dat.temp)
summary(mod.bias)
plot(mod.bias)

# Saving the mean predicted & residuals
dat.temp$pred  <- predict(mod.bias)
dat.temp$resid <- resid(mod.bias)
summary(dat.temp)

# Checkign the residuals to see if we can assume normality
plot(resid ~ pred, data=dat.temp); abline(h=0, col="red")
plot(resid ~ doy, data=dat.temp); abline(h=0, col="red")
hist(dat.temp$resid)

# Modeling the anomalies with a similar framework
# Calculating the anomalies from the current means
lm.anom <- lm(X ~ as.factor(doy)-1, data=dat.temp)
dat.pred$anomaly <- dat.pred$X - predict(lm.anom, newdata=dat.pred)
summary(dat.pred)

dat.pred$pred <- predict(mod.bias, newdata=dat.pred)
summary(dat.pred)

mod.anom <- gam(anomaly ~ s(doy) + pred, data=dat.pred)
summary(mod.anom)
plot(mod.anom)

dat.pred$resid2 <- resid(mod.anom)
plot(resid2 ~ pred, data=dat.pred); abline(h=0, col="red")
plot(resid2 ~ doy, data=dat.pred); abline(h=0, col="red")
hist(dat.pred$resid)

# --------
# Predicting a bunch of potential posteriors over the full dataset
# --------
# Get the model coefficients
coef.gam <- coef(mod.bias)
coef.anom <- coef(mod.anom)

# Generate a random distribution of betas using the covariance matrix
Rbeta <- mvrnorm(n=n, coef(mod.bias), vcov(mod.bias))
Rbeta.anom <- mvrnorm(n=n, coef(mod.anom), vcov(mod.anom))

# Create the prediction matrix
Xp <- predict(mod.bias, newdata=dat.pred, type="lpmatrix")
Xp.anom <- predict(mod.anom, newdata=dat.pred, type="lpmatrix")

# -----
# Simulate predicted met variables
# -----
# Adding in some residual error 
# ** QUESTION -- do I add a constant error per time series or random error to each day?
# -----
sim1 <- Xp %*% t(Rbeta) + Xp.anom %*% t(Rbeta.anom)

# # Option 1: Adding a constant error per time series
# res <- rnorm(n, mean(dat.temp$resid), sd(dat.temp$resid))
# sim1 <- sweep(sim1, 2, res, FUN="+")

# # Option 2: Adding a random error to each observation
# sim1 <- sim1 + array(rnorm(length(sim1), mean(dat.temp$resid), sd(dat.temp$resid)), dim=dim(sim1))

# # Option 3: Adding a day-specific error (although this should be part of the smoother)
# -----

# NOTE: Don't need to add in anomalies because training & fitting data are both observed annual deviations!
dat.pred$mean <- apply(sim1, 1, mean)
dat.pred$lwr  <- apply(sim1, 1, quantile, 0.025)
dat.pred$upr  <- apply(sim1, 1, quantile, 0.975)

# Doing som exploratory graphing
dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))

# Plotting the observed and the bias-corrected 95% CI
ggplot(data=dat.pred[dat.pred$year>=2005 & dat.pred$year<=2006,]) +
  geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=time, y=mean), color="red", size=0.5) +
  geom_line(aes(x=time, y=X), color='black', size=0.5) +
  theme_bw()

dat.yr <- aggregate(dat.pred[,c("X", "mean", "lwr", "upr")],
                    by=list(dat.pred$year),
                    FUN=mean)
names(dat.yr)[1] <- "year"

ggplot(data=dat.yr[,]) +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=year, y=mean), color="red", size=0.5) +
  geom_line(aes(x=year, y=X), color='black', size=0.5) +
  theme_bw()
# --------

# --------
# Storing the output
# --------
dat.sims <- data.frame(dataset=dat.bias, met=var.now, dat.pred[,c("year", "doy", "time")])
dat.sims <- cbind(dat.sims, sim1)

dat.out$ci   <- rbind(dat.out$ci, data.frame(dataset=dat.bias, met=var.now, dat.pred[,names(dat.out$ci)[3:ncol(dat.out$ci)]]))
dat.out$sims <- rbind(dat.out$sims, data.frame(dat.sims))
# --------
# --------------------


# --------------------
# Working with GCM historical
# Trick: Propagating uncertainty from first correction!
# --------------------
# Define the training dataset
dat.train = "CRUNCEP"
dat.bias  = "MIROC-ESM.hist"
var.now   = "tmax"
yr.min    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.min"]
yr.max    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.max"]

# establish a dummy data frame
dat.temp0 <- data.frame(year = dat.out$sims[dat.out$sims$dataset==dat.train  & dat.out$sims$year>=yr.min & dat.out$sims$year<=yr.max, "year" ],
                        doy  = dat.out$sims[dat.out$sims$dataset==dat.train  & dat.out$sims$year>=yr.min & dat.out$sims$year<=yr.max, "doy"  ],
                        stack(dat.out$sims[dat.out$sims$dataset==dat.train & dat.out$sims$year>=yr.min & dat.out$sims$year<=yr.max, paste0("X", 1:n)])
                        )

dat.temp <- aggregate(dat.temp0$values, by=dat.temp0[,c("doy", "ind")], FUN=mean)
names(dat.temp)[3] <- "Y"
summary(dat.temp)

dat.temp2 <- aggregate(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, var.now],
                       by=list(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"]),
                       FUN=mean)
names(dat.temp2) <- c("doy", "X")
summary(dat.temp2)

# merging together the two sets of daily means 
# Note: working off of the climatological means rather than trying to pair years gets us borader CIs
dat.temp <- merge(dat.temp[,], dat.temp2)
summary(dat.temp)


# The prediction Data
dat.pred <- data.frame(year = met.bias[met.bias$dataset==dat.bias, "year" ],
                       doy  = met.bias[met.bias$dataset==dat.bias, "doy"  ],
                       X    = met.bias[met.bias$dataset==dat.bias, var.now]
                       )


# Plotting the correlation between dat.train & dat.bias
plot(Y ~ X, data=dat.temp) 
abline(a=0, b=1, col="red")
abline(lm(Y ~ X, data=dat.temp), col="red", lty="dashed")

# Running a model
mod.bias <- gam(Y ~ s(doy) + X, data=dat.temp)
summary(mod.bias)
plot(mod.bias)

# Saving the mean predicted & residuals
dat.temp$pred  <- predict(mod.bias)
dat.temp$resid <- resid(mod.bias)
summary(dat.temp)

# Checkign the residuals to see if we can assume normality
plot(resid ~ pred, data=dat.temp); abline(h=0, col="red")
plot(resid ~ doy, data=dat.temp); abline(h=0, col="red")
hist(dat.temp$resid)

# Modeling the anomalies with a similar framework
# Calculating the anomalies from the current means
lm.anom <- lm(X ~ as.factor(doy)-1, data=dat.temp)
dat.pred$anomaly <- dat.pred$X - predict(lm.anom, newdata=dat.pred)
summary(dat.pred)

dat.pred$pred <- predict(mod.bias, newdata=dat.pred)
summary(dat.pred)

mod.anom <- gam(anomaly ~ s(doy) + pred, data=dat.pred)
summary(mod.anom)
plot(mod.anom)

dat.pred$resid2 <- resid(mod.anom)
plot(resid2 ~ pred, data=dat.pred); abline(h=0, col="red")
plot(resid2 ~ doy, data=dat.pred); abline(h=0, col="red")
hist(dat.pred$resid)
# --------
# Predicting a bunch of potential posteriors over the full dataset
# --------
# Get the model coefficients
coef.gam <- coef(mod.bias)
coef.anom <- coef(mod.anom)

# Generate a random distribution of betas using the covariance matrix
Rbeta <- mvrnorm(n=n, coef(mod.bias), vcov(mod.bias))
Rbeta.anom <- mvrnorm(n=n, coef(mod.anom), vcov(mod.anom))

# Create the prediction matrix
Xp <- predict(mod.bias, newdata=dat.pred, type="lpmatrix")
Xp.anom <- predict(mod.anom, newdata=dat.pred, type="lpmatrix")

# -----
# Simulate predicted met variables
# -----
# Adding in some residual error 
# ** QUESTION -- do I add a constant error per time series or random error to each day?
# -----
sim1 <- Xp %*% t(Rbeta) + Xp.anom %*% t(Rbeta.anom)

# # Option 1: Adding a constant error per time series
# res <- rnorm(n, mean(dat.temp$resid), sd(dat.temp$resid))
# sim1 <- sweep(sim1, 2, res, FUN="+")

# # Option 2: Adding a random error to each observation
# sim1 <- sim1 + array(rnorm(length(sim1), mean(dat.temp$resid), sd(dat.temp$resid)), dim=dim(sim1))

# # Option 3: Adding a day-specific error (although this should be part of the smoother)
# -----

dat.pred$mean <- apply(sim1, 1, mean)
dat.pred$lwr  <- apply(sim1, 1, quantile, 0.025)
dat.pred$upr  <- apply(sim1, 1, quantile, 0.975)

# Doing some exploratory graphing
dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))

# Plotting the observed and the bias-corrected 95% CI
ggplot(data=dat.pred[dat.pred$year>=1854 & dat.pred$year<=1856,]) +
  geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=time, y=mean), color="red", size=0.5) +
  geom_line(aes(x=time, y=X), color='black', size=0.5) +
  theme_bw()

dat.yr <- aggregate(dat.pred[,c("X", "mean", "lwr", "upr")],
                    by=list(dat.pred$year),
                    FUN=mean)
names(dat.yr)[1] <- "year"

ggplot(data=dat.yr[,]) +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=year, y=mean), color="red", size=0.5) +
  geom_line(aes(x=year, y=X), color='black', size=0.5) +
  theme_bw()
# --------

# --------
# Storing the output
# --------
dat.sims <- data.frame(dataset=dat.bias, met=var.now, dat.pred[,c("year", "doy", "time")])
dat.sims <- cbind(dat.sims, sim1)

dat.out$ci   <- rbind(dat.out$ci, data.frame(dataset=dat.bias, met=var.now, dat.pred[,names(dat.out$ci)[3:ncol(dat.out$ci)]]))
dat.out$sims <- rbind(dat.out$sims, data.frame(dat.sims))
# --------
# --------------------


# --------------------
# Working with GCM p1000
# Trick: Propagating uncertainty from first correction!
# NOTE: Need to match climatic means for to different periods
# --------------------
# Define the training dataset
dat.train = "MIROC-ESM.hist"
dat.bias  = "MIROC-ESM.p1000"
var.now   = "tmax"
yr.min    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.min"]
yr.max    = yrs.cal[yrs.cal$dataset==dat.bias, "yr.max"]

# Need to specify a separate range for the training data -- use the same number of years!
yr.min.train = min(dat.out$sims[dat.out$sims$dataset==dat.train, "year"])
yr.max.train = yr.min.train + (yr.max - yr.min)

# establish a dummy data frame
dat.temp0 <- data.frame(year = dat.out$sims[dat.out$sims$dataset==dat.train  & dat.out$sims$year>=yr.min.train & dat.out$sims$year<=yr.max.train, "year" ],
                        doy  = dat.out$sims[dat.out$sims$dataset==dat.train  & dat.out$sims$year>=yr.min.train & dat.out$sims$year<=yr.max.train, "doy"  ],
                        stack(dat.out$sims[dat.out$sims$dataset==dat.train & dat.out$sims$year>=yr.min.train & dat.out$sims$year<=yr.max.train, paste0("X", 1:n)])
                        )

dat.temp <- aggregate(dat.temp0$values, by=dat.temp0[,c("doy", "ind")], FUN=mean)
names(dat.temp)[3] <- "Y"
summary(dat.temp)

dat.temp2 <- aggregate(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, var.now],
                       by=list(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"]),
                       FUN=mean)
names(dat.temp2) <- c("doy", "X")
summary(dat.temp2)

# merging together the two sets of daily means 
# Note: working off of the climatological means rather than trying to pair years gets us borader CIs
dat.temp <- merge(dat.temp[,], dat.temp2)
summary(dat.temp)

# The prediction Data
dat.pred <- data.frame(year = met.bias[met.bias$dataset==dat.bias, "year" ],
                       doy  = met.bias[met.bias$dataset==dat.bias, "doy"  ],
                       X    = met.bias[met.bias$dataset==dat.bias, var.now]
                       )


# Plotting the correlation between dat.train & dat.bias
plot(Y ~ X, data=dat.temp) 
abline(a=0, b=1, col="red")
abline(lm(Y ~ X, data=dat.temp), col="red", lty="dashed")

# Running a model
mod.bias <- gam(Y ~ s(doy) + X, data=dat.temp)
summary(mod.bias)
plot(mod.bias)

# Saving the mean predicted & residuals
dat.temp$pred  <- predict(mod.bias)
dat.temp$resid <- resid(mod.bias)
summary(dat.temp)

# Checkign the residuals to see if we can assume normality
plot(resid ~ pred, data=dat.temp); abline(h=0, col="red")
plot(resid ~ doy, data=dat.temp); abline(h=0, col="red")
hist(dat.temp$resid)

# Modeling the anomalies with a similar framework
# Calculating the anomalies from the current means
lm.anom <- lm(X ~ as.factor(doy)-1, data=dat.temp)
dat.pred$anomaly <- dat.pred$X - predict(lm.anom, newdata=dat.pred)
summary(dat.pred)

dat.pred$pred <- predict(mod.bias, newdata=dat.pred)
summary(dat.pred)

mod.anom <- gam(anomaly ~ s(doy) + pred, data=dat.pred)
summary(mod.anom)
plot(mod.anom)

dat.pred$resid2 <- resid(mod.anom)
plot(resid2 ~ pred, data=dat.pred); abline(h=0, col="red")
plot(resid2 ~ doy, data=dat.pred); abline(h=0, col="red")
hist(dat.pred$resid)
# --------
# Predicting a bunch of potential posteriors over the full dataset
# --------
# Get the model coefficients
coef.gam <- coef(mod.bias)
coef.anom <- coef(mod.anom)

# Generate a random distribution of betas using the covariance matrix
Rbeta <- mvrnorm(n=n, coef(mod.bias), vcov(mod.bias))
Rbeta.anom <- mvrnorm(n=n, coef(mod.anom), vcov(mod.anom))

# Create the prediction matrix
Xp <- predict(mod.bias, newdata=dat.pred, type="lpmatrix")
Xp.anom <- predict(mod.anom, newdata=dat.pred, type="lpmatrix")

# -----
# Simulate predicted met variables
# -----
# Adding in some residual error 
# ** QUESTION -- do I add a constant error per time series or random error to each day?
# -----
sim1 <- Xp %*% t(Rbeta) + Xp.anom %*% t(Rbeta.anom)

# # Option 1: Adding a constant error per time series
# res <- rnorm(n, mean(dat.temp$resid), sd(dat.temp$resid))
# sim1 <- sweep(sim1, 2, res, FUN="+")

# # Option 2: Adding a random error to each observation
# sim1 <- sim1 + array(rnorm(length(sim1), mean(dat.temp$resid), sd(dat.temp$resid)), dim=dim(sim1))

# # Option 3: Adding a day-specific error (although this should be part of the smoother)
# -----

dat.pred$mean <- apply(sim1, 1, mean)
dat.pred$lwr  <- apply(sim1, 1, quantile, 0.025)
dat.pred$upr  <- apply(sim1, 1, quantile, 0.975)

# Doing some exploratory graphing
dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))

# Plotting the observed and the bias-corrected 95% CI
ggplot(data=dat.pred[dat.pred$year>=1844 & dat.pred$year<=1846,]) +
  geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=time, y=mean), color="red", size=0.5) +
  geom_line(aes(x=time, y=X), color='black', size=0.5) +
  theme_bw()

dat.yr <- aggregate(dat.pred[,c("X", "mean", "lwr", "upr")],
                    by=list(dat.pred$year),
                    FUN=mean)
names(dat.yr)[1] <- "year"

ggplot(data=dat.yr[,]) +
  geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  geom_line(aes(x=year, y=mean), color="red", size=0.5) +
  geom_line(aes(x=year, y=X), color='black', size=0.5) +
  theme_bw()
# --------

# --------
# Storing the output
# --------
dat.sims <- data.frame(dataset=dat.bias, met=var.now, dat.pred[,c("year", "doy", "time")])
dat.sims <- cbind(dat.sims, sim1)

dat.out$ci   <- rbind(dat.out$ci, data.frame(dataset=dat.bias, met=var.now, dat.pred[,names(dat.out$ci)[3:ncol(dat.out$ci)]]))
dat.out$sims <- rbind(dat.out$sims, data.frame(dat.sims))
# --------
# --------------------

# Looking at the new bias-corrected time series! 
dat.yr.bias <- aggregate(dat.out$ci[,c("X", "anomaly", "mean", "lwr", "upr")],
                         by=dat.out$ci[,c("dataset", "met", "year")],
                         FUN=mean)
dat.yr.bias$dataset <- factor(dat.yr.bias$dataset, levels=c("MIROC-ESM.p1000", "MIROC-ESM.hist", "CRUNCEP", "NLDAS"))
# dat.yr.bias$met     <- factor(dat.yr.bias$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
summary(dat.yr.bias)
summary(met.year[met.year$met==var.now,]) #raw annual data
# ry.noresid <- dat.yr.bias

ggplot(data=dat.yr.bias) +
  facet_grid(met~., scales="free_y") +
  #   geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=dataset)) +
  geom_line(aes(x=year, y=mean, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")

LDAS.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="NLDAS",]
CRU.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="CRUNCEP" & dat.yr.bias$year<min(LDAS.use$year),]
GCM.hist.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="MIROC-ESM.hist" & dat.yr.bias$year<min(CRU.use$year),]
GCM.p1k.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="MIROC-ESM.p1000" & dat.yr.bias$year<min(GCM.hist.use$year),]
met.final.yr <- rbind(LDAS.use.yr, CRU.use.yr, GCM.hist.use.yr, GCM.p1k.use.yr)

ggplot(data=met.final.yr) +
  facet_grid(met~., scales="free_y") +
  #   geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=dataset)) +
  geom_line(aes(x=year, y=mean, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")


LDAS.use <- dat.out$ci[dat.out$ci$dataset=="NLDAS",]
CRU.use <- dat.out$ci[dat.out$ci$dataset=="CRUNCEP" & dat.out$ci$year<min(LDAS.use$year),]
GCM.hist.use <- dat.out$ci[dat.out$ci$dataset=="MIROC-ESM.hist" & dat.out$ci$year<min(CRU.use$year),]
GCM.p1k.use <- dat.out$ci[dat.out$ci$dataset=="MIROC-ESM.p1000" & dat.out$ci$year<min(GCM.hist.use$year),]

met.final.day <- rbind(LDAS.use, CRU.use, GCM.hist.use, GCM.p1k.use)
times.spice <- data.frame(Splice=c("GCM.p1000-GCM.hist", "GCM.hist-CRUNCEP", "CRUNCEP-LDAS"),
                          yr.start=c(1847, 1898, 1977),
                          yr.end  =c(1852, 1903, 1982))

met.final.day$splice <- NA
met.final.day[met.final.day$year>=1847 & met.final.day$year<=1852,"splice"] <- "GCM.p1000-GCM.hist"
met.final.day[met.final.day$year>=1899 & met.final.day$year<=1904,"splice"] <- "GCM.hist-CRUNCEP"
met.final.day[met.final.day$year>=1977 & met.final.day$year<=1982,"splice"] <- "CRUNCEP-LDAS"
met.final.day$splice <- as.factor(met.final.day$splice)
met.final.day$year.frac <- met.final.day$year + met.final.day$doy/366
summary(met.final.day)

ggplot(data=met.final.day[!is.na(met.final.day$splice),]) +
  facet_wrap(~splice, scales="free_x", ncol=1) +
  #   geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=dataset)) +
  geom_line(aes(x=year.frac, y=mean, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")


# Original 
summary(met.year)
ggplot(data=met.year[met.year$met==var.now,]) +
  facet_grid(met~., scales="free_y") +
  geom_line(aes(x=year, y=value, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
# -----------------------------------

# -----------------------------------
# 4. Generate an ensemble of the predicted mean day of year
# -----------------------------------
#     - use covariance matrix and make array of n ensemble members & run raw data
# -----------------------------------

# -----------------------------------
# 5. Add in anomaly
# -----------------------------------
#    - use anomaly from original climate mean rather than residual from predicted model
# -----------------------------------

# -----------------------------------
# 6. Save an .Rdata file with an ensemble of bias-corrected daily data
# -----------------------------------
#    - this ensemble file can either be read into the temporal downscaling or fed into
#      a script that writes to netcdf format for models that take daily drivers
# -----------------------------------
