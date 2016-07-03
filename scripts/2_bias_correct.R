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
source("scripts/bias_correct_day.R")

# The met vars we need (in the order we want to do them)
vars.met <- c("tmax", "tmin", "precipf", "swdown", "qair", "lwdown", "press", "wind")
dat.cal = "NLDAS"

# Note This dataframe is only for the datasets to be bias-corrected!
yrs.cal = data.frame(dataset = c("CRUNCEP", "MIROC-ESM.hist", "MIROC-ESM.p1000"),
                     cal.min = c(1980,             1901,              1829),
                     cal.max = c(2010,             1921,              1849)
                     )

# -----------------------------------

# -----------------------------------
# 3. Doing the DOY bias-corrections : 
# -----------------------------------
# bias.correct <- function(met.day, met.var, dat.train, yrs.cal, n){
# yrs.cal = yrs.cal[1,]

# Make precipf in total kg/m2/day rather than per s
met.day$precipf <- met.day$precipf*60*60*24
summary(met.day)

source("scripts/bias_correct_day.R")
met.bias <- met.day
dat.out.full <- bias.correct(met.bias=met.bias, vars.met=vars.met, dat.train="NLDAS", yrs.cal=yrs.cal, n=100)
# --------------------


# --------------------
# Looking at the new bias-corrected time series! 
# --------------------
for(met.var in vars.met){
  print(met.var)
  dat.out <- dat.out.full[[met.var]]
  dat.yr.bias <- aggregate(dat.out$ci[,c("X", "anom.raw", "mean", "lwr", "upr")],
                           by=dat.out$ci[,c("dataset", "met", "year")],
                           FUN=mean)
  dat.yr.bias$dataset <- factor(dat.yr.bias$dataset, levels=c("MIROC-ESM.p1000", "MIROC-ESM.hist", "CRUNCEP", "NLDAS"))
  # dat.yr.bias$met     <- factor(dat.yr.bias$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
  summary(dat.yr.bias)
  summary(met.year[met.year$met==met.var,]) #raw annual data
  
  LDAS.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="NLDAS",]
  CRU.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="CRUNCEP" & dat.yr.bias$year<min(LDAS.use.yr$year),]
  GCM.hist.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="MIROC-ESM.hist" & dat.yr.bias$year<min(CRU.use.yr$year),]
  GCM.p1k.use.yr <- dat.yr.bias[dat.yr.bias$dataset=="MIROC-ESM.p1000" & dat.yr.bias$year<min(GCM.hist.use.yr$year),]
  met.final.yr <- rbind(LDAS.use.yr, CRU.use.yr, GCM.hist.use.yr, GCM.p1k.use.yr)
  
  # Original 
  summary(met.year)
  png(file.path(path.out, paste0("Met_", met.var, "_original_year.png")), height=8.5, width=11, "in", res=180)
  print(
  ggplot(data=met.final.yr[met.final.yr$met==met.var,]) +
    facet_grid(met~., scales="free_y") +
    geom_line(aes(x=year, y=X, color=dataset)) +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(name="Annual Mean Value") +
    theme_bw() +
    theme(legend.position="top",
          legend.direction="horizontal")
  )
  dev.off()
  
  # Final bias-corrected time series
  png(file.path(path.out, paste0("Met_", met.var, "_bias-corrected_year.png")), height=8.5, width=11, "in", res=180)
  print(
  ggplot(data=met.final.yr[met.final.yr$met==met.var,]) +
    facet_grid(met~., scales="free_y") +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr, fill=dataset), alpha=0.3) +
    geom_line(aes(x=year, y=mean, color=dataset)) +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(name="Annual Mean Value") +
    theme_bw() +
    theme(legend.position="top",
          legend.direction="horizontal")
  )
  dev.off()
  
  LDAS.use <- dat.out$ci[dat.out$ci$dataset=="NLDAS",]
  CRU.use <- dat.out$ci[dat.out$ci$dataset=="CRUNCEP" & dat.out$ci$year<min(LDAS.use$year),]
  GCM.hist.use <- dat.out$ci[dat.out$ci$dataset=="MIROC-ESM.hist" & dat.out$ci$year<min(CRU.use$year),]
  GCM.p1k.use <- dat.out$ci[dat.out$ci$dataset=="MIROC-ESM.p1000" & dat.out$ci$year<min(GCM.hist.use$year),]
  
  met.final.day <- rbind(LDAS.use, CRU.use, GCM.hist.use, GCM.p1k.use)
  met.final.day$splice <- NA
  met.final.day[met.final.day$year>=1847 & met.final.day$year<=1852,"splice"] <- "GCM.p1000-GCM.hist"
  met.final.day[met.final.day$year>=1898 & met.final.day$year<=1903,"splice"] <- "GCM.hist-CRUNCEP"
  met.final.day[met.final.day$year>=1977 & met.final.day$year<=1982,"splice"] <- "CRUNCEP-LDAS"
  met.final.day$splice <- as.factor(met.final.day$splice)
  met.final.day$year.frac <- met.final.day$year + met.final.day$doy/366
  summary(met.final.day)
  met.final.day$dataset <- factor(met.final.day$dataset, levels=c("MIROC-ESM.p1000", "MIROC-ESM.hist", "CRUNCEP", "NLDAS"))
  
  png(file.path(path.out, paste0("Met_", met.var, "_bias-corrected_splices.png")), height=8.5, width=11, "in", res=180)
  print(
  ggplot(data=met.final.day[!is.na(met.final.day$splice),]) +
    facet_wrap(~splice, scales="free_x", ncol=1) +
    geom_ribbon(aes(x=year.frac, ymin=lwr, ymax=upr, fill=dataset), alpha=0.3) +
    geom_line(aes(x=year.frac, y=mean, color=dataset)) +
    scale_x_continuous(expand=c(0,0), name="Year") +
    scale_y_continuous(expand=c(0,0), name="Annual Mean Value") +
    theme_bw() +
    theme(legend.position="top",
          legend.direction="horizontal"))
  dev.off()
}
# -----------------------------------

save(dat.out.full, file=file.path(path.out, "MetAll_Daily_Ensemble_MIROC-ESM.Rdata"))
