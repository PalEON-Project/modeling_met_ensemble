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
# wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble"
# wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

# Defining a site name -- this can go into a function later
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18
GCM.list=c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
# GCM.list=c("CCSM4")
LDAS="NLDAS"
n=25 # Number of ensemble members
# -----------------------------------

for(GCM in GCM.list){
  
print(GCM)
path.dat <- file.path(wd.base, "data/paleon_sites", site.name)
# path.out <- file.path(wd.base, "data/met_ensembles", site.name, GCM, "day")
path.out <- file.path("~/Desktop/met_bias_day", site.name, GCM)
if(!dir.exists(path.out)) dir.create(path.out, recursive=T)  

met.done <- dir(path.out, ".csv")


# -----------------------------------
# 1. Read in & format the different datasets
# -----------------------------------
# Find the appropriate file name for each
file.ldas <- dir(path.dat, LDAS)
file.cru <- dir(path.dat, "CRUNCEP")
file.hist <- dir(path.dat, paste0(GCM, "_historical"))
file.p1k <- dir(path.dat, paste0(GCM, "_p1000"))

ldas     <- read.csv(file.path(path.dat, file.ldas))
cruncep  <- read.csv(file.path(path.dat, file.cru))
gcm.hist <- read.csv(file.path(path.dat, file.hist))
gcm.p1k  <- read.csv(file.path(path.dat, file.p1k))

# Adding an hour field to the gcm; setting as noon (middle of window) for simplicity
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
# ldas$doy    <- ldas$doy - 1
# cruncep$doy <- cruncep$doy -1
# ldas$precipf <- ldas$precipf*.1

summary(ldas)
summary(cruncep)
summary(gcm.p1k)
summary(gcm.hist)



# vars.met <- c("tair", "tmax", "tmin", "precipf", "press", "qair", "wind", "swdown", "lwdown")
vars.met <- c("tair", "tmax", "tmin", "qair", "precipf", "swdown", "press", "lwdown", "wind")
cols.bind <- c("dataset", "year", "doy", "hour", vars.met)
met.all <- rbind(ldas[,cols.bind], cruncep[,cols.bind], gcm.p1k[,cols.bind], gcm.hist[,cols.bind])
met.all$Date <- as.Date(met.all$doy, origin=as.Date(paste(met.all$year, "01", "01", sep="-")))
summary(met.all)
# ----------------
# Creating a daily met record
# ----------------
met.day <- aggregate(met.all[,vars.met], by=met.all[,c("dataset", "year", "doy")], FUN=mean)

# getting tmax & tmin for ldas & cru
met.day2 <- aggregate(met.all[met.all$dataset %in% c(LDAS, "CRUNCEP"),"tair"], 
                      by=met.all[met.all$dataset %in% c(LDAS, "CRUNCEP"),c("dataset", "year", "doy")], 
                      FUN=max)
names(met.day2)[4] <- "tmax"
met.day2$tmin <- aggregate(met.all[met.all$dataset %in% c(LDAS, "CRUNCEP"),"tair"], 
                           by=met.all[met.all$dataset %in% c(LDAS, "CRUNCEP"),c("dataset", "year", "doy")], 
                           FUN=min)[,4]
summary(met.day2)

# merging tmax & tmin back into met.day 
met.day[met.day$dataset %in% c(LDAS, "CRUNCEP"), "tmax"] <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                                                                         by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
                                                                         FUN=max)[,4]
met.day[met.day$dataset %in% c(LDAS, "CRUNCEP"), "tmin"] <- aggregate(met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),"tair"], 
                                                                         by=met.all[met.all$dataset %in% c("NLDAS", "CRUNCEP"),c("dataset", "year", "doy")], 
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

met.year$dataset <- factor(met.year$dataset, levels=c(paste0(GCM, ".p1000"), paste0(GCM, ".hist"), "CRUNCEP", LDAS))
met.doy $dataset <- factor(met.doy $dataset, levels=c(paste0(GCM, ".p1000"), paste0(GCM, ".hist"), "CRUNCEP", LDAS))
met.year$met     <- factor(met.year$met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
met.doy $met     <- factor(met.doy $met,  levels=c("tair", "tmax", "tmin", "precipf", "swdown", "lwdown", "press", "qair", "wind"))

# Putting in a dummy upper bound for precip to deal with some of the real oddballs
precip.cutoff <- quantile(met.doy[met.doy$met=="precipf","upr"], 0.95)
met.doy$upr2 <- ifelse(met.doy$met=="precipf" & met.doy$upr>=precip.cutoff, precip.cutoff, met.doy$upr)

png(file.path(path.out, paste0(GCM, "_Raw_Year_0850-2015.png")), height=11, width=8.5, "in", res=180)
ggplot(data=met.year) +
  facet_grid(met~., scales="free_y") +
  geom_line(aes(x=year, y=value, color=dataset)) +
  scale_x_continuous(expand=c(0,0), name="Year") +
  scale_y_continuous(name="Annual Mean Value") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
dev.off()

png(file.path(path.out, paste0(GCM, "_Raw_DOY_All.png")), height=11, width=8.5, "in", res=180)
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
# 3. Doing the DOY bias-corrections : 
# -----------------------------------
# bias.correct <- function(met.day, met.var, dat.train, yrs.cal, n){
# yrs.cal = yrs.cal[1,]

# The met vars we need (in the order we want to do them)
# vars.met <- c("tmax", "tmin", "swdown", "lwdown", "precipf", "qair", "press", "wind")
vars.met <- c("tmax", "tmin", "qair", "precipf", "swdown", "press", "lwdown", "wind")
dat.cal = LDAS

# Note This dataframe is only for the datasets to be bias-corrected!
yrs.cal = data.frame(dataset = c("CRUNCEP", paste0(GCM, ".hist"), paste0(GCM, ".p1000")),
                     cal.min = c(1980,             1901,              1829),
                     cal.max = c(2010,             1921,              1849)
                     )


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

source("scripts/bias_correct_day.R")
met.bias <- met.day
dat.out.full <- bias.correct(met.bias=met.bias, vars.met=vars.met, dat.train=LDAS, GCM=GCM, yrs.cal=yrs.cal, n=n, path.out=path.out)
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
  dat.yr.bias$dataset <- factor(dat.yr.bias$dataset, levels=c(paste0(GCM, ".p1000"), paste0(GCM, ".hist"), "CRUNCEP", LDAS))
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

  png(file.path(path.out, paste0(GCM,"_", met.var, "_bias-correction.png")), height=8.5, width=14, "in", res=180)
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

png(file.path(path.out, paste0(GCM,"_Ensembles_Smoothed.png")), height=8.5, width=14, "in", res=180)
print(
ggplot(data=dat.final) +
  facet_wrap(~Met, scales="free_y") +
  geom_line(aes(x=year, y=annual, color=EnsMem), size=0.2, alpha=0.2) +
  geom_line(aes(x=year, y=decadal, color=EnsMem)) +
  scale_x_continuous(name="Year (A.D.)", expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  guides(color=F) +
  geom_vline(xintercept=c(1850, 1901, 1980), linetype="dashed") +
  theme_bw()
)
dev.off()
# --------------

save(dat.out.full, file=file.path(path.out, paste0(GCM, "_day_alldata.Rdata")))
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


yr.bins <- c(min(dat.out.full$met.bias$year), seq(min(dat.out.full$met.bias$year)+50, round(max(dat.out.full$met.bias$year),-2), by=100))
for(i in 1:n){
  # Make a directory for each ensemble member
  out.name <- paste0(site.name, "_", GCM, "_day_", str_pad(i, 3, pad=0))
  new.dir <- file.path(wd.base, "data/met_ensembles", site.name, GCM, "day", out.name)
  if(!dir.exists(new.dir)) dir.create(new.dir, recursive=T)  
  
  for(j in 1:length(yr.bins)){
    date.start <- as.Date(paste0(yr.bins[j], "-01-01"))
    if(j < length(yr.bins)){
      date.end   <- as.Date(paste0(yr.bins[j+1]-1, "-12-31"))
    } else {
      date.end   <- as.Date(paste0(max(dat.out.full$met.bias$year), "-12-31"))
    }
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
    
    # Naming convention: [SITE]_[GCM]_day_[member]_[YEAR].nc
    nc <- nc_create(file.path(new.dir, paste0(out.name, "_", str_pad(year(date.start), 4, pad=0), ".nc")), var.list)
    for(v in 1:length(var.list)) {
      ncvar_put(nc, var.list[[v]], dat.list[[v]])
    }
    nc_close(nc)    
  }
  # Compress the ensemble file  setwd(dir.out) # Go to the output directory so we don't get annoying file paths
  setwd(path.out)
  system(paste0("tar -jcvf ", out.name, ".tar.bz2 ", out.name)) # Compress the folder
  system(paste0("rm -rf ", out.name)) # remove the uncompressed folder
  setwd(wd.base) # Go back to our base directory
  

} # end splitting ensemble members
# -----------------------------------

} # End GCM loop
