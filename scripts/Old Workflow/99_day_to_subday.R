# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Temporal downscaling of daily met data to hourly or 3-hourly data
# Creator: Christy Rollinson, 1 July 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Take the daily, bias-corrected met files that come out of step 2 and temporally 
# downscale to subdaily (hourly or 3-hourly) using the raw LDAS data.
#
# Default output will be hourly data (NLDAS timestep), but that can be changed to 
# 3-hourly as needed
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# Key things to keep in mind:
#  1. Need to use a filter or rollapply-like function so that we smooth over any 
#     gaps between days
#  2. think about setting up as a function with set seed and ensemble member numbers 
#     so that we can either exactly replicate our work and/or make this a function 
#     that can get shipped to participants so that can make the files themselves rather 
#     than download large files

# 0. Set up file structure, etc.
# 1. Read in & format data sets
# -----------------------------------

rm(list=ls())

# -----------------------------------
# 0. Set up file structure, load packages, etc
# -----------------------------------
# Load libraries
library(mgcv); library(ggplot2)

# Set the working directory
# wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble"
wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

# Defining a site name -- this can go into a function later
site.name="HARVARD"

path.dat <- file.path(wd.base, "data/met_ensembles", site.name)
path.out <- file.path(wd.base, "data/met_ensembles", site.name, "subdaily")
if(!dir.exists(path.out)) dir.create(path.out, recursive=T)  

met.done <- dir(path.out, ".csv")
# -----------------------------------


# -----------------------------------
# 1. Read in & format data sets
# -----------------------------------
# Just reading in the LDAS training data and doing a couple quick 
# exploratory figures for now
ldas <- read.csv(file.path("data/paleon_sites", site.name, "NLDAS_1980-2015.csv"))

# ******* TYPO CORRECTIONS!! *******
# **** NOTE: WILL NEED TO CANCEL THIS OUT UPON NEW EXTRACTION
# make first day of year start on 0
ldas$doy    <- ldas$doy - 1
ldas$precipf <- ldas$precipf*.1

summary(ldas)

# Aggregating LDAS to look at the hourly distribution by day of year
vars.met <- c("tair", "precipf", "press", "qair", "wind", "swdown", "lwdown")

ldas.doy0 <- aggregate(ldas[,vars.met], by=ldas[,c("dataset", "doy", "hour")], FUN=mean)
ldas.doy <- stack(ldas.doy0[,vars.met])
names(ldas.doy) <- c("value", "met")
ldas.doy[,c("dataset", "doy", "hour")] <- ldas.doy0[,c("dataset", "doy", "hour")]
ldas.doy$lwr <- stack(aggregate(ldas[,vars.met], by=ldas[,c("dataset", "doy", "hour")], FUN=quantile, 0.025)[,vars.met])[,1]
ldas.doy$upr <- stack(aggregate(ldas[,vars.met], by=ldas[,c("dataset", "doy", "hour")], FUN=quantile, 0.975)[,vars.met])[,1]
summary(ldas.doy)

# Doing some factor ordering
ldas.doy$met <- factor(ldas.doy$met,  levels=c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind"))

# Putting in a dummy upper bound for precip to deal with some of the real oddballs
precip.cutoff <- quantile(ldas.doy[ldas.doy$met=="precipf","upr"], 0.95)
ldas.doy$upr2 <- ifelse(ldas.doy$met=="precipf" & ldas.doy$upr>=precip.cutoff, precip.cutoff, ldas.doy$upr)


days.plot <- round(seq((365/4)/2, 365, by=365/4),0)
1 <- as.Date(days.plot, origin="2000-01-01")

png(file.path(path.out, "NLDAS_SubDaily_Cycle_Means.png"), height=11, width=8.5, "in", res=180)
ggplot(data=ldas.doy[ldas.doy$doy %in% days.plot,]) +
  facet_grid(met~., scales="free_y") +
  geom_ribbon(aes(x=hour, ymin=lwr, ymax=upr2, fill=factor(doy)), alpha=0.3) +
  geom_line(aes(x=hour, y=value, color=factor(doy))) +
  scale_x_continuous(expand=c(0,0), name="Hour of Day") +
  scale_y_continuous(name="Hourly Mean") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
dev.off()

# Stacking ldas to make it easier to graph
ldas.stack <- stack(ldas[,vars.met])
names(ldas.stack) <- c("value", "met")
ldas.stack[,c("dataset", "year", "doy", "hour")] <- ldas[,c("dataset", "year", "doy", "hour")]
ldas.stack$met <- factor(ldas.stack$met,  levels=c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind"))
summary(ldas.stack)

png(file.path(path.out, "NLDAS_SubDaily_Cycle_Examples.png"), height=11, width=8.5, "in", res=180)
ggplot(data=ldas.stack[ldas.stack$doy %in% days.plot & ldas.stack$year %in% c(1990, 1995, 2000, 2005, 2010),]) +
  facet_grid(met~doy, scales="free_y") +
#   geom_ribbon(aes(x=hour, ymin=lwr, ymax=upr2, fill=factor(doy)), alpha=0.3) +
  geom_line(aes(x=hour, y=value, color=factor(year))) +
  scale_x_continuous(expand=c(0,0), name="Hour of Day") +
  scale_y_continuous(name="Hourly Mean") +
  theme_bw() +
  theme(legend.position="top",
        legend.direction="horizontal")
dev.off()

# -----------------------------------

