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
library(ncdf4)
library(mgcv); library(ggplot2)
library(stringr)
library(lubridate)

# Set the working directory
wd.base <- "~/Dropbox/PalEON_CR/met_ensemble/"
out.base <- "~/Desktop/met_ensembles/"
setwd(wd.base)

# Setting some important file paths
path.pecan <- "~/Desktop/Research/pecan"

# Defining a site name -- this can go into a function later
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18
GCM.list=c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
# GCM.list=NULL
# GCM.list=c("bcc-csm1-1")
LDAS="NLDAS"
# ens=1:50
ens=1:10
n.ens=length(ens)
# Set up the appropriate seed
set.seed(1159)
seed.vec <- sample.int(1e6, size=500, replace=F)
seed <- seed.vec[min(ens)] # This makes sure that if we add ensemble members, it gets a new, but reproducible seed

# -----------------------------------


# -----------------------------------
# Working through the datasets to be bias-corrected
# General Workflow
# 1. Align Data:
#    1.1. Common temporal resolution 
#         - match training dataset (daily)
#    1.2. Pair Ensemble members
#         A. n train = n out --> pair by member ID
#         B. n train < n out -->
#            --> if n out is multiple of n train --> save multiple out for each in
#            --> if n out not multiple of n train --> randomly pick which members to save multiples for
#         C. n train > n out --> randomly subset members
# 2. Debias Met 
#    2.1. Climatology Debias: Annual means + seasonal cycle
#    2.2. Anomaly Debias: adjust anomaly variance/distribution to try to maintain variance
# 3. Save Met
#    - write years of debiased (or training) met to keep with new ensemble & member ID
# -----------------------------------
GCM.list
n.ens = 10
ens.ID = "TEST"

# 1. Align CRU 6-hourly with LDAS daily
source(file.path(path.pecan, "modules/data.atmosphere/R", "align_met.R"))
train.path <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/NLDAS_day"
source.path <- "~/Desktop/Research/met_ensembles/data/paleon_sites/HARVARD/CRUNCEP"

# For first round, we only want a single in & out
met.out <- align.met(train.path, source.path, n.ens=1, seed=201708, pair.mems = FALSE)

# Calculate wind speed if it's not already there
if(!"wind_speed" %in% names(met.out$dat.source)){
	met.out$dat.source$wind_speed <- sqrt(met.out$dat.source$eastward_wind^2 + met.out$dat.source$northward_wind^2)
}

# 2. Pass the training & source met data into the bias-correction functions; this will get written to the ensemble
source(file.path(path.pecan, "modules/data.atmosphere/R", "debias_met_regression.R"))
debias.met.regression(train.data=met.out$dat.train, source.data=met.out$dat.source, n.ens=10, vars.debias=NULL, CRUNCEP=TRUE,
                      pair.anoms = TRUE, pair.ens = FALSE, uncert.prop="mean", resids = FALSE, seed=Sys.Date(),
                      outfolder="~/Desktop/Research/met_ensembles/data/met_ensembles/TEST", 
                      yrs.save=1901:1979, ens.name="TEST2", ens.mems=NULL, lat.in=site.lat, lon.in=site.lon,
                      save.diagnostics=TRUE, path.diagnostics="~/Desktop/Research/met_ensembles/data/met_ensembles/TEST/bias_correct_qaqc",
                      parallel = FALSE, n.cores = NULL, overwrite = TRUE, verbose = FALSE) 



# -----------------------------------






