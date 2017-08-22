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
library(lubridate)
library(ggplot2)
# library(tictoc)
rm(list=ls())

wd.base <- "~/Dropbox/PalEON_CR/met_ensemble/"
setwd(wd.base)

path.out <- "~/Desktop/Research/met_ensembles/data/met_ensembles/"
path.pecan <- "~/Desktop/Research/pecan/modules/data.atmosphere/R/"

fig.dir <- file.path(path.out, "model_qaqc")

if(!dir.exists(path.out)) dir.create(path.out, recursive = T)
if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)
# ------------------------------------------

# ------------------------------------------
# 1. Load and format training data
# ------------------------------------------
# ----------
# 1.0 Convert the training data to an nc file
# ----------
{
  var.table <- data.frame(paleon=c("tair", "precipf",  "swdown", "lwdown", "press", "qair", "uas", "vas", "wind"),
                          CF.name=c("air_temperature", "precipitation_flux", "surface_downwelling_shortwave_flux_in_air", 
                                    "surface_downwelling_longwave_flux_in_air", "air_pressure", "specific_humidity", 
                                    "eastward_wind", "northward_wind", "wind_speed"))
  
  dat.train <- read.csv("data/paleon_sites/HARVARD/NLDAS_1980-2015.csv")

  # Renaming our variables to CF convention
  for(i in 1:ncol(dat.train)){
    if(!names(dat.train)[i] %in% var.table$paleon) next
    
    names(dat.train)[i] <- paste(var.table[var.table$paleon==names(dat.train)[i], "CF.name"])
    
  }
  summary(dat.train)
  
  

}
# ----------
# ------------------------------------------

# ------------------------------------------
# 2. Generate the sub-daily models
# ------------------------------------------
# Name of dat.train file in netcdf format meeting CF standards
# dat.trian.nc <- ()
scripts.tdm <- dir(path.pecan, "tdm")
source(file.path(path.pecan, "tdm_generate_subdaily_models.R"))
source(file.path(path.pecan, "tdm_temporal_downscale_functions.R"))
source(file.path(path.pecan, "tdm_model_train.R"))

gen.subdaily.models(outfolder=path.out, dat.train_file=dat.train.nc, in.prefix="HARVARD", 
                                n.beta=500, day.window=2, resids = FALSE, parallel = FALSE, n.cores = NULL, overwrite = TRUE, 
                                verbose = FALSE) 
# ------------------------------------------
