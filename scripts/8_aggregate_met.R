# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: 
# Creator: 
# Contact: crollinson@mortonarb.org
# -----------------------------------
# Description
# -----------------------------------
# 
# -----------------------------------
# General Workflow Components
# -----------------------------------
# 0. Set up file structure, load packages, etc
# 1. 
# -----------------------------------
# rm(list=ls())

path.in = "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/1hr/ensembles/CCSM4/CCSM4_001.01/"
years.agg=NULL
save.day=T
save.month=T
out.base = "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/"
day.dir="day/TEST1" 
mo.dir="month/TEST1"
add.vars=c("daylength", "temp.max", "temp.min")
parallel=F
n.cores=NULL
