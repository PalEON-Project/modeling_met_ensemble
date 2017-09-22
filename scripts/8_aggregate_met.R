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
source("aggregate_met.R")
source("aggregate_file.R")

in.base = "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/1hr/ensembles/"
out.base = "~/Desktop/Research/met_ensembles/data/met_ensembles/HARVARD/aggregated"

GCM.list <- dir(in.base)
for(GCM in GCM.list){
  gcm.ens <- dir(file.path(in.base, GCM))
  pb <- txtProgressBar(min=0, max=length(gcm.ens), style=3)
  pb.ind=1
  for(ens in gcm.ens){
    aggregate.met(path.in=file.path(in.base, GCM, ens), 
                  years.agg=NULL, save.day=T, save.month=T, 
                  out.base=out.base, day.dir=file.path("day", GCM, ens), mo.dir=file.path("month", GCM, ens), 
                  add.vars=c("daylength", "air_temperature_maximum", "air_temperature_minimum"),
                  parallel=T, n.cores=8, 
                  print.progress=F, verbose=FALSE)
    
    setTxtProgressBar(pb, pb.ind)
    pb.ind=pb.ind+1
  }
}
