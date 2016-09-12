# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Comparison of low-frequency trends in the daily GCM output
# Creator: Christy Rollinson, 19 August 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Graph the low-frequency trends in the GCM ensembles against each other
#  -- 10-, 100- year runnign means
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
library(zoo); library(ggplot2)

# Set the working directory
wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble"
# wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

# Defining a site name -- this can go into a function later
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18
GCM.list=c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
n=10 # Number of ensemble members
# -----------------------------------

# -----------------------------------
# Load & Smooth the data
# -----------------------------------
dat.out <- NULL
for(GCM in GCM.list){
  # load the data
  load(file.path(wd.base, "data/met_ensembles", site.name, GCM, "day", paste0(GCM, "_day_alldata.Rdata")))
  
  # Extract the sims for each variable
  for(v in names(dat.out.full)[!names(dat.out.full)=="met.bias"]){
    # Aggregate the data to annual
    dat.ann <- aggregate(dat.out.full[[v]]$sims[,paste0("X",1:n)],
                         by=dat.out.full[[v]]$sims[,c("dataset", "met", "year")],
                         FUN=mean)
    
    # Order the data
    dat.ann <- dat.ann[order(dat.ann$year),]
    
    # Calculate 10- and 100-year running averages
    smooth.10  <- rollapply(dat.ann[,paste0("X",1:n)], width=10 , fill=NA, FUN=mean)
    smooth.100 <- rollapply(dat.ann[,paste0("X",1:n)], width=100, fill=NA, FUN=mean)

    # Calculate the mean & CI *from the smoothed ensemble*
    dat.smooth <- data.frame(GCM=GCM,
                             met=rep(dat.ann$met, 3),
                             year=rep(dat.ann$year, 3),
                             res=c(rep("annual", nrow(dat.ann)),
                                   rep("decadal", nrow(smooth.10)),
                                   rep("centennial", nrow(smooth.100))),
                             mean=c(apply(dat.ann[,paste0("X",1:n)],1,mean, na.rm=T),
                                    apply(smooth.10[,paste0("X",1:n)],1,mean, na.rm=T),
                                    apply(smooth.100[,paste0("X",1:n)],1,mean, na.rm=T)),
                             lwr =c(apply(dat.ann[,paste0("X",1:n)],1,quantile, 0.025, na.rm=T),
                                    apply(smooth.10[,paste0("X",1:n)],1,quantile, 0.025, na.rm=T),
                                    apply(smooth.100[,paste0("X",1:n)],1,quantile, 0.025, na.rm=T)),
                             upr =c(apply(dat.ann[,paste0("X",1:n)],1,quantile, 0.975, na.rm=T),
                                    apply(smooth.10[,paste0("X",1:n)],1,quantile, 0.975, na.rm=T),
                                    apply(smooth.100[,paste0("X",1:n)],1,quantile, 0.975, na.rm=T))
                             )
    
    
    # Append variables and GCMs into one document for graphing
    if(is.null(dat.smooth)){
      dat.out <- dat.smooth
    } else {
      dat.out <- rbind(dat.out, dat.smooth)
    }
  }
  
} # End GCM loop

summary(dat.out)
# -----------------------------------

# -----------------------------------
# Graph the output
# -----------------------------------
dat.out$res <- factor(dat.out$res, levels=c("annual", "decadal", "centennial"))

pdf(file.path(wd.base, "data/met_ensembles", site.name, "LowFrequencyTrends_BiasCorrected.pdf"), height=8, width=11)
for(v in unique(dat.out$met)){
  print(
  ggplot(data=dat.out[dat.out$met==v,]) +
    facet_grid(res~met, scale="free_y") +
    geom_ribbon(data=dat.out[dat.out$met==v & dat.out$res=="annual",], aes(x=year, ymin=lwr, ymax=upr, fill=GCM), alpha=0.2) +
    geom_line(data=dat.out[dat.out$met==v & dat.out$res=="annual",], aes(x=year, y=mean, color=GCM), size=0.5, alpha=0.5) +
    geom_ribbon(data=dat.out[dat.out$met==v & dat.out$res=="decadal",], aes(x=year, ymin=lwr, ymax=upr, fill=GCM), alpha=0.3) +
    geom_line(data=dat.out[dat.out$met==v & dat.out$res=="decadal",], aes(x=year, y=mean, color=GCM), size=1) +
    geom_ribbon(data=dat.out[dat.out$met==v & dat.out$res=="centennial",], aes(x=year, ymin=lwr, ymax=upr, fill=GCM), alpha=0.3) +
    geom_line(data=dat.out[dat.out$met==v & dat.out$res=="centennial",], aes(x=year, y=mean, color=GCM), size=1.5) +
    geom_vline(xintercept=c(1850, 1901, 1980), linetype="dashed") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(name=v) +
    theme_bw() +
    theme(legend.position="top")    
  )
}
dev.off()
# -----------------------------------
