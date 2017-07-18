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
library(zoo); library(ggplot2); library(stringr); library(ncdf4)

# Set the working directory
# wd.base <- "~/Desktop/Research/PalEON_CR/met_ensemble"
wd.base <- "~/met_ensemble"
# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)
dir.dat <- "~/met_ensemble/data/met_ensembles"
# Defining a site name -- this can go into a function later
site.name="HARVARD"
site.lat=42.54
site.lon=-72.18
GCM.list=c("MIROC-ESM", "MPI-ESM-P", "bcc-csm1-1", "CCSM4")
ens.type <- "1hr" # set for day or sub-day to check on
# n=25 # Number of ensemble members
# -----------------------------------

# -----------------------------------
# Extract the data
# -----------------------------------
vars <- c("tair", "precipf", "swdown", "lwdown", "press", "qair", "wind")
dat.out <- list()
for(v in vars){
  dat.out[[v]] <- data.frame(Year = 850:2015)
}


for(GCM in GCM.list){
  print(GCM)
  # Get list of ensemble members
  ens.list <- dir(file.path(dir.dat, site.name, GCM, ens.type), paste0(site.name, "_", GCM, "_", ens.type))

  pb.index=1
  pb <- txtProgressBar(min=1, max=length(850:2015)*length(ens.list), style=3)
  # loop through ensemble members and calculate annual means
  for(ens in ens.list){
    files.ens <- dir(file.path(dir.dat, site.name, GCM, ens.type, ens))
    
    # loop through individual files
    for(fnow in files.ens){
      yr.now <- as.numeric(str_split(str_split(fnow, "_")[[1]][[5]], "[.]")[[1]][1])
      
      # Open the nc file
      ncT <- nc_open(file.path(file.path(dir.dat, site.name, GCM, ens.type, ens, fnow)))
      
      for(v in vars){
        dat.out[[v]][dat.out[[v]]$Year == yr.now, ens] <- mean(ncvar_get(ncT, v))
      }
      
      nc_close(ncT)
      
      # Update our Progress Bar
      setTxtProgressBar(pb, pb.index) 
      pb.index = pb.index+1
    }
  }
  print("")
} # End GCM loop

summary(dat.out)
# -----------------------------------

# -----------------------------------
# Smooth & summarize the data
# -----------------------------------
dat.out.10 <- dat.out.100 <- list()
for(v in names(dat.out)){
  dat.out.10[[v]]  <- data.frame(Year = 850:2015)
  dat.out.100[[v]] <- data.frame(Year = 850:2015)
}

# Getting the running means
for(v in names(dat.out)){
  dat.out.10[[v]][,2:ncol(dat.out[[v]])]  <- rollapply(dat.out[[v]][,2:ncol(dat.out[[v]])], width=10 , fill=NA, FUN=mean, by.column=T, align="center")
  dat.out.100[[v]][,2:ncol(dat.out[[v]])] <- rollapply(dat.out[[v]][,2:ncol(dat.out[[v]])], width=100, fill=NA, FUN=mean, by.column=T, align="center")
}

# Creating an index for GCMs
gcm.names <- vector(length=ncol(dat.out[[1]]))
for(i in 2:ncol(dat.out[[1]])){
  gcm.names[i] <- str_split(names(dat.out[[1]])[i], "_")[[1]][[2]]
}

# Packaging the data for graphing
years <- 850:2015
dat.smooth <- NULL
for(GCM in GCM.list){
  cols.gcm <- which(gcm.names == GCM)
  for(v in names(dat.out)){
    dat.gcm <- data.frame(GCM=GCM,
                          met.var = v,
                          res=rep(c("annual", "decadal", "centennial"), each=length(years)),
                          mean = c(apply(dat.out    [[v]][,cols.gcm], 1, mean , na.rm=T),
                                   apply(dat.out.10 [[v]][,cols.gcm], 1, mean , na.rm=T),
                                   apply(dat.out.100[[v]][,cols.gcm], 1, mean , na.rm=T)),
                          lwr  = c(apply(dat.out    [[v]][,cols.gcm], 1, quantile, 0.025, na.rm=T),
                                   apply(dat.out.10 [[v]][,cols.gcm], 1, quantile, 0.025, na.rm=T),
                                   apply(dat.out.100[[v]][,cols.gcm], 1, quantile, 0.025, na.rm=T)),
                          upr  = c(apply(dat.out    [[v]][,cols.gcm], 1, quantile, 0.975, na.rm=T),
                                   apply(dat.out.10 [[v]][,cols.gcm], 1, quantile, 0.975, na.rm=T),
                                   apply(dat.out.100[[v]][,cols.gcm], 1, quantile, 0.975, na.rm=T))
                          )
    if(is.null(dat.smooth)){
      dat.smooth <- dat.gcm
    } else {
      dat.smooth <- rbind(dat.smooth, dat.gcm)
    }
  } # end met var loop
} # end GCM loop

# -----------------------------------



# -----------------------------------
# Graph the output
# -----------------------------------
dat.out$res <- factor(dat.out$res, levels=c("annual", "decadal", "centennial"))

pdf(file.path(dir.dat, site.name, paste0("LowFrequencyTrends_", ens.type, "_Ensemble.pdf")), height=8, width=11)
for(v in unique(dat.out$met)){
  print(
  ggplot(data=dat.out[dat.out$met==v,]) +
    facet_grid(res~met, scales="free_y") +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr, fill=GCM, alpha=res)) +
    geom_line(aes(x=year, y=mean, color=GCM, size=res)) +
    geom_vline(xintercept=c(1850, 1901, 1980), linetype="dashed") +
    scale_alpha_discrete(values=c(0.2, 0.3, 0.3)) +
    scale_size_discrete(values=c(0.3, 1, 1.5)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(name=v) +
    guides(size=F, alpha=F) +
    theme_bw() +
    theme(legend.position="top")    
  )
}
dev.off()
# -----------------------------------
