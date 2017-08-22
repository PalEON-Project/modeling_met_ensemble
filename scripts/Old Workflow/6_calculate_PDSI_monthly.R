# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Calculate PDSI from monthly output
# Creator: Christy Rollinson, 12 June 2018
# Contact: crollinson@mortonarb.org
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Calculating Palmer Drought Severity Index from the met ensemble using
# monthly values gained through aggregating the hourly data in script 
# 5_aggregate_timesteps.R.  For now we'll calibrate on the 1931-1990 
# time period that should be comparable to what B. Cook did and NADA.
#
# AWC is being calculated from our soil driver rather than extracted 
# from SSURGO because SSURGO is a pain and right now we're planning to 
# use the old environmental drivers.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Load libraries, set up file paths, etc
# 1. Load data, calculate AWC for the site
# 3. Calculate & Save PDSI
# -----------------------------------

# -----------------------------------
# 0. Load libraries, set up file paths, etc
# -----------------------------------
library(ncdf4)

path.met <- "~/met_ensemble/data/met_ensembles/HARVARD/monthly_all"
path.soil <- "~/ED_PalEON/MIP2_Region/phase2_env_drivers_v2/soil/"

site.lat <- 42.54
site.lon <- -72.18
# -----------------------------------


# -----------------------------------
# 1. Load data & format, calculate AWC for the site
# -----------------------------------
# -----------------
# Calculate AWC
# -----------------
sand.t <- nc_open(file.path(path.soil, "paleon_soil_t_sand.nc"))
sand.s <- nc_open(file.path(path.soil, "paleon_soil_s_sand.nc"))
clay.t <- nc_open(file.path(path.soil, "paleon_soil_t_clay.nc"))
clay.s <- nc_open(file.path(path.soil, "paleon_soil_s_clay.nc"))
depth  <- nc_open(file.path(path.soil, "paleon_soil_soil_depth.nc"))

lon <- ncvar_get(sand.t, "longitude")
lat <- ncvar_get(sand.t, "latitude")

x.ind <- which(lon-0.25<=site.lon & lon+0.25>=site.lon)
y.ind <- which(lat-0.25<=site.lat & lat+0.25>=site.lat)

sand1 <- ncvar_get(sand.t, "t_sand", c(x.ind, y.ind), c(1,1))
sand2 <- ncvar_get(sand.s, "s_sand", c(x.ind, y.ind), c(1,1))
clay1 <- ncvar_get(clay.t, "t_clay", c(x.ind, y.ind), c(1,1))
clay2 <- ncvar_get(clay.s, "s_clay", c(x.ind, y.ind), c(1,1))
depth2 <- ncvar_get(depth, "soil_depth", c(x.ind, y.ind), c(1,1))


source("calc.awc.R")
awc1 <- calc.awc(sand1, clay1)
awc2 <- calc.awc(sand2, clay2)

wcap1 <- awc1*ifelse(depth2>30, 30, depth2-1) * 1/2.54 # 30 cm top depth * 1 in / 2.54 cm
wcap2 <- awc2*ifelse(depth2>30, depth2-30, 1) * 1/2.54 # remaining depth * 1 in / 2.54 cm
# -----------------

# -----------------
# Load met data, put in Fahrenheit & inches
# -----------------
tair   <- read.csv(file.path(path.met, "tair.csv"  ), row.names=1)
precip <- read.csv(file.path(path.met, "precip.csv"), row.names=1)
daylen <- read.csv(file.path(path.met, "daylength.csv"), row.names=1)

# Formatting the climate data and calculating PDSI
source("pdsi1.R")
pdsi.final <- matrix(NA, nrow=nrow(tair), ncol=ncol(tair)) # A place holder matrix
pdsi2.final <- matrix(NA, nrow=nrow(tair), ncol=ncol(tair)) # A place holder matrix
pdsi3.final <- matrix(NA, nrow=nrow(tair), ncol=ncol(tair)) # A place holder matrix
# pdsi.all <- array(dim=c(1166, 12, ncol(tair))) # Save it all so we can make some easier graphs
# pdsi2.all <- array(dim=c(1166, 12, ncol(tair))) # Save it all so we can make some easier graphs
# pdsi3.all <- array(dim=c(1166, 12, ncol(tair))) # Save it all so we can make some easier graphs
# pdsiM.all <- array(dim=c(1166, 12, ncol(tair))) # Save it all so we can make some easier graphs
# pdsiH.all <- array(dim=c(1166, 12, ncol(tair))) # Save it all so we can make some easier graphs

pb <- txtProgressBar(min = 0, max = ncol(tair), style = 3)

for(i in 1:ncol(tair)){
  if(is.na(mean(tair[,i]))) next
  # Format the climate data into monthly columns
  TEMP1   <- matrix(tair[,i], nrow=nrow(tair)/12, ncol=12, byrow=T)
  PRECIP1 <- matrix(precip[,i], nrow=nrow(precip)/12, ncol=12, byrow=T)
  dayfact <- matrix(daylen[,i], nrow=nrow(daylen)/12, ncol=12, byrow=T)
  row.names(TEMP1) <- 850:2015
  colnames (TEMP1) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  row.names(PRECIP1) <- 850:2015
  colnames (PRECIP1) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  row.names(dayfact) <- 850:2015
  colnames (dayfact) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  # Convert Temp
  # TEMP native units: K
  C2F <- function(x){x*9/5 + 32}
  TEMP1 <- C2F(TEMP1-273.15)
  
  # Convert Precip: 
  # PRECIP native units: kg/m2/mo = mm/mo 
  PRECIP1 <- PRECIP1/25.4 # mm to in
  
  datmet <- list(Temp=TEMP1, Precip=PRECIP1)
  
  # Calculate the proper daylength adjustment
  library(lubridate)
  dpm <- days_in_month(1:12)
  dayfact <- dayfact/12*dpm/30
  
  
  # Adding in everything else we need
  datother <- list()
  datother$pdsi.fun <- "."
  datother$metric <- F
  datother$lat <- site.lat
  datother$watcap <- list(awcs=wcap1, awcu=wcap2)
  datother$yrs.calib <- c(1931, 1990) # copied from B. Cook
  datother$dayz      <- NULL
  datother$dayfact   <- dayfact
  
  siteID <- "HARVARD"
  
  # Run the actual PDSI calculation
  pdsi.out <- pdsi1(datmet, datother, metric=F, siteID, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL)
  
  setTxtProgressBar(pb, i)
  
  # pdsi.all[,,i] <- pdsi.out$X
  # pdsiM.all[,,i] <- pdsi.out$XM
  # pdsiH.all[,,i] <- pdsi.out$XH
  pdsi.final[,i] <- as.vector(t(pdsi.out$X))
  
  # Doing a second calibration
  datother$yrs.calib <- c(850, 950) # modified
  # Run the actual PDSI calculation
  pdsi2.out <- pdsi1(datmet, datother, metric=F, siteID, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL)
  # pdsi2.all[,,i] <- pdsi2.out$X
  pdsi2.final[,i] <- as.vector(t(pdsi2.out$X))
  
  # Doing a third calibration
  datother$yrs.calib <- c(1800, 1850) # modified
  # Run the actual PDSI calculation
  pdsi3.out <- pdsi1(datmet, datother, metric=F, siteID, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL)
  
  # pdsi3.all[,,i] <- pdsi3.out$X
  pdsi3.final[,i] <- as.vector(t(pdsi3.out$X))
}

pdsi.final <- data.frame(pdsi.final)
pdsi2.final <- data.frame(pdsi2.final)
pdsi3.final <- data.frame(pdsi3.final)

names(pdsi.final) <- names(pdsi2.final) <- names(pdsi3.final) <- names(tair)
row.names(pdsi.final) <- row.names(pdsi2.final) <- row.names(pdsi3.final) <- row.names(tair)

# Save the files as .csv
write.csv(pdsi.final , file.path(path.met, "pdsi_calib_1931-1990.csv"), row.names=T)
write.csv(pdsi2.final, file.path(path.met, "pdsi_calib_1890-1850.csv"), row.names=T)
write.csv(pdsi3.final, file.path(path.met, "pdsi_calib_0850-0869.csv"), row.names=T)

pdsi.final  <- read.csv(file.path(path.met, "pdsi_calib_1931-1990.csv"), row.names=1)
pdsi2.final <- read.csv(file.path(path.met, "pdsi_calib_1890-1850.csv"), row.names=1)
pdsi3.final <- read.csv(file.path(path.met, "pdsi_calib_0850-0869.csv"), row.names=1)


pdsi.all <- array(dim=c(1166, 12, ncol(pdsi.final))) # Save it all so we can make some easier graphs
pdsi2.all <- array(dim=c(1166, 12, ncol(pdsi2.final))) # Save it all so we can make some easier graphs
pdsi3.all <- array(dim=c(1166, 12, ncol(pdsi3.final))) # Save it all so we can make some easier graphs
tair.all <- array(dim=c(1166, 12, ncol(pdsi3.final))) # Save it all so we can make some easier graphs
precip.all <- array(dim=c(1166, 12, ncol(pdsi3.final))) # Save it all so we can make some easier graphs
daylen.all <- array(dim=c(1166, 12, ncol(pdsi3.final))) # Save it all so we can make some easier graphs

for(i in 1:ncol(tair)){
  pdsi.all[,,i] <- matrix(pdsi.final[,i], ncol=12, nrow=1166, byrow = T)
  pdsi2.all[,,i] <- matrix(pdsi2.final[,i], ncol=12, nrow=1166, byrow = T)
  pdsi3.all[,,i] <- matrix(pdsi3.final[,i], ncol=12, nrow=1166, byrow = T)
  
  tair.all[,,i] <- matrix(tair[,i], ncol=12, nrow=1166, byrow = T)
  precip.all[,,i] <- matrix(precip[,i], ncol=12, nrow=1166, byrow = T)
  daylen.all[,,i] <- matrix(daylen[,i], ncol=12, nrow=1166, byrow = T)
}

# Some graphing of output
pdsi.ann <- data.frame(apply(pdsi.all, c(1,3), mean, na.rm=T))
pdsi2.ann <- data.frame(apply(pdsi2.all, c(1,3), mean, na.rm=T))
pdsi3.ann <- data.frame(apply(pdsi3.all, c(1,3), mean, na.rm=T))
tair.ann <- data.frame(apply(tair.all, c(1,3), mean, na.rm=T))
precip.ann <- data.frame(apply(precip.all, c(1,3), sum, na.rm=T))
dayl.ann <- data.frame(apply(daylen.all, c(1,3), mean, na.rm=T))
names(pdsi.ann) <- names(tair)
names(pdsi2.ann) <- names(tair)
names(pdsi3.ann) <- names(tair)
names(tair.ann) <- names(tair)
names(precip.ann) <- names(tair)
# names(dayl.ann) <- names(tair)


pdsi <- data.frame(year=850:2015, 
                   median=apply(pdsi.ann, 1, median, na.rm=T),
                   lwr =apply(pdsi.ann, 1, quantile, 0.025, na.rm=T),
                   upr =apply(pdsi.ann, 1, quantile, 0.975, na.rm=T))
pdsi2 <- data.frame(year=850:2015, 
                    median=apply(pdsi2.ann, 1, median, na.rm=T),
                    lwr =apply(pdsi2.ann, 1, quantile, 0.025, na.rm=T),
                    upr =apply(pdsi2.ann, 1, quantile, 0.975, na.rm=T))
pdsi3 <- data.frame(year=850:2015, 
                    median=apply(pdsi3.ann, 1, median, na.rm=T),
                    lwr =apply(pdsi3.ann, 1, quantile, 0.025, na.rm=T),
                    upr =apply(pdsi3.ann, 1, quantile, 0.975, na.rm=T))

tair.summ <- data.frame(year=850:2015, 
                        median=apply(tair.ann, 1, median, na.rm=T),
                        lwr =apply(tair.ann, 1, quantile, 0.025, na.rm=T),
                        upr =apply(tair.ann, 1, quantile, 0.975, na.rm=T))
precip.summ <- data.frame(year=850:2015, 
                          median=apply(precip.ann, 1, median, na.rm=T),
                          lwr =apply(precip.ann, 1, quantile, 0.025, na.rm=T),
                          upr =apply(precip.ann, 1, quantile, 0.975, na.rm=T))
daylen.summ <- data.frame(year=850:2015,
                     median=apply(dayl.ann, 1, median, na.rm=T),
                     lwr =apply(dayl.ann, 1, quantile, 0.025, na.rm=T),
                     upr =apply(dayl.ann, 1, quantile, 0.975, na.rm=T))


library(ggplot2)

pdf(file.path(path.met, "Met_ensemble_Tair_Precip.pdf"))
print(
  ggplot(data=tair.summ) +
    geom_ribbon(aes(x=year, ymin=lwr-273.15, ymax=upr-273.15), alpha=0.5) +
    geom_line(aes(x=year, y=median-273.15), size=1) +
    ylab("Mean Annual Temperature (degrees C)") +
    theme_bw()
)
print(
  ggplot(data=precip.summ) +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
    geom_line(aes(x=year, y=median), size=1) +
    ylab("Total Annual Precipitation (mm)") +
    theme_bw()
)
print(
  ggplot(data=daylen.summ) +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
    geom_line(aes(x=year, y=median), size=1) +
    theme_bw()
)

dev.off()


pdf(file.path(path.met, "PDSI_ensemble.pdf"))
print(
  ggplot(data=pdsi) +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
    geom_line(aes(x=year, y=median), size=1) +
    ggtitle("NADA Calibration (1931-1990)") +
    coord_cartesian(ylim=c(-15, 25)) +
    theme_bw()
)
print(
  ggplot(data=pdsi2) +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
    geom_line(aes(x=year, y=median), size=1) +
    ggtitle("Pre-settlement Calibration 1800-1850") +
    coord_cartesian(ylim=c(-15, 25)) +
    theme_bw()
)
print(
  ggplot(data=pdsi3) +
    geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), alpha=0.5) +
    geom_line(aes(x=year, y=median), size=1) +
    ggtitle("Spinup Calibration (850-869)") +
    coord_cartesian(ylim=c(-25, 25)) +
    theme_bw()
)
dev.off()

# -----------------

# -----------------------------------

