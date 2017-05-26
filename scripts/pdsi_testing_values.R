# Values for testing PDSI

# For test: awc for Harvard Forest
awcs <- 0.32
awcu <- 0.14


# Climate pulling one site from my dissertation
dat.clim <- read.csv("~/Dropbox/carca/Data - csv/ClimateData-BlueKnob.csv")
dat.clim$tavg <- rowMeans(dat.clim[,c("tmin", "tmax")])
dat.clim <- dat.clim[dat.clim$plotID=="BLU1A",]
dat.clim$year <- as.factor(paste(dat.clim$year))
dat.clim$month <- as.ordered(paste(dat.clim$month))
dat.clim$month <- factor(dat.clim$month, levels=paste(1:12))
summary(dat.clim)

library(reshape2)
temp <- recast(dat.clim[,c("year", "month", "tavg")], year ~ month)
row.names(temp) <- temp$year
temp <- temp[complete.cases(temp), 2:13]
summary(temp)

precip <- recast(dat.clim[,c("year", "month", "precip")], year ~ month)
row.names(precip) <- precip$year
precip <- precip[complete.cases(precip), 2:13]
summary(precip)

awcs <- 0.32 # cm3/cm3 = in3/in3
awcu <- 0.14 # cm3/cm3 = in3/in3

Temp <- temp
Precip <- precip
awcs <- awcs*10*30 # in mm w/ 30 cm depth
awcu <- awcu*10*30 # in mm w/ 30 cm depth

# English unit conversions
# C2F <- function(x){x*9/5 + 32}
# Temp <- C2F(temp)
# Precip <- precip/25.4
# lat <- mean(dat.clim$lat.plot)

# awcs*(1/2.54^3)
# awcs <- awcs*11.81 # assuming top depth of 30 cm; x10^3 to put in mm3
# awcu <- awcu*11.81 # assuming under depth of 30 cm; x10 to put in mm


# Converting 
library(R.matlab)
dayz <- readMat("PDSI_fromBenCook/PDSICODE/daylennh.mat")$dayz


pdsi.inches <- pdsi$x
