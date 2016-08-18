# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Extract raw met for an individual site
# Creator: Christy Rollinson, 1 July 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Function to extract raw met for a given site for model and save as .csv files for easy handling
# Note: this script (at the moment) only deals with the extraction and saves the 
#       splicing & downscaling for a different script
#
# Note: This script requires CRUNCEP & GCM data to be already downloaded
#
# Vars we need for model input:
# 1. tair
# 2. precipf
# 3. swdown
# 4. lwdown
# 5. qair
# 6. psurf
# 7. wind (ugrd = ws*cos(theta))
# -----------------------------------

# -----------------------------------
# Workflow Overview
# -----------------------------------
# 0) Set up file structure
# Extract Data:
# 1) Extract NLDAS (1980-present) or GLDAS data
#     -- if not already downloaded, will take 36+ hours 
#     -- Use for spatial downscaling
# 2) CRUNCEP: 1901 - 2009
#     -- use for extending empirical/subday back to 1901
# 3) GCM - historical: 1850-2010 - mostly day, (radiation = monthly)
#     - use to fill 1850-1901
# 4) GCM - past millenium: 850-1849 - mostly day, (radiation = monthly)
# -----------------------------------

# Making this into a function
get.raw <- function(wd.base, site.name, lat, lon, ldas.type, GCM.list){
  # Description of declared variables
  # wd.base   = the base directory for this github repository.  
  # site.name = character name of the site; preference for all caps
  # lat       = latitude of the site of interest
  # lon       = longitude of the site of interest
  # ldas.type = type of LDAS data to pull: NLDAS (0.25 degree, hourly)
  # GCM
  # 
  
  
  
# -----------------------------------
# Set up file structure, etc
# -----------------------------------
# Load libraries
library(ncdf4)

# Set the working directory
# wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
# wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

path.out <- file.path(wd.base, "data/paleon_sites", site.name)
if(!dir.exists(path.out)) dir.create(path.out)  

met.done <- dir(path.out, ".csv")
# -----------------------------------

# -----------------------------------
# Get NLDAS / GLDAS data for site
#
# This is our hourly dat that will be used raw for the most recent years
#  and will be used to drive the models from the point at which they start
# ------------
setwd(wd.base)
# Figure out if we already have the LDAS we need
if(!ldas.type %in% substr(met.done, 1, 5)) {
  
  if(!ldas.type %in% c("NLDAS", "GLDAS")) stop("Invalid ldas.type!  Must be either 'NLDAS' or 'GLDAS'")
  
  # Get the LDAS training data for each site if we don't already have it downloaded.  
  path.ldas <- "data/paleon_sites/LDAS" 
  sites.met <- dir(path.ldas)

  # Getting the right capitalization for the site
  dir.site <- sites.met[which(toupper(sites.met) == toupper(site.name))]
  
  # See if the data has already been downloaded
  # Note: making everything uppercase to make case not matter
  if(!(toupper(site.name) %in% toupper(sites.met)) | 
     !length(dir(file.path(path.ldas, dir.site), ldas.type))>0){ 
    print("No LDAS training data found. Downloading. This will take a while")  
    source(paste0("scripts/GetMet_Scripts/download.", ldas.type, ".R"))
    # dir.out <- file.path(path.ldas, site.name)
    if(ldas.type=="NLDAS"){
      download.NLDAS(outfolder=path.ldas, start_date="1980-01-01", end_date="2015-12-31", site_id=site.name, lat.in=lat, lon.in=lon)
    } else {
      download.GLDAS(outfolder=path.ldas, start_date="1980-01-01", end_date="2015-12-31", site_id=site.name, lat.in=lat, lon.in=lon)
    }
  }
  
  # Now that we have training data, move to the next step
  files.site <- dir(file.path(path.ldas, dir.site), ldas.type)
  files.remove <- which(substr(files.site,nchar(files.site)-3, nchar(files.site)) %in% c(".csv", "data")) # make sure we don't try to read an already post-processed .csv or .Rdata file
  if(length(files.remove)>=1){
    files.site <- files.site[!files.remove] 
  }
  
  ldas <- list() # create an empty list to store the data
  for(i in files.site){
    ncT <- nc_open(file.path(path.ldas, dir.site, i))
    
    # Loop through and get all of the variables into 1 list
    for(v in names(ncT$var)){
      if(is.null(ldas[[v]])){ # first time through need to initialize the vectors
        ldas[[v]] <- ncvar_get(ncT, v)
      } else {
        ldas[[v]] <- c(ldas[[v]], ncvar_get(ncT, v))
      }
    }
    
    # write in the year, day, & hour stamps
    nsteps <- length(ncvar_get(ncT, names(ncT$var[1])))
    year=as.numeric(strsplit(i, "[.]")[[1]][2])
    nday = ifelse(lubridate:: leap_year(year), 366, 365)
    time.step = 24/((nsteps)/nday)
    hours = seq(time.step, 24, by=time.step)-1 # Note: Hour marker at the END of the hour!!
  
    time.df <- data.frame(year=year, doy=rep(0:(nday-1), each=24/time.step), hour=rep(hours, nday))
    if(is.null(ldas$time)){
      ldas[["time"]] <- time.df
    } else {
      ldas[["time"]] <- rbind(ldas[["time"]], time.df)
    }
  
    nc_close(ncT)
    
  }
  
  # Converting the list to a dataframe -- not sure why I made it a list in the first place, but whatever
  ldas.df <- data.frame(dataset=ldas.type,
                        ldas$time, 
                        tair=ldas$air_temperature, 
                        precipf=ldas$precipitation_flux, 
                        swdown=ldas$surface_downwelling_shortwave_flux_in_air,
                        lwdown=ldas$surface_downwelling_longwave_flux_in_air,
                        press=ldas$air_pressure,
                        qair=ldas$specific_humidity,
                        uas=ldas$eastward_wind,  # y = x = east  = zonal     = east-west 
                        vas=ldas$northward_wind # v = y = north = meridonal = north=south
                        )
  # calculate single wind: Pathagorus or abs(vwind)/sin(atan(abs(vwind/uwind)))?
  ldas.df$wind <- sqrt(ldas.df$uas^2 + ldas.df$vas^2)
  summary(ldas.df)
  
  write.csv(ldas.df, file.path(path.out, paste0(ldas.type, "_", min(as.numeric(paste(ldas.df$year))), "-", max(as.numeric(paste(ldas.df$year))), ".csv")), row.names = F)
}
# -----------------------------------

# -----------------------------------
# Get CRUNCEPdata
# -----------------------------------
setwd(wd.base)

# Figure out if we've already processed the necessary CRUNCEP data
if(!"CRUNCEP" %in% substr(met.done, 1, 7)) {
  # Get CRUNCEP data for each site if we don't already have it downloaded.  
  #  -- Uses Pecan scrip download.CRUNCEP_Global.R
  path.cruncep <- "data/paleon_sites/CRUNCEP" 
  if(!dir.exists(path.cruncep)) dir.create(path.cruncep)
  sites.met <- dir(path.cruncep)
  
  # Getting the right capitalization for the site
  dir.site <- sites.met[which(toupper(sites.met) == toupper(site.name))]
  
  # See if the data has already been downloaded
  # Note: making everything uppercase to make case not matter
  if(!(toupper(site.name) %in% toupper(sites.met)) | !length(dir(file.path(path.cruncep, dir.site)))>0){ 
    print("No CRUNCEP data found. Downloading. This could take a while")  
    source(paste0("scripts/GetMet_Scripts/download.CRUNCEP_Global.R"))
    # dir.out <- file.path(path.cruncep, site.name)
    download.CRUNCEP(outfolder=path.cruncep, start_date="1901-01-01", end_date="2010-12-31", site_id=site.name, lat.in=lat, lon.in=lon)
  }
  
  # Now that we have training data, move to the next step
  files.site <- dir(file.path(path.cruncep, dir.site))
  files.remove <- which(substr(files.site,nchar(files.site)-3, nchar(files.site)) %in% c(".csv", "data")) # make sure we don't try to read an already post-processed .csv or .Rdata file
  if(length(files.remove)>=1){
    files.site <- files.site[!files.remove] 
  }
  
  cruncep <- list() # create an empty list to store the data
  for(i in files.site){
    ncT <- nc_open(file.path(path.cruncep, dir.site, i))
    
    # Loop through and get all of the variables into 1 list
    for(v in names(ncT$var)){
      if(is.null(cruncep[[v]])){ # first time through need to initialize the vectors
        cruncep[[v]] <- ncvar_get(ncT, v)
      } else {
        cruncep[[v]] <- c(cruncep[[v]], ncvar_get(ncT, v))
      }
    }
    
    # write in the year, day, & hour stamps
    nsteps <- length(ncvar_get(ncT, names(ncT$var[1])))
    year=as.numeric(strsplit(i, "[.]")[[1]][2])
    nday = ifelse(lubridate:: leap_year(year), 366, 365)
    time.step = round(24/((nsteps)/nday),0)
    hours = seq(time.step, 24, by=time.step)-1 # Note: Hour marker at the END of the hour!!
    
    time.df <- data.frame(year=year, doy=rep(0:(nday-1), each=24/time.step), hour=rep(hours, nday))
    if(is.null(cruncep$time)){
      cruncep[["time"]] <- time.df
    } else {
      cruncep[["time"]] <- rbind(cruncep[["time"]], time.df)
    }
    
    nc_close(ncT)
    
  }
  
  # Converting the list to a dataframe -- not sure why I made it a list in the first place, but whatever
  cruncep.df <- data.frame(dataset="CRUNCEP",
                           cruncep$time, 
                           tair=cruncep$air_temperature, 
                           precipf=cruncep$precipitation_flux, 
                           swdown=cruncep$surface_downwelling_shortwave_flux_in_air,
                           lwdown=cruncep$surface_downwelling_longwave_flux_in_air,
                           press=cruncep$air_pressure,
                           qair=cruncep$specific_humidity,
                           uas=cruncep$eastward_wind,  # y = x = east  = zonal     = east-west 
                           vas=cruncep$northward_wind # v = y = north = meridonal = north=south
                           )

  # calculate single wind: Pathagorus or abs(vwind)/sin(atan(abs(vwind/uwind)))?
  cruncep.df$wind <- sqrt(cruncep.df$uas^2 + cruncep.df$vas^2)
  summary(cruncep.df)
  
  write.csv(cruncep.df, file.path(path.out, paste0("CRUNCEP_", min(as.numeric(paste(cruncep.df$year))), "-", max(as.numeric(paste(cruncep.df$year))), ".csv")), row.names = F)
}
# -----------------------------------


# -----------------------------------
# Get GCM data
# -----------------------------------
setwd(wd.base)

for(GCM in GCM.list){

print(paste0(" ------ ", GCM, " ------ "))  
# Right now haven't extracted the paleon domain & are going to work with just MIROC-ESM
dir.gcm <- file.path("data/full_raw", GCM)

if(!dir.exists(dir.gcm)) stop("Invalid GCM!  Are you sure the raw data is downloaded? Check file structure.")

gcm.recode <- data.frame(gcm = c("tas", "tasmax", "tasmin", "pr"     , "psl"  , "huss", "hus" , "sfcWind", "rsds"  , "rlds"  ),
                         met = c("tair", "tmax"  , "tmin"  , "precipf", "press", "qair", "qair", "wind"   , "swdown", "lwdown")
                         )

# --------------
# start with p1000
# --------------
setwd(wd.base)
# Figure out if we need to extract p1000 GCM data
if(!paste0(GCM, "_p1000") %in% substr(met.done, 1, nchar(GCM)+6)) {
  # Start by getting daily for what we can
  vars.gcm.day <- dir(file.path(dir.gcm, "p1000", "day"))
  vars.gcm.mo <- dir(file.path(dir.gcm, "p1000", "month"))
  vars.gcm.mo <- vars.gcm.mo[!vars.gcm.mo %in% vars.gcm.day]
  
  vars.gcm <- c(vars.gcm.day, vars.gcm.mo)
  dat.gcm.p1k=NULL # Giving a value to evaluate to determine if we need to make a new dataframe or not
  for(v in vars.gcm){
    # Getting the daily files first
    freq=ifelse(v %in% vars.gcm.day, "day", "month")
    setwd(file.path(dir.gcm, "p1000", freq, v))
    
    print(paste0("** processing: ", v, " (p1000, ", freq, ")"))
    
    files.v <- dir()
    for(i in 1:length(files.v)){
      fnow=files.v[i]

      # Open the file
      ncT <- nc_open(fnow)
      lat_bnd <- ncvar_get(ncT, "lat_bnds")
      lon_bnd <- ncvar_get(ncT, "lon_bnds")
      nc.time <- ncvar_get(ncT, "time")
      
      # Find the closest grid cell for our site (using harvard as a protoype)
      ind.lat <- which(lat_bnd[1,]<=lat & lat_bnd[2,]>=lat)
      if(max(lon)>=180){
        ind.lon <- which(lon_bnd[1,]>=lon & lon_bnd[2,]<=lon)
      } else {
        ind.lon <- which(lon_bnd[1,]<=180+lon & lon_bnd[2,]>=180+lon)
      }
      
      # extracting the met variable
      if(v %in% c("hus", "ua", "va")){
        plev <- ncvar_get(ncT, "plev")
        puse <- which(plev==max(plev)) # Get humidity at the place of highest pressure (closest to surface)
        dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,puse,] 
        
        # If dat.now has missing values, try the next layer
        puse.orig <- puse
        while(is.na(mean(dat.now))){
          if(puse.orig==1) { puse = puse + 1 } else { puse = puse -1 }
          dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,puse,] 
        }
      } else {
        dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,]
      }
      
      # creating a vector of the time stamp (daily)
      #fday <- (nc.time-0.5)/365 # Note: this is the midpoint
      if(freq=="day"){
        date.start <- as.Date(substr(strsplit(fnow, "_")[[1]][6], 1,8),"%Y%m%d") 
      } else {
        date.start <- as.Date(paste0(substr(strsplit(fnow, "_")[[1]][6], 1,6), "01"),"%Y%m%d") 
      }
      
      
      
      # Some models do NOT do leap year 
      #  -- there's probably an better way to soft code this, but right now I'm going to be 
      #     lazy and hard-code it
      # No leap-year models: bcc-csm1-1
      no.leap <- c("bcc-csm1-1")
      if(freq=="day" & GCM %in% no.leap){
        dates <- seq(from=date.start, length.out=length(dat.now), by=freq)
        leaps <- which(substr(dates,nchar(paste(dates))-4,nchar(paste(dates)))=="02-29")
        dates <- dates[!(1:length(dates) %in% leaps)] # Exclude leap years
        dates <- c(dates, seq(max(dates)+1, length.out=length(leaps), by=freq)) # add extra dates
      } else {
        dates <- seq(from=date.start, length.out=length(dat.now), by=freq)
      }
  
      # Combine everything into a dataframe
      dat.var <- data.frame(year=format(dates, "%Y"), month=format(dates, "%m"), day=format(dates, "%d"), doy=as.POSIXlt(dates, "%Y%b%d")$yday, value=dat.now)
      
      # Add in leap year as duplicate of Feb-28
      if(freq=="day" & GCM %in% no.leap){        
        for(y in unique(dat.var$year)){
          if(lubridate:: leap_year(as.numeric(paste(y)))){
            row.copy <- which(dat.var$year==y & dat.var$month=="02" & dat.var$day=="28")
            
            date.fill <- as.Date(paste0(y, "-02-29"))
            fill <- data.frame(year=format(date.fill, "%Y"), month=format(date.fill, "%m"), day=format(date.fill, "%d"), doy=as.POSIXlt(date.fill, "%Y%b%d")$yday, value=dat.now[row.copy])
            
            dat.var <- rbind(dat.var, fill)
          }
        }
        # Sort the data so it's ordered properly
        dat.var <- dat.var[order(dat.var$year, dat.var$doy),]
      }
      
      # Sticking stuff together
      if(i==1){
        dat.gcm.p1k[[v]] <- dat.var
      } else {
        dat.gcm.p1k[[v]] <- rbind(dat.gcm.p1k[[v]], dat.var)
      }
      nc_close(ncT)
    } # End file loop
    
    setwd(wd.base)
  } # End variable
  
  p1k.day <- data.frame(dataset=paste0(GCM, ".p1000"),
                        dat.gcm.p1k[[vars.gcm.day[1]]][,1:(ncol(dat.gcm.p1k[[vars.gcm.day[1]]])-1)]
                        )
  p1k.mo <- data.frame(dataset=paste0(GCM, ".p1000"),
                       dat.gcm.p1k[[vars.gcm.mo[1]]][,1:(ncol(dat.gcm.p1k[[vars.gcm.mo[1]]])-1)]
                       )

  for(v in vars.gcm.day){
    if(v %in% c("uas", "vas")){
      wind <- sqrt(dat.gcm.p1k[["uas"]]$value^2 + dat.gcm.p1k[["vas"]]$value^2)
      p1k.day[,"wind"] <- wind
    } else if(v %in% c("ua", "va")){
      wind <- sqrt(dat.gcm.p1k[["ua"]]$value^2 + dat.gcm.p1k[["va"]]$value^2)
      p1k.day[,"wind"] <- wind
    } else {
      var2 <- paste(gcm.recode[gcm.recode$gcm==v,"met"])
      p1k.day[,var2] <- dat.gcm.p1k[[v]]$value
    }
  }
  for(v in vars.gcm.mo){
    if(v %in% c("uas", "vas")){
      wind <- sqrt(dat.gcm.p1k[["uas"]]$value^2 + dat.gcm.p1k[["vas"]]$value^2)
      p1k.mo[,"wind"] <- wind
    } else if(v %in% c("ua", "va")){
      wind <- sqrt(dat.gcm.p1k[["ua"]]$value^2 + dat.gcm.p1k[["va"]]$value^2)
      p1k.mo[,"wind"] <- wind
    } else {
      var2 <- paste(gcm.recode[gcm.recode$gcm==v,"met"])
      p1k.mo[,var2] <- dat.gcm.p1k[[v]]$value
    }
  }
  
  p1k.df <- merge(p1k.day, p1k.mo[,!names(p1k.mo) %in% c("day", "doy")], all=T)
  summary(p1k.df)
  
  write.csv(p1k.df, file.path(path.out, paste0(GCM, "_p1000_", min(as.numeric(paste(p1k.df$year))), "-", max(as.numeric(paste(p1k.df$year))), ".csv")), row.names = F)
}
# --------------


# --------------
# Getting historical output
# --------------
setwd(wd.base)

# Figure out if we need to get historical GCM output
if(!paste0(GCM, "_historical") %in% substr(met.done, 1, nchar(GCM)+11)){ 
  # Start by getting daily for what we can
  vars.gcm.day <- dir(file.path(dir.gcm, "historical", "day"))
  vars.gcm.mo <- dir(file.path(dir.gcm, "historical", "month"))
  vars.gcm.mo <- vars.gcm.mo[!vars.gcm.mo %in% vars.gcm.day] # get rid of anything already in day
  
  vars.gcm <- c(vars.gcm.day, vars.gcm.mo)
  dat.gcm.hist=NULL # Giving a value to evaluate to determine if we need to make a new dataframe or not
  for(v in vars.gcm){
    freq=ifelse(v %in% vars.gcm.mo, "month", "day")
    setwd(file.path(dir.gcm, "historical", freq, v))
    
    print(paste0("** processing: ", v, " (historical, ", freq, ")"))
    
    files.v <- dir()
    for(i in 1:length(files.v)){
      fnow=files.v[i]

      # Open the file
      ncT <- nc_open(fnow)
      lat_bnd <- ncvar_get(ncT, "lat_bnds")
      lon_bnd <- ncvar_get(ncT, "lon_bnds")
      nc.time <- ncvar_get(ncT, "time")
      
      # Find the closest grid cell for our site (using harvard as a protoype)
      ind.lat <- which(lat_bnd[1,]<=lat & lat_bnd[2,]>=lat)
      if(lon>=180){
        ind.lon <- which(lon_bnd[1,]>=lon & lon_bnd[2,]<=lon)
      } else {
        ind.lon <- which(lon_bnd[1,]<=180+lon & lon_bnd[2,]>=180+lon)
      }
      
      # extracting the met variable
      if(v %in% c("hus", "ua", "va")){
        plev <- ncvar_get(ncT, "plev")
        puse <- which(plev==max(plev)) # Get humidity at the place of highest pressure (closest to surface)
        dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,puse,] 

        # If dat.now has missing values, try the next layer
        puse.orig <- puse
        while(is.na(mean(dat.now))){
          if(puse.orig==1) { puse = puse + 1 } else { puse = puse -1 }
          dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,puse,] 
        }
        
      } else {
        dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,]
      }
      
      # creating a vector of the time stamp (daily)
      #fday <- (nc.time-0.5)/365 # Note: this is the midpoint
      if(freq=="day"){
        date.start <- as.Date(substr(strsplit(fnow, "_")[[1]][6], 1,8),"%Y%m%d") 
      } else {
        date.start <- as.Date(paste0(substr(strsplit(fnow, "_")[[1]][6], 1,6), "01"),"%Y%m%d") 
      }
      
      
      dates <- seq(from=date.start, length.out=length(dat.now), by=freq)
      
      # Combine everything into a dataframe
      dat.var <- data.frame(year=format(dates, "%Y"), month=format(dates, "%m"), day=format(dates, "%d"), doy=as.POSIXlt(dates, "%Y%b%d")$yday, value=dat.now)
      if(i==1){
        dat.gcm.hist[[v]] <- dat.var
      } else {
        dat.gcm.hist[[v]] <- rbind(dat.gcm.hist[[v]], dat.var)
      }
      nc_close(ncT)
    } # End file loop
    
    setwd(wd.base)
  } # End variable
  
  # Double checking to make sure things are in the right place
  for(v in names(dat.gcm.hist)){
    if(!v %in% vars.gcm.mo){ # if the variable isn't already in month, double check
      if(nrow(dat.gcm.hist[[v]])<nrow(dat.gcm.hist[["pr"]])){ 
        # If this variable has fewer observations than temp, add it to monthly
        vars.gcm.mo <- c(vars.gcm.mo, v)
        vars.gcm.day <- vars.gcm.day[!vars.gcm.day==v]
      }
    }
  }
  
  hist.day <- data.frame(dataset=paste0(GCM, ".hist"),
                         dat.gcm.hist[[vars.gcm.day[1]]][,1:(ncol(dat.gcm.hist[[vars.gcm.day[1]]])-1)]
                         )

  for(v in vars.gcm.day){
    if(v %in% c("uas", "vas")){
      wind <- sqrt(dat.gcm.hist[["uas"]]$value^2 + dat.gcm.hist[["vas"]]$value^2)
      hist.day[,"wind"] <- wind
    } else if(v %in% c("ua", "va")){
      wind <- sqrt(dat.gcm.hist[["ua"]]$value^2 + dat.gcm.hist[["va"]]$value^2)
      hist.day[,"wind"] <- wind
    } else {
      var2 <- paste(gcm.recode[gcm.recode$gcm==v,"met"])
      hist.day[,var2] <- dat.gcm.hist[[v]]$value
    }
  }

  if(length(vars.gcm.mo)>0){
    hist.mo <- data.frame(dataset=paste0(GCM, ".hist"),
                          dat.gcm.hist[[vars.gcm.mo[1]]][,1:(ncol(dat.gcm.hist[[vars.gcm.mo[1]]])-1)]
                          )
    for(v in vars.gcm.mo){
      if(v %in% c("uas", "vas")){
        wind <- sqrt(dat.gcm.hist[["uas"]]$value^2 + dat.gcm.hist[["vas"]]$value^2)
        hist.mo[,"wind"] <- wind
      } else if(v %in% c("ua", "va")){
        wind <- sqrt(dat.gcm.hist[["ua"]]$value^2 + dat.gcm.hist[["va"]]$value^2)
        hist.mo[,"wind"] <- wind
      } else {
        var2 <- paste(gcm.recode[gcm.recode$gcm==v,"met"])
        hist.mo[,var2] <- dat.gcm.hist[[v]]$value
      }
    }

    hist.df <- merge(hist.day, hist.mo[,!names(hist.mo) %in% c("day", "doy")], all=T)
    summary(hist.df)    
  } else {
    hist.df <- hist.day
  }

  
  
  write.csv(hist.df, file.path(path.out, paste0(GCM, "_historical_", min(as.numeric(paste(hist.df$year))), "-", max(as.numeric(paste(hist.df$year))), ".csv")), row.names = F)
}
# --------------
}
# -----------------------------------
}