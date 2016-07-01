# Code snippet from 1_get_raw_site.R to get CRUNCEP data so that I can play with formatting things

# -----------------------------------
# Set up file structure, etc
# -----------------------------------
# Load libraries
library(ncdf4)

# Set the working directory
# wd.base <- "~/Dropbox/PalEON_CR/met_ensemble"
wd.base <- "/projectnb/dietzelab/paleon/met_ensemble"
setwd(wd.base)

sites <- data.frame(name = c("HARVARD"),
                    lat  = c( 42.54),
                    lon  = c(-72.18))
# -----------------------------------

# -----------------------------------
# Get CRUNCEPdata
# -----------------------------------
{
  setwd(wd.base)
  dir.cruncep <- "data/paleon_domain/paleon1/cruncep"
  # Get a list of what variables we have to loop through
  vars.cru <- dir(dir.cruncep)
  vars.cru <- vars.cru[!substr(vars.cru, nchar(vars.cru)-3, nchar(vars.cru)) %in% c(".txt", ".htm") & vars.cru!="wind"]
  
  dat.cru=list() # Giving a value to evaluate to determine if we need to make a new dataframe or not
  # NOTE: This will also need to be fixed for multiple sites!
  for(v in vars.cru){
    print(paste0("** processing: ", v))
    
    setwd(file.path(dir.cruncep, v))
    
    files.v <- dir()
    for(i in 1:length(files.v)){
      fnow=files.v[i]
      yr <- as.numeric(substr(fnow, nchar(fnow)-6, nchar(fnow)-3) )
      
      print(paste0("     ", yr))
      
      # Open the file
      ncT <- nc_open(fnow)
      lat <- ncvar_get(ncT, "lat")
      lon <- ncvar_get(ncT, "lon")
      nc.time <- ncvar_get(ncT, "time")
      
      # Turning lat & lon into a bands
      # hard-coding the 0.5-degree resolution of CRUNCEP
      lat <- data.frame(min=lat-0.25, max=lat+0.25) 
      lon <- data.frame(min=lon-0.25, max=lon+0.25)
      
      # Find the closest grid cell for our site (using harvard as a protoype)
      ind.lat <- which(lat[,1]<=sites[1,"lat"] & lat[,2]>=sites[1,"lat"])
      ind.lon <- which(lon[,2]>=sites[1,"lon"] & lon[,1]<=sites[1,"lon"])
      
      # extracting the met variable
      dat.now <- ncvar_get(ncT, v)[ind.lon,ind.lat,]
      
      # creating a vector of the time stamp
      nsteps <- length(dat.now)
      nday = ifelse(lubridate:: leap_year(yr), 366, 365)
      time.step = 24/((nsteps)/nday)
      # time.step = ifelse(time.step==24, 1, time.step) # change 24 hrs into 1 hr
      hours = seq(time.step, 24, by=time.step)-1 # Note: Hour marker at the END of the hour!!
      
      # Combine everything into a dataframe
      dat.var <- data.frame(year=yr, doy=rep(1:nday, each=24/time.step), hour=rep(hours, nday), value=dat.now)
      
      if(i==1){
        dat.cru[[v]] <- dat.var
      } else {
        dat.cru[[v]] <- rbind(dat.cru[[v]], dat.var)
      }
    } # End file loop
    
    setwd(wd.base)
  } # End variable
  
  cru.df <- data.frame(dataset="CRUNCEP",
                       dat.cru[[1]][,1:(ncol(dat.cru[[1]])-1)], 
                       tair=dat.cru$tair$value, 
                       precipf=dat.cru$precipf$value,
                       swdown=dat.cru$swdown$value,
                       lwdown=dat.cru$lwdown$value,
                       pres=dat.cru$psurf$value,
                       qair=dat.cru$qair$value,
                       uas=dat.cru$uwind$value,
                       vas=dat.cru$vwind$value
  )
  # calculate single wind: Pathagorus or abs(vwind)/sin(atan(abs(vwind/uwind)))?
  cru.df$wind <- sqrt(cru.df$uas^2 + cru.df$vas^2)
  summary(cru.df)
  
  write.csv(cru.df, file.path("data/paleon_sites", sites$name, "cruncep_1901-2010.csv"), row.names=F)
  
}
# -----------------------------------
