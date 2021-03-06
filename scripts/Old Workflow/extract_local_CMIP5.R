##' Extract NLDAS from local download
##' Extract NLDAS meteorology for a poimt from a local download of the full grid
# ----------------------------------- 
# Description
# -----------------------------------
##' @title extract.local.CMIP5
##' @family 
##' @author Christy Rollinson, 
##' @description This function extracts CMIP5 data from grids that have been downloaded and stored locally.
##'              Files are saved as a netCDF file in CF conventions at *DAILY* resolution.  Note: At this point
##'              in time, variables that are only available at a native monthly resolution will be repeated to
##'              give a pseudo-daily record (and can get dealt with in the downscaling workflow).  These files 
##'              are ready to be used in the general PEcAn workflow or fed into the downscaling workflow.
# ----------------------------------- 
# Parameters
# -----------------------------------
##' @param outfolder - directory where output files will be stored
##' @param in.path - path to the raw full grids
##' @param start_date - first day for which you want to extract met (yyyy-mm-dd)
##' @param end_date - last day for which you want to extract met (yyyy-mm-dd)
##' @param site_id name to associate with extracted files
##' @param lat.in site latitude in decimal degrees
##' @param lon.in site longitude in decimal degrees
##' @param model which GCM to extract data from
##' @param scenario which experiment to pull (p1000, historical, ...)
##' @param ensemble_member which CMIP5 experiment ensemble member
##' @param overwrite logical. Download a fresh version even if a local file with the same name already exists?
##' @param verbose logical. Passed on to \code{\link[ncdf4]{ncvar_def}} and \code{\link[ncdf4]{nc_create}}
##'   to control printing of debug info
##' @param ... Other arguments, currently ignored##' @export
##' @examples
# -----------------------------------
extract.local.CMIP5 <- function(outfolder, in.path, start_date, end_date, site_id, lat.in, lon.in, 
                                model , scenario , ensemble_member = "r1i1p1",
                                overwrite = FALSE, verbose = FALSE, ...){
  library(lubridate)
  library(ncdf4)
  library(stringr)
  
  # Some GCMs don't do leap year; we'll have to deal with this separately
  no.leap <- c("bcc-csm1-1", "CCSM4")
  
  # Date stuff
  start_date <- as.POSIXlt(start_date, tz = "GMT")
  end_date <- as.POSIXlt(end_date, tz = "GMT")
  start_year <- year(start_date)
  end_year   <- year(end_date)
  
  lat.in = as.numeric(lat.in)
  lon.in = as.numeric(lon.in)
  # dir.nldas="http://hydro1.sci.gsfc.nasa.gov/thredds/dodsC/NLDAS_FORA0125_H.002"
  dir.create(outfolder, showWarnings=FALSE, recursive=TRUE)
  
  ylist <- seq(start_year,end_year,by=1)
  rows = length(ylist)
  results <- data.frame(file=character(rows), host=character(rows),
                        mimetype=character(rows), formatname=character(rows),
                        startdate=character(rows), enddate=character(rows),
                        dbfile.name = "NLDAS",
                        stringsAsFactors = FALSE
  )
  
  # The table of var name conversion
  # psl; sfcWind; tasmax; tasmin; huss
  var <- data.frame(DAP.name = c("tas", "tasmax", "tasmin", "rlds", "ps", "rsds", "uas", "vas", "sfcWind", "huss", "pr"), 
                    CF.name = c("air_temperature", "air_temperature_maximum", "air_temperature_minimum", 
                                "surface_downwelling_longwave_flux_in_air",
                                "air_pressure", "surface_downwelling_shortwave_flux_in_air", 
                                "eastward_wind", "northward_wind", "wind_speed", "specific_humidity", "precipitation_flux"), 
                    units = c("Kelvin", "Kelvin", "Kelvin", "W/m2", "Pascal", "W/m2", "m/s", "m/s", "m/s", "g/g", "kg/m2/s"))
  
  # Figuring out what we have daily for and what we only have monthly for
  vars.gcm.day <- dir(file.path(in.path, "day"))
  vars.gcm.mo <- dir(file.path(in.path, "month"))
  vars.gcm.mo <- vars.gcm.mo[!vars.gcm.mo %in% vars.gcm.day]
  
  vars.gcm <- c(vars.gcm.day, vars.gcm.mo)
  
  # Rewriting the dap name to get the closest variable that we have for the GCM (some only give uss stuff at sea level)
  library(car) # having trouble gettins stuff to work otherwise
  if(!("huss" %in% vars.gcm)) var$DAP.name <- recode(var$DAP.name, "'huss'='hus'")
  if(!("ps" %in% vars.gcm  )) var$DAP.name <- recode(var$DAP.name, "'ps'='psl'")
  
  # Making sure we're only trying to grab the variables we have (i.e. don't try sfcWind if we don't have it)
  var <- var[var$DAP.name %in% vars.gcm,]
  
  # Native CMIP5 file structure is organized by variable and then with multiple years per file
  # this means we need to do some funky things to get all variables for one year into a single file
  var$DAP.name <- as.character(var$DAP.name)
  
  files.var <- list()
  n.file=0
  for(v in var$DAP.name){
  	files.var[[v]] <- list()
    if(v %in% vars.gcm.day){
	  # Get a list of file names
      files.var[[v]][["files"]] <- dir(file.path(in.path, "day", v))  		
  	} else {
  	  files.var[[v]][["files"]] <- dir(file.path(in.path, "month", v))
  	}
  	
	# Set up an index to help us find out which file we'll need
    files.var[[v]][["years"]] <- data.frame(first.year=NA, last.year=NA)
    for(i in 1:length(files.var[[v]][["files"]])){
    	yr.str <- str_split(str_split(files.var[[v]][["files"]][[i]], "_")[[1]][6], "-")[[1]]
  		
    	# Don't bother storing this file if we don't want those years
    	if(as.numeric(substr(yr.str[1], 1, 4)) > end_year | as.numeric(substr(yr.str[2], 1, 4))< start_year) next
    	files.var[[v]][["years"]][i, "first.year"] <- as.numeric(substr(yr.str[1], 1, 4))
  		files.var[[v]][["years"]][i, "last.year" ] <- as.numeric(substr(yr.str[2], 1, 4))

  		n.file=n.file+1
  	 } # End file loop
  } # end variable loop
  
  
  # Querying large netcdf files 1,000 times is slow.  So lets open the connection once and 
  # pull the full time series
  # Loop through using the files using the first variable; shoudl be tair & should be highest res avail
  # This will require quite a bit of memory, but it's doable
  dat.all <- list()
  dat.time <- seq(start_date, end_date, by="day")  # Everything shoudl end up being a day
  
  print("- Extracting files: ")
  pb <- txtProgressBar(min=1, max=n.file, style=3)
  pb.ind=1
  # Loop through each variable so that we don't have to open files more than once
  for(v in 1:nrow(var)){
    
    var.now <- var[v,"DAP.name"]
    # print(var.now)
    
    dat.all[[v]] <- vector() # initialize the layer
    # Figure out the temporal resolution of the variable
    v.res <- ifelse(var.now %in% vars.gcm.day, "day", "month")

    # Figure out what file we need
    # file.ind <- which(files.var[[var.now]][i])
    for(i in 1:length(files.var[[var.now]]$files)){
      setTxtProgressBar(pb, pb.ind)
      pb.ind=pb.ind+1
      f.now <- files.var[[var.now]]$files[i]
      # print(f.now)
      
      # Open up the file
      ncT <- nc_open(file.path(in.path, v.res, var.now, f.now))
      
      # Extract our dimensions
      lat_bnd <- ncvar_get(ncT, "lat_bnds")
      lon_bnd <- ncvar_get(ncT, "lon_bnds")
      nc.time <- ncvar_get(ncT, "time")
      
      # splt.ind <- ifelse(GCM %in% c("MPI-ESM-P"), 4, 3)
      # date.origin <- as.Date(str_split(ncT$dim$time$units, " ")[[1]][splt.ind])
 
      # Find the closest grid cell for our site (using harvard as a protoype)
      ind.lat <- which(lat_bnd[1,]<=lat.in & lat_bnd[2,]>=lat.in)
      if(max(lon.in)>=180){
        ind.lon <- which(lon_bnd[1,]>=lon.in & lon_bnd[2,]<=lon.in)
      } else {
        ind.lon <- which(lon_bnd[1,]<=180+lon.in & lon_bnd[2,]>=180+lon.in)
      }
      
      # Extract all of the available data
      if(var.now %in% c("hus", "ua", "va")){ # These have multiple strata; we only want 1
        plev <- ncvar_get(ncT, "plev")
        puse <- which(plev==max(plev)) # Get humidity at the place of highest pressure (closest to surface)
        dat.temp <- ncvar_get(ncT, var.now, c(ind.lon, ind.lat, puse, 1), c(1,1,1,length(nc.time)))
        # If dat.list has missing values, try the next layer
        puse.orig <- puse
        while(is.na(mean(dat.temp))){
          if(puse.orig==1) { puse = puse + 1 } else { puse = puse -1 }
          dat.temp <- ncvar_get(ncT, var.now, c(ind.lon, ind.lat, puse, 1), c(1,1,1,length(nc.time)))
        }
      } else {
        dat.temp <- ncvar_get(ncT, var.now, c(ind.lon, ind.lat, 1), c(1,1,length(nc.time)))
      }
      
      # If we have daily data and we're dealing with a model that skips leap year, add it in
      if(GCM %in% no.leap & v.res == "day" & leap_year(y.now) & length(dat.temp[[v]]) < nday){
        dat.temp <- append(dat.temp, dat.temp[sum(dpm[1:2])], sum(dpm[1:2]))
      }
      
      # If we have monthly data, lets trick it into being daily
      if(v.res == "month"){
        mo.ind <- rep(1:12, length.out=length(dat.temp))
        dat.trick <- vector()
        for(j in 1:length(dat.temp)){
          dat.trick <- c(dat.trick, rep(dat.temp[j], dpm[mo.ind[j]]))
        }
        dat.temp <- dat.trick
      } # End leap day trick
      
      dat.all[[v]] <- append(dat.all[[v]], dat.temp, length(dat.all[[v]]))
      nc_close(ncT)    	
    } # End file loop
  } # End variable loop
    

  print("")
  print("- Writing to NetCDF: ")
  pb <- txtProgressBar(min=1, max=rows, style=3)
  for (i in 1:rows){
    setTxtProgressBar(pb, i)
    
    y.now = ylist[i]    
    yr.ind <- which(year(dat.time)==y.now)
    
    
    dpm <- days_in_month(1:12)
    if(leap_year(y.now)) dpm[2] <- dpm[2] + 1 # make sure Feb has 29 days if we're dealing with a leap year
    
    # figure out how many days we're working with
    if(rows>1 & i!=1 & i!=rows){ # If we have multiple years and we're not in the first or last year, we're taking a whole year
      nday  = ifelse(lubridate:: leap_year(y.now), 366, 365) # leap year or not; days per year
      day1 = 1
      day2 = nday
      days.use = day1:day2
    } else if(rows==1){
      # if we're working with only 1 year, lets only pull what we need to
      nday  = ifelse(lubridate:: leap_year(y.now), 366, 365) # leap year or not; days per year
      day1 <- yday(start_date)
      # Now we need to check whether we're ending on the right day
      day2 <- yday(end_date)
      days.use = day1:day2
      nday=length(days.use) # Update nday
    } else if(i==1) {
      # If this is the first of many years, we only need to worry about the start date
      nday  = ifelse(lubridate:: leap_year(y.now), 366, 365) # leap year or not; days per year
      day1 <- yday(start_date)
      day2 = nday
      days.use = day1:day2
      nday=length(days.use) # Update nday
    } else if(i==rows) {
      # If this is the last of many years, we only need to worry about the start date
      nday  = ifelse(lubridate:: leap_year(y.now), 366, 365) # leap year or not; days per year
      day1 = 1
      day2 <- yday(end_date)
      days.use = day1:day2
      nday=length(days.use) # Update nday
    }
    ntime = nday # leap year or not; time slice (coerce to daily)
    
    loc.file <- file.path(outfolder, paste(model, scenario, ensemble_member, str_pad(y.now, width=4, side="left",  pad="0"), "nc", sep = "."))
    
    
    ## Create dimensions
    dim.lat <- ncdim_def(name='latitude', units='degree_north', vals=lat.in, create_dimvar=TRUE)
    dim.lon <- ncdim_def(name='longitude', units='degree_east', vals=lon.in, create_dimvar=TRUE)
    dim.time <- ncdim_def(name='time', units="sec", vals=seq((min(days.use)+1-1/24)*24*360, (max(days.use)+1-1/24)*24*360, length.out=ntime), create_dimvar=TRUE, unlim=TRUE)
    nc.dim=list(dim.lat,dim.lon,dim.time)
    
    
    # Defining our dimensions up front
    var.list = list()
    dat.list = list()

    for(j in 1:nrow(var)){
      var.list[[j]] = ncvar_def(name=as.character(var$CF.name[j]), units=as.character(var$units[j]), dim=nc.dim, missval=-999, verbose=verbose)
      dat.list[[j]] <- array(NA, dim=c(length(lat.in), length(lon.in), ntime)) # Go ahead and make the arrays
    }
    names(var.list) <- names(dat.list) <- var$CF.name
    
    # Loop through each variable in the order of everything else
    for(v in 1:nrow(var)){
	    	dat.list[[v]] <- dat.all[[v]][yr.ind]	
    } # End variable loop
        
    ## put data in new file
    loc <- nc_create(filename=loc.file, vars=var.list, verbose=verbose)
    for(j in 1:nrow(var)){
      ncvar_put(nc=loc, varid=as.character(var$CF.name[j]), vals=dat.list[[j]])
    }
    nc_close(loc)
    
    results$file[i] <- loc.file
    # results$host[i] <- fqdn()
    results$startdate[i]  <- paste0(as.Date(paste(y.now, day1, sep="-"), format = "%Y-%j"), " 00:00:00")
    results$enddate[i]    <- paste0(as.Date(paste(y.now, day2, sep="-"), format = "%Y-%j"), " 00:00:00")
    results$mimetype[i]   <- 'application/x-netcdf'
    results$formatname[i] <- 'CF Meteorology'
    
  } # End i loop (rows/years)
  
} # End function

