##' Download and conver to CF CRUNCEP single grid point from MSTIMIP server using OPENDAP interface
##' @name download.NLDAS
##' @title download.NLDAS
##' @export
##' @param outfolder
##' @param start_date
##' @param end_date
##' @param lat
##' @param lon
##'
##' @author Christy Rollinson (based on downlad.CRUNCEP)

download.NLDAS <- function(outfolder, start_date, end_date, site_id, lat.in, lon.in, overwrite=FALSE, verbose=FALSE, ...){  
  # require(PEcAn.utils)
  require(lubridate)
  require(ncdf4)
  require(stringr)
  require(abind)

  # Because the format is grib1, the handy lat/lon extraction of grib2 doesn't work
  # To get a single point, it looks like we're going to have to hard-code the indices
  nlds.lon <- seq(-124.938000, -67.063000, by=0.125)
  nlds.lat <- seq(  25.063000,  52.938000, by=0.125)
  lon.vec <- rep(nlds.lon, times=length(nlds.lat))
  lat.vec <- rep(nlds.lat, length.out=length(lon.vec), each=length(nlds.lon))
  
  lon.lim <- data.frame(min=lon.vec-0.125/2, max=lon.vec+0.125/2)
  lat.lim <- data.frame(min=lat.vec-0.125/2, max=lat.vec+0.125/2)
  
  # Find the cells we want
  ind.rec <- which(lon.lim$min<=lon.in & lon.lim$max>=lon.in & lat.lim$min<=lat.in & lat.lim$max>=lat.in)
  
  # Date stuff
  start_date <- as.POSIXlt(start_date, tz = "GMT")
  end_date <- as.POSIXlt(end_date, tz = "GMT")
  start_year <- year(start_date)
  end_year   <- year(end_date)
  # site_id = as.numeric(site_id)
  outfolder = paste0(outfolder,"_site_",paste0(site_id %/% 1000000000, "-", site_id %% 1000000000))
  
  lat.in = as.numeric(lat.in)
  lon.in = as.numeric(lon.in)
#   lat_trunc = floor(2*(90-as.numeric(lat.in))) # Rounds to nearest degree
#   lon_trunc = floor(2*(as.numeric(lon.in)+180)) # Rounds to nearest degree, adds 180; not sure why
  dap_base="http://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002"
  dir.create(outfolder, showWarnings=FALSE, recursive=TRUE)
  
  ylist <- seq(start_year,end_year,by=1)
  rows = length(ylist)
  results <- data.frame(file=character(rows), host=character(rows),
                        mimetype=character(rows), formatname=character(rows),
                        startdate=character(rows), enddate=character(rows),
                        dbfile.name = "NLDAS",
                        stringsAsFactors = FALSE)
  
  var = data.frame(DAP.name = c("TMP","DLWRF","PRES","DSWRF","UGRD","VGRD","SPFH","APCP"),
                   CF.name = c("air_temperature","surface_downwelling_longwave_flux_in_air","air_pressure","surface_downwelling_shortwave_flux_in_air","eastward_wind","northward_wind","specific_humidity","precipitation_flux"),
                   units = c('Kelvin',"W/m2","Pascal","W/m2","m/s","m/s","g/g","kg/m2/s")
  )
  
  for (i in 1:rows){
    year = ylist[i]    
    
    # figure out how many days we're working with
    if(rows>1 & !i==1 & !i==rows){ # If we have multiple years and we're not in the first or last year, we're taking a whole year
      nday  = ifelse(year%%4 == 0,366,365) # leap year or not; days per year
      days.use = 1:nday
    } else if(rows==1){
      # if we're working with only 1 year, lets only pull what we need to
      nday  = ifelse(year%%4 == 0,366,365) # leap year or not; days per year
      day1 <- yday(start_date)
      # Now we need to check whether we're ending on the right day
      day2 <- yday(end_date)
      days.use = day1:day2
      nday=length(days.use) # Update nday
    } else if(i==1) {
      # If this is the first of many years, we only need to worry about the start date
      nday  = ifelse(year%%4 == 0,366,365) # leap year or not; days per year
      day1 <- yday(start_date)
      days.use = day1:nday
      nday=length(days.use) # Update nday
    } else if(i==1) {
      # If this is the last of many years, we only need to worry about the start date
      nday  = ifelse(year%%4 == 0,366,365) # leap year or not; days per year
      day2 <- yday(end_date)
      days.use = 1:day2
      nday=length(days.use) # Update nday
    }
    ntime = nday*24 # leap year or not;time slice (hourly)

    loc.file = file.path(outfolder,paste("NLDAS",year,"nc",sep="."))
    
    ## Create dimensions
    lat <- ncdim_def(name='latitude', units='degree_north', vals=lat.in, create_dimvar=TRUE)
    lon <- ncdim_def(name='longitude', units='degree_east', vals=lon.in, create_dimvar=TRUE)
    time <- ncdim_def(name='time', units="sec", vals=seq((min(days.use)*24*360), (max(days.use)+1-1/24)*24*360, length.out=ntime), create_dimvar=TRUE, unlim=TRUE)
    dim=list(lat,lon,time)
    
    var.list = list()
    dat.list = list()

    # Defining our dimensions up front
    for(j in 1:nrow(var)){
      var.list[[j]] = ncvar_def(name=as.character(var$CF.name[j]), units=as.character(var$units[j]), dim=dim, missval=-999, verbose=verbose)
    }

    ## get data off OpenDAP
    for(j in days.use){
      date.now <- as.Date(j, origin=as.Date(paste0(year-1,"-12-31")))
      mo.now <- str_pad(month(date.now), 2, pad="0")
      day.mo <- str_pad(day(date.now), 2, pad="0")
      doy <- str_pad(j, 3, pad="0")
      for(h in seq(0000, 2300, by=100)){
        hr <- str_pad(h,4,pad="0")
        dap_file = paste0(dap_base, "/",year, "/", doy, "/","NLDAS_FORA0125_H.A",year,mo.now,day.mo, ".",hr, ".002.grb")
        
        # Temporarily download the file using curl
        system(paste0("curl -f -v -s ", dap_file, " -o tmp.grb"), ignore.stderr=T)
        for (v in unique(var$DAP.name)) {
          # create the wgrib command
          wg.str <- paste0("wgrib -s tmp.grb | grep \":", 
                           v, ":\" | wgrib -V -i -text tmp.grb -o var_out.txt")
          system(wg.str, ignore.stderr = TRUE) # execute the wgrib command, creating a temporary text file
          vals <- scan("var_out.txt", skip = 1,quiet = TRUE)[ind.rec] # grab the cells we want
          vals <- array(vals, dim=c(length(lat.in), length(lon.in), length(vals)))
          if(v %in% names(dat.list)){
            # If the variable allready exists, use abind to add the new data
            dat.list[[v]] <- abind(dat.list[[v]], vals, along=3)
          } else {
            # If this is the first time through, make a new 3-D array for that variable
            dat.list[[v]] <- vals
          }
          system("rm -f var_out.txt", ignore.stderr=T) # clean up after ourselves!
        } # end variable loop
        system("rm -f tmp.grb", ignore.stderr=T) # clean up after ourselves!
      } # end hour
    } # end day
    ## change units of precip to kg/m2/s instead of hour accumulated precip
    dat.list[["APCP"]] = dat.list[["APCP"]]/360
    
    ## put data in new file
    loc <- nc_create(filename=loc.file, vars=var.list, verbose=verbose)
    for(j in 1:nrow(var)){
      ncvar_put(nc=loc, varid=as.character(var$CF.name[j]), vals=dat.list[[j]])
    }
    nc_close(loc)
    
    results$file[i] <- loc.file
#     results$host[i] <- fqdn()
    results$startdate[i] <- paste0(year,"-01-01 00:00:00")
    results$enddate[i] <- paste0(year,"-12-31 23:59:59")
    results$mimetype[i] <- 'application/x-netcdf'
    results$formatname[i] <- 'CF Meteorology'
    
  }
  
  invisible(results)
}


