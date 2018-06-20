##' Aggregate File 
# ----------------------------------- 
# Description
# -----------------------------------
##' @title aggregate.file
##' @family Meteorology
##' @author Christy Rollinson
##' @description Take a single sub-daily met file and turn it into daily &/or monthly data; saving by year as we go
##'              This gets called by aggregate_met and is a separate function so that it can be easily parallelized
# ----------------------------------- 
# Parameters
# -----------------------------------
##' @param f.agg - file path passed from function above
##' @param path.in File path to where the .nc files should be housed
##' @param day [TRUE]/FALSE, whether to save daily data
##' @param month [TRUE]/FALSE whether to calculate/save daily data
##' @param out.base root file path for where you want to write data
##' @param day.dir subdirectory you want to save daily output in (defaults to "day")
##' @param mo.dir subdirectory you want to save monthly output in (defaults to "month")
##' @param add.vars c("daylength", "temp.max", "temp.min") additional handy variables that we can compute other than just mean
##' @return Automatically saved files withthe same naming structure as previously, just in a different location (this is dangerous!)
##' @export
# -----------------------------------
<<<<<<< HEAD
aggregate.file <- function(f.agg, path.in, save.day, save.month, out.base, day.dir, mo.dir, add.vars, ...){
=======
aggregate.file <- function(f.agg, path.in, save.day, save.month, out.base, day.dir, mo.dir, add.vars, verbose=FALSE, ...){
>>>>>>> Ulmus
  yr.now <- strsplit(f.agg, "[.]")[[1]]
  yr.now <- as.numeric(yr.now[length(yr.now)-1])
  nday <- ifelse(lubridate::leap_year(yr.now), 366, 365)
  
  # Open the file
  ncT <- ncdf4::nc_open(file.path(path.in, f.agg))
  
  # Extract some useful info
  vars.file <- names(ncT$var)
  nc.time <- ncdf4::ncvar_get(ncT, "time")
<<<<<<< HEAD
=======
  
  # Calculate the timestep
  dt <- nday/length(nc.time) # Time step in days
  nc.time <- seq(dt/2, length.out=length(nc.time), by=dt)
  
>>>>>>> Ulmus
  lat.in <- ncdf4::ncvar_get(ncT, "latitude")
  lon.in <- ncdf4::ncvar_get(ncT, "longitude")
  
  dim.lat <- ncT$dim$latitude
  dim.lon <- ncT$dim$latitude
  
  # Creating a handy time dataframe
  dat.all <- data.frame(nc.time=nc.time, date=as.Date(paste0(yr.now, "-01-01"))+nc.time)
  dat.all$yday <- lubridate::yday(dat.all$date)
  dat.all$month <- lubridate::month(dat.all$date)
  
  # Loop through each variable and add it to the data frame
  for(v in vars.file){
    dat.all[,v] <- ncdf4::ncvar_get(ncT, v)
  }
  
  # Aggregate to daily; note: we need to do this regardless of whether or not we're saving daily data
  dat.day <- aggregate(dat.all[,c("nc.time", vars.file)], dat.all[,c("yday", "month")], FUN=mean)
  
  # Now lets do our special cases
  if("daylength" %in% add.vars){
    for(j in 1:nrow(dat.day)){
      d.now <- dat.day[j, "yday"]
      dat.day[j, "daylength"] <- nrow(dat.all[dat.all$yday==d.now & dat.all$surface_downwelling_shortwave_flux_in_air>0, ])
    }
  }
  if("air_temperature_maximum" %in% add.vars){
    dat.day[,"air_temperature_maximum"] <- aggregate(dat.all$air_temperature, list(dat.all$yday), FUN=max)[,2]
  }
  if("air_temperature_minimum" %in% add.vars){
    dat.day[,"air_temperature_minimum"] <- aggregate(dat.all$air_temperature, list(dat.all$yday), FUN=min)[,2]
  }
  # summary(dat.day)
  
  if(save.day==TRUE){
    dim.day <- ncdf4::ncdim_def(name = "time", units = paste0("days since ", yr.now, "-01-01T00:00:00Z"), 
                                vals = dat.day$nc.time, create_dimvar = TRUE, unlim = TRUE)
    
    nc.dim <- list(dim.lat, dim.lon, dim.day)
    
    # We already have most of the info stored, so lets 
    var.list <- list()
    for(j in vars.file){
      var.list[[j]] <- ncdf4::ncvar_def(name = j,
                                        units = ncT$var[[v]]$units, dim = nc.dim, missval = -9999,
                                        verbose = verbose)
    }
    
    # add our additional variables
    if("daylength" %in% add.vars){
      var.list[["daylength"]] <- ncdf4::ncvar_def(name = "daylength", units = "hours", dim = nc.dim, missval = -9999,
                                                  verbose = verbose)
    }
    if("air_temperature_maximum" %in% add.vars){
      var.list[["air_temperature_maximum"]] <- ncdf4::ncvar_def(name = "air_temperature_maximum", 
                                                                units = ncT$var$air_temperature$units, 
                                                                dim = nc.dim, missval = -9999, verbose = verbose)
    }
    if("air_temperature_minimum" %in% add.vars){
      var.list[["air_temperature_minimum"]] <- ncdf4::ncvar_def(name = "air_temperature_minimum", 
                                                                units = ncT$var$air_temperature$units, 
                                                                dim = nc.dim, missval = -9999, verbose = verbose)
    }
    
    
    loc <- ncdf4::nc_create(filename = file.path(out.base, day.dir, f.agg), vars = var.list, verbose = verbose)
    
    for (j in names(var.list)) {
      ncdf4::ncvar_put(nc = loc, varid = as.character(j), vals = dat.day[,j])
    }
    ncdf4::nc_close(loc)
    
  } # End daving days
  
  if(save.month==TRUE){
    dat.mo <- aggregate(dat.day[,c("nc.time", vars.file)], by=list(dat.day[,"month"]), FUN=mean)
    if(length(add.vars)>0){
      dat.mo[,add.vars] <- aggregate(dat.day[,add.vars], by=list(dat.day[,"month"]), FUN=mean)[,2:length(add.vars)]
    }
    
    dim.mo <- ncdf4::ncdim_def(name = "time", units = paste0("days since ", yr.now, "-01-01T00:00:00Z"), vals = dat.mo$nc.time, create_dimvar = TRUE, unlim = TRUE)
    
    nc.dim <- list(dim.lat, dim.lon, dim.mo)
    
    # We already have most of the info stored, so lets 
    var.list <- list()
    for(j in vars.file){
      var.list[[j]] <- ncdf4::ncvar_def(name = j,
                                        units = ncT$var[[v]]$units, dim = nc.dim, missval = -9999,
                                        verbose = verbose)
    }
    if("daylength" %in% add.vars){
      var.list[["daylength"]] <- ncdf4::ncvar_def(name = "daylength", units = "hours", dim = nc.dim, missval = -9999,
                                                  verbose = verbose)
    }
    if("air_temperature_maximum" %in% add.vars){
      var.list[["air_temperature_maximum"]] <- ncdf4::ncvar_def(name = "air_temperature_maximum", 
                                                                units = ncT$var$air_temperature$units, 
                                                                dim = nc.dim, missval = -9999, verbose = verbose)
    }
    if("air_temperature_minimum" %in% add.vars){
      var.list[["air_temperature_minimum"]] <- ncdf4::ncvar_def(name = "air_temperature_minimum", 
                                                                units = ncT$var$air_temperature$units, 
                                                                dim = nc.dim, missval = -9999, verbose = verbose)
    }
    
    
    loc <- ncdf4::nc_create(filename = file.path(out.base, mo.dir, f.agg), vars = var.list, verbose = verbose)
    
    for (j in names(var.list)) {
      ncdf4::ncvar_put(nc = loc, varid = as.character(j), vals = dat.mo[,j])
    }
    ncdf4::nc_close(loc)
    
  } # End month
  
  # Close file!
  ncdf4::nc_close(ncT)
  
  return(paste0("Processed: ", f.agg, " (Day=", save.day, "; Month=", save.month, ")"))
  
}

