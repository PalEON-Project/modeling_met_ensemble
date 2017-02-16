Ameriflux2traindata <- function(filepath, Ameriflux, site.name, outfile, start_date, end_date,
                                upscale=FALSE, CF.names=TRUE,
                                overwrite = FALSE, verbose = FALSE, ...) {
  # Upscale can either be set for FALSE (leave alone) or to the temporal resolution you want to aggregate to
  # options are: year, doy (day of year), or hour
  library(car)
  
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date <- as.POSIXlt(end_date, tz = "UTC")
  
  start_year <- lubridate::year(start_date)
  end_year <- lubridate::year(end_date)
  
  yr_seq <- seq(start_year,end_year)
  
  # make sure output folder exists
  # if (!file.exists(outfolder)) {
  #   dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
  # }
  #------------ Now let's read in the data
  
  input_met = list()
  for (j in seq_along(yr_seq)){
    input_met[[j]] = file.path(filepath,paste0(Ameriflux,".",yr_seq[j],".nc"))
  }
  
  
  var <- data.frame(CF.name = c("air_temperature", "air_temperature_max", "air_temperature_min", 
                                "surface_downwelling_longwave_flux_in_air", "air_pressure", "surface_downwelling_shortwave_flux_in_air", 
                                "eastward_wind", "northward_wind", "specific_humidity", "precipitation_flux"), 
                    units = c("Kelvin", "Kelvin", "Kelvin", "W/m2", "Pascal", "W/m2", "m/s", 
                              "m/s", "g/g", "kg/m2/s"))
  
  step = list() # list of time steps just in case they vary
  
  #Have to do 1 go first to initialize the dataframe 
  raw_train_data <- list()
  tem <- ncdf4::nc_open(input_met[[1]])
  dim <- tem$dim
  for (j in seq_along(var$CF.name)) {
    if (exists(as.character(var$CF.name[j]), tem$var) == FALSE) {
      raw_train_data[[j]] = NA
    } else {
      raw_train_data[[j]] = ncdf4::ncvar_get(tem, as.character(var$CF.name[j]))
    }
  }
  names(raw_train_data)<- var$CF.name 
  train_df = data.frame(raw_train_data)
  # colnames(train_df) <- var$CF.name 
  # train_df$year <- yr_seq[1]
  
  raw_train_data = data.frame(raw_train_data)
  if ((nrow(raw_train_data) == 17520) | (nrow(raw_train_data) == 17568)) {
    step[1] = 2
  }
  if ((nrow(raw_train_data) == 8760 ) | (nrow(raw_train_data) == 8784 ) ) {
    step[1] = 1
  }
  
  # Add a time stamp
  start_time <- as.POSIXlt(paste0(yr_seq[1], "-01-01"), tz = "UTC")
  end_time <- as.POSIXlt(paste0(yr_seq[1]+1, "-01-01"), tz = "UTC")
  train_df$date <- seq.POSIXt(from=start_time, by=60*60/step[[1]], length.out=nrow(train_df))
  summary(train_df)
  
  lat_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "latitude"))
  lon_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "longitude"))
  ncdf4::nc_close(tem)
  
  # Loop through remaining files and format
  if(length(input_met)>1){ # Just a safe guard in case we run this for a single year
    for (i in 2:length(input_met)) {
      if (file.exists(input_met[[i]])) {
        raw_train_data <- list()
        tem <- ncdf4::nc_open(input_met[[i]])
        dim <- tem$dim
        for (j in seq_along(var$CF.name)) {
          if (exists(as.character(var$CF.name[j]), tem$var) == FALSE) {
            raw_train_data[[j]] = NA
          } else {
            raw_train_data[[j]] = ncdf4::ncvar_get(tem, as.character(var$CF.name[j]))
          }
        }
        names(raw_train_data)<- var$CF.name 
        raw_train_data = data.frame(raw_train_data)
        # colnames(raw_train_data) <- var$CF.name 
        # raw_train_data$year = yr_seq[i]
        
        #assign(paste0("train_df", i),raw_train_data)
        if ((nrow(raw_train_data) == 17520) | (nrow(raw_train_data) == 17568)) {
          step[i] = 2
        }
        if ((nrow(raw_train_data) == 8760 ) | (nrow(raw_train_data) == 8784 ) ) {
          step[i] = 1
        }
        
        # Add a time stamp
        start_time <- as.POSIXlt(paste0(yr_seq[i], "-01-01"), tz = "UTC")
        end_time <- as.POSIXlt(paste0(yr_seq[i]+1, "-01-01"), tz = "UTC")
        raw_train_data$date <- seq.POSIXt(from=start_time, by=60*60/step[[i]], length.out=nrow(raw_train_data))
        summary(raw_train_data)
        
        train_df = rbind(train_df,raw_train_data)
      } else {
        step[i] = NA
      }
      ncdf4::nc_close(tem)
    } # End year loop
  }
  
  # Quick & dirty way of gap-filling any lingering NAs
  for(i in 1:ncol(train_df)) {
    train_df[is.na(train_df[,i]), i] <- mean(train_df[,i], na.rm = TRUE)
  }

  # Getting additional time stamps
  train_df$year <- year(train_df$date)
  train_df$doy <- yday(train_df$date)
  train_df$hour <- hour(train_df$date)
  
  # ---------------------------------
  # Temporal upscaling
  # ---------------------------------
  # year = unique(train_df$year)
  # reso_len = list()
  # sp = list()
  # for (i in 1:length(year)) {
  #   if (lubridate::leap_year(year[i]) == TRUE) {
  #     sp[i] <- 366
  #   } else {
  #     sp[i] <- 365
  #   }
  #   reso_len[i] <- as.numeric(sp[i]) * 24
  # }
  
  # upscale_data <- data.frame()
  # for (n in seq_along(colnames(train_df))) {
  #   for (i in seq_along(step)) {
  #     if (is.na(step[[i]]) == FALSE) {
  #       upscale_data[1:sum(as.numeric(reso_len)),n] = colMeans(matrix(train_df[[n]], nrow=step[[i]]))
  #     }
  #   }
  # }
  if(upscale==FALSE){
    dat.train=train_df
  } else {
    # Figure out which temporal variables we're aggregating over
    time.vars <- c("year", "doy", "hour")
    agg.ind <- which(time.vars==upscale)
    time.vars <- time.vars[1:agg.ind]
    
    dat.train <- aggregate(train_df[,names(train_df)[!names(train_df) %in% c("year", "doy", "hour")]],
                              by=train_df[,time.vars],
                              FUN=mean, na.rm=F)
    dat.train <- dat.train[order(dat.train$date),]
  }
  # ---------------------------------
  
  
  # Add in wind speed & dataset name
  dat.train$wind_speed <- sqrt(dat.train$eastward_wind^2+dat.train$northward_wind^2)
  dat.train$dataset <- paste0("Ameriflux-", site.name)
  
  # creating the columns required for downscaling functions
  if(CF.names==TRUE){
    # Just reorder the columns
    dat.train <- dat.train[,c("date", "year", "doy","hour", "air_temperature", "precipitation_flux", "air_temperature_max", "air_temperature_min",
                              "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air","air_pressure",
                              "specific_humidity", "eastward_wind", "northward_wind", "wind_speed")]
  } else {
    # Re-order the column, plus use the codes Christy uses
    dat.train <- dat.train[,c("date", "year", "doy","hour", "air_temperature", "precipitation_flux", "air_temperature_max", "air_temperature_min",
                              "surface_downwelling_shortwave_flux_in_air", "surface_downwelling_longwave_flux_in_air","air_pressure",
                              "specific_humidity", "eastward_wind", "northward_wind", "wind_speed")]
    names(dat.train) <- recode(names(dat.train), "'air_temperature'='tair';
                                                  'precipitation_flux'='precipf';
                                                  'air_temperature_max'='tmax';
                                                  'air_temperature_min'='tmin';
                                                  'surface_downwelling_shortwave_flux_in_air'='swdown';
                                                  'surface_downwelling_longwave_flux_in_air'='lwdown';
                                                  'air_pressure'='press';
                                                  'specific_humidity'='qair';
                                                  'eastward_wind'='uas';
                                                  'northward_wind'='vas';
                                                  'wind_speed'='wind'"
                               )
    
  }
  
  rm(raw_train_data)
  rm(train_df)
  write.csv(x = dat.train,file = outfile, row.names = FALSE)
  
}
