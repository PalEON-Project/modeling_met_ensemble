Ameriflux2traindata <- function(filepath,sitename, outfolder, start_date, end_date,
                                overwrite = FALSE, verbose = FALSE, ...) {
  
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date <- as.POSIXlt(end_date, tz = "UTC")
  
  start_year <- lubridate::year(start_date)
  end_year <- lubridate::year(end_date)
  
  yr_seq <- seq(start_year,end_year)
  
  # make sure output folder exists
  if (!file.exists(outfolder)) {
    dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
  }
  #------------ Now let's read in the data
  
  input_met = list()
  for (j in seq_along(yr_seq)){
    input_met[[j]] = paste0(filepath,sitename,".",yr_seq[j],".nc")
  }
  
  
  var <- data.frame(CF.name = c("air_temperature", "air_temperature_max", "air_temperature_min", 
                                "surface_downwelling_longwave_flux_in_air", "air_pressure", "surface_downwelling_shortwave_flux_in_air", 
                                "eastward_wind", "northward_wind", "specific_humidity", "precipitation_flux"), 
                    units = c("Kelvin", "Kelvin", "Kelvin", "W/m2", "Pascal", "W/m2", "m/s", 
                              "m/s", "g/g", "kg/m2/s"))
  
  step = list()
  raw_train_data <- list()
  #Have to do 1 go first to initialize the dataframe 
  tem <- ncdf4::nc_open(input_met[[1]])
  dim <- tem$dim
  for (j in seq_along(var$CF.name)) {
    if (exists(as.character(var$CF.name[j]), tem$var) == FALSE) {
      raw_train_data[[j]] = NA
    } else {
      raw_train_data[[j]] = ncdf4::ncvar_get(tem, as.character(var$CF.name[j]))
    }
  }
  train_df = data.frame(raw_train_data)
  colnames(train_df) <- var$CF.name 
  train_df$year <- yr_seq[1]
  raw_train_data = data.frame(raw_train_data)
  if ((nrow(raw_train_data) == 17520) | (nrow(raw_train_data) == 17568)) {
    step[1] = 2
  }
  if ((nrow(raw_train_data) == 8760 ) | (nrow(raw_train_data) == 8784 ) ) {
    step[1] = 1
  }
  lat_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "latitude"))
  lon_raw_train_data <- as.numeric(ncdf4::ncvar_get(tem, "longitude"))
  ncdf4::nc_close(tem)
  
  
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
      raw_train_data = data.frame(raw_train_data)
      colnames(raw_train_data) <- var$CF.name 
      raw_train_data$year = yr_seq[i]
      train_df = rbind(train_df,raw_train_data)
      #assign(paste0("train_df", i),raw_train_data)
      if ((nrow(raw_train_data) == 17520) | (nrow(raw_train_data) == 17568)) {
        step[i] = 2
      }
      if ((nrow(raw_train_data) == 8760 ) | (nrow(raw_train_data) == 8784 ) ) {
        step[i] = 1
      }
    } else {
      step[i] = NA
    }
  }
  ncdf4::nc_close(tem)
  
  year = unique(train_df$year)
  reso_len = list()
  sp = list()
  for (i in 1:length(year)) {
    if (lubridate::leap_year(year[i]) == TRUE) {
      sp[i] <- 366
    } else {
      sp[i] <- 365
    }
    reso_len[i] <- as.numeric(sp[i]) * 24
  }
  
  for(i in 1:ncol(train_df)) {
    train_df[is.na(train_df[,i]), i] <- mean(train_df[,i], na.rm = TRUE)
  }
  
  
  upscale_data <- data.frame()
  for (n in seq_along(colnames(train_df))) {
    for (i in seq_along(step)) {
      if (is.na(step[[i]]) == FALSE) {
        upscale_data[1:sum(as.numeric(reso_len)),n] = colMeans(matrix(train_df[[n]], nrow=step[[i]]))
      }
    }
  }
  
  colnames(upscale_data) = colnames(train_df)
  wind_speed = sqrt(upscale_data$eastward_wind^2+upscale_data$northward_wind^2)
  
  # creating the columns required for downscaling functions
  dat.train = data.frame(rep("US-PFa", nrow(upscale_data)))
  colnames(dat.train) = "dataset"
  dat.train$year = upscale_data$year
  
  doy = vector()
  for (i in seq_along(sp)) {
    a = vector()
    a = floor(seq(0,as.numeric(sp[i])-1/24,1/24))
    doy = append(doy,a)
  }
  
  dat.train$doy = floor(doy)
  
  hour = seq(0,23,1)
  dt_hr = vector()
  for ( i in seq_along(sp)) {
    hr_tem = vector()
    hr_tem = rep(hour,sp[[i]])
    dt_hr = append(dt_hr,hr_tem)
  }
  dat.train$hour = dt_hr
  dat.train$tair = upscale_data$air_temperature
  dat.train$precipf = upscale_data$precipitation_flux
  dat.train$swdown = upscale_data$surface_downwelling_shortwave_flux_in_air
  dat.train$lwdown = upscale_data$surface_downwelling_longwave_flux_in_air
  dat.train$press = upscale_data$air_pressure
  dat.train$qair = upscale_data$specific_humidity
  dat.train$uas = upscale_data$eastward_wind
  dat.train$vas = upscale_data$northward_wind
  dat.train$wind = wind_speed
  
  rm(raw_train_data)
  rm(upscale_data)
  write.csv(x = dat.train,file = file.path(outfolder,paste0(sitename,"_training_data")), row.names = FALSE)
}
