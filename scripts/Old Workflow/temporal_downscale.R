predict.subdaily <- function(dat.mod, n.ens, path.model, lags.list=NULL, lags.init=NULL, dat.train){
  # dat.mod    = data to be predicted at the time step of the training data
  # n.ens      = number of hourly ensemble members to generate
  # path.model = path to where the training model & betas is stored
  # lags.list  = a list with layers of same name dat.sim and n=n.ens that provide the initial lags; 
  #              used if entering the function from a parallel apply function
  # lags.init  = a data frame of initialization paramters to match the data in dat.mod
  # dat.train  = the training data used to fit the model; needed for night/day in swdown
  
  # --------------------------------
  # Models each variable separately at the moment rather than a generalized eqauiton that could potentially be parallilizeable
  # --------------------------------
  # #    2.1 tair (air temperature)
  # #    2.2 precipf (preciptiation, water equivalent)
  # #    2.3 swdown (downwelling shortwave radiation)
  # #    2.4 lwdown (downwelling longwave radiation)
  # #    2.5 press (surface pressure)
  # #    2.6 qair (specific humidity)
  # #    2.7 wind (surface wind speed)
  # --------------------------------# 
  # Load libraries
  library(ncdf4)
  library(mgcv)
  library(MASS)
  # library(lubridate)
  library(ggplot2)
  # library(tictoc)

  # Figure out if we need to extract the approrpiate 
  if(is.null(lags.init)){
    lags.init <- lags.list[[unique(dat.mod$ens.day)]]
  }
  
  # Set up the ensemble members in a list so the uncertainty can be propogated
  dat.sim <- list() 

  # DOY indexing is now off from the original; fix by subtracting 1
  # dat.mod$doy = dat.mod$doy-1
  # ------------------------------------------
  # Modeling SWDOWN 
  # Note: this can be generalized to just run by DOY for all years at once since there's no memory in the system
  # ------------------------------------------
  {
    dat.sim[["swdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      # For SWDOWN, we only want to model daylight hours -- make sure this matches what's in the swdown function
      day.now = unique(dat.mod[dat.mod$time.day==i, "doy"])-1
      
      # Use the training data to figure out night/day
      hrs.day = unique(dat.train[dat.train$doy==day.now & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
      
      rows.now = which(dat.mod$time.day==i)
      rows.mod = which(dat.mod$time.day==i & dat.mod$hour %in% hrs.day)
      
      dat.temp <- dat.mod[rows.mod,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      
      dat.temp$swdown = 99999 # Dummy value so there's a column
      # day.now = unique(dat.temp$doy)
      
      # Load the saved model
      load(file.path(path.model, "swdown", paste0("model_swdown_", day.now, ".Rdata")))
      mod.swdown.doy <- mod.save
      rm(mod.save)
    
      # Pull coefficients (betas) from our saved matrix
      betas.swdown <- nc_open(file.path(path.model, "swdown", paste0("betas_swdown_", day.now, ".nc")))
      n.beta <- betas.swdown$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.swdown$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.swdown, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.swdown)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.swdown.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      dat.pred[dat.pred<0] <- 0
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$swdown), replace=T)
      for(j in 1:ncol(dat.sim$swdown)){
        dat.sim[["swdown"]][rows.mod,j] <- dat.pred[,cols.prop[j]]
      }
      
      # For night time hours, value shoudl be 0
      dat.sim[["swdown"]][rows.now[!rows.now %in% rows.mod],] <- 0
    rm(mod.swdown.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling Temperature 
  # ------------------------------------------
  {
    dat.sim[["tair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$tair = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$tair)
        names(sim.lag) <- c("lag.tair", "ens")
        
        sim.lag$lag.tmin <- stack(lags.init$tmin)[,1]
        sim.lag$lag.tmax <- stack(lags.init$tmax)[,1]
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["tair"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$tair)))))
        names(sim.lag) <- c("lag.tair", "ens")
        sim.lag$lag.tmin <- stack(apply(dat.sim[["tair"]][dat.mod$time.day==(i+1),], 2, min))[,1]
        sim.lag$lag.tmax <- stack(apply(dat.sim[["tair"]][dat.mod$time.day==(i+1),], 2, max))[,1]
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      
      # Loading the saved model file
      load(file.path(path.model, "tair", paste0("model_tair_", day.now, ".Rdata")))
      mod.tair.doy <- mod.save
      rm(mod.save)

      # Pull coefficients (betas) from our saved matrix
      betas.tair <- nc_open(file.path(path.model, "tair", paste0("betas_tair_", day.now, ".nc")))
      n.beta <- betas.tair$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.tair$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.tair, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.tair)
     
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.tair.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$tair), replace=T)
      
      for(j in 1:ncol(dat.sim$tair)){
        dat.prop <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
        tmax.ens <- max(dat.temp[dat.temp$ens==paste0("X",j), "tmax.day"])
        tmin.ens <- min(dat.temp[dat.temp$ens==paste0("X",j), "tmin.day"])
        
        # Hard-coding in some bounds so we don't drift too far away from our given maxes & mins
        # Not going to worry for the moment about what happens if we undershoot out max/min since
        # it looks like we're normally decently close
        dat.prop[dat.prop>tmax.ens+2] <- tmax.ens+2
        dat.prop[dat.prop<tmin.ens-2] <- tmin.ens-2
        
        dat.sim[["tair"]][rows.now,j] <- dat.prop
      }
    rm(mod.tair.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  
  # ------------------------------------------
  # Modeling PRECIPF 
  # NOTE: For Precipf, we're doing this differently:
  #       We're basicallly modeling the proportion of that day's precipitation as a function of 
  #       the hour of day.  Because our beta fitting method leverages the covariance among betas,
  #       we'll end up with a daily sum that's pretty close to our daily total (90-110%)
  # ------------------------------------------
  {
    dat.sim[["precipf"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$precipf = 99999 # Dummy value so there's a column
      dat.temp$rain.prop = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$precipf)
        names(sim.lag) <- c("lag.precipf", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["precipf"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$precipf)))))
        names(sim.lag) <- c("lag.precipf", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      
      # Load the saved model
      load(file.path(path.model, "precipf", paste0("model_precipf_", day.now, ".Rdata")))
      mod.precipf.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      betas.precipf <- nc_open(file.path(path.model, "precipf", paste0("betas_precipf_", day.now, ".nc")))
      n.beta <- betas.precipf$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.precipf$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.precipf, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.precipf)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.precipf.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Re-distribute negative probabilities -- add randomly to make more peaky
      if(max(dat.pred)>0){ # If there's no rain on this day, skip the re-proportioning
        tmp <- 1:nrow(dat.pred) # A dummy vector of the 
        for(j in 1:ncol(dat.pred)){
          if(min(dat.pred[,j])>=0) next
          rows.neg <- which(dat.pred[,j]<0)
          rows.add <- sample(tmp[!tmp %in% rows.neg],length(rows.neg), replace=T)
          
          for(z in 1:length(rows.neg)){
            dat.pred[rows.add[z],j] <- dat.pred[rows.add[z],j] - dat.pred[rows.neg[z],j]
            dat.pred[rows.neg[z],j] <- 0  
          }
        }
        dat.pred <- dat.pred/rowSums(dat.pred)
        dat.pred[is.na(dat.pred)] <- 0
      }
      # Convert precip into real units
      dat.pred <- dat.pred*as.vector((dat.temp$precipf.day*length(unique(dat.temp$hour))))
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$precipf), replace=T)
      
      for(j in 1:ncol(dat.sim$precipf)){
        dat.sim[["precipf"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    rm(mod.precipf.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  
  # ------------------------------------------
  # Modeling LWDOWN 
  # ------------------------------------------
  {
    dat.sim[["lwdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$lwdown = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$lwdown)
        names(sim.lag) <- c("lag.lwdown", "ens")
        
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["lwdown"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$lwdown)))))
        names(sim.lag) <- c("lag.lwdown", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      dat.temp[,c("sqrt(lwdown)")] <- sqrt(dat.temp$lwdown)
      # dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      
      
      # Load the saved model
      load(file.path(path.model, "lwdown", paste0("model_lwdown_", day.now, ".Rdata")))
      mod.lwdown.doy <- mod.save
      rm(mod.save)
      
      
      # Pull the coefficients (betas) from our saved matrix
      betas.lwdown <- nc_open(file.path(path.model, "lwdown", paste0("betas_lwdown_", day.now, ".nc")))
      n.beta <- betas.lwdown$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.lwdown$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.lwdown, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.lwdown)

      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.lwdown.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      dat.pred <- dat.pred^2 # because squared to prevent negative numbers
      
      # Hard-coding some sanity bounds by ball-parking things from NLDAS & CRUNCEP
      # This is necessary if you have poorly constrained training models
      dat.pred[dat.pred<100] <- 100
      dat.pred[dat.pred>600] <- 600
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$lwdown), replace=T)
      
      for(j in 1:ncol(dat.sim$lwdown)){
        # test <- which(dat.temp$ens==paste0("X", j))
        dat.sim[["lwdown"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
    rm(mod.lwdown.doy) # Clear out the model to save memory
    }
    
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling PRESS 
  # ------------------------------------------
  {
    dat.sim[["press"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$press = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$press)
        names(sim.lag) <- c("lag.press", "ens")
        
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["press"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$press)))))
        names(sim.lag) <- c("lag.press", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      # Load the saved model
      load(file.path(path.model, "press", paste0("model_press_", day.now, ".Rdata")))
      mod.press.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      betas.press <- nc_open(file.path(path.model, "press", paste0("betas_press_", day.now, ".nc")))
      n.beta <- betas.press$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.press$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.press, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.press)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.press.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$press), replace=T)
      
      for(j in 1:ncol(dat.sim$press)){
        dat.sim[["press"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    rm(mod.press.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling QAIR 
  # ------------------------------------------
  {
    dat.sim[["qair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$qair = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$qair)
        names(sim.lag) <- c("lag.qair", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["qair"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$qair)))))
        names(sim.lag) <- c("lag.qair", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      dat.temp[,"log(qair)"] <- log(dat.temp$qair)
      
      # Load the saved model
      load(file.path(path.model, "qair", paste0("model_qair_", day.now, ".Rdata")))
      mod.qair.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      betas.qair <- nc_open(file.path(path.model, "qair", paste0("betas_qair_", day.now, ".nc")))
      n.beta <- betas.qair$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.qair$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.qair, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.qair)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.qair.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      dat.pred <- exp(dat.pred) # because log-transformed
      
      # qair sometimes ends up with high or infinite values, so lets make sure those get brought down a bit
      if(max(dat.pred)>0.03)  {
        qair.fix <- ifelse(quantile(dat.pred, 0.99)<0.03, quantile(dat.pred, 0.99), 0.03)
        dat.pred[dat.pred>qair.fix] <- qair.fix
      }

      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$qair), replace=T)
      
      for(j in 1:ncol(dat.sim$qair)){
        dat.sim[["qair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    rm(mod.qair.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling WIND 
  # ------------------------------------------
  {
    dat.sim[["wind"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$wind = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)-1
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$wind)
        names(sim.lag) <- c("lag.wind", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["wind"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$wind)))))
        names(sim.lag) <- c("lag.wind", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      dat.temp[,"sqrt(wind)"] <- sqrt(dat.temp$wind)
      
      # Load the saved model
      load(file.path(path.model, "wind", paste0("model_wind_", day.now, ".Rdata")))
      mod.wind.doy <- mod.save
      rm(mod.save)
      
      # Pull the coefficients (betas) from our saved matrix
      betas.wind <- nc_open(file.path(path.model, "wind", paste0("betas_wind_", day.now, ".nc")))
      n.beta <- betas.wind$var[[1]]$dim[[1]]$len # Number of rows
      col.beta <- betas.wind$var[[1]]$dim[[2]]$len # number of coefficients
      rows.beta <- sample(1:(n.beta-n.ens), 1, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.wind, paste(day.now), c(rows.beta,1), c(n.ens,col.beta)), nrow=length(rows.beta), ncol=col.beta)
      nc_close(betas.wind)
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.wind.doy, 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      dat.pred <- dat.pred^2 # because square-rooted to prevent negative numbers
      
      # Hard code an upper-level sanity check on the wind
      # 20 m/s = 45 mph; not hurricane strength, but plenty strong enough for most models
      # dat.pred[dat.pred>20] <- 20
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$wind), replace=T)
      
      for(j in 1:ncol(dat.sim$wind)){
        dat.sim[["wind"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    rm(mod.wind.doy) # Clear out the model to save memory
    }
  }
  # ------------------------------------------
  
  return(dat.sim)

}