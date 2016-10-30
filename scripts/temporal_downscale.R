predict.subdaily <- function(dat.mod, n.ens, path.model, lags.init, dat.train){
  # dat.mod    = data to be predicted at the time step of the training data
  # n.ens      = number of hourly ensemble members to generate
  # path.model = path to where the training model & betas is stored
  # lags.init  = a list with layers of same name dat.sim and n=n.ens that provide the initial lags
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

  # Set up the ensemble members in a list so the uncertainty can be propogated
  dat.sim <- list() 

  # ------------------------------------------
  # Modeling SWDOWN 
  # Note: this can be generalized to just run by DOY for all years at once since there's no memory in the system
  # ------------------------------------------
  {
    # Load the saved model
    load(file.path(path.model,"model_swdown.Rdata"))
    mod.swdown.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.swdown <- nc_open(file.path(path.model,"betas_swdown.nc"))
    n.betas <- nrow(ncvar_get(betas.swdown, "1"))
    
    dat.sim[["swdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      # For SWDOWN, we only want to model daylight hours -- make sure this matches what's in the swdown function
      day.now = unique(dat.mod[dat.mod$time.day==i, "doy"])
      
      # Use the training data to figure out night/day
      hrs.day = unique(dat.train[dat.train$doy==day.now & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
      
      rows.now = which(dat.mod$time.day==i)
      rows.mod = which(dat.mod$time.day==i & dat.mod$hour %in% hrs.day)
      
      dat.temp <- dat.mod[rows.mod,c("time.day", "year", "doy", "hour", "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")]
      dat.temp$swdown = 99999 # Dummy value so there's a column

      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.swdown, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      dimnames(Rbeta)[[2]] <- names(coef(mod.swdown.doy[[paste(day.now)]]))
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.swdown.doy[[paste(day.now)]], 
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
    }

    nc_close(betas.swdown)
    rm(mod.swdown.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling Temperature 
  # ------------------------------------------
  {
    # Load the saved model
    load(file.path(path.model,"model_tair.Rdata"))
    mod.tair.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.tair <- nc_open(file.path(path.model,"betas_tair.nc"))
    n.betas <- nrow(ncvar_get(betas.tair, "1"))
    
    dat.sim[["tair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$tair = -99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
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
      
      dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.tair, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      dimnames(Rbeta)[[2]] <- names(coef(mod.tair.doy[[paste(day.now)]]))
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.tair.doy[[paste(day.now)]], 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$tair), replace=T)

      for(j in 1:ncol(dat.sim$tair)){
        dat.sim[["tair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    }

    nc_close(betas.tair)
    rm(mod.tair.doy) # Clear out the model to save memory
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
    # Load the saved model
    load(file.path(path.model,"model_precipf.Rdata"))
    mod.precipf.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.precipf <- nc_open(file.path(path.model,"betas_precipf.nc"))
    n.betas <- nrow(ncvar_get(betas.precipf, "1"))
    
    dat.sim[["precipf"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))
    
    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      r
      ows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press")]
      dat.temp$precipf = 99999 # Dummy value so there's a column
      dat.temp$rain.prop = 99999 # Dummy value so there's a column
      
      day.now = unique(dat.temp$doy)
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.precipf, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
      # dimnames(Rbeta)[[2]] <- names(coef(mod.precipf.doy[[paste(day.now)]]))[1:ncol(Rbeta)]
      
      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.precipf.doy[[paste(day.now)]], 
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
      dat.pred <- dat.pred*(dat.temp$precipf.day)
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$precipf), replace=T)
      for(j in 1:ncol(dat.sim$precipf)){
        dat.sim[["precipf"]][rows.now,j] <- dat.pred[,cols.prop[j]]
      }
    }
    
    nc_close(betas.precipf)
    rm(mod.precipf.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  
  # ------------------------------------------
  # Modeling LWDOWN 
  # ------------------------------------------
  {
    # Load the saved model
    load(file.path(path.model,"model_lwdown.Rdata"))
    mod.lwdown.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.lwdown <- nc_open(file.path(path.model,"betas_lwdown.nc"))
    n.betas <- nrow(ncvar_get(betas.lwdown, "1"))
    
    dat.sim[["lwdown"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press")]
      dat.temp$lwdown = -99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$lwdown)
        names(sim.lag) <- c("lag.lwdown", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["lwdown"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$lwdown)))))
        names(sim.lag) <- c("lag.lwdown", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      dat.temp$swdown <- stack(dat.sim$swdown[rows.now,])[,1]
      dat.temp$tair <- stack(dat.sim$tair[rows.now,])[,1]
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.lwdown, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))

      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.lwdown.doy[[paste(day.now)]], 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)

      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$lwdown), replace=T)
      for(j in 1:ncol(dat.sim$lwdown)){
        dat.sim[["lwdown"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
      dat.mod[rows.now,"mod.lwdown"] <- apply(dat.sim[["lwdown"]][rows.now,],1, mean)
      if(i>min(dat.mod$time.day)){ 
        dat.mod[dat.mod$time.day==i-1,"lag.lwdown" ] <- dat.mod[dat.mod$time.day==i & dat.mod$hour==0,"mod.lwdown"] 
      }
    }
    
    nc_close(betas.lwdown)
    rm(mod.lwdown.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling PRESS 
  # ------------------------------------------
  {
    # Load the saved model
    load(file.path(path.model,"model_press.Rdata"))
    mod.press.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.press <- nc_open(file.path(path.model,"betas_press.nc"))
    n.betas <- nrow(ncvar_get(betas.press, "1"))
    
    dat.sim[["press"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press")]
      dat.temp$press = -99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$press)
        names(sim.lag) <- c("lag.press", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["press"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$press)))))
        names(sim.lag) <- c("lag.press", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.press, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))

      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.press.doy[[paste(day.now)]], 
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
      
    }

    nc_close(betas.press)
    rm(mod.press.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling QAIR 
  # ------------------------------------------
  {
    
    # Load the saved model
    load(file.path(path.model,"model_qair.Rdata"))
    mod.qair.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.qair <- nc_open(file.path(path.model,"betas_qair.nc"))
    n.betas <- nrow(ncvar_get(betas.qair, "1"))
    
    dat.sim[["qair"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$qair = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$qair)
        names(sim.lag) <- c("lag.qair", "ens")
        
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["qair"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$qair)))))
        names(sim.lag) <- c("lag.qair", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.qair, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))

      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.qair.doy[[paste(day.now)]], 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
      
      dat.pred <- exp(dat.pred) # because log-transformed
      
      # having a *very* hard time keeping humidity reasonable, even with transformations so truncating the max
      dat.pred[dat.pred>max(dat.train$qair)] <- max(dat.train$qair) 
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$qair), replace=T)
      for(j in 1:ncol(dat.sim$qair)){
        dat.sim[["qair"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
      
    }

    nc_close(betas.qair)
    rm(mod.qair.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  # ------------------------------------------
  # Modeling WIND 
  # ------------------------------------------
  {
    # Load the saved model
    load(file.path(path.model,"model_wind.Rdata"))
    mod.wind.doy <- mod.list
    rm(mod.list)
    
    # Load the meta info for the betas
    betas.wind <- nc_open(file.path(path.model,"betas_wind.nc"))
    n.betas <- nrow(ncvar_get(betas.wind, "1"))
    
    dat.sim[["wind"]] <- data.frame(array(dim=c(nrow(dat.mod), n.ens)))

    for(i in max(dat.mod$time.day):min(dat.mod$time.day)){
      rows.now = which(dat.mod$time.day==i)
      dat.temp <- dat.mod[rows.now,c("time.day", "year", "doy", "hour", 
                                     "tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day",
                                     "next.tmax", "next.tmin", "next.precipf", "next.swdown", "next.lwdown", "next.press", "next.qair", "next.wind")]
      dat.temp$wind = 99999 # Dummy value so there's a column
      day.now = unique(dat.temp$doy)
      
      # Set up the lags
      if(i==max(dat.mod$time.day)){
        sim.lag <- stack(lags.init$wind)
        names(sim.lag) <- c("lag.wind", "ens")
      } else {
        sim.lag <- stack(data.frame(array(dat.sim[["wind"]][dat.mod$time.day==(i+1)  & dat.mod$hour==0,], dim=c(1, ncol(dat.sim$wind)))))
        names(sim.lag) <- c("lag.wind", "ens")
      }
      dat.temp <- merge(dat.temp, sim.lag, all.x=T)
      
      rows.beta <- sample(1:n.betas, n.ens, replace=T)
      Rbeta <- as.matrix(ncvar_get(betas.wind, paste(day.now))[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))

      dat.pred <- predict.met(newdata=dat.temp, 
                              model.predict=mod.wind.doy[[paste(day.now)]], 
                              Rbeta=Rbeta, 
                              resid.err=F,
                              model.resid=NULL, 
                              Rbeta.resid=NULL, 
                              n.ens=n.ens)
       dat.pred <- exp(dat.pred) # because log-transformed
      
      # Randomly pick which values to save & propogate
      cols.prop <- sample(1:n.ens, ncol(dat.sim$wind), replace=T)
      for(j in 1:ncol(dat.sim$wind)){
        dat.sim[["wind"]][rows.now,j] <- dat.pred[dat.temp$ens==paste0("X", j),cols.prop[j]]
      }
    }

    nc_close(betas.wind)
    rm(mod.wind.doy) # Clear out the model to save memory
  }
  # ------------------------------------------
  
  return(dat.sim)

}