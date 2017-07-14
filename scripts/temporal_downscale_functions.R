model.tair <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=1237){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)

  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    dat.subset$year <- as.ordered(dat.subset$year)  
    # day model works pretty good
      # Note: Tried fitting with hourly swdown & it didn't improve things (visually or with AIC), so for the sake of simplicity, using swdown.day
      # mod.doy <- lm(tair ~ as.ordered(hour)*tmax.day*(lag.tair + lag.tmin + tmin.day) +  as.ordered(hour)*tmin.day*next.tmax + as.ordered(hour)*swdown.day*(tmax.day + tmin.day) - 1 - as.ordered(hour) - swdown.day - lag.tair - lag.tmin - next.tmax - tmax.day - tmin.day - tmin.day*tmax.day - swdown.day*tmax.day*tmin.day, data=dat.subset) #
      mod.doy <- lm(tair ~ as.ordered(hour)*tmax.day*(lag.tair + lag.tmin + tmin.day) +  as.ordered(hour)*tmin.day*next.tmax - 1 - as.ordered(hour) - lag.tair - lag.tmin - next.tmax - tmax.day - tmin.day, data=dat.subset) #

      # If we can't estimate the covariance matrix, double our data and try again
      # NOTE: THIS IS NOT A GOOD PERMANENT FIX!!
      if(is.na(summary(mod.doy)$adj.r.squared)){
        warning(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy)))
        dat.subset <- rbind(dat.subset, dat.subset)
        mod.doy <- lm(tair ~ as.ordered(hour)*tmax.day*(lag.tair + lag.tmin + tmin.day) +  as.ordered(hour)*tmin.day*next.tmax - 1 - as.ordered(hour) - lag.tair - lag.tmin - next.tmax - tmax.day - tmin.day, data=dat.subset) #
      }
      
      # Generate a bunch of random coefficients that we can pull from 
      # without needing to do this step every day
      mod.coef <- coef(mod.doy)
      mod.cov  <- vcov(mod.doy)
      piv <- as.numeric(which(!is.na(mod.coef)))
      Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
      
      list.out <- list(model=mod.doy, 
                       betas=Rbeta)

      
      # Model residuals as a function of hour so we can increase our uncertainty
      if(resids==T){
        dat.subset[!is.na(dat.subset$lag.tair) & !is.na(dat.subset$next.tmax),"resid"] <- resid(mod.doy)
        resid.model <- lm(resid ~ as.ordered(hour)*(tmax.day*tmin.day)-1, data=dat.subset[!is.na(dat.subset$lag.tair),])
        res.coef <- coef(resid.model)
        res.cov  <- vcov(resid.model)
        res.piv <- as.numeric(which(!is.na(res.coef)))
        
        beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
        
        list.out[["model.resid"]] <- resid.model
        list.out[["betas.resid"]] <- beta.resid
      }
      return(list.out)
      
    }

  dat.list <- list()
  mod.out <- list()

  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_tair_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      # mod.save <- mod.out[[i]][["model"]]
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_tair_", i, ".Rdata")))
    }
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_tair_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      # mod.save <- mod.out$model
      # save(mod.save, file=file.path(path.out, paste0("model_tair_", i, ".Rdata")))
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      # mod.save$factors[mod.save$factors=="as.ordered(hour)"] <- "hour"
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_tair_", i, ".Rdata")))
   }
  }
  # return(mod.out)
}

model.swdown <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=1341){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  # The model we're going to use
  model.train <- function(dat.subset, threshold=NULL, n.beta, resids=resids){ 
    
    # Don't bother trying to fit hours that are completely or pretty darn close to dark
    hrs.day <- unique(dat.subset[dat.subset$swdown>threshold, "hour"])
    
    # Note: played around with a log-transformation of swdown to prevent negative values, but that resulted in bias at upper range
    # Solution was to just say anything <0 = 0
    # mod.doy <- lm(swdown ~ as.ordered(hour)*swdown.day, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    mod.doy <- lm(swdown ~ as.ordered(hour)*swdown.day-1 - swdown.day - as.ordered(hour), data=dat.subset[dat.subset$hour %in% hrs.day,]) ###

    # If we can't estimate the covariance matrix, double our data and try again
    # NOTE: THIS IS NOT A GOOD PERMANENT FIX!!
    if(is.na(summary(mod.doy)$adj.r.squared)){
      warning(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy)))
      dat.subset <- rbind(dat.subset, dat.subset)
      mod.doy <- lm(swdown ~ as.ordered(hour)*swdown.day-1 - swdown.day - as.ordered(hour), data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    }
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      dat.subset[dat.subset$hour %in% hrs.day,"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*swdown.day-1, data=dat.subset[dat.subset$hour %in% hrs.day,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
      
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids, threshold=quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05))
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_swdown_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_swdown_", i, ".Rdata")))
    }
    
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], threshold=quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_swdown_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_swdown_", i, ".Rdata")))
    }
  }
  
  # return(mod.out)
}

model.lwdown <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=341){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(lwdown ~ as.ordered(hour)*lwdown.day*(lag.lwdown + next.lwdown + swdown.day + tmax.day + tmin.day) - as.ordered(hour) - tmax.day - tmin.day - swdown.day - 1, data=dat.subset) ###
    mod.doy <- lm(sqrt(lwdown) ~ as.ordered(hour)*lwdown.day*(lag.lwdown + next.lwdown) - as.ordered(hour) - 1 - lag.lwdown - next.lwdown - lwdown.day - lwdown.day*lag.lwdown - lwdown.day*next.lwdown, data=dat.subset) ###
    
    # If we can't estimate the covariance matrix, stop & increase the moving window
    if(is.na(summary(mod.doy)$adj.r.squared)){
      stop(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy), ";  Increase day.window and try again"))
    }
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      dat.subset[!is.na(dat.subset$lag.lwdown) & !is.na(dat.subset$next.lwdown),"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*lwdown.day-1, data=dat.subset[,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_lwdown_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_lwdown_", i, ".Rdata")))
    }
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_lwdown_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      # mod.save$factors[mod.save$factors=="as.ordered(hour)"] <- "hour"
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_lwdown_", i, ".Rdata")))
      
    }
  }
  
  # return(mod.out)
}

model.press <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=1347){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(press ~ as.ordered(hour)*(press.day + lag.press + next.press)-as.ordered(hour)-1, data=dat.subset) ###
    mod.doy <- lm(press ~ as.ordered(hour)*(press.day + lag.press + next.press)-as.ordered(hour)-1-press.day - lag.press - next.press, data=dat.subset) ###

    # If we can't estimate the covariance matrix, double our data and try again
    # NOTE: THIS IS NOT A GOOD PERMANENT FIX!!
    if(is.na(summary(mod.doy)$adj.r.squared)){
      stop(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy), ";  Increase day.window and try again"))
      # dat.subset <- rbind(dat.subset, dat.subset)
      # mod.doy <- lm(press ~ as.ordered(hour)*(press.day + lag.press + next.press)-as.ordered(hour)-1-press.day - lag.press - next.press, data=dat.subset) ###
    }
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      dat.subset[!is.na(dat.subset$lag.press) & !is.na(dat.subset$next.press),"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*press.day-1, data=dat.subset[,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_press_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_press_", i, ".Rdata")))
    }
    
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_press_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_press_", i, ".Rdata")))
      
    }
  }
  
  # return(mod.out)
}

model.wind <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=708){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(log(wind) ~ as.ordered(hour)*log(wind.day)*(log(lag.wind) + log(next.wind) + press.day + tmin.day + tmax.day)-as.ordered(hour)-1 - press.day - tmin.day - tmax.day - log(wind.day)*press.day - log(wind.day)*tmin.day- log(wind.day)*tmax.day - as.ordered(hour)*tmin.day - as.ordered(hour)*tmax.day - as.ordered(hour)*press.day, data=dat.subset) ###
    # mod.doy <- lm(log(wind) ~ as.ordered(hour)*wind.day*(lag.wind + next.wind)-as.ordered(hour)-1 - wind.day - lag.wind - next.wind - wind.day*lag.wind - wind.day*next.wind, data=dat.subset) ###
    mod.doy <- lm(sqrt(wind) ~ as.ordered(hour)*wind.day*(lag.wind + next.wind)-as.ordered(hour)-1 - wind.day - lag.wind - next.wind - wind.day*lag.wind - wind.day*next.wind, data=dat.subset) ###
    
    # If we can't estimate the covariance matrix, stop & increase the moving window
    if(is.na(summary(mod.doy)$adj.r.squared)){
      stop(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy), ";  Increase day.window and try again"))
    }
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      dat.subset[!is.na(dat.subset$lag.wind) & !is.na(dat.subset$next.wind),"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*wind.day-1, data=dat.subset[,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_wind_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_wind_", i, ".Rdata")))
    }
    
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_wind_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_wind_", i, ".Rdata")))
      
    }
  }
  
  # return(mod.out)
}

model.precipf <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=1562){
  library(MASS)
  # library(fitdistrplus)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # Precip needs to be a bit different.  We're going to calculate the fraction of precip occuring in each hour
    # we're going to estimate the probability distribution of rain occuring in a given hour
    dat.subset$rain.prop <- dat.subset$precipf/(dat.subset$precipf.day*length(unique(dat.subset$hour)))
    mod.doy <- lm(rain.prop ~ as.ordered(hour)*precipf.day-1 - as.ordered(hour)-precipf.day, data=dat.subset)
    
    # If we can't estimate the covariance matrix, increase the moving window
    if(is.na(summary(mod.doy)$adj.r.squared)){
      stop(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy), ";  Increase day.window and try again"))
      # dat.subset <- rbind(dat.subset, dat.subset)
      # mod.doy <- lm(rain.prop ~ as.ordered(hour)*precipf.day-1 - as.ordered(hour)-precipf.day, data=dat.subset)
    }
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      # dat.subset[!is.na(dat.subset$lag.precipf) & !is.na(dat.subset$next.precipf),"resid"] <- resid(mod.doy)
      dat.subset[,"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*precipf.day-1, data=dat.subset[,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_precipf_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_precipf_", i, ".Rdata")))
    }
    
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_precipf_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_precipf_", i, ".Rdata")))
      
    }
  }
  
  # return(mod.out)
}

model.qair <- function(dat.train, n.beta=1000, path.out, resids=F, parallel=F, n.cores=NULL, day.window=5, seed=1009){
  library(MASS)
  set.seed(seed)
  if(!dir.exists(path.out)) dir.create(path.out, recursive=T)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(log(qair) ~ as.ordered(hour)*log(qair.day)*(log(lag.qair) + log(next.qair) + precipf.day + tmin.day + tmax.day)-as.ordered(hour)-1 - precipf.day - tmin.day - tmax.day - log(qair.day)*precipf.day - log(qair.day)*tmin.day- log(qair.day)*tmax.day - as.ordered(hour)*tmin.day - as.ordered(hour)*tmax.day - as.ordered(hour)*precipf.day, data=dat.subset) ###
    # mod.doy <- lm(log(qair) ~ as.ordered(hour)*qair.day*(lag.qair + next.qair)-as.ordered(hour)-1 - qair.day - lag.qair - next.qair - qair.day*lag.qair - qair.day*next.qair, data=dat.subset) ###
    mod.doy <- lm(log(qair) ~ as.ordered(hour)*qair.day*(lag.qair + next.qair + tmax.day)-as.ordered(hour)-1 - tmax.day, data=dat.subset) ###
    # mod.doy <- glm(qair ~ as.ordered(hour)*qair.day*(lag.qair + next.qair)-as.ordered(hour)-1 - qair.day - lag.qair - next.qair - qair.day*lag.qair - qair.day*next.qair, data=dat.subset, family="quasibinomial") ###
    
    # If we can't estimate the covariance matrix, stop & increasing the moving window
    if(is.na(summary(mod.doy)$adj.r.squared)){
      stop(paste0("Can not estimate covariance matrix for day of year: ", unique(dat.subset$doy), ";  Increase day.window and try again"))
      # dat.subset <- rbind(dat.subset, dat.subset)
      # mod.doy <- lm(log(qair) ~ as.ordered(hour)*qair.day*(lag.qair + next.qair + tmax.day)-as.ordered(hour)-1 - tmax.day, data=dat.subset) ###
    }
    
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    list.out <- list(model=mod.doy, 
                     betas=Rbeta)
    # Model residuals as a function of hour so we can increase our uncertainty
    if(resids==T){
      dat.subset[!is.na(dat.subset$lag.qair) & !is.na(dat.subset$next.qair),"resid"] <- resid(mod.doy)
      resid.model <- lm(resid ~ as.ordered(hour)*qair.day-1, data=dat.subset[,])
      res.coef <- coef(resid.model)
      res.cov  <- vcov(resid.model)
      res.piv <- as.numeric(which(!is.na(res.coef)))
      
      beta.resid <- mvrnorm(n=n.beta, res.coef[res.piv], res.cov)
      
      list.out[["model.resid"]] <- resid.model
      list.out[["betas.resid"]] <- beta.resid
    }
    
    return(list.out)
  }
  
  dat.list <- list()
  mod.out <- list()
  
  # Make the data into a list
  # Training the model on ax X-day window around the actual DOY we're trying to model
  # this helps avoid problems with lack of data in small datasets like Ameriflux
  # Default window is 5 days (+/- 2)
  for(i in unique(dat.train$doy)){
    if(i >= 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=365-day.window/2 | dat.train$doy<=day.window/2,]
    } else if(i == 1){
      dat.list[[paste(i)]] <- dat.train[dat.train$doy<=i+day.window/2 | dat.train$doy>=365-day.window/2,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=i-day.window/2 & dat.train$doy<=i+day.window/2,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    warning("Running model calculation in parallel.  This WILL crash if you do not have access to a LOT of memory!")
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
    
    # Use a loop to sace each day of year independently
    for(i in names(mod.out)){
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_qair_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[[i]][["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[[i]][["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[[i]][["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out[[i]]$model$call
      mod.save$coef  <- coef(mod.out[[i]]$model)
      mod.save$formula <- parse(text=mod.out[[i]]$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out[[i]]$model$terms, "factors"))
      mod.save$xlev  <- mod.out[[i]]$model$xlevels
      mod.save$contr <- mod.out[[i]]$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_qair_", i, ".Rdata")))
    }
    
  } else {
    for(i in names(dat.list)){
      mod.out <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
      
      # Save the betas as .nc
      outfile=file.path(path.out, paste0("betas_qair_", i, ".nc"))
      dimY <- ncdim_def( paste0("coeffs_", i), units="unitless", longname="model.out coefficients", vals=1:ncol(mod.out[["betas"]]))
      dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(mod.out[["betas"]]))
      var.list <- ncvar_def(i, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", i, " model.out coefficients"))
      nc <- nc_create(outfile, var.list)
      ncvar_put(nc, var.list, mod.out[["betas"]])
      nc_close(nc)
      
      # Save the model as a .Rdata
      mod.save <- list()
      mod.save$call  <- mod.out$model$call
      mod.save$coef  <- coef(mod.out$model)
      mod.save$formula <- parse(text=mod.out$model$call[[2]][c(1,3)])
      mod.save$factors  <- rownames(attr(mod.out$model$terms, "factors"))
      mod.save$xlev  <- mod.out$model$xlevels
      mod.save$contr <- mod.out$model$contrasts
      save(mod.save, file=file.path(path.out, paste0("model_qair_", i, ".Rdata")))
      
    }
  }
  
  # return(mod.out)
}

predict.met <- function(newdata, model.predict, Rbeta, resid.err=F, model.resid=NULL, Rbeta.resid=NULL, n.ens){
  # set.seed(9321)
  err.resid = 0 # dummy residual error term; if we want to add residual error, we're modeling it by hour

  # mod.terms <- terms(model.predict)
  # mod.coef <- coef(model.predict)
  # mod.cov  <- vcov(model.predict)
  # mod.resid <- resid(model.predict)
  piv <- as.numeric(which(!is.na(model.predict$coef)))
  
  # m <- model.frame(mod.terms, newdata, xlev = model.predict$xlevels)
  # Xp <- model.matrix(mod.terms, m, contrasts.arg = model.predict$contrasts)
  df.hr <- data.frame(hour = model.predict$xlev[[1]])
  df.hr[,"as.ordered(hour)"] <- df.hr$hour
  
  model.predict$factors[model.predict$factors=="as.ordered(hour)"] <- "hour"
  m  <- newdata[complete.cases(newdata[,model.predict$factors]),model.predict$factors]
  m[,"as.ordered(hour)"] <- m$hour
  if(length(df.hr$hour)!= length(m$hour)) m <- merge(m, df.hr, all=T)
  
  Xp <-  model.matrix(eval(model.predict$formula), m, contrasts.arg=model.predict$contr)
  
  if(resid.err==T){
    newdata$resid <- 99999
    resid.terms <- terms(model.resid)
    resid.coef <- coef(model.resid)
    resid.cov  <- vcov(model.resid)
    resid.resid <- resid(model.resid)
    resid.piv <- as.numeric(which(!is.na(resid.coef)))
    
    m2 <- model.frame(resid.terms, newdata, xlev = model.resid$xlevels)
    Xp.res <- model.matrix(resid.terms, m2, contrasts.arg = model.resid$contrasts)
    
    err.resid <- Xp.res[,resid.piv] %*% t(Rbeta.resid)
  }
  
  dat.sim <- Xp[,piv] %*% t(Rbeta) + err.resid
  
  return(dat.sim)
  
}

save.betas <- function(model.out, betas, outfile){
  # Function to save betas as a .nc file
  # model.out = the model output list
  # betas = the name of the layer of betas to save (e.g. "betas", or "betas.resid")
  
  # Save a vector of the dimnames
  # names.coefs <- dimnames(model.out[[1]][[betas]])[[2]]
  
  var.list <- list()
  for(v in names(model.out)){
    # Note: Need a separate list of coefficients for each variable to make my life easier if the is swdown which has varying
    #       predictors by day
    dimY <- ncdim_def( paste0("coeffs_", v), units="unitless", longname="model.out coefficients", vals=1:ncol(model.out[[v]][[betas]]))
    dimX <- ncdim_def( "random", units="unitless", longname="random betas", vals=1:nrow(model.out[[v]][[betas]]))
    
    var.list[[v]] <- ncvar_def(v, units="coefficients", dim=list(dimX, dimY), longname=paste0("day ", v, " model.out coefficients"))
  }
  
  # dim.string <- ncdim_def("names", "", 1:max(nchar(names.coefs)), create_dimvar=FALSE)
  # var.list[["names.coefs"]] <- ncvar_def("names.coefs", units="", dim=list(dim.string, dimY), longname="model coefficient names", prec="char")
  
  nc <- nc_create(outfile, var.list)
  # ncvar_put(nc, var.list$names.coefs, names.coefs)
  for(v in names(model.out)) {
    ncvar_put(nc, var.list[[v]], model.out[[v]][[betas]])
  }
  nc_close(nc)    
  
  
} 

save.model <- function(model.out, model, outfile){
  # Function to save linear models as a .Rdata file
  # model.out = the model output list
  # model = the name of the layer of betas to save (e.g. "model", or "model.resid")
  
  mod.list <- list()
  for(v in names(model.out)){
    mod.list[[v]] <- model.out[[v]][[model]]
  }
  
  save(mod.list, file=outfile)
} 

graph.resids <- function(var, dat.train, model.var, fig.dir){
  dat.train$var     <- dat.train[,var]
  dat.train$lag.var <- dat.train[,paste0("lag.", var)]
  if(var=="tair"){
    dat.train$next.var <- dat.train[,paste0("next.", "tmax")]
  } else {
    dat.train$var.day <- dat.train[,paste0(var, ".day")]
    dat.train$next.var <- dat.train[,paste0("next.", var)]
  }
  
  if(var == "precipf"){
    for(i in names(model.var)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      if(length(dat.train[dat.train$doy==as.numeric(i) & dat.train$var.day>0, "resid"])==0) next # Skip the days with no rain!
      dat.train[dat.train$doy==as.numeric(i) & dat.train$var.day>0, "predict"] <- predict(model.var[[i]]$model, newdata=dat.train[dat.train$doy==as.numeric(i) & dat.train$var.day>0, ])
      # dat.train[dat.train$doy==as.numeric(i) & dat.train$var.day>0, "resid"] <- resid(model.var[[i]]$model)
    }
  } else if(var=="swdown"){ 
    for(i in names(model.var)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      hrs.day = unique(dat.train[dat.train$doy==as.numeric(i) & dat.train$swdown>quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), "hour"])
      dat.train[dat.train$doy==as.numeric(i) & dat.train$hour %in% hrs.day, "predict"] <- predict(mod.swdown.doy[[i]]$model, newdata=dat.train[dat.train$doy==as.numeric(i) & dat.train$hour %in% hrs.day, ])
      # dat.train[dat.train$doy==as.numeric(i) & dat.train$hour %in% hrs.day, "resid"] <- resid(mod.swdown.doy[[i]]$model)
      
      dat.train[dat.train$doy==as.numeric(i) & !(dat.train$hour %in% hrs.day), "predict"] <- 0
    }
  } else {
    for(i in names(model.var)){
      if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
      dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.var) & !is.na(dat.train$next.var), "predict"] <- predict(model.var[[i]]$model, newdata=dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.var) & !is.na(dat.train$next.var), ])
      # dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.var) & !is.na(dat.train$next.var), "resid"] <- resid(model.var[[i]]$model)
    }
  }
  if(var %in% c("qair")){
    dat.train$predict <- exp(dat.train$predict)
    if(var=="qair"){
      dat.train$predict[dat.train$predict>max(dat.train[,var])] <- max(dat.train[,var])
    }
  }
  if(var %in% c("wind", "lwdown")){
    dat.train$predict <- (dat.train$predict)^2
  }
  dat.train$resid <- dat.train[,var] - dat.train$predict
  
  # summary(dat.train)
  
  png(file.path(fig.dir, paste0(var, "_Resid_vs_Hour.png")), height=8, width=8, units="in", res=180)
  plot(resid ~ hour, data=dat.train, cex=0.5); abline(h=0, col="red")
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_Resid_vs_DOY.png")), height=8, width=8, units="in", res=180)
  plot(resid ~ doy, data=dat.train, cex=0.5); abline(h=0, col="red") # slightly better in summer, but no clear temporal over-dispersion
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_Resid_vs_Predict.png")), height=8, width=8, units="in", res=180)
  plot(resid ~ predict, data=dat.train); abline(h=0, col="red")
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_Resid_vs_Obs.png")), height=8, width=8, units="in", res=180)
  plot(resid ~ var, data=dat.train, cex=0.5); abline(h=0, col="red")
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_Predict_vs_Obs.png")), height=8, width=8, units="in", res=180)
  plot(predict ~ var, data=dat.train, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()
  
  # Looking at the daily maxes & mins
  day.stats <- aggregate(dat.train[,c("var", "resid", "predict")],
                         by=dat.train[,c("time.day2", "year", "doy")],
                         FUN=mean, na.rm=T)
  day.stats$mod.max <- aggregate(dat.train[,c("predict")],
                                 by=dat.train[,c("time.day2", "year", "doy")],
                                 FUN=max)[,"x"]
  day.stats$mod.min <- aggregate(dat.train[,c("predict")],
                                 by=dat.train[,c("time.day2", "year", "doy")],
                                 FUN=min)[,"x"]
  
  png(file.path(fig.dir, paste0(var, "_Predict_vs_Mean.png")), height=8, width=8, units="in", res=180)
  plot(predict~ var, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_ResidualHistograms_Day.png")), height=6, width=8, units="in", res=180)
  # par(mfrow=c(3,1))
  par(mfrow=c(1,1))
  hist(day.stats$resid, main="Daily Mean Residuals")
  dev.off()
}

graph.predict <- function(dat.mod, dat.ens, var, yr, fig.dir){
  # ---------
  # Graph the output
  # ---------
  {

    date.vec <- dat.mod[dat.mod$ens.day==unique(dat.mod$ens.day)[1],"date"]
    dat.mod <- aggregate(dat.mod[,c("tmax.day", "tmin.day", "precipf.day", "swdown.day", "lwdown.day", "press.day", "qair.day", "wind.day")],
                         by=dat.mod[,c("year", "doy", "time.day", "time.hr")],
                         FUN=mean)
    dat.mod$date <- date.vec
    head(dat.mod)
    
    dat.mod$var.pred <- apply(dat.ens[[var]],1,mean)
    dat.mod$var.025 <- apply(dat.ens[[var]],1,quantile, 0.025)
    dat.mod$var.975 <- apply(dat.ens[[var]],1,quantile, 0.975)
    
    
    if(var=="tair"){
      
      # If this is temperature, add lines fo the daily max and min
      png(file.path(fig.dir, paste0(var, yr, "_year.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==yr,]) +
          geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=var.pred), color="blue") +
          geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
          # geom_point(aes(x=date, y=tmax.day), color="black", size=0.1, alpha=0.5) +
          # geom_point(aes(x=date, y=tmin.day), color="black", size=0.1, alpha=0.5) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(paste0(var, ": ", yr)) +
          theme_bw()
      )
      dev.off()
      
      dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
      dat.graph1$season <- as.factor("winter")
      dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
      dat.graph2$season <- as.factor("spring")
      dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
      dat.graph3$season <- as.factor("summer")
      dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
      dat.graph4$season <- as.factor("fall")
      
      dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
      
      png(file.path(fig.dir, paste(var, yr,"examples.png", sep="_")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==yr,]) +
          facet_wrap(~season, scales="free") +
          # geom_point(aes(x=date, y=tmax.day), color="black", size=0.1, alpha=0.5) +
          # geom_point(aes(x=date, y=tmin.day), color="black", size=0.1, alpha=0.5) +
          geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=var.pred), color="blue") +
          geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
          scale_y_continuous(name=var) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(paste0(var, ": ", yr)) +
          theme_bw()
      )
      dev.off()
      
    } else {
      png(file.path(fig.dir, paste(var, yr, "year.png", sep="_")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==yr,]) +
          geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=var.pred), color="blue") +
          geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
          # geom_line(aes(x=date, y=swdown), color="black", alpha=0.5) +
          # geom_point(aes(x=date, y=swdown), color="black", size=0.3, alpha=0.5) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(paste0(var, ": ", yr)) +
          theme_bw()
      )
      dev.off()
      
      dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
      dat.graph1$season <- as.factor("winter")
      dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
      dat.graph2$season <- as.factor("spring")
      dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
      dat.graph3$season <- as.factor("summer")
      dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
      dat.graph4$season <- as.factor("fall")
      
      dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
      
      png(file.path(fig.dir, paste0(var, yr,"_examples.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==yr,]) +
          facet_wrap(~season, scales="free") +
          # geom_line(aes(x=date, y=swdown), color="black") +
          # geom_point(aes(x=date, y=swdown), color="black", size=0.5) +
          geom_ribbon(aes(x=date, ymin=var.025, ymax=var.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=var.pred), color="blue") +
          geom_point(aes(x=date, y=var.pred), color="blue", size=0.5) +
          scale_y_continuous(name=var) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(paste0(var, ": ", yr)) +
          theme_bw()
      )
      dev.off()
    }
  }
  # ---------
  
}