model.tair <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=1237){
  library(MASS)
  set.seed(seed)

  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    dat.subset$year <- as.ordered(dat.subset$year)  
    # day model works pretty good
      # Note: Tried fitting with hourly swdown & it didn't improve things (visually or with AIC), so for the sake of simplicity, using swdown.day
      mod.doy <- lm(tair ~ as.ordered(hour)*tmax.day*(lag.tair + lag.tmin + tmin.day) +  as.ordered(hour)*tmin.day*next.tmax + as.ordered(hour)*swdown.day*(tmax.day + tmin.day) - 1 - as.ordered(hour) - swdown.day - lag.tair - lag.tmin - next.tmax - tmax.day - tmin.day - tmin.day*tmax.day - swdown.day*tmax.day*tmin.day, data=dat.subset) #

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
        resid.model <- lm(resid ~ as.factor(hour)*(tmax.day*tmin.day)-1, data=dat.subset[!is.na(dat.subset$lag.tair),])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.swdown <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=1341){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, threshold=NULL, n.beta, resids=resids){ 
    
    # Don't bother trying to fit hours that are completely or pretty darn close to dark
    hrs.day <- unique(dat.subset[dat.subset$swdown>threshold, "hour"])
    
    # Note: played around with a log-transformation of swdown to prevent negative values, but that resulted in bias at upper range
    # Solution was to just say anything <0 = 0
    mod.doy <- lm(swdown ~ as.factor(hour)*swdown.day, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###

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
      resid.model <- lm(resid ~ as.factor(hour)*swdown.day-1, data=dat.subset[dat.subset$hour %in% hrs.day,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], threshold=quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.lwdown <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=341){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    mod.doy <- lm(lwdown ~ as.factor(hour)*lwdown.day*(lag.lwdown + next.lwdown + swdown.day + tmax.day + tmin.day) - as.factor(hour) - tmax.day - tmin.day - swdown.day - 1, data=dat.subset) ###
    # mod.doy <- lm(lwdown ~ as.factor(hour)*lwdown.day*(lag.lwdown + tair + swdown), data=dat.subset) ###
    # mod.doy0 <- lm(lwdown ~ as.factor(hour)*lwdown.day, data=dat.subset) ###
    # mod.doy1 <- lm(lwdown ~ as.factor(hour)*lwdown.day*(tmean.day + swdown.day), data=dat.subset) ###
    # mod.doy2 <- lm(lwdown ~ as.factor(hour)*lwdown.day*(tair + swdown), data=dat.subset) ###
    # mod.doy3 <- lm(lwdown ~ as.factor(hour)*lwdown.day*(tmean.day + swdown.day + lag.lwdown), data=dat.subset) ###
    # mod.doy4 <- lm(lwdown ~ as.factor(hour)*lwdown.day*(lag.lwdown + tmean.day + swdown.day + tmax.day), data=dat.subset) ###
    # mod.doy5 <- lm(lwdown ~ as.factor(hour)*lwdown.day*(lag.lwdown + swdown.day + tair), data=dat.subset) ###
    # AIC(mod.doy5)
    
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
      resid.model <- lm(resid ~ as.factor(hour)*lwdown.day-1, data=dat.subset[,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.press <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=1347){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy0 <- lm(press ~ as.factor(hour)*press.day-1-as.factor(hour), data=dat.subset) ###
    # mod.doy1 <- lm(press ~ as.factor(hour)*press.day*(lag.press)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy2 <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy3 <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press)-1-as.factor(hour)- next.press - lag.press, data=dat.subset) ###
    # mod.doy4 <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press + tmean.day)-1-as.factor(hour)- next.press - lag.press, data=dat.subset) ###
    # mod.doy5 <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press + tmax.day)-1-as.factor(hour)- next.press - lag.press, data=dat.subset) ###
    # mod.doy6 <- lm(press ~ as.factor(hour)*press.day*(lag.press + next.press + tmin.day)-1-as.factor(hour)- next.press - lag.press - press.day - tmin.day, data=dat.subset) ###
    mod.doy <- lm(press ~ as.factor(hour)*(press.day + lag.press + next.press)-as.factor(hour)-1, data=dat.subset) ###

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
      resid.model <- lm(resid ~ as.factor(hour)*press.day-1, data=dat.subset[,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.wind <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=708){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy0 <- lm(wind ~ as.factor(hour)*wind.day-1-as.factor(hour), data=dat.subset) ###
    # mod.doy1 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy2 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy3 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind)-1-as.factor(hour)- next.wind - lag.wind, data=dat.subset) ###
    # mod.doy4 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind + tmean.day)-1-as.factor(hour)- next.wind - lag.wind, data=dat.subset) ###
    # mod.doy5 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind + tmax.day)-1-as.factor(hour)- next.wind - lag.wind, data=dat.subset) ###
    # mod.doy6 <- lm(wind ~ as.factor(hour)*wind.day*(lag.wind + next.wind + tmin.day)-1-as.factor(hour)- next.wind - lag.wind - wind.day - tmin.day, data=dat.subset) ###
    # mod.doy <- lm(log(wind) ~ as.factor(hour)*log(wind.day)*(log(lag.wind) + log(next.wind) + press.day + tmin.day + tmax.day)-as.factor(hour)-1 - press.day - tmin.day - tmax.day - log(wind.day)*press.day - log(wind.day)*tmax.day - log(wind.day), data=dat.subset) ###
    # mod.doy <- lm(log(wind) ~ as.factor(hour)*wind.day*(log(lag.wind) + next.wind + press.day + tmin.day + tmax.day)-as.factor(hour)-1 - press.day - tmin.day - tmax.day - wind.day*press.day - wind.day*tmax.day - wind.day*tmin.day, data=dat.subset) ###
    # mod.doy <- lm(log(wind) ~ as.factor(hour)*log(wind.day)*(log(lag.wind) + log(next.wind) + press.day + tmin.day)-as.factor(hour)-1 - press.day - tmin.day - log(wind.day)*press.day - log(wind.day)*tmin.day, data=dat.subset) ###
    mod.doy <- lm(log(wind) ~ as.factor(hour)*log(wind.day)*(log(lag.wind) + log(next.wind) + press.day + tmin.day + tmax.day)-as.factor(hour)-1 - press.day - tmin.day - tmax.day - log(wind.day)*press.day - log(wind.day)*tmin.day- log(wind.day)*tmax.day - as.factor(hour)*tmin.day - as.factor(hour)*tmax.day - as.factor(hour)*press.day, data=dat.subset) ###
    
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
      resid.model <- lm(resid ~ as.factor(hour)*wind.day-1, data=dat.subset[,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.precipf <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=1562){
  library(MASS)
  # library(fitdistrplus)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # Precip needs to be a bit different.  We're going to calculate the fraction of precip occuring in each hour
    # we're going to estimate the probability distribution of rain occuring in a given hour
    dat.subset$rain.prop <- dat.subset$precipf/(dat.subset$precipf.day*24)
    mod.doy <- lm(rain.prop ~ as.factor(hour)*precipf.day-1 - as.factor(hour), data=dat.subset)
    
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
      resid.model <- lm(resid ~ as.factor(hour)*precipf.day-1, data=dat.subset[,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

model.qair <- function(dat.train, n.beta=1000, resids=F, parallel=F, n.cores=NULL, seed=1009){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, n.beta, resids=resids){ 
    
    # mod.doy <- lm(qair ~ as.factor(hour)*qair.day*(lag.qair + next.qair)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy0 <- lm(qair ~ as.factor(hour)*qair.day-1-as.factor(hour), data=dat.subset) ###
    # mod.doy1 <- lm(qair ~ as.factor(hour)*qair.day*(lag.qair)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy2 <- lm(qair ~ as.factor(hour)*qair.day*(lag.qair + next.qair)-1-as.factor(hour), data=dat.subset) ###
    # mod.doy3 <- lm(qair ~ as.factor(hour)*qair.day*(lag.qair + next.qair)-1-as.factor(hour)- next.qair - lag.qair, data=dat.subset) ###
    # mod.doy4 <- lm(qair ~ as.factor(hour)*qair.day*(lag.qair + next.qair + tmean.day + tmax.day + tmin.day + precipf.day)-1-as.factor(hour)- next.qair - lag.qair, data=dat.subset) ###
    mod.doy <- lm(log(qair) ~ as.factor(hour)*log(qair.day)*(log(lag.qair) + log(next.qair) + precipf.day + tmin.day + tmax.day)-as.factor(hour)-1 - precipf.day - tmin.day - tmax.day - log(qair.day)*precipf.day - log(qair.day)*tmin.day- log(qair.day)*tmax.day - as.factor(hour)*tmin.day - as.factor(hour)*tmax.day - as.factor(hour)*precipf.day, data=dat.subset) ###
    
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
      resid.model <- lm(resid ~ as.factor(hour)*qair.day-1, data=dat.subset[,])
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
  for(i in unique(dat.train$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.train[dat.train$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.train[dat.train$doy==i,]
    }
  }
  
  # Do the computation and save a list
  # Final list will have 2 layers per DOY: the model, and a bunch of simulated betas
  if(parallel==T){
    library(parallel)
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta, resids=resids)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta, resids=resids)
    }
  }
  
  return(mod.out)
}

predict.met <- function(newdata, model.predict, betas, resid.err=F, model.resid=NULL, betas.resid=NULL, n.ens, seed=9321){
  set.seed(9321)
  err.resid = 0 # dummy residual error term; if we want to add residual error, we're modeling it by hour

  mod.terms <- terms(model.predict)
  mod.coef <- coef(model.predict)
  mod.cov  <- vcov(model.predict)
  mod.resid <- resid(model.predict)
  piv <- as.numeric(which(!is.na(mod.coef)))
  
  m <- model.frame(mod.terms, newdata, xlev = model.predict$xlevels)
  Xp <- model.matrix(mod.terms, m, contrasts.arg = model.predict$contrasts)
  
  rows.beta <- sample(1:nrow(betas), n.ens, replace=T)

  Rbeta <- as.matrix(betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
  
  if(resid.err==T){
    newdata$resid <- 99999
    resid.terms <- terms(model.resid)
    resid.coef <- coef(model.resid)
    resid.cov  <- vcov(model.resid)
    resid.resid <- resid(model.resid)
    resid.piv <- as.numeric(which(!is.na(resid.coef)))
    
    m2 <- model.frame(resid.terms, newdata, xlev = model.resid$xlevels)
    Xp.res <- model.matrix(resid.terms, m2, contrasts.arg = model.resid$contrasts)
    
    rows.beta2 <- sample(1:nrow(betas.resid), n.ens, replace=T)
    
    Rbeta.res <- as.matrix(betas.resid[rows.beta2,], nrow=length(rows.beta2), ncol=ncol(betas.resid))
    
    err.resid <- Xp.res[,resid.piv] %*% t(Rbeta.res)
  }
  
  dat.sim <- Xp[,piv] %*% t(Rbeta) + err.resid
  
  return(dat.sim)
  
}

graph.resids <- function(var, dat.train, model.var, fig.dir){
  dat.train$var <- dat.train[,var]
  dat.train$lag.var <- dat.train[,paste0("lag.", var)]
  dat.train$next.var <- dat.train[,paste0("next.", var)]
  for(i in names(model.var)){
    if(as.numeric(i) == 365) next # 365 is weird, so lets skip it
    dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.var) & !is.na(dat.train$next.var), "resid"] <- resid(model.var[[i]]$model)
    dat.train[dat.train$doy==as.numeric(i) & !is.na(dat.train$lag.var) & !is.na(dat.train$next.var), "predict"] <- predict(model.var[[i]]$model)
  }
  summary(dat.train)
  
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
  plot(predict~ tair, data=day.stats, cex=0.5); abline(a=0, b=1, col="red")
  dev.off()
  
  png(file.path(fig.dir, paste0(var, "_ResidualHistograms_Day.png")), height=6, width=8, units="in", res=180)
  # par(mfrow=c(3,1))
  par(mfrow=c(1,1))
  hist(day.stats$resid, main="Daily Mean Residuals")
  dev.off()
}

graph.predict <- function(){
  # ---------
  # Graph the output
  # ---------
  {
    for(y in unique(dat.mod$year)){
      png(file.path(fig.dir, paste0("swdown_", y, "_year.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_ribbon(aes(x=date, ymin=mod.swdown.025, ymax=mod.swdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.swdown), color="blue") +
          geom_point(aes(x=date, y=mod.swdown), color="blue", size=0.5) +
          geom_line(aes(x=date, y=swdown), color="black", alpha=0.5) +
          geom_point(aes(x=date, y=swdown), color="black", size=0.3, alpha=0.5) +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("swdown_", y, "_year_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.mod[dat.mod$year==y,]) +
          geom_point(aes(x=swdown, y=mod.swdown), color="black", size=0.5) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
    }  
    
    
    dat.graph1 <- dat.mod[dat.mod$doy>=32 & dat.mod$doy<=(32+14),]
    dat.graph1$season <- as.factor("winter")
    dat.graph2 <- dat.mod[dat.mod$doy>=123 & dat.mod$doy<=(123+14),]
    dat.graph2$season <- as.factor("spring")
    dat.graph3 <- dat.mod[dat.mod$doy>=214 & dat.mod$doy<=(213+14),]
    dat.graph3$season <- as.factor("summer")
    dat.graph4 <- dat.mod[dat.mod$doy>=305 & dat.mod$doy<=(305+14),]
    dat.graph4$season <- as.factor("fall")
    
    dat.graph <- rbind(dat.graph1, dat.graph2, dat.graph3, dat.graph4)
    
    for(y in unique(dat.graph$year)){
      png(file.path(fig.dir, paste0("swdown_",y,"_examples.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_line(aes(x=date, y=swdown), color="black") +
          geom_point(aes(x=date, y=swdown), color="black", size=0.5) +
          geom_ribbon(aes(x=date, ymin=mod.swdown.025, ymax=mod.swdown.975), alpha=0.5, fill="blue") +
          geom_line(aes(x=date, y=mod.swdown), color="blue") +
          geom_point(aes(x=date, y=mod.swdown), color="blue", size=0.5) +
          scale_y_continuous(name="shortwave radiation") +
          scale_x_datetime(expand=c(0,0)) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
      
      png(file.path(fig.dir, paste0("swdown_", y, "_examp.es_scatter.png")), height=8, width=10, units="in", res=220)
      print(
        ggplot(data=dat.graph[dat.graph$year==y,]) +
          facet_wrap(~season, scales="free") +
          geom_point(aes(x=swdown, y=mod.swdown), color="black", size=0.5) +
          ggtitle(y) +
          theme_bw()
      )
      dev.off()
    }
  }
  # ---------
  
}