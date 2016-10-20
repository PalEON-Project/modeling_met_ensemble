model.tair <- function(dat.train, n.cores=4, n.beta=1000, resids=F, parallel=F, ncores=NULL, seed=1237){
  library(MASS)
  set.seed(seed)

  # The model we're going to use
  model.train <- function(dat.subset, n.beta){ 
      # day model works pretty good
      mod.doy <- lm(tair ~ as.factor(hour)*tmean.day*(lag.tair + lag.tmin + max.dep*min.dep) + as.factor(hour)*swdown.day*(max.dep + min.dep) - 1 - as.factor(hour) - swdown.day - tmean.day - lag.tair - lag.tmin - min.dep - max.dep - min.dep*max.dep - tmean.day*max.dep*min.dep - swdown.day*max.dep*min.dep, data=dat.subset) #

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
        mod.resid <- resid(mod.doy)
        resid.model <- lm(mod.resid ~ as.factor(dat.subset$hour)*(dat.subset$tmean + dat.subset$max.dep + dat.subset$min.dep)-1)
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
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], n.beta=n.beta)
    }
  }
  
  return(mod.out)
}

model.swdown <- function(dat.train, n.cores=4, n.beta=1000, resids=F, parallel=F, ncores=NULL, seed=1341){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, threshold=NULL, n.beta){ 
    
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
      mod.resid <- resid(mod.doy)
      resid.model <- lm(mod.resid ~ as.factor(dat.subset[dat.subset$hour %in% hrs.day,"hour"])*dat.subset[dat.subset$hour %in% hrs.day,"swdown.day"]-1)
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
    mod.out <- mclapply(dat.list, model.train, mc.cores=n.cores, n.beta=n.beta)
  } else {
    for(i in names(dat.list)){
      mod.out[[i]] <- model.train(dat.subset=dat.list[[i]], threshold=quantile(dat.train[dat.train$swdown>0,"swdown"], 0.05), n.beta=n.beta)
    }
  }
  
  return(mod.out)
}




predict.met <- function(newdata, mod.predict, betas, resids=F, mod.resid=NULL, betas.resid=NULL, n.ens, seed=9321){
  set.seed(9321)
  err.resid = 0 # dummy residual error term; if we want to add residual error, we're modeling it by hour

  mod.terms <- terms(mod.predict)
  mod.coef <- coef(mod.predict)
  mod.cov  <- vcov(mod.predict)
  mod.resid <- resid(mod.predict)
  piv <- as.numeric(which(!is.na(mod.coef)))
  
  m <- model.frame(mod.terms, newdata, xlev = mod.predict$xlevels)
  Xp <- model.matrix(mod.terms, m, contrasts.arg = mod.predict$contrasts)
  
  rows.beta <- sample(1:nrow(betas), n.ens, replace=T)

  Rbeta <- as.matrix(betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
  
  if(resids==T){
    resid.terms <- terms(mod.resid)
    resid.coef <- coef(mod.resid)
    resid.cov  <- vcov(mod.resid)
    resid.resid <- resid(mod.resid)
    resid.piv <- as.numeric(which(!is.na(resid.coef)))
    
    m2 <- model.frame(resid.terms, newdata, xlev = mod.resid$xlevels)
    Xp.res <- model.matrix(resid.terms, m2, contrasts.arg = mod.resid$contrasts)
    
    rows.beta <- sample(1:nrow(betas), n.ens, replace=T)
    
    Rbeta.res <- as.matrix(betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
    
    err.resid <- Xp.res[,resid.piv] %*% t(Rbeta.res)
  }
  
  dat.sim <- Xp[,piv] %*% t(Rbeta) + err.resid
  
  return(dat.sim)
  
}