model.tair <- function(dat.train, n.cores=4, n.beta=1000, parallel=F, ncores=NULL, seed=1237){
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
      
      return(list(model=mod.doy, betas=Rbeta))
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

model.swdown <- function(dat.train, n.cores=4, n.beta=1000, parallel=F, ncores=NULL, seed=1341){
  library(MASS)
  set.seed(seed)
  
  # The model we're going to use
  model.train <- function(dat.subset, threshold=NULL, n.beta){ 
    # day model works pretty good
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day + as.factor(hour)*tair*swdown.day - as.factor(hour) - swdown.day - tair -1, data=dat.subset) #
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*tmax.day  - as.factor(hour) - tmax.day - swdown.day - swdown.day*tmax.day -1, data=dat.subset) #
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*max.dep*tmean.day  - as.factor(hour) - max.dep - swdown.day - tmean.day - tmean.day*max.dep -1, data=dat.subset) # Pretty good, but max swdown too high
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*max.dep*min.dep*tmean.day  - as.factor(hour) - max.dep - min.dep - swdown.day - tmean.day - tmean.day*max.dep*min.dep - max.dep*min.dep - as.factor(hour)*max.dep - as.factor(hour)*min.dep - as.factor(hour)*tmean -1, data=dat.subset) # Pretty good, but max swdown too high
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*max.dep*tmean.day  - as.factor(hour) - max.dep - swdown.day - tmean.day - swdown*tmean.day*max.dep -1, data=dat.subset) # Pretty good, but max swdown too high
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*tmax.day - as.factor(hour) - tmax.day - swdown.day - swdown.day*tmax.day -1, data=dat.subset) ### Pretty good
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day + press.day) - as.factor(hour) - tmax.day - press.day - swdown.day -1, data=dat.subset) ### BEST SO FAR
    # mod.doy1 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day + press.day + lwdown.day) - as.factor(hour) - swdown.day - tmax.day - press.day - lwdown.day -1, data=dat.subset) ###
    # mod.doy2 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset) ###
    # mod.doy3 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(press.day) - as.factor(hour) - swdown.day - press.day -1, data=dat.subset) ###
    # mod.doy4 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(lwdown.day) - as.factor(hour) - swdown.day - lwdown.day -1, data=dat.subset) ###
    # mod.doy5 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day + press.day) - as.factor(hour) - swdown.day - tmax.day - press.day -1, data=dat.subset) ###
    # mod.doy6 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day - swdown.day*tmax.day -1, data=dat.subset) ###
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day - swdown.day*tmax.day -1, data=dat.subset) ### 
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - swdown.day - tmax.day - swdown.day*tmax.day -1, data=dat.subset) ### BEST
    # mod.doyA <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - tmax.day -1, data=dat.subset) ### 
    # mod.doyB <- lm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - tmax.day -1, data=dat.subset) ### 
    # mod.doyC <- lm(log(swdown) ~ as.factor(hour)*log(swdown.day)*(tmax.day) - tmax.day -1, data=dat.subset) ### 
    # mod.doyB <- glm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) -1, data=dat.subset, family="Gamma") ### 
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day)  - swdown.day - tmax.day -1, data=dat.subset) ### 
    # mod.doy <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset) ### 


    # Trying different distributions
    # mod.doy0 <- glm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset) ### 
    # mod.doy1 <- glm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset) ### 
    # mod.doy2 <- glm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset, family=Gamma(link="log")) ### 
    # mod.doy3 <- glm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset, family=quasipoisson(link="log")) ### 
    # mod.doy4 <- glm(sqrt(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day - tmax.day -1, data=dat.subset) ### 
    
    hrs.day <- unique(dat.subset[dat.subset$swdown>threshold, "hour"])
    
    # mod.doy <- lm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - as.factor(hour) - swdown.day*tmax.day - swdown.day - tmax.day -1, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    mod.doy <- lm(swdown ~ as.factor(hour)*swdown.day - as.factor(hour) -1, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    # mod.doy1 <- lm(swdown ~ as.factor(hour)*swdown.day*(tmax.day) - swdown.day - tmax.day -1, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    # mod.doy2 <- lm(swdown ~ as.factor(hour)*swdown.day*(tmax.day + lag.swdown) - swdown.day - lag.swdown - tmax.day -1, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    # mod.doy2 <- lm(log(swdown) ~ as.factor(hour)*swdown.day*(tmax.day) - swdown.day - tmax.day -1, data=dat.subset[dat.subset$hour %in% hrs.day,]) ###
    
    # AIC(mod.doyA)
    # AIC(mod.doyB)
    # AIC(mod.doyC)
    # AIC(mod.doy0)
    # AIC(mod.doy1)
    # AIC(mod.doy2)
    # AIC(mod.doy3)
    # AIC(mod.doy4)
    # AIC(mod.doy5)
    # AIC(mod.doy6)

    # plot(predict(mod.doy1)~dat.subset[dat.subset$hour %in% hrs.day,"swdown"]); abline(a=0,b=1, col="red")
    # plot(exp(predict(mod.doy2))~dat.subset[dat.subset$hour %in% hrs.day,"swdown"]); abline(a=0,b=1, col="red")
    # # plot(predict(mod.doy0)~log(dat.subset$swdown)); abline(a=0,b=1, col="red")
    # plot(exp(predict(mod.doy0))~dat.subset$swdown); abline(a=0,b=1, col="red")
    # plot(predict(mod.doy1)~dat.subset$swdown); abline(a=0,b=1, col="red")
    # plot(exp(predict(mod.doy2))~dat.subset$swdown); abline(a=0,b=1, col="red")
    # plot(exp(predict(mod.doy3))~dat.subset$swdown); abline(a=0,b=1, col="red")
    # plot((predict(mod.doy4)^2)~dat.subset$swdown); abline(a=0,b=1, col="red")
    
    # Generate a bunch of random coefficients that we can pull from 
    # without needing to do this step every day
    mod.coef <- coef(mod.doy)
    mod.cov  <- vcov(mod.doy)
    piv <- as.numeric(which(!is.na(mod.coef)))
    Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
    
    return(list(model=mod.doy, betas=Rbeta))
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




predict.met <- function(newdata, mod.predict, betas, n.ens, seed=9321){
  set.seed(9321)

  mod.terms <- terms(mod.predict)
  mod.coef <- coef(mod.predict)
  mod.cov  <- vcov(mod.predict)
  mod.resid <- resid(mod.predict)
  piv <- as.numeric(which(!is.na(mod.coef)))
  
  m <- model.frame(mod.terms, newdata, xlev = mod.predict$xlevels)
  Xp <- model.matrix(mod.terms, m, contrasts.arg = mod.predict$contrasts)
  
  rows.beta <- sample(1:nrow(betas), n.ens, replace=T)

  Rbeta <- as.matrix(betas[rows.beta,], nrow=length(rows.beta), ncol=ncol(betas))
  
  dat.sim <- Xp[,piv] %*% t(Rbeta) #+ rnorm(n.ens, mean=mean(mod.resid), sd=sd(mod.resid))
  
  return(dat.sim)
  
}