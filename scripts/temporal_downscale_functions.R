model.tair <- function(dat.tair, n.cores=4, n.beta=100, parallel=F, ncores=NULL, seed=1237){
  library(MASS)
  set.seed(seed)

  # The model we're going to use
  model.train <- function(dat.subset, n.beta){ 
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr  + lag.hr*tmean + as.factor(hour)*tmean*max.dep*min.dep - as.factor(hour) - 1 - tmean - lag.hr - max.dep - min.dep - max.dep*min.dep-tmean*max.dep*min.dep, data=dat.subset) # REALLY REALLY GOOD
      
      # Tryign with a single met covariate
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr + as.factor(hour)*swdown + lag.hr*tmean + as.factor(hour)*tmean*max.dep*min.dep + as.factor(hour)*swdown*max.dep*min.dep - as.factor(hour) - 1 - swdown*tmean - tmean - lag.hr - max.dep - min.dep - max.dep*min.dep, data=dat.subset) # REALLY REALLY GOOD
      
      # Tryign with two covariates
      mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr + as.factor(hour)*swdown + lag.hr*tmean + as.factor(hour)*tmean*max.dep*min.dep + as.factor(hour)*swdown*max.dep*min.dep + as.factor(hour)*press - as.factor(hour) - 1 - swdown - press - tmean - lag.hr - max.dep - min.dep - max.dep*min.dep, data=dat.subset) # REALLY REALLY GOOD
      

      # Rejected options
      # mod.doy <- lm(tair ~ lag.hr*tmean + as.factor(hour)*max.dep*min.dep - as.factor(hour), data=dat.subset) # Not bad
      # mod.doy <- lm(tair ~ as.factor(hour)*lag.hr  + lag.hr*tmean + as.factor(hour)*max.dep*min.dep - as.factor(hour), data=dat.subset) # Good
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean  + lag.hr*tmean + as.factor(hour)*max.dep*min.dep - as.factor(hour), data=dat.subset) # better
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr  + lag.hr*tmean + as.factor(hour)*max.dep*min.dep - as.factor(hour) - 1 - tmean - lag.hr, data=dat.subset) # REALLY GOOD
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr  + lag.hr*tmean + tmean*max.dep*min.dep + as.factor(hour)*max.dep*min.dep - as.factor(hour) - 1 - tmean - lag.hr, data=dat.subset) # REALLY REALLY GOOD
      # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr  + lag.hr*tmean + as.factor(hour)*tmean*max.dep*min.dep - as.factor(hour) - 1 - tmean - lag.hr - max.dep - min.dep - max.dep*min.dep-tmean*max.dep*min.dep, data=dat.subset) # REALLY REALLY GOOD
      # mod.doy <- lm(tair ~ as.factor(hour)*(tmean + tmin + tmax) + as.factor(hour)*lag.hr  + lag.hr*tmean + as.factor(hour)*tmean*max.dep*min.dep - as.factor(hour) - 1 - tmean - tmin - tmax - lag.hr - max.dep - min.dep - max.dep*min.dep-tmean*max.dep*min.dep, data=dat.subset) # REALLY REALLY GOOD
    # # mod.doy <- lm(tair ~ as.factor(hour)*tmean + as.factor(hour)*lag.hr  + as.factor(hour)*tmean*max.dep*min.dep + as.factor(hour)*swdown + swdown*max.dep*min.dep - as.factor(hour) - 1 - tmean - lag.hr -swdown, data=dat.subset) #
    
    
      # # Okay in short runs, but wonky values when run over a year
      # mod.doy <- lm(tair ~ as.factor(hour)*lag.hr*(tmax*tmin +tmean) - as.factor(hour)-1, data=dat.subset) # Original; not bad
      # mod.doy <- lm(tair ~ tmax*tmin*tmean*(as.factor(hour) + lag.hr) - as.factor(hour)-1, data=dat.subset) # Original; not bad
      # mod.doy <- lm(tair ~ as.factor(hour)*lag.hr*(tmax*tmin +tmean) + tmean*max.dep*min.dep - as.factor(hour)-1, data=dat.subset) # good
      # mod.doy <- lm(tair ~ as.factor(hour)*lag.hr*(tmax*tmin + tmean)  + tmean*max.dep*min.dep + lag.hr*max.dep*min.dep - as.factor(hour)-1, data=dat.subset) #  best
      mod.coef <- coef(mod.doy)
      mod.cov  <- vcov(mod.doy)
      piv <- as.numeric(which(!is.na(mod.coef)))
      Rbeta <- mvrnorm(n=n.beta, mod.coef[piv], mod.cov)
      
      return(list(model=mod.doy, betas=Rbeta))
    }

  dat.list <- list()
  mod.out <- list()

  # Make the data into a list
  for(i in unique(dat.tair$doy)){
    if(i == 365){ # Lump leap day in with non-leap Dec 31
      dat.list[[paste(i)]] <- dat.tair[dat.tair$doy>=364,]
    } else {
      dat.list[[paste(i)]] <- dat.tair[dat.tair$doy==i,]
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
  
}