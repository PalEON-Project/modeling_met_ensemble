bias.correct <- function(met.bias, vars.met, dat.train, yrs.cal, n){
# Inputs:
# met.bias   = data frame of all daily met data
# vars.met  = vector of variables to downscale, **in the order they should be done**
# dat.train = name of the training dataset
# yrs.cal - data frame of the series to be bias-corrected with three columns:
#             1. dataset name, 
#             2. first year for calibration
#             3. last year for calibration
# n         = number of ensemble members

dat.train.orig <- dat.train
# Order yrs.cal so that we work from most recent to furthest back, propogating uncertainty as we go
yrs.cal <- yrs.cal[order(yrs.cal$cal.min, decreasing=T),]

vars.pred <- vector()
for(v in 1:length(vars.met)){
  # define what met variable we're working with
  met.var = vars.met[v]
  dat.train <- dat.train.orig
  
  print(met.var)

  # Getting building covariance matrices leveraging what we've already gapfilled
  vars.pred <- vector()
  if(v>1) vars.pred <- c(vars.pred, vars.met[v-1])

  # -------------
  # Creating the output list with LDAS as our base before we forget
  # -------------
  # dat.train = yrs.cal[1,"dataset"] # declared
  yr.min    = min(met.bias[met.bias$dataset==dat.train, "year"])
  yr.max    = max(met.bias[met.bias$dataset==dat.train, "year"])
  
  # Calculating the anomalies from the current means
  dat.temp <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min & met.bias$year<=yr.max, ]
  dat.temp$Y <- dat.temp[,met.var]
  lm.train <- lm(Y ~ as.factor(doy)-1, data=dat.temp)
  train.anom <- met.bias[met.bias$dataset==dat.train, met.var] - predict(lm.train, newdata=met.bias[met.bias$dataset==dat.train, c("year", "doy")])
  summary(train.anom)
  
  train.ci <- data.frame(dataset=dat.train, met=met.var, met.bias[met.bias$dataset==dat.train, c(c("year", "doy", vars.met))],
                         X=met.bias[met.bias$dataset==dat.train, met.var], 
                         anom.raw=train.anom,
                         mean=met.bias[met.bias$dataset==dat.train, met.var], 
                         lwr=met.bias[met.bias$dataset==dat.train, met.var], 
                         upr=met.bias[met.bias$dataset==dat.train, met.var],
                         time=as.Date(met.bias[met.bias$dataset==dat.train, "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset==dat.train, "year"], "-01-01")))
                         )
  train.sims <- data.frame(dataset=dat.train, met=met.var, met.bias[met.bias$dataset==dat.train, c("year", "doy", vars.met)], 
                           time=as.Date(met.bias[met.bias$dataset==dat.train, "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset==dat.train, "year"], "-01-01")))
                           )
  train.sims[,paste0("X", 1:n)] <- met.bias[met.bias$dataset==dat.train, met.var]
  
  if(v==1) dat.out <- list()
  dat.out[[met.var]] <- list()
  dat.out[[met.var]]$ci <- train.ci
  dat.out[[met.var]]$sims <- train.sims
  # -------------
  
  for(i in 1:nrow(yrs.cal)){
    # --------------------
    # Working through one dataset at a time
    # --------------------
    # Define the training dataset
    dat.bias  = paste(yrs.cal[i, "dataset"])
    if(i!=1) dat.train = paste(yrs.cal[i-1, "dataset"]  )
    yr.min    = yrs.cal[yrs.cal$dataset==dat.bias, "cal.min"]
    yr.max    = yrs.cal[yrs.cal$dataset==dat.bias, "cal.max"]
    
    # If there is insufficient overlap in the calibration years between the training data and the bias-corrected data,
    # then grab the closest time period with a similar number of years
    # Need to specify a separate range for the training data -- use the same number of years!
    if(yr.max < min(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train, "year"]) + (yr.max - yr.min)){
      yr.min.train = min(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train, "year"])
      yr.max.train = yr.min.train + (yr.max - yr.min)
    } else {
      yr.min.train = yr.min
      yr.max.train = yr.max
    }
    
    
    # ---------
    # Set up the calibration & prediction data frames
    # ---------
    # A. Calibration -- Aggregate so we're looking at the climatic mean for the first part
    dat.temp0 <- data.frame(year = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "year" ],
                            doy  = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "doy"  ],
                            stack(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, paste0("X", 1:n)])
                            )
        
    dat.temp <- aggregate(dat.temp0$values, by=dat.temp0[,c("doy", "ind")], FUN=mean)
    names(dat.temp)[3] <- "Y"
    summary(dat.temp)
    
    dat.temp2 <- aggregate(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, c(met.var, vars.met)],
                           by=list(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"]),
                           FUN=mean)
    names(dat.temp2) <- c("doy", c("X", vars.met))
    summary(dat.temp2)
    
    # merging together the two sets of daily means 
    # Note: working off of the climatological means because we'll model the anomalies separately    
    dat.temp <- merge(dat.temp[,], dat.temp2)
    summary(dat.temp)
        
    # getting the raw data for the calibration periods so we can model the anomalies
    raw.train <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min.train & met.bias$year<=yr.max.train, c("dataset", "year", "doy", met.var)]
    raw.bias  <- met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min       & met.bias$year<=yr.max      , c("dataset", "year", "doy", met.var, vars.met)]
    names(raw.train) <- c("dataset", "year", "doy", "X")
    names(raw.bias) <- c("dataset", "year", "doy", "X", vars.met)
    
    # The prediction data
    dat.pred <- data.frame(met.bias[met.bias$dataset==dat.bias, c("year", "doy",  vars.met)],
                           X    = met.bias[met.bias$dataset==dat.bias, met.var]
                           )
    
    # 0 out unnecessary predictors in all data frames
    dat.temp [,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    raw.train[,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    raw.bias [,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    dat.pred [,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    # ---------
    
    
    # ---------
    # Doing the climatological bias correction
    # ---------
  #   # Plotting the correlation between dat.train & dat.bias
  #   plot(Y ~ X, data=dat.temp) 
  #   abline(a=0, b=1, col="red")
  #   abline(lm(Y ~ X, data=dat.temp), col="red", lty="dashed")
    
    # Running a model
    mod.bias <- gam(Y ~ s(doy) + X + tmax + tmin + precipf + swdown + qair + lwdown + press + wind, data=dat.temp)
    summary(mod.bias)
  #   plot(mod.bias)
    
    # Saving the mean predicted & residuals
    dat.temp$pred  <- predict(mod.bias)
    dat.temp$resid <- resid(mod.bias)
    summary(dat.temp)
  
  #   # Checkign the residuals to see if we can assume normality
  #   plot(resid ~ pred, data=dat.temp); abline(h=0, col="red")
  #   plot(resid ~ doy, data=dat.temp); abline(h=0, col="red")
  #   hist(dat.temp$resid)
    
    dat.pred$pred <- predict(mod.bias, newdata=dat.pred)
    # ---------
    
    # ---------
    # Modeling the anomalies with a similar framework
    # ---------
    # Calculating the anomalies from the current means
    anom.train <- gam(X ~ as.factor(doy)-1, data=raw.train)
    anom.bias  <- gam(X ~ as.factor(doy)-1, data=raw.bias)
    raw.train$anom.train <- resid(anom.train)
    raw.bias $anom.raw   <- resid(anom.bias)
    
    dat.pred$anom.raw <- dat.pred$X - predict(anom.bias, newdata=dat.pred)
  
    # Modeling the anomalies of the other predictors
    
    for(j in vars.met){
      raw.bias$Q <- raw.bias[,j]
      dat.pred$Q <- dat.pred[,j]
      
      if(j %in% vars.pred){
        anom.bias2  <- gam(Q ~ as.factor(doy)-1, data=raw.bias)
        raw.bias[,paste0(j, ".anom")] <- resid(anom.bias2)
        dat.pred[,paste0(j, ".anom")] <- dat.pred$Q - predict(anom.bias2, newdata=dat.pred)
      } else {
        raw.bias[,paste0(j, ".anom")] <- 0
        dat.pred[,paste0(j, ".anom")] <- 0
      }
    }
  
  
    # Modeling the anomalies
    empirical = c("NLDAS", "GLDAS", "CRUNCEP", "PRISM", "DAYMET")
    # Figure out if we can line up the anomalies for better accuracy
    if(dat.bias %in% empirical & dat.train %in% empirical){ # if it's empirical we can, pair the anomalies for best estimation
      dat.anom <- merge(raw.bias[,c("year", "doy", "X", "anom.raw", vars.met, paste0(vars.met, ".anom"))], raw.train[,c("year", "doy", "anom.train")])
    } else { # If we can't, just cbind them and we'll be simulating a bunch of anomalies
      dat.anom <- cbind(raw.bias[,c("year", "doy", "X", "anom.raw", vars.met, paste0(vars.met, ".anom"))], anom.train=raw.train$anom.train)
    }
    
  #   # The anomaly model
  #   plot(anom.train ~ anom.raw, data=dat.anom)
  #   abline(a=0, b=1, col="red")
  #   abline(lm(anom.train ~ anom.raw, data=dat.anom), col="red", lty="dashed")
  
    # First choice for modeling the anomalies: 
    #  1. Use a correlation between anom.raw & anom.train to adjust variance
    mod.anom <- gam(anom.train ~ s(doy) + anom.raw + tmax.anom + tmin.anom + precipf.anom + swdown.anom + qair.anom + lwdown.anom + press.anom + wind.anom, data=dat.anom)
    summary(mod.anom)
  #   plot(mod.anom)
    
  #   dat.anom$resid.anom <- resid(mod.anom)
  #   plot(resid.anom ~ doy, data=dat.anom); abline(h=0, col="red")
  #   hist(dat.anom$resid.anom)
    
    # If there's no correlation between the anomalies (either anom.raw or one of our other anoms), 
    #  then just model the inerannual variability in a way that
    #  tries to preserve low- and mid-frequence events in the GCM
    #  by fitting the anomalies of the prediction data
    anom.pv <- summary(mod.anom)$p.pv
    if(!anom.pv["anom.raw"]<0.05){
      # Modeling in the predicted value from mod.bias
      dat.anom$pred <- predict(mod.bias, newdata=dat.anom)
      
      k=round(length(dat.pred$year)/(25*366),0)
      mod.anom <- gam(anom.raw ~ s(doy) + s(year, k=k) + tmax.anom + tmin.anom + precipf.anom + swdown.anom + qair.anom + lwdown.anom + press.anom + wind.anom, data=dat.pred)
  
      summary(mod.anom)
  #     plot(mod.anom, pages=1)
      
      dat.pred$resid.anom <- resid(mod.anom)
  #     plot(resid.anom ~ doy, data=dat.pred); abline(h=0, col="red")
  #     plot(resid.anom ~ year, data=dat.pred); abline(h=0, col="red")
  #     hist(dat.pred$resid.anom)  
    }
  
    # ---------
    
    # --------
    # Predicting a bunch of potential posteriors over the full dataset
    # --------
    # Get the model coefficients
    coef.gam <- coef(mod.bias)
    coef.anom <- coef(mod.anom)
    
    # Generate a random distribution of betas using the covariance matrix
    Rbeta <- mvrnorm(n=n, coef(mod.bias), vcov(mod.bias))
    Rbeta.anom <- mvrnorm(n=n, coef(mod.anom), vcov(mod.anom))
    
    # Create the prediction matrix
    Xp <- predict(mod.bias, newdata=dat.pred, type="lpmatrix")
    Xp.anom <- predict(mod.anom, newdata=dat.pred, type="lpmatrix")
    
    # -----
    # Simulate predicted met variables
    # -----
    # Adding in some residual error 
    # ** QUESTION -- do I add a constant error per time series or random error to each day?
    # -----
    sim1 <- Xp %*% t(Rbeta) + Xp.anom %*% t(Rbeta.anom)
    
    # # Option 1: Adding a constant error per time series
    # res <- rnorm(n, mean(dat.temp$resid), sd(dat.temp$resid))
    # sim1 <- sweep(sim1, 2, res, FUN="+")
    
    # # Option 2: Adding a random error to each observation
    # sim1 <- sim1 + array(rnorm(length(sim1), mean(dat.temp$resid), sd(dat.temp$resid)), dim=dim(sim1))
    
    # # Option 3: Adding a day-specific error (although this should be part of the smoother)
    # -----
    if(met.var=="precipf"){
  #     # because of how the log transformation works, our big values get super exagerated & we don't want to 
  #     # propagate unreasonable values
  #     # Solution: hard code in a max where anything >= max(training)*1.25 will get that hard boundary 
  #     # For Harvard, this is 13% of the CRUNCEP values, plus 10 from p1000
  #     max.empirical <- max(exp(met.bias$precipf))*1.25 # doing this in real units just to make life easier
  #     sim1[sim1>log(max.empirical)] <- log(max.empirical)
      # just in case, get rid of any negative values in precip
      sim1[sim1<0] = 0
    }
  
    # NOTE: Don't need to add in anomalies because training & fitting data are both observed annual deviations!
    dat.pred$mean <- apply(sim1, 1, mean)
    dat.pred$lwr  <- apply(sim1, 1, quantile, 0.025)
    dat.pred$upr  <- apply(sim1, 1, quantile, 0.975)
    summary(dat.pred)
    
    # Doing som exploratory graphing
    dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))
    
  #   # Plotting the observed and the bias-corrected 95% CI
  #   ggplot(data=dat.pred[dat.pred$year>=(yr.max-5) & dat.pred$year<=(yr.max-3),]) +
  #     geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  #     geom_line(aes(x=time, y=mean), color="red", size=0.5) +
  #     geom_line(aes(x=time, y=X), color='black', size=0.5) +
  #     theme_bw()
  #   
  #   dat.yr <- aggregate(dat.pred[,c("X", "mean", "lwr", "upr")],
  #                       by=list(dat.pred$year),
  #                       FUN=mean)
  #   names(dat.yr)[1] <- "year"
  #   
  #   ggplot(data=dat.yr[,]) +
  #     geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
  #     geom_line(aes(x=year, y=mean), color="red", size=0.5) +
  #     geom_line(aes(x=year, y=X), color='black', size=0.5) +
  #     theme_bw()
    # --------
    
    # --------
    # Storing the output
    # --------
    dat.sims <- data.frame(dataset=dat.bias, met=met.var, dat.pred[,c("year", "doy", vars.met, "time")])
    dat.sims <- cbind(dat.sims, sim1)
    
    dat.out[[met.var]]$ci   <- rbind(dat.out[[met.var]]$ci, data.frame(dataset=dat.bias, met=met.var, dat.pred[,names(dat.out[[met.var]]$ci)[3:ncol(dat.out[[met.var]]$ci)]]))
    dat.out[[met.var]]$sims <- rbind(dat.out[[met.var]]$sims, data.frame(dat.sims))
    
    # overwriting the raw with the bias-corrected output 
    met.bias[met.bias$dataset==dat.bias, met.var] <- dat.pred$mean
    # --------
    # --------------------
  }
}
dat.out$met.bias <- met.bias
return(dat.out)
}