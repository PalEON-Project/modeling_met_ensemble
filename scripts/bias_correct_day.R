bias.correct <- function(met.bias, vars.met, dat.train, GCM, yrs.cal, n, path.out, seed=2509){
# Inputs:
# met.bias   = data frame of all daily met data
# vars.met  = vector of variables to downscale, **in the order they should be done**
# dat.train = name of the training dataset
# yrs.cal - data frame of the series to be bias-corrected with three columns:
#             1. dataset name, 
#             2. first year for calibration
#             3. last year for calibration
# n         = number of ensemble members
# path.out  = file path for where the automatically generated QAQC graphs should go

# Load required libraries
library(MASS)
library(mgcv)
library(ggplot2)
set.seed(seed)
dir.create(file.path(path.out, "BiasCorrect_QAQC"), recursive=T, showWarnings=F)  
empirical = c("NLDAS", "GLDAS", "CRUNCEP", "PRISM", "DAYMET") # Empirical datasets

dat.train.orig <- dat.train

# Order yrs.cal so that we work from most recent to furthest back, propogating uncertainty as we go
yrs.cal <- yrs.cal[order(yrs.cal$cal.min, decreasing=T),]

# # Smoothing over the Jan 1 transition by adding in the previous December
# for(i in unique(met.bias$dataset)){
#   
# }

vars.pred <- vector()
dat.out <- list()
for(v in 1:length(vars.met)){
# for(v in 1:7){
  # -------------
  # Doing some initial setup to prepare for looping for a variable for the first time
  # -------------
  met.var = vars.met[v]
  dat.train <- dat.train.orig
  
  print(met.var)

  # -------------
  # Creating the output list with LDAS as our base before we forget
  # -------------
  # dat.train = yrs.cal[1,"dataset"] # declared
  yr.min    = min(met.bias[met.bias$dataset==dat.train, "year"])
  yr.max    = max(met.bias[met.bias$dataset==dat.train, "year"])
  
  # Calculating the anomalies from the current means
  dat.temp <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min & met.bias$year<=yr.max, ]
  dat.temp$Y <- dat.temp[,met.var]
  lm.train <- gam(Y ~ s(doy), data=dat.temp)
  train.anom <- met.bias[met.bias$dataset==dat.train, met.var] - predict(lm.train, newdata=met.bias[met.bias$dataset==dat.train, c("year", "doy")])
  summary(train.anom)
  
  train.ci <- data.frame(dataset=dat.train, met=met.var, met.bias[met.bias$dataset==dat.train, c("year", "doy")],
                         X=met.bias[met.bias$dataset==dat.train, met.var], 
                         anom.raw=train.anom,
                         mean=met.bias[met.bias$dataset==dat.train, met.var], 
                         lwr=met.bias[met.bias$dataset==dat.train, met.var], 
                         upr=met.bias[met.bias$dataset==dat.train, met.var],
                         time=as.Date(met.bias[met.bias$dataset==dat.train, "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset==dat.train, "year"], "-01-01")))
                         )
  train.sims <- data.frame(dataset=dat.train, met=met.var, met.bias[met.bias$dataset==dat.train, c("year", "doy")], 
                           time=as.Date(met.bias[met.bias$dataset==dat.train, "doy"], origin=as.Date(paste0(met.bias[met.bias$dataset==dat.train, "year"], "-01-01")))
                           )
  train.sims[,paste0("X", 1:n)] <- met.bias[met.bias$dataset==dat.train, met.var]
  
  dat.out[[met.var]] <- list()
  dat.out[[met.var]]$ci <- train.ci
  dat.out[[met.var]]$sims <- train.sims
  # -------------
  
  # -------------
  # If we're dealing with precip, lets keep the training data handy &
  # calculate the number of rainless days in each year to make sure we 
  # don't get a constant drizzle
  # -------------
  if(met.var=="precipf"){
    rain.train <- met.bias[met.bias$dataset==dat.train.orig,]
    rainless <- vector()
    for(y in min(rain.train$year):max(rain.train$year)){
      n.rainless <- nrow(rain.train[rain.train$year==y & rain.train$precipf==0,])
      rainless <- c(rainless, n.rainless)
    }
    mean(rainless); sd(rainless)
    rain.max <- max(rain.train$precipf)*1.5
  }
  # -------------
  
  for(i in 1:nrow(yrs.cal)){ # Working through one dataset at a time
    # ---------
    # Define the raw & training datasets and some accompanying info
    # ---------
    dat.bias  = paste(yrs.cal[i, "dataset"])
    if(i!=1) dat.train = paste(yrs.cal[i-1, "dataset"]  ) 
    yr.min    = yrs.cal[yrs.cal$dataset==dat.bias, "cal.min"]
    yr.max    = yrs.cal[yrs.cal$dataset==dat.bias, "cal.max"]
    
    print(paste0(" -- ", dat.bias))
    
    # Ideally we'll have overlap between our datasets, but if there is little to no overlap 
    # between the calibration years in the raw & training data, grab the closest window 
    # period with a similar number of years
    if(yr.max < min(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train, "year"]) + (yr.max - yr.min)){
      yr.min.train = min(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train, "year"])
      yr.max.train = yr.min.train + (yr.max - yr.min)
    } else {
      yr.min.train = yr.min
      yr.max.train = yr.max
    }
    # ---------
    
    
    # ---------
    # Set up the calibration & prediction data frames
    # ---------
    # 1. Grab the training data -- this will be called "Y" in our bias correction equations
    #     -- preserving the different simulations so we can look at a distribution of potential values
    #     -- This will get aggregated right off the bat so we so we're looking at the climatic means 
    #        for the first part of bias-correction
    dat.stack <- stack(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, paste0("X", 1:n)])
    dat.temp0 <- data.frame(year    = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "year" ],
                            doy     = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "doy"  ],
                            values  = dat.stack$values,
                            ind     = dat.stack$ind
                            )
    
    # For precip, we want to adjust the total annual precipitation, and then calculate day of year adjustment &
    # anomaly as fraction of total annual precipitation
    if(met.var=="precipf"){
        temp.ann <- aggregate(dat.temp0$values, by=dat.temp0[,c("year", "ind")], FUN=sum) 
    	names(temp.ann)[3] <- "Y.tot"
    	
    	dat.temp0 <- merge(dat.temp0, temp.ann)
    	dat.temp0$values <- dat.temp0$values/dat.temp0$Y.tot # our "Y" is now faction of annual precip in each day
    }
    
    dat.temp <- aggregate(dat.temp0$values, by=dat.temp0[,c("doy", "ind")], FUN=mean) # Note: this gets rid of years
    names(dat.temp)[3] <- "Y"
    summary(dat.temp)
    
    
    
    # 2. Pull the raw data that needs to be bias-corrected -- this will be called "X"
    #    -- this gets aggregated to the climatological mean right off the bat
    if(met.var == "precipf"){
		met.bias2 <- met.bias[met.bias$dataset==dat.bias, ]
		bias.ann <- aggregate(met.bias2[, c(met.var)],
							   by=list(met.bias2[, "year"]),
							   FUN=sum)
		names(bias.ann) <- c("year", "X.tot")
		summary(bias.ann)
		
		met.bias2 <- merge(met.bias2, bias.ann, all.x=T)
		
		met.bias2[,met.var] <- met.bias2[,met.var]/met.bias2$X.tot # Putting precip as fraction of the year again
		
		dat.temp2 <- aggregate(met.bias2[met.bias2$year>=yr.min & met.bias2$year<=yr.max, c(met.var, vars.met)],
                           by=list(met.bias2[met.bias2$year>=yr.min & met.bias2$year<=yr.max, "doy"]),
                           FUN=mean)
	    names(dat.temp2) <- c("doy", "X", vars.met)
	    summary(dat.temp2)
    } else {
        dat.temp2 <- aggregate(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, c(met.var, vars.met)],
                           by=list(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"]),
                           FUN=mean)
		names(dat.temp2) <- c("doy", "X", vars.met)
		summary(dat.temp2)
    }

    
    
    # 3. Merge the training & raw data together the two sets of daily means 
    #    -- this ends up pairing each daily climatological mean of the raw data with each simulation from the training data
    dat.temp <- merge(dat.temp[,], dat.temp2)
    summary(dat.temp)
    
    if(met.var == "precipf"){
      # If we don't have overlapping time periods for the annual model, need to just pari things together
      if(max(bias.ann$year)<=min(temp.ann$year)){
        bias.ann2 <- bias.ann[bias.ann$year >= yr.min & bias.ann$year <= yr.max,]
        dat.ann <- temp.ann
        dat.ann$X.tot <- bias.ann2$X.tot
      } else {
        dat.ann <- merge(temp.ann, bias.ann, all=T)
      }
    }
        
    # 4. Getting the raw ("bias") & training data for the calibration periods so we can model the anomalies and 
    # 5. Pulling all of the raw data to predict the full bias-corrected time series
    #    -- NOTE: When possible, we're pulling covariates from the ensembles we've already done.

    # The data to be bias-corrected
    if(met.var == "precipf"){
        raw.bias  <- met.bias2[met.bias2$year>=yr.min & met.bias2$year<=yr.max, c("dataset", "year", "doy", met.var)]
    	names(raw.bias) <- c("dataset", "year", "doy", "X")
    	raw.bias <- merge(raw.bias, data.frame(ind=paste0("X", 1:n)), all=T)
    } else {
        raw.bias  <- met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min       & met.bias$year<=yr.max      , c("dataset", "year", "doy", met.var)]
	    names(raw.bias) <- c("dataset", "year", "doy", "X")
	    raw.bias <- merge(raw.bias, data.frame(ind=paste0("X", 1:n)), all=T)
    }
    
    # The training data
    raw.train <- data.frame(dataset=dat.train,
                            year = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "year" ],
                            doy  = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "doy"  ],
                            X    = dat.stack$values,
                            ind  = dat.stack$ind
                            )
    
    if(met.var == "precipf"){
    	raw.ann <- aggregate(raw.train[,"X"],
    						 by=raw.train[,c("year", "ind")],
    						 FUN=sum)
    	names(raw.ann)[3] <- "X.tot"
    	
    	raw.train <- merge(raw.train, raw.ann, all.x=T)
    	raw.train$X <- raw.train$X/raw.train$X.tot
    }
    
    # The prediction data
    if(met.var == "precipf"){
        dat.pred <- data.frame(met.bias2[, c("year", "doy")],
                           X    = met.bias2[, met.var]
                           )
        dat.pred <- merge(dat.pred, bias.ann, all.x=T)

    } else {
        dat.pred <- data.frame(met.bias[met.bias$dataset==dat.bias, c("year", "doy")],
                           X    = met.bias[met.bias$dataset==dat.bias, met.var]
                           )
    }
    dat.pred <- merge(dat.pred, data.frame(ind=paste0("X", 1:n)), all=T)
    

    # Adding in covariates, propogating from what we've already downscaled as we go
    for(v.pred in vars.met[!vars.met==met.var]){
      if(v.pred %in% names(dat.out)){
        raw.train[,v.pred] <- stack(dat.out[[v.pred]]$sims[dat.out[[v.pred]]$sims$dataset==dat.train & dat.out[[v.pred]]$sims$year>=yr.min.train & dat.out[[v.pred]]$sims$year<=yr.max.train, paste0("X", 1:n)])[,1]
        raw.bias [,v.pred] <- stack(dat.out[[v.pred]]$sims[dat.out[[v.pred]]$sims$dataset==dat.bias  & dat.out[[v.pred]]$sims$year>=yr.min & dat.out[[v.pred]]$sims$year<=yr.max, paste0("X", 1:n)])[,1]
        dat.pred [,v.pred] <- stack(dat.out[[v.pred]]$sims[dat.out[[v.pred]]$sims$dataset==dat.bias, paste0("X", 1:n)])[,1]
      } else {
        raw.train[,v.pred] <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min.train & met.bias$year<=yr.max.train, v.pred]
        raw.bias [,v.pred] <- met.bias[met.bias$dataset==dat.bias & met.bias$year>=yr.min & met.bias$year<=yr.max, v.pred]
        dat.pred [,v.pred] <- met.bias[met.bias$dataset==dat.bias, v.pred]
      }
    }

    # We have several variables where we run into zero-truncation problems
    # where negative values are not possible (and many for which zero is unlikely)
    #  -- this was first noticed with humidity, but was also a problem for radiation
    #     and wind
    if(met.var %in% c("swdown", "qair", "lwdown", "wind")){
      # After some playing a square-root transformation seems to work best
      dat.temp$X  <- sqrt(dat.temp$X)
      dat.temp$Y  <- sqrt(dat.temp$Y)
      raw.train$X <- sqrt(raw.train$X)
      raw.bias$X  <- sqrt(raw.bias$X)
      dat.pred$X  <- sqrt(dat.pred$X)
    }
    # ---------
    
    
    # ---------
    # Doing the climatological bias correction
    # In all variables except precip, this adjusts the climatological means closest to the splice point
    # -- because precip is relatively stochastic without a clear seasonal pattern, a zero-inflated distribution,  
    #    and low correlation with other met variables, we'll instead model potential low-frequency patterns in
    #    the data that is to be bias-corrected.  In this instance we essentially consider any daily precip to be 
    #    an anomaly
    # ---------
    # mod.bias2 <- gam(Y ~ s(doy, by=ind, k=12) + X + ind, data=dat.temp)
    mod.bias <- gam(Y ~ s(doy, by=ind, k=6) + X + ind, data=dat.temp)
    summary(mod.bias)
    
    # Saving the mean predicted & residuals
    dat.temp$pred  <- predict(mod.bias)
    dat.temp$resid <- resid(mod.bias)
    summary(dat.temp)
    
    # Storing the model residuals to add in some extra error
    resid.bias <- resid(mod.bias)
    
    # # Checking the residuals to see if we can assume normality
    # plot(resid ~ pred, data=dat.temp); abline(h=0, col="red")
    # plot(resid ~ doy, data=dat.temp); abline(h=0, col="red")
    # hist(dat.temp$resid)
    dat.pred$pred <- predict(mod.bias, newdata=dat.pred)
    
    # For Precip we need to bias-correct the total annual preciptiation + seasonal distribution
    if(met.var == "precipf"){
    	mod.ann <- lm(Y.tot ~ X.tot + ind, data=dat.ann)
    	summary(mod.ann)
    	
    	dat.ann$pred.ann <- predict(mod.ann)
    	dat.ann$resid.ann <- resid(mod.ann)
    	
    	dat.pred$pred.ann <- predict(mod.ann, newdata=dat.pred)
    }


    # ---------
    
    # ---------
    # Modeling the anomalies
    # In most cases, this is the deviation of each observation from the climatic mean for that day (estimated using a smoother)
    #  -- This is done to adjust for differences in the anomaly distributions between data products as well as adjust for seasonal
    #     biases in products such as the GCMs (which show exaggerated seasonal trends relative to CRUNCEP & LDAS products)
    #  -- Again, precipf is handled differently because to get distributions right, we consider any precipitation event to be anomalous
    #     -- One big challenge in precip was that the GCMs essentially had a summer monsoon season, even for Harvard, which is totally bogus
    #        and would cause major problems with snow effects
    # ---------
    # Calculating the anomalies from the current means
    # The spline is a lot faster than fitting 366 means    
    raw.train$pred <- predict(mod.bias, newdata=raw.train)
    raw.bias$pred <- predict(mod.bias, newdata=raw.bias)

    # We want to look at anomalies relative to the raw expected seasonal pattern, so we need to fit training and data to be debiased separately
    anom.train <- gam(X ~ s(doy, k=6, by=ind) + ind, data=raw.train) # Need to account for the climatic differences in the simulations 
    anom.bias  <- gam(X ~ s(doy, k=6) + ind, data=raw.bias) # Note: the "bias" dataset has not corrected yet, so there should only be no variation around the values

    raw.train$anom.train <- resid(anom.train)
    raw.bias $anom.raw   <- resid(anom.bias)      
    dat.pred$anom.raw <- dat.pred$X - predict(anom.bias, newdata=dat.pred)
    # par(mfrow=c(2,1))
    # plot(anom.train~doy, data=raw.train)
    # plot(anom.raw~doy, data=raw.bias)
    # par(mfrow=c(1,1))


    # Modeling the anomalies of the other predictors 
    #  -- note: the downscaling & bias-correction of precip should have removed the monsoonal trend if there is no empirical basis for it
    #     so this should be pretty straight-forward now
    for(j in vars.met[vars.met!=met.var]){
      # print(j)
      raw.train$Q <- raw.train[,j]
      raw.bias$Q <- raw.bias[,j]
      dat.pred$Q <- dat.pred[,j]
      
      anom.train2  <- gam(Q ~ s(doy, k=6, by=ind) + ind, data=raw.train) # Need to account for the climatic differences in the simulations 
      anom.bias2  <- gam(Q ~ s(doy, k=6, by=ind) + ind, data=raw.bias) # Note: the "bias" dataset has not corrected yet, so there should only be 1 value
      
      raw.train[,paste0(j, ".anom")] <- resid(anom.train2)
      raw.bias[,paste0(j, ".anom")] <- resid(anom.bias2)
      dat.pred[,paste0(j, ".anom")] <- dat.pred$Q - predict(anom.bias2, newdata=dat.pred)
      
      rm(anom.train2, anom.bias2)
    }

    # CRUNCEP has a few variables that assume a constant pattern from 1901-1950; 
    # so we don't want to use their anomaly as a predictor otherwise we will perpetuate that less than ideal situation
    if(dat.bias=="CRUNCEP" & met.var %in% c("lwdown", "press", "wind")) raw.bias$anom.raw <- 0
  
    # Actually Modeling the anomalies
    #  -- If we have empirical data, we can pair the anomalies to find a way to bias-correct those
    #  -- If one of our datasets is a GCM, the patterns observed are just what underly the climate signal and no actual
    #     event is "real".  In this case we just want to leverage use the covariance our other met drivers to try and get
    #     the right distribution of anomalies
    if(dat.bias %in% empirical & dat.train %in% empirical){ 
      # if it's empirical we can, pair the anomalies for best estimation
      # Note: Pull the covariates from the training data to get any uncertainty &/or try to correct covariances
      #        -- this makes it mroe consistent with the GCM calculations
      dat.anom <- merge(raw.bias[,c("year", "doy", "ind", "X", "anom.raw")], raw.train[,c("year", "doy", "anom.train", "ind", vars.met[vars.met!=met.var], paste0(vars.met[vars.met!=met.var], ".anom"))])

      k=round(length(unique(dat.pred$year))/50,0)
      k=max(k, 4) # we can't have less than 4 knots
      
      # plot(anom.train ~ anom.raw, data=dat.anom)
      # abline(a=0, b=1, col="red")
      # abline(lm(anom.train ~ anom.raw, data=dat.anom), col="red", lty="dashed")
            
      # Modeling in the predicted value from mod.bias
      dat.anom$pred <- predict(mod.bias, newdata=dat.anom)
      
      if (met.var %in% c("tair", "tmax", "tmin")){
        # ** We want to make sure we do these first **
        # These are the variables that have quasi-observed values for their whole time period, 
        # so we can use the the seasonsal trend, and the observed anaomalies
        # Note: because we can directly model the anomalies, the inherent long-term trend should be preserved
        mod.anom <- gam(anom.train ~ s(doy, k=6, by=ind) + ind*anom.raw -1, data=dat.anom)
        # mod.anom <- gam(anom.raw ~ s(doy, k=4) + s(year, k=k) + ind*tmax.anom*tmin.anom -1 - ind, data=dat.pred)
      } else if(met.var %in% c("swdown", "qair")){
        # CRUNCEP swdown and qair have been vary hard to fit to NLDAS because it has a different variance for some reason, 
        # and the only way I've been able to fix it is to model the temporal pattern seen in the dataset based on 
        # its own anomalies (not ideal, but it works)
        mod.anom <- gam(anom.raw ~ s(doy, k=6, by=ind) + s(year, k=k, by=ind) + ind*tmax.anom*tmin.anom -1 , data=dat.pred)
      } else if(met.var=="precipf"){
        # Precip is really only different from the others in that I deliberately chose a more rigid seasonal pattern and we need to force the intercept
        # through 0 so we can try and reduce the likelihood of evenly distributed precipitation events
        # k=round(length(dat.pred$year)/(25*366),0)
        # k=max(k, 4) # we can't have less than 4 knots
        
        mod.anom <- gam(anom.raw ~ s(year, k=k, by=ind) + ind*(tmax.anom + tmin.anom + swdown.anom + lwdown.anom + qair.anom) -1, data=dat.pred)
      } else if(met.var %in% c("wind", "press", "lwdown")) {
        # These variables are constant in CRU pre-1950.  
        # This means that we can not use information about the long term trend OR the actual annomalies 
        # -- they must be inferred from the other met we have
        # mod.anom <- gam(anom.train ~ s(doy) + tmax.anom*tmin.anom + swdown.anom + qair.anom -1, data=dat.anom)
        # mod.anom <- gam(anom.train ~ s(doy) + ind*swdown.anom*qair -1, data=dat.anom)
        mod.anom <- gam(anom.raw ~ s(doy, k=6, by=ind) + ind*(tmin*tmax.anom + swdown.anom + qair.anom) -1, data=dat.pred)
      }      
    } else { 
      # If we're dealing with non-empirical datasets, we can't pair anomalies to come up with a direct adjustment 
      # In this case we have 2 options:
      #   1) If we've already done at least one variable, we can leverage the covariance of the met drivers we've already downscaled 
      #      to come up with a relationship that we an use to predict the new set of anomalies
      #   2) If we don't have any other variables to leverage (i.e. this is our first met variable), we incorporate both the seasonal
      #      trend (doy spline) and potential low-frequency trends in the data (year spline)
      k=round(length(unique(dat.pred$year))/50,0)
      k=max(k, 4) # we can't have less than 4 knots
      
      # vars.met <- c("tair", "tmax", "tmin", "qair", "precipf", "swdown", "press", "lwdown", "wind")
      # Vars that are at daily and we just need to adjust the variance
      # We have some other anomaly to use! that helps a lot. -- use that to try and get low-frequency trends in the past
      if(met.var %in% c("tmax", "tmin")){
        mod.anom <- gam(anom.raw ~ s(year, k=k) -1, data=dat.pred[dat.pred$ind=="X1",]) # Because all will be the same
      } else if(met.var=="precipf"){ 
        # If we're working with precipf, need to make the intercept 0 so that we have plenty of days with little/no rain
        mod.anom <- gam(anom.raw ~  s(year, k=k, by=ind) + ind*(tmax.anom*tmin.anom + swdown.anom + lwdown.anom + qair.anom) -1, data=dat.pred)  
      } else if(met.var %in% c("swdown", "lwdown")){
        # See if we have some other anomaly that we can use to get the anomaly covariance & temporal trends right
        # This relies on the assumption that the low-frequency trends are in proportion to the other met variables
        # (this doesn't seem unreasonable, but that doesn't mean it's right)
        mod.anom <- gam(anom.train ~ s(doy, k=4, by=ind) + ind*(tmax.anom*tmin.anom + qair.anom + press.anom + wind.anom) -1, data=raw.train)
      } else {
        # If we haven't already done another met product, our best shot is to just model the existing variance 
        # and preserve as much of the low-frequency cylce as possible
        # THis should be tair, tmax, tmin, qair, press, wind
        mod.anom <- gam(anom.raw ~ s(doy, k=6, by=ind) + s(year, k=k, by=ind) +ind*(tmax.anom*tmin.anom)-1, data=dat.pred)
      }
    }
    summary(mod.anom)
    # plot(mod.anom, pages=1)
 
    resid.anom <- resid(mod.anom)
    # ---------
    
    # --------
    # Predicting a bunch of potential posteriors over the full dataset
    # --------
    # Get the model coefficients
    coef.gam <- coef(mod.bias)
    coef.anom <- coef(mod.anom)
    if(met.var == "precipf") coef.ann <- coef(mod.ann)
    
    # Generate a random distribution of betas using the covariance matrix
    Rbeta <- mvrnorm(n=n, coef(mod.bias), vcov(mod.bias))
    Rbeta.anom <- mvrnorm(n=n, coef(mod.anom), vcov(mod.anom))
    if(met.var == "precipf") Rbeta.ann <- mvrnorm(n=n, coef(mod.ann), vcov(mod.ann))
    
    # Create the prediction matrix
    Xp <- predict(mod.bias, newdata=dat.pred, type="lpmatrix")
    Xp.anom <- predict(mod.anom, newdata=dat.pred, type="lpmatrix")
    if(met.var == "precipf"){
      # Linear models have a bit of a difference
      # Xp.ann <- predict(mod.ann, newdata=dat.pred, type="lpmatrix")
      
      dat.pred$Y.tot <- dat.pred$pred.ann
      mod.terms <- terms(mod.ann)
      m <- model.frame(mod.terms, dat.pred, xlev=mod.ann$xlevels)
      Xp.ann <- model.matrix(mod.terms, m, constrasts.arg <- mod.ann$contrasts)
    } 
    
    # -----
    # Simulate predicted met variables & add in some residual error
    # NOTE: Here we're assuming normal distribution of the errors, which looked pretty valid
    #       in the tests I ran when doing the intial code development
    # We do have a couple options for how to add the residual error/uncertainty back in
    # -----
    # Options for adding in residual error
    # # Option 1: Adding a constant error per time series
    #    -- This is currently used for the climatological bias-correction because we're going to assume
    #       that we've biased the mean offset in the climatology (the seasonal bias is encorporated in the 
    #       spline estimation)
    #    -- Note: Precipitation doesn't get residual error added here because that sort of bias is funneled into
    #             the anomaly model.  The error in the Rbetas should adequately represent the uncertainty in the 
    #             low-frequency trends in the data
    # # Option 2: Adding a random error to each observation
    #    -- This is used for the anomalies because they are by definition stochastic, highly unpredictable
    #    -- Note: this option currently ignores potential autocorrelation in anomalies (i.e. if 1 Jan was 
    #             unseasonably warm, odds are that the days around it weren't record-breaking cold)
    #              -- I'm rolling with this for now and will smooth some of these over in the downscaling to
    #                 subdaily data
    # # Option 3: explicitly modeling the errors in some way
    #    -- I tried this and it made my brain hurt
    # -----
    # Default option: no residual error; all error from the downscaling parameters
    sim1a <- Xp %*% t(Rbeta)  # Seasonal Climate component with uncertainty
    sim1b <- Xp.anom %*% t(Rbeta.anom) # Weather component with uncertainty
    if(met.var == "precipf") sim1c <- Xp.ann %*% t(Rbeta.ann) # Mean annual precip uncertainty

    # If we're dealing with the temperatures where there's basically no anomaly, 
    # we'll get the uncertainty subtract the multi-decadal trend out of the anomalies; not a perfect solution, but it will increase teh variability
    if(!(dat.bias %in% empirical) & (met.var %in% c("tmax", "tmin"))){
      sim1b.norm <- apply(sim1b, 1, mean) 
      sim1b <- sim1b - sim1b.norm + as.vector(dat.pred$anom.raw) # Get the range around that medium-frequency trend 
    }
    
    
    
    # Option 1: Adding a constant error per time series for the cliamte correction 
    #             (otherwise we're just doubling anomalies)
    # sim1a <- sweep(sim1a, 2, rnorm(n, mean(resid.bias), sd(resid.bias)), FUN="+")
    # if(met.var!="precipf") sim1a <- sweep(sim1a, 2, rnorm(n, mean(resid.bias), sd(resid.bias)), FUN="+") # Only apply if not working with precipf
    # sim1b <- sweep(sim1b, 2, rnorm(n, mean(resid.anom), sd(resid.anom)), FUN="+")

    # # # Option 2: Adding a random error to each observation (anomaly error)
    # if(met.var!="precipf") sim1a <- sim1a + rnorm(length(sim1a), mean(resid.bias), sd(resid.bias))
    # sim1b <- sim1b + rnorm(length(sim1b), mean(resid.anom), sd(resid.anom))
    
    # # Option 3: explicitly modeling the errors in some way
    # -----
    
    # Adding climate and anomaly together
    sim1 <- sim1a + sim1b # climate + weather = met driver!!
    
    # If we're dealing with precip, sim1 is proportion of rain
    if(met.var == "precipf") sim1 <- sim1*sim1c
    
    # Un-transform variables where we encounter zero-truncation issues
    # NOTE: Need to do this *before* we sum the components!! 
    if(met.var %in% c("swdown", "qair", "lwdown", "wind")){
      sim1 <- sim1^2
      dat.pred$X <- dat.pred$X^2
    } 
    
    
    # For preciptiation, we need to make sure we don't have constant drizzel and have 
    # at least some dry days.  To deal with this, I make the assumption that there hasn't
    # been a trend in number of rainless days over the past 1000 years and use the mean & 
    # sd of rainless days in the training data to randomly distribute the rain in the past
    if(met.var=="precipf"){
      for(j in 1:ncol(sim1)){
        for(y in min(dat.pred$year):max(dat.pred$year)){
          # Figure out which rows belong to this particular year
          rows.yr <- which(dat.pred$year==y)
          
          # Before adjusting rainless days, make sure we get rid of our negative days first
          dry <- rows.yr[which(sim1[rows.yr,j] < 0)]
          while(length(dry)>0){ # until we have our water year balanced
            for(r in 1:length(dry)){
              # Pick a year with some rain and take the rain from it
              #  -- this *should* make sure we don't get an infinite loop by making one rainless day have negative rain
              row.steal <- sample(rows.yr[which(sim1[rows.yr,j]>0)], 1) # The row we're stealing precip out of to balance the budget
              sim1[row.steal,j] <- sim1[row.steal,j] + sim1[dry[r],j]
              sim1[dry[r],j] <- 0
            }
            dry <- rows.yr[which(sim1[rows.yr,j] < 0)] # update our dry days
          }
          
          n.now <- round(rnorm(1, mean(rainless), sd(rainless)), 0) 
          cutoff <- quantile(sim1[rows.yr, j], n.now/366)
          
          # Figure out which days are currently below our cutoff and randomly distribute 
          # their precip to days that are not below the cutoff (this causes a more bi-modal 
          # distribution hwere dry days get drier), but other options ended up with either 
          # too few rainless days because of only slight redistribution (r+1) or buildup 
          # towards the end of the year (random day that hasn't happened)
          dry <- rows.yr[which(sim1[rows.yr,j] < cutoff)]
          
          # Figure out where to put the extra rain; allow replacement for good measure
          wet <- sample(rows.yr[!rows.yr %in% dry], length(dry), replace=T)
          
          # Go through and randomly redistribute the precipitation to days we're not designating as rainless
          # Note, if we don't loop through, we might lose some of our precip
          for(r in 1:length(dry)){
            sim1[wet[r],j] <- sim1[wet[r],j] + sim1[dry[r],j]
            sim1[dry[r],j] <- 0
          }
        }
      }
    }
  
    # Randomly pick one from this meta-ensemble to save
    # this *should* be propogating uncertainty because we have the ind effects in all of the models and we're randomly adding as we go
    sim.final <- data.frame(array(dim=c(nrow(met.bias[met.bias$dataset==dat.bias,]), n)))
    names(sim.final) <- paste0("X", 1:n)
    for(ens in 1:n){
      sim.final[,ens] <- sim1[which(dat.pred$ind==paste0("X", ens)),sample(1:n,1)]
    }
  
    # Aggregate the data to something that's easier to deal with in graphing etc.
    # If we're working with precip, lets turn everything back into numbers rather than probabilities
    if(met.var == "precipf"){
      dat.pred$X <- dat.pred$X * dat.pred$X.tot
      dat.pred$pred <- dat.pred$pred * dat.pred$pred.ann
      dat.pred$anom.raw <- dat.pred$anom.raw * dat.pred$X.tot
    }
    
    dat.pred <- aggregate(dat.pred[,c("X", "pred", "anom.raw")],
                          by=dat.pred[,c("year", "doy")],
                          FUN=mean)
      
    dat.pred$mean <- apply(sim.final, 1, mean)
    dat.pred$lwr  <- apply(sim.final, 1, quantile, 0.025)
    dat.pred$upr  <- apply(sim.final, 1, quantile, 0.975)
    dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))
    summary(dat.pred)
    
    # Formatting the output
    dat.sims <- data.frame(dataset=dat.bias, met=met.var, dat.pred[,c("year", "doy", "time")])
    dat.sims <- cbind(dat.sims, sim.final)
    # ---------

    # ---------
    # Doing some exploratory graphing    
    # ---------
    # Plotting the observed and the bias-corrected 95% CI
    png(file.path(path.out, "BiasCorrect_QAQC", paste0("BiasCorr_", met.var, "_", dat.bias, "_day.png")))
    print(
    ggplot(data=dat.pred[dat.pred$year>=(yr.max-5) & dat.pred$year<=(yr.max-3),]) +
      geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
      geom_line(aes(x=time, y=mean), color="red", size=0.5) +
      geom_line(aes(x=time, y=X), color='black', size=0.5) +
      theme_bw()
    )
    dev.off()
    
    # Plotting a few random series to get an idea for what an individual pattern looks liek
    stack.sims <- stack(data.frame(sim.final))
    stack.sims[,c("year", "doy", "time")] <- dat.pred[,c("year", "doy", "time")]

    png(file.path(path.out, "BiasCorrect_QAQC", paste0("BiasCorr_", met.var, "_", dat.bias, "_day2.png")))
    print(
      ggplot(data=stack.sims[stack.sims$ind %in% paste0("X", sample(1:n, 3, replace=F)) & stack.sims$year>=(yr.max-10) & stack.sims$year<=(yr.max-3),]) +
        geom_line(aes(x=time, y=values, color=ind), size=0.2, alpha=0.8) +
        theme_bw()
    )
    dev.off()

    # Looking tat the annual means over the whole time series to make sure we're getting decent interannual variability
    dat.yr <- aggregate(dat.pred[,c("X", "mean", "lwr", "upr")],
                        by=list(dat.pred$year),
                        FUN=mean)
    names(dat.yr)[1] <- "year"
    
    png(file.path(path.out, "BiasCorrect_QAQC", paste0("BiasCorr_", met.var, "_", dat.bias, "_annual.png")))
    print(
    ggplot(data=dat.yr[,]) +
      geom_ribbon(aes(x=year, ymin=lwr, ymax=upr), fill="red", alpha=0.5) +
      geom_line(aes(x=year, y=mean), color="red", size=0.5) +
      geom_line(aes(x=year, y=X), color='black', size=0.5) +
      theme_bw()
    )
    dev.off()
    # --------
    
    rm(mod.bias, anom.train, anom.bias, mod.anom)
    
    # --------
    # Storing the output
    # --------
    cols.bind <- c("year", "doy", "X", "anom.raw", "mean", "lwr", "upr", "time")
    dat.out[[met.var]]$ci   <- rbind(dat.out[[met.var]]$ci, data.frame(dataset=dat.bias, met=met.var, dat.pred[,c("year", "doy", "X", "anom.raw", "mean", "lwr", "upr", "time")]))
    dat.out[[met.var]]$sims <- rbind(dat.out[[met.var]]$sims, data.frame(dat.sims))
    # --------
    # --------------------
  } # End dataset loop
} # end variable loop
dat.out$met.bias <- met.bias
return(dat.out)
} # end function