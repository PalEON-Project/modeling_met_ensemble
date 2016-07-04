bias.correct <- function(met.bias, vars.met, dat.train, GCM, yrs.cal, n, path.out){
# Inputs:
# met.bias   = data frame of all daily met data
# vars.met  = vector of variables to downscale, **in the order they should be done**
# dat.train = name of the training dataset
# yrs.cal - data frame of the series to be bias-corrected with three columns:
#             1. dataset name, 
#             2. first year for calibration
#             3. last year for calibration
# n         = number of ensemble members

# Load required libraries
library(MASS)
library(mgcv)
library(ggplot2)
set.seed(2509)
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
  # define what met variable we're working with
  met.var = vars.met[v]
  dat.train <- dat.train.orig
  
  print(met.var)

  # Getting building covariance matrices leveraging what we've already gapfilled
  if(v>1) vars.pred <- vars.met[1:(v-1)]
  # because I don't feel good about qair as a predictor
  vars.pred <- vars.pred[!vars.pred=="qair"]

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
    
    print(paste0(" -- ", dat.bias))
    
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
#     # Set up a data frame with the met that we've already processed by index
#     if(length(vars.pred)>0){
#       for(j in vars.pred){
#         met.df0 <- data.frame(year = dat.out[[j]]$sims[dat.out[[j]]$sims$dataset==dat.bias, "year" ],
#                               doy  = dat.out[[j]]$sims[dat.out[[j]]$sims$dataset==dat.bias, "doy"  ],
#                               stack(dat.out[[j]]$sims[dat.out[[j]]$sims$dataset==dat.bias, paste0("X", 1:n)])
#                               )
#         names(met.df0)[which(names(met.df0)=="values")] <- j
#         
#         if(j == vars.pred[1]){
#           met.df <- met.df0
#         } else {
#           met.df[,j] <- met.df0[,j]
#         }      
#       }
#     } else {
#       met.df <- data.frame(ind="X1")
#     }   
    
#     # insert 0s into vars that haven't been done yet
#     for(j in vars.met[!vars.met %in% names(met.df)]){
#       met.df[,j] <- 0
#     }
    
    # A. Calibration -- Aggregate so we're looking at the climatic mean for the first part
    dat.temp0 <- data.frame(year = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "year" ],
                            doy  = dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train  & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, "doy"  ],
                            stack(dat.out[[met.var]]$sims[dat.out[[met.var]]$sims$dataset==dat.train & dat.out[[met.var]]$sims$year>=yr.min.train & dat.out[[met.var]]$sims$year<=yr.max.train, paste0("X", 1:n)])
                            )
        
    dat.temp <- aggregate(dat.temp0$values, by=dat.temp0[,c("doy", "ind")], FUN=mean)
    names(dat.temp)[3] <- "Y"
    summary(dat.temp)
    
    # Note: because we're dealing with climatological means in this step, we don't need to pull the full distribution
    dat.temp2 <- aggregate(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, c(met.var, vars.met)],
                           by=list(met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min & met.bias$year<=yr.max, "doy"]),
                           FUN=mean)
    names(dat.temp2) <- c("doy", "X", vars.met)
    summary(dat.temp2)
    
    # merging together the two sets of daily means 
    # Note: working off of the climatological means because we'll model the anomalies separately    
    dat.temp <- merge(dat.temp[,], dat.temp2)
    summary(dat.temp)
        
    # getting the raw data for the calibration periods so we can model the anomalies
    raw.train <- met.bias[met.bias$dataset==dat.train & met.bias$year>=yr.min.train & met.bias$year<=yr.max.train, c("dataset", "year", "doy", met.var, vars.met)]
    raw.bias  <- met.bias[met.bias$dataset==dat.bias  & met.bias$year>=yr.min       & met.bias$year<=yr.max      , c("dataset", "year", "doy", met.var, vars.met)]
    names(raw.train) <- c("dataset", "year", "doy", "X", vars.met)
    names(raw.bias) <- c("dataset", "year", "doy", "X", vars.met)
    
    #  # Merge in the full permutation of anomalies
    #  raw.bias <- merge(raw.bias, met.df)
    
    # The prediction data
    dat.pred <- data.frame(met.bias[met.bias$dataset==dat.bias, c("year", "doy", vars.met)],
                           X    = met.bias[met.bias$dataset==dat.bias, met.var]
                           )
    # dat.pred <- merge(dat.pred, met.df)
    
    # 0 out unnecessary predictors in all data frames
    dat.temp [,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    raw.train[,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    raw.bias [,vars.met[which(!vars.met %in% vars.pred)]] <- 0
    dat.pred [,vars.met[which(!vars.met %in% vars.pred)]] <- 0

    # If we're working with qair, we need to transform to prevent 0 humidity 
    # (which isn't possible otherwise everything would die)
    # lwdown is also having some distribution problems
    if(met.var %in% c("qair", "lwdown", "wind")){
      # 0 out unnecessary predictors in all data frames
      dat.temp [,"X"] <- sqrt(dat.temp [,"X"])
      dat.temp [,"Y"] <- sqrt(dat.temp [,"Y"])
      raw.train[,"X"] <- sqrt(raw.train[,"X"])
      raw.bias [,"X"] <- sqrt(raw.bias [,"X"])
      dat.pred [,"X"] <- sqrt(dat.pred [,"X"])
    }
    # ---------
    
    
    # ---------
    # Doing the climatological bias correction
    # ---------
    # # Plotting the correlation between dat.train & dat.bias
    # plot(Y ~ X, data=dat.temp) 
    # abline(a=0, b=1, col="red")
    # abline(lm(Y ~ X, data=dat.temp), col="red", lty="dashed")
    
    # Running a model
    # If we're looking at precip, make sure there's no seasonal trend unless it's clear in everything else
    if(met.var == "precipf"){
      mod.bias <- gam(Y ~ X + tmax + tmin + precipf + swdown + qair + lwdown + press + wind -1, data=dat.temp)
    } else {
      mod.bias <- gam(Y ~ s(doy) + X + tmax + tmin + precipf + swdown + qair + lwdown + press + wind, data=dat.temp)
    }
    summary(mod.bias)
    # plot(mod.bias, pages=1)
    
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
    # ---------
    
    # ---------
    # Modeling the anomalies with a similar framework
    # ---------
    # Calculating the anomalies from the current means
    # The spline is a lot faster than fitting 366 means
    
    anom.train <- gam(X ~ s(doy), data=raw.train)
    anom.bias  <- gam(X ~ s(doy), data=raw.bias)
    raw.train$anom.train <- resid(anom.train)
    raw.bias $anom.raw   <- resid(anom.bias)
    
    dat.pred$anom.raw <- dat.pred$X - predict(anom.bias, newdata=dat.pred)
      
    # par(mfrow=c(2,1))
    # plot(anom.train~doy, data=raw.train)
    # plot(anom.raw~doy, data=raw.bias)
    # par(mfrow=c(1,1))


    # Modeling the anomalies of the other predictors
    for(j in vars.met){
      raw.train$Q <- raw.train[,j]
      raw.bias$Q <- raw.bias[,j]
      dat.pred$Q <- dat.pred[,j]
      
      if(j %in% vars.pred){ 
        anom.train2  <- gam(Q ~ s(doy), data=raw.train)
        raw.train[,paste0(j, ".anom")] <- resid(anom.train2)
        
        anom.bias2  <- gam(Q ~ s(doy), data=raw.bias)
        raw.bias[,paste0(j, ".anom")] <- resid(anom.bias2)
        dat.pred[,paste0(j, ".anom")] <- dat.pred$Q - predict(anom.bias2, newdata=dat.pred)
      } else {
        raw.train[,paste0(j, ".anom")] <- 0
        raw.bias[,paste0(j, ".anom")] <- 0
        dat.pred[,paste0(j, ".anom")] <- 0
      }
    }

    # CRUNCEP has a few variables that assume a constant pattern from 1901-1950; try modeling those purely from the covariates
    if(dat.bias=="CRUNCEP" & met.var %in% c("lwdown", "press", "wind")) raw.bias$anom.raw <- 0
  
    # Modeling the anomalies
    # Figure out if we can line up the anomalies for better accuracy
    if(dat.bias %in% empirical & dat.train %in% empirical){ 
      # if it's empirical we can, pair the anomalies for best estimation
      dat.anom <- merge(raw.bias[,c("year", "doy", "X", "anom.raw", vars.met, paste0(vars.met, ".anom"))], raw.train[,c("year", "doy", "anom.train")])

      # plot(anom.train ~ anom.raw, data=dat.anom)
      # abline(a=0, b=1, col="red")
      # abline(lm(anom.train ~ anom.raw, data=dat.anom), col="red", lty="dashed")
            
      # Modeling in the predicted value from mod.bias
      dat.anom$pred <- predict(mod.bias, newdata=dat.anom)
      
      # Because anom.train and anom.raw shoudl be correlated, we can use that as a predictor
      mod.anom <- gam(anom.train ~ s(doy) + anom.raw + tmax.anom + tmin.anom + precipf.anom + swdown.anom + qair.anom + lwdown.anom + press.anom + wind.anom, data=dat.anom)
      
      # CRUNCEP swdown has been vary hard to fit to NLDAS because it has a different variance, so here we're just 
      # going to model the pattern seen in the dataset based on its own anomalies
      if(met.var=="swdown"){
        mod.anom <- gam(anom.raw ~ s(doy) + s(year, k=3) + tmax.anom + tmin.anom, data=dat.pred)
      }
      
    } else { 
      # If we can't, then we'll be using the two datasets separately and don't need to bind them together      

      # See if we have some other anomaly that we can use to get the anomaly covariance right
      if(max(raw.bias[,vars.met])>0){
        # We have some other anomaly to use! that helps a lot. -- use that to try and get low-frequency trends in the past
        if(met.var=="precipf"){ 
          # If we're working with precipf, need to make the intercept 0 so that we have plenty of days without rain
          mod.anom <- gam(anom.train ~ tmax.anom + tmin.anom + precipf.anom + swdown.anom + qair.anom + lwdown.anom + press.anom + wind.anom - 1, data=raw.train)  
        } else {
          mod.anom <- gam(anom.train ~ s(doy) + tmax.anom + tmin.anom + precipf.anom + swdown.anom + qair.anom + lwdown.anom + press.anom + wind.anom, data=raw.train)  
        }
      } else {
        # If we haven't already done another met product, our best shot is to just model the variance adding in as much of the low-frequency cylce as possible
        k=round(length(dat.pred$year)/(25*366),0)
        k=max(k, 4) # we can't have less than 4 knots
        mod.anom <- gam(anom.raw ~ s(doy) + s(year, k=k) + pred, data=dat.pred)
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
    sim1a <- Xp %*% t(Rbeta)  # Climate component with uncertainty
    sim1b <- Xp.anom %*% t(Rbeta.anom) # Weather component with uncertainty
    
    sim1 <- sim1a + sim1b # climate + weather
    
    # Un-transform qair
    if(met.var %in% c("qair", "lwdown", "wind")){
      sim1 <- sim1^2
      dat.pred$X <- dat.pred$X^2
    }


    # # Option 1: Adding a constant error per time series
    # res <- rnorm(n, mean(dat.temp$resid), sd(dat.temp$resid))
    # sim1 <- sweep(sim1, 2, res, FUN="+")
    #  sim1 <- sweep(Xp %*% t(Rbeta), 2, rnorm(n, mean(resid.bias), sd(resid.bias)), FUN="+") + 
    #            sweep(Xp.anom %*% t(Rbeta.anom),2, rnorm(n, mean(resid.bias), sd(resid.bias)), FUN="+")

    # # Option 2: Adding a random error to each observation
    # sim1 <- sim1 + array(rnorm(length(sim1), mean(dat.temp$resid), sd(dat.temp$resid)), dim=dim(sim1))
    
    # # Option 3: Adding a day-specific error (although this should be part of the smoother)
    # -----
    # Get rid of any negative values
    sim1[sim1<0] = 0
  
    # NOTE: Don't need to add in anomalies because training & fitting data are both observed annual deviations!
    dat.pred$mean <- apply(sim1, 1, mean)
    dat.pred$lwr  <- apply(sim1, 1, quantile, 0.025)
    dat.pred$upr  <- apply(sim1, 1, quantile, 0.975)
    summary(dat.pred)
    
    # CRUNCEP swdown is producing impossiblly low values

    # Aggregating across iterations in dat.pred
    
    # Doing som exploratory graphing
    dat.pred$time <- as.Date(dat.pred$doy, origin=as.Date(paste0(dat.pred$year, "-01-01")))
    
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
    
    # --------
    # Storing the output
    # --------
    dat.sims <- data.frame(dataset=dat.bias, met=met.var, dat.pred[,c("year", "doy", vars.met, "time")])
    dat.sims <- cbind(dat.sims, sim1)
  
    # aggregate dat.pred and dat.sims to get mean of the n iterations of the pervious data
    cols.pred <- names(dat.out[[met.var]]$ci)[3:ncol(dat.out[[met.var]]$ci)] # The columns we're going to bind on
    dat.pred <- aggregate(dat.pred[,cols.pred[3:length(cols.pred)]], by=dat.pred[,c("year", "doy")], FUN=mean, na.rm=T)
    
    dat.sims <- aggregate(dat.sims[,5:ncol(dat.sims)], by=dat.sims[,c("dataset", "met", "year", "doy")], FUN=mean, na.rm=T)
    
    dat.out[[met.var]]$ci   <- rbind(dat.out[[met.var]]$ci, data.frame(dataset=dat.bias, met=met.var, dat.pred[,cols.pred]))
    dat.out[[met.var]]$sims <- rbind(dat.out[[met.var]]$sims, data.frame(dat.sims))
    
    # overwriting the raw with the mean bias-corrected output 
    met.bias[met.bias$dataset==dat.bias, met.var] <- dat.pred$mean
    # --------
    # --------------------
  }
}
dat.out$met.bias <- met.bias
return(dat.out)
}