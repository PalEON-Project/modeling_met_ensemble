# Function to calculate PDSI 
# Translated from Ben Cook's matlab code
# Author: Christy Rollinson, crollinson@mortonarb.org
#
# Description & Notes
# This calculates PDSI using from Temp, Precip, and PE
#  - current implementation is monthly, but shouldn't be hard to tweak to 
#    get daily
# PE Options:
# 1. Thornthwaite: 
#     - uses precip, awc, and day length adjustment factor
#     - awc can be provided or calculated based on soil texture
#     - day length adjustment currently based on hard-coded tables, but 
#       can be adjusted to be provided (based on SWdown or other dataset)
# 2. Penman-Montheith (not implemented)
#     - PM PE is more robust, but it requires Rnet, which we can't get from
#       met drivers alone
# 
# Unit Notes:
#  - Outputs DO need to be in INCHES because of coefficients; patterns are similar with metric, 
#    but not identical; Temperature is only used in PE, which requires celcius, so that doesn't 
#    need to be converted; 
#  - The only place that needs to be adjusted to work with daily data is PE.thorntwaite
#
# --------------
# Inputs:
# --------------
#  1. datmet: list, length=2  
#     ** redo to work with daily?
#     ** redo units to metric??
#     1. Precip = precipitation; data frame; units=total units per observation; 
#                 preferred units: "mm" (per month)
#                 dim=c(nyr, ntime); row names = years
#     2. Temp   = temperature; data frame; 
#                 units="Fahrenheit"; 
#                 dim=c(nyr, ntime); row names = years
#  2. datother: list, length=4
#     1. pdsi.fun = file path to where the necessary PDSI helper functions are
#     2. metric = T/F; is moisture in metric (mm) or inches?
#     3. lat = latitude; numeric
#              units="decimal degrees"
#              length=1
#     4. watcap = water capacity; list, length=2; units="mm"  (volumetric * depth)
#         ** Note: to get to units multiply volumetric awc by depth
#            for PalEON drivers, topsoil = 1-30 cm; subsoil = 30-depth
#         awcs = awc surface layer (standard = 1"; paleon drivers: 30 cm (or depth-1))
#         awcu = awc underlying layer (standard = 5"; paleon drivers: depth-30 cm)
#     5. yrs.calib = window for normals & calibrations
#     6. dayz = lookup table for percentage of possible sunshine
#     7. daylength = provided values of day length (in hours); alternative to using dayz
#  5. method.PE = method of potential evapotranspiration
#     - "Thornthwaite" (default)
#     - "Penman-Monteith" (not implemented)
#     - units = precip units per obs time (e.g. mm/day or mm/mo; written for in/mo)
#  6. snow = what to do with snow option
#     - NULL = all ppt as rain regardless of temp (null)
#     - "redistribute" = convert rain to snow & re-distribute based on 
#        temperature & Tcrit (see snowopts)
#  3. snowopts: list, length=3
#     1. Tcrit = temperature to convert ppt to snow if kopt[[1]] = 1; default = 32F
#     2. mons1 = months to convert precip to snow (typically oct to may)
#     3. melttbl = lookup table to convert snow accumulation & snow melt
#  6. penopts = list; length = 2 (options for using Penman-Montheith PE)
#     1. wind
#        1 = monthly timeseries input
#        2 = monthly means only
#     2. relative humidity
#        1 = monthly time series
#        2 = monthly means only
#  7. datpen = list, length=2; (input for Penman-Montheith PE)
#     1. wind speed (see penopts for dimensions)
#     2. relative humidity (see penopts for dimensions)
# --------------

# --------------
# Outputs (datout)
# --------------
# datout: list; length = 8
#  1. Z  = Z index (unitless); data frame; dim=c(nyr, 13)
#  2. X  = PDSI value (unitless); data frame; dim=c(nyr, 13)
#  3. XM = modified PDSI (unitless); data frame; dim=c(nyr, 13)
#  4. W  = avg soil moist (unit="in"); data frame; dim=c(nyr, 13)
#  5. RO = monthly runoff (unit=??); dataframe; dim=c(nyr, 13)
#  6. S1 = effective ppt (unit="in"); array; dim=c(nyr, 13, 10)
#        = max(0, p - f*pe)
#  7. P  = precip input; data frame; dim=c(nyr, 13)
#  8. T  = temperature input; data frame; dim=c(nyr, 13)
# --------------

# --------------
# Functions Called
# --------------
# 1. pdsix - calculates X values rom
# 2. pethorn - calculate Thornthwaite PE
# 3. soilmoi1 - does soil moisture accounting
# --------------

# --------------
# Workflow
# --------------
#  1. Check input, calculate AWC if needed, unit conversions, make pointers
#  2. Compute PE (potential evapotranspiration)
#  3. Redistribute monthly P if needed (not implemented)
#  4. Calculate soil moisture using saturated initial conditions
#     4.1. Get calib preiod PE, P, T, & calculate 'normals'
#     4.2. build synthetic 10-year time series for P & PE
#     4.3. Calculate soil moisture with synthetic time series
#     4.4. Extract initial conditions
#  5. Calculate water balance normals & coefficients for PDSI calculations
#     5.1. Calculate soil moisture using the calibration period & initial soil moisture
#     5.2. Calculate means
#     5.3. Calculate coefficients: 
#          5.3.a. alpha 
#          5.3.b. beta
#          5.3.c. gamma
#          5.3.d. delta
#  6. Soil moisture accoutning on full time series
#     6.1 run moisture calculation
#     6.2 extract values
#     6.3 calculate effective preciptiation (precip - fraction of PE)
#     6.4 Calculate CAFEC precip (pg. 14, eqns 10-14)
#     6.5 Calculate excesses and deficits
#  7. Calculate K' and and weighted-K values
#     7.1 Calculate K' for each time step (p25, eqn 26)
#     7.2 Calculate weighted departures for each month and its sum
#     7.3 Caluclate final K values
#  8. Calculate PDSI 
#     8.1 Calculate Z-values from K & d
#     8.2 Calculate PDSI using pdsix function
#  9. Format & return output
# --------------


pdsi1 <- function(datmet, datother, metric=F, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL){

  # ------------------------------------------
  #  1. Check input, calculate AWC if needed, unit conversions, make pointers
  # ------------------------------------------
  # Do more robust input checks when I have time; in the meanwhile just a couple warnings
  # about things I haven't implemented yet
  if(method.PE!="Thornthwaite") stop("invalid method.PE! Only Thornthwaite available currently")
  if(!is.null(snow)) stop("Snow redistribution not implemented yet. Please set to NULL for now")
  if(any(dim(datmet$Temp) != dim(datmet$Precip))) stop("Temperature & Precipitation time series do not align!")
  
  # Extract some constants from datother
  Temp <- as.matrix(datmet$Temp)
  Precip <- as.matrix(datmet$Precip)
  
  yrs <- as.numeric(row.names(Temp))
  nyrs <- length(yrs)
  yrs.calib <- datother$yrs.calib
  rows.calib <- which(yrs>=yrs.calib[1] & yrs<=yrs.calib[2])
  lat     <- datother$lat
  awcs    <- datother$watcap$awcs
  awcu    <- datother$watcap$awcu
  dayz    <- datother$dayz
  dayfact <- datother$dayfact

  # Convert Precip from mm to inches
  # Do unit conversions on moisture if necessary
  #  Assumes that if not in inches, it's in mm
  if(metric==T){
    Precip <- Precip/25.4
    awcs   <- awcs/25.4
    awcu   <- awcu/25.4
    
    C2F <- function(x){x*9/5 + 32}
    Temp   <- C2F(Temp)
  }

  # # Daylength stuff
  # library(R.matlab)
  # ------------------------------------------
  
  # ------------------------------------------
  #  2. Compute PE (potential evapotranspiration) for full time series
  #     Output Units: mm/time
  # ------------------------------------------
  if(method.PE=="Thornthwaite"){ 
    if(is.null(dayfact)) dayz <- R.matlab::readMat("PDSI_fromBenCook/PDSICODE/daylennh.mat")$dayz
    # dayl <- NULL
    
    # source(file.path(pdsi.fun, "PE.thornthwaite.R"))
    if(ncol(Temp)==12) timestep = "monthly"
    # dayfact = calc.dayfact(timestep=timestep, daylength=dayl, lat=lat, dayz=dayz)
    # dayfact = calc.dayfact(timestep="daily", daylength=dayl, lat=lat, dayz=dayz)
    PE <- PE.thorn(Temp, yrs.calib, lat, dayz=dayz, dayfact=dayfact, celcius=F)
  }
  row.names(PE) <- row.names(Temp)
  # ------------------------------------------

  # ------------------------------------------
  #  3. Redistribute monthly P if needed (not implemented)
  # ------------------------------------------
  if(!is.null(snow)) warning("ignoring snow redistribution; not yet implemented")
  # ------------------------------------------
  
  # ------------------------------------------
  #  4. Calculate soil moisture using saturated initial conditions
  #     4.1. Get calib preiod PE, P, T, & calculate 'normals'
  #     4.2. build synthetic 10-year time series for P & PE
  #     4.3. Calculate soil moisture with synthetic time series
  #          output units: unit/time; default PE is mm/time
  #     4.4. Extract initial conditions
  # ------------------------------------------
  # 4.1 Getting the mean Temp, Precip, & PE for the calibration period
  Pbar  = colMeans(Precip[rows.calib,])
  Tbar  = colMeans(Temp  [rows.calib,])
  PEbar = colMeans(PE    [rows.calib,])
  
  # 4.2 Building a synthetic n-yr time series for P & PE norms fo soil moisture
  # NOTE: We want these to be vectors!
  nsyn = 10
  Psyn  <- rep(Pbar , nsyn)
  PEsyn <- rep(PEbar, nsyn)
  
  # 4.3 Calculate soil moisture from synthetic time series
  moist1 <- calc.soilmoist(p=Psyn, pe=PEsyn, awcs=awcs, awcu=awcu, ssgo=awcs, sugo=awcu)
  
  # 4.4 Extract initial conditions
  ssgo <- matrix(moist1$ss1, nrow=nsyn, byrow=T)[nsyn,1]
  sugo <- matrix(moist1$su1, nrow=nsyn, byrow=T)[nsyn,1]
  
  rm(moist1)
  # ------------------------------------------
  
  # ------------------------------------------
  #  5. Calculate water balance normals & coefficients for PDSI calculations
  #      5.1. Calculate soil moisture using the calibration period & initial soil moisture
  #      5.2. Calculate means
  #      5.3. Calculate coefficients: 
  #           5.3.a. alpha 
  #           5.3.b. beta
  #           5.3.c. gamma
  #           5.3.d. delta
  # ------------------------------------------
  # creating a vector of conditions from the normals period
  # The default is to go by columns when creating a vector, so the easiest hack
  # I could figure out is to just transpose the matrix
  P.norm  <- as.vector(as.matrix(t(Precip[rows.calib,])))
  PE.norm <- as.vector(as.matrix(t(PE    [rows.calib,])))
  
  # 5.1 Calculate soil moisture for the calibration period
  moist.norm <- calc.soilmoist(p=P.norm, pe=PE.norm, 
                               awcs=awcs, awcu=awcu, 
                               ssgo=ssgo, sugo=sugo)
  
  rm(P.norm, PE.norm) # cleaning up to reduce our memory needs
  
  # reformatting output back into arrays
  ET    <- matrix(moist.norm$et   , nrow=length(rows.calib), byrow=T)
  R     <- matrix(moist.norm$r    , nrow=length(rows.calib), byrow=T)
  PR    <- matrix(moist.norm$pr   , nrow=length(rows.calib), byrow=T)
  RO    <- matrix(moist.norm$ro   , nrow=length(rows.calib), byrow=T)
  PRO   <- matrix(moist.norm$pro  , nrow=length(rows.calib), byrow=T)
  LOSS  <- matrix(moist.norm$loss , nrow=length(rows.calib), byrow=T)
  PLOSS <- matrix(moist.norm$ploss, nrow=length(rows.calib), byrow=T)
  W     <- matrix(moist.norm$smean, nrow=length(rows.calib), byrow=T) # mean soil moisture
  # SS2   <- matrix(moist.norm$ss2, nrow=length(rows.calib), byrow=T) # surface ending soil moist
  # SU2   <- matrix(moist.norm$su2, nrow=length(rows.calib), byrow=T) # under ending soil moist
  
  # 5.2 Finding the means for the calibration period
  ETbar    <- apply(ET   , 2, mean)
  Rbar     <- apply(R    , 2, mean)
  PRbar    <- apply(PR   , 2, mean)
  RObar    <- apply(RO   , 2, mean)
  PRObar   <- apply(PRO  , 2, mean)
  LOSSbar  <- apply(LOSS , 2, mean)
  PLOSSbar <- apply(PLOSS, 2, mean)
  
  # 5.3 Calculating the coefficients
  # 5.3.a. alpha
  alpha <- rep(NA, length=length(ETbar))
  L1 <- which(PEbar==0 & ETbar==0)
  L2 <- which(PEbar==0 & ETbar!=0)
  L3 <- which(!1:length(alpha) %in% c(L1, L2))
  
  if(length(L1)>0) alpha[L1] = 1
  if(length(L2)>0) alpha[L2] = 0
  if(length(L3)>0) alpha[L3] = ETbar[L3] / PEbar[L3]
  
  # 5.3.b. beta
  beta <- rep(NA, length=length(Rbar))
  L1 <- which(PRbar==0 & Rbar==0)
  L2 <- which(PRbar==0 & Rbar!=0)
  L3 <- which(!1:length(beta) %in% c(L1, L2))
  
  if(length(L1)>0) beta[L1] = 1
  if(length(L2)>0) beta[L2] = 0
  if(length(L3)>0) beta[L3] = Rbar[L3] / PRbar[L3]
  
  # 5.3.c. gamma
  gamma <- rep(NA, length=length(RObar))
  L1 <- which(PRObar==0 & RObar==0)
  L2 <- which(PRObar==0 & RObar!=0)
  L3 <- which(!1:length(gamma) %in% c(L1, L2))
  
  if(length(L1)>0) gamma[L1] = 1
  if(length(L2)>0) gamma[L2] = 0
  if(length(L3)>0) gamma[L3] = RObar[L3] / PRObar[L3]
  
  # 5.3.d. delta
  delta <- rep(NA, length=length(PLOSSbar))
  L1 <- which(PLOSSbar==0)
  L3 <- which(!1:length(delta) %in% L1)
  
  if(length(L1)>0) delta[L1] = 0
  if(length(L3)>0) delta[L3] = LOSSbar[L3] / PLOSSbar[L3]
  # ------------------------------------------

  # ------------------------------------------
  #  6. Soil moisture accoutning on full time series
  #     6.1 run moisture calculation
  #     6.2 extract values
  #     6.3 calculate effective preciptiation (precip - fraction of PE)
  #     6.4 Calculate CAFEC precip (pg. 14, eqns 10-14)
  #     6.5 Calculate excesses and deficits
  # ------------------------------------------
  P.temp  <- as.vector(as.matrix(t(Precip[,])))
  PE.temp <- as.vector(as.matrix(t(PE    [,])))
  
  # 6.1 run moisture calculation
  soilmoist <- calc.soilmoist(p=P.temp, pe=PE.temp, 
                              awcs=awcs, awcu=awcu, 
                              ssgo=ssgo, sugo=sugo)
  
  # 6.2 extract values
  ET    <- matrix(soilmoist$et   , nrow=nrow(Precip), byrow=T) # evapotranspiration
  R     <- matrix(soilmoist$r    , nrow=nrow(Precip), byrow=T) # recharge
  PR    <- matrix(soilmoist$pr   , nrow=nrow(Precip), byrow=T) # potential recharge
  RO    <- matrix(soilmoist$ro   , nrow=nrow(Precip), byrow=T) # runoff
  PRO   <- matrix(soilmoist$pro  , nrow=nrow(Precip), byrow=T) # potential runoff
  LOSS  <- matrix(soilmoist$loss , nrow=nrow(Precip), byrow=T) # loss
  PLOSS <- matrix(soilmoist$ploss, nrow=nrow(Precip), byrow=T) # potential loss
  W     <- matrix(soilmoist$smean, nrow=nrow(Precip), byrow=T) # mean soil moisture
  # SS2     <- matrix(soilmoist$ss2, nrow=nrow(Precip), byrow=T) # mean soil moisture
  # SU2     <- matrix(soilmoist$su2, nrow=nrow(Precip), byrow=T) # mean soil moisture
  dimnames(W)[[1]] <- row.names(Precip)
 
  # # Calculating a running deficit for diagnostics
  # d <- P.temp - PE.temp
  # d.run <- rep(NaN, length(d))
  # d.run[1] <- d[1]
  # for(i in 2:length(d)){
  #   d.run[i] <- d[i] + d.run[i-1]
  # }
  # d.mat <- matrix(d   , nrow=nrow(Precip), byrow=T)
  # d.run.mat <- matrix(d.run   , nrow=nrow(Precip), byrow=T)
  # plot(rowSums(d.run.mat), type="l")
  # 
  # plot(rowSums(Precip), type="l")
  # plot(rowSums(PE), type="l")
  # plot(rowMeans(W), type="l")
  
   
  # 6.3 calculate effective preciptiation (precip - fraction of PE)
  S1 <- array(dim=c(dim(Precip), 10))
  for(i in 1:10){
    nfract = i/10
    dat.temp <- as.matrix(Precip-nfract*PE)
    dat.temp[dat.temp<0] <- 0 # Has a min of 0 (no effective precip; no negative)
    S1[,,i]  <- dat.temp
  }
  
  # 6.4 Calculate CAFEC precip (pg. 14, eqns 10-14)
  EThat   = t(apply(PE   , 1, FUN=function(x){x * alpha}))
  Rhat    = t(apply(PR   , 1, FUN=function(x){x * beta}))
  ROhat   = t(apply(PRO  , 1, FUN=function(x){x * gamma}))
  LOSShat = t(apply(PLOSS, 1, FUN=function(x){x * delta}))
  Phat    = EThat + Rhat + ROhat - LOSShat
  
  # 6.5 Calculate excesses and deficits
  d= Precip - Phat
  # ------------------------------------------

  # ------------------------------------------
  #  7. Calculate K' and and weighted-K values
  #     7.1 Calculate K' for each time step (p25, eqn 26)
  #     7.2 Calculate weighted departures for each month and its sum
  #     7.3 Caluclate final K values
  # ------------------------------------------
  # Get the mean absolute departures for each timestep of the calibration period
  Dbar <- apply(abs(d[rows.calib,]), 2, mean)
  Dbar[Dbar==0] <- 1e-6 # replace 0 with a tiny, tiny number
  
  # 7.1 Compute K' for each time step (eqn 26)
  # IFF For some reason PE, R, and RO are all 0, we need to make it a tiny number
  bar.denom <- (Pbar+LOSSbar)
  bar.denom[bar.denom==0] <- 1e-6
  Kprime <- 1.5*log10((((PEbar+Rbar+RObar)/bar.denom) + 2.80) / Dbar) + 0.50
  
  # 7.2 Calculate weighted departures for each month and its sum
  DKprime <- Dbar * Kprime
  denom <- sum(DKprime) # denominator of eqn 27, pg. 26
  
  # 7.3 Caluclate final K values
  K <- (17.67/denom) * Kprime
  # ------------------------------------------
  
  # ------------------------------------------
  #  8. Calculate PDSI 
  #     8.1 Calculate Z-values from K & d
  #     8.2 Calculate PDSI using pdsix function
  # ------------------------------------------
  # d is a time series matrix (something x [ntime]) and K is [ntime] monthly constants
  # just do the same thing we've done everywhere else
  Z = t(apply(d   , 1, FUN=function(x){x * K}))
  
  ztemp <- as.vector(t(Z))
  pdsi  <- pdsix(z=ztemp)
  # ------------------------------------------

  # ------------------------------------------
  #  9. Format & return output
  # ------------------------------------------
  # datout: list; length = 8
  #  1. Z  = Z index (unitless); data frame; dim=c(nyr, 12)
  #  2. X  = PDSI value (unitless); data frame; dim=c(nyr, 12)
  #  3. XM = modified PDSI (unitless); data frame; dim=c(nyr, 12)
  #  4. W  = avg soil moist (unit="in"); data frame; dim=c(nyr, 12)
  #     R  = recharge
  #  5. RO = monthly runoff (unit=??); dataframe; dim=c(nyr, 12)
  #  6. S1 = effective ppt (unit="in"); array; dim=c(nyr, 12, 10)
  #        = max(0, p - f*pe)
  #  7. P  = precip input; data frame; dim=c(nyr, 12)
  #  8. T  = temperature input; data frame; dim=c(nyr, 12)

  datout <- list()
  datout$Z  <- Z
  datout$X  <- matrix(pdsi$x, nrow=nrow(Precip), byrow=T)
  datout$XM <- matrix(pdsi$xm, nrow=nrow(Precip), byrow=T)
  datout$XH <- matrix(pdsi$x4, nrow=nrow(Precip), byrow=T)
  datout$PE <- PE
  datout$W  <- W
  datout$R  <- R
  datout$RO <- RO
  datout$S1 <- S1
  datout$P  <- Precip
  datout$T  <- Temp
  
  return(datout)
  # ------------------------------------------
 
  
}

