# Function to calculate potential evapotranspiration by Thornthwaite method
# Translated from Ben Cook's matlab code
# Author: Christy Rollinson, crollinson@mortonarb.org
# Note: Deleted Ben's hemisphere argument since this will 
#       only work for N hemisphere righ tnow

# -----------------------
# Inputs:
# -----------------------
#  1. Temp = matrix/dataframe/array of and 12 monthly temperatures (mean values) 
#            dim = nyear X 12; with years as rownames
#            units = degrees Celcius (Ben's script was Fahrenheit) 
#            * Note: this deviates from Ben's script because he has year in column 1
#  2. yrs.calib = start & end years of 'calibration' period to compute heat index
#           length = 2
#           units = years
#  3. lat = decimal latitude of the locaiton of interest
#           length = 1
#           units = decimal degrees
#  4. dayz = table of mean possible monthly duration of sunlght in 
#            dim=51 x 12
#            units = ?????
#            this is a matlab file from Ben! 
#            Original source = Thornthwaite & Mather, p 228, table 6
# -----------------------

# -----------------------
# Outputs
# -----------------------
# 1. PE (Potential Evapotranspiration) = maxtrix with col 1 = year, and 12 monthly values
#      dim= nyear  x 12; rownames = years
#      units = mm/time (**NOTE** Ben's Script returns in/mo; we're going to stick with 
#              metric whever possible)
# -----------------------

# -----------------------
# Notes & Citations
# -----------------------
# Global tables needed: daylen -- these must also be global in calling pgm
#
# Sources:  
#  Sellers, 1960, Physical Climatology, 
#        Heat index equations for I and a from p. 171
#  Pelton, King, and Tanner, 1960, An evaluation of the Thornthwaite and Mean
#        temperature methods for determining potential evapotranspiration. 
#        Argonomy Journal, 387-395 --  eqns 3, 4 for computing undadjusted
#        PE in cm/month and for adjusting for deviation from 12-hr day and
#        30-day month
#  Thornthwaite and Mather 1957, Instructions ...
#        Table 5, p. 226 for unadjusted PE when mean temp above 26.5C
#        Table 6, 7, p. 228,229   mean poss monthly duration of sunlight
# -----------------------

# -----------------------
# Workflow
# -----------------------
# 1. Check inputs, build some tables for reference
# 2. Calculate temperature normals
# 3. Apply heat index equation (if T>0C)
# 4. Calculate *un-adjusted PE) (30-day month, 12 hrs sun)
#     ** NOTE: This could probably be adjusted for different time steps
# 5. Adjust for sunshine duration (non-30 day months; sun >/< 12 hrs)
#     ** NOTE: This could probably be adjusted for different time steps
# 6. Convert PE from mm to month to inches per month
# 7. Format & return output
# -----------------------

PE.thorn <- function(Temp, yrs.calib, lat, dayfact, celcius=T){

  # ------------------------------------------
  # 1. Check inputs, build some tables for reference
  # ------------------------------------------
  library(lubridate)
  # library(R.matlab)
  # dayz <- readMat("PDSI_fromBenCook/PDSICODE/daylennh.mat")$dayz
  
  # Set up some constants
  # daysmon=c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) # number of days in month
  dpm <- days_in_month(1:12)
  anan=NaN
  nyr=nrow(Temp)
  yrs <- as.numeric(row.names(Temp))
  yrs.leap <- which(leap_year(yrs))
  
  # Build Table for unadjusted PE for t greater than 26.5 oC, or 80 oF. Table
  # values are in mm/day, and you specify T in deg C in Thornthwaite.  
  Thot = c(4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.6, 4.6, 4.6, 
           4.6, 4.7, 4.7, 4.7, 4.8, 4.8, 4.8, 4.8, 4.9, 4.9,
           4.9, 5.0, 5.0, 5.0, 5.0, 5.1, 5.1, 5.1, 5.1, 5.2, 
           5.2, 5.2, 5.2, 5.2, 5.3, 5.3, 5.3, 5.3, 5.4, 5.4,
           5.4, 5.4, 5.4, 5.5, 5.5, 5.5, 5.5, 5.5, 5.6, 5.6,
           5.6, 5.6, 5.6, 5.6, 5.7, 5.7, 5.7, 5.7, 5.7, 5.8, 
           5.8, 5.8, 5.8, 5.8, 5.8, 5.8, 5.9, 5.9, 5.9, 5.9, 
           5.9, 5.9, 5.9, 5.9, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0,
           6.0, 6.0, 6.0, 6.0, 6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 
           6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 6.1, 
           6.1, 6.1, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 
           6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 
           6.2)
            
   xThot <- seq(26.0, 38.0, by=0.1)
   
   
   # handy functions to convert to/from Fahrenheit
   F2C <- function(x){5/9*(x-32)}
   C2F <- function(x){x*9/5 + 32}
   mm2in <- 1/25.4

   if(celcius==F){ 
     Temp = F2C(Temp)
   }
   # ------------------------------------------
   
   # ------------------------------------------
   # 2. Calculate temperature normals
   # ------------------------------------------
   yrs.norm <- which(yrs>=yrs.calib[1] & yrs<=yrs.calib[2])
   norms <- colMeans(Temp[yrs.norm,])
   # norms <- F2C(norms)
   # ------------------------------------------

   
   # ------------------------------------------
   # 3. Apply heat index equation (if T>0C)
   # Apply eqn for heat index I in Sellers, P. 171
   # Note that heat index is computed on months whos mean T is above freezing,
   # and that for purposes of the computation, any monthly normal above 26.5C is
   # set to 26.5C
   # ------------------------------------------
   # Selecting only greater than 0 and capping at 26.5
   Tmean1 <- norms[norms>0]
   Tmean1[Tmean1>26.5] <- 26.5
   I <- sum((Tmean1/5)^ 1.514)
   
   # Heat index equation to get exponent "a"
   a=1E-6 * (0.675*I^3 - 77.1*I^2 + 17920*I + 492390)
   # ------------------------------------------

   # ------------------------------------------
   # 4. Calculate *un-adjusted PE (30-day month, 12 hrs sun)
   #     ** NOTE: This could probably be adjusted for different time steps
   # ------------------------------------------
   # Temp = F2C(Temp)
   
   # Identifying cells that go out of boundary conditions
   Lwarm <- which(Temp>=26.5)
   Lhot  <- which(Temp> 38.0)
   Lcold <- which(Temp<= 0.0)
   
   # Replace anything above 38C with 38 C
   if(length(Lhot)>0) Temp[Lhot] <- 38.0
   
   # Compute unadjusted PE; in mm/day; assumes 30 days per month
   # PE = 16 * ((10.0 * Temp / I)^a)/30
   PE = 16 * ((10.0 * Temp / I)^a) # mm/mo
   
   # Replace anything with Temp <=0, as 0
   if(length(Lcold)>0) PE[Lcold] <- 0
   
   # Dealing with high temperature cases using the table "Thot" defined above
   if(length(Lwarm)>0){
     toast <- Temp[Lwarm]
     S <- matrix(0, ncol=ncol(Temp), nrow=nyr)
     
     # Extract PE from Table using interpolation/approximation
     # in matlab this was interp1; in R, it looks like approx() works
     # This returns values in mm/day
     S[Lwarm] <- approx(xThot, Thot, toast)$y
     
     # # convert mm/day to mm/mo
     S <- t(apply(S, 1, FUN=function(x){x * dpm}))
     
     # # Add an extra day to leap year
     S[yrs.leap,2] <- S[yrs.leap,2]*29/28
     
     # Putting our warm-adjusted values in our PE matrix
     # Units = mm/day
     PE[Lwarm] <- S[Lwarm]
   }
   
   # If we have monthly data, convert mm/day to mm/mo
   # if(ncol(PE)==12) {
   #   PE <- t(apply(PE, 1, FUN=function(x){x * dpm}))
   #   PE[yrs.leap,2] <- PE[yrs.leap,2]*29/28
   # }
   # ------------------------------------------

   # ------------------------------------------
   # 5. Adjust for sunshine duration (non-30 day months; sun >/< 12 hrs)
   #     ** NOTE: This could probably be adjusted to work daily
   # ------------------------------------------
   # Use linear interpolation with the dayz table to pull the adjustment factor
   # Should return a vector with lenght 12 (b/c apply to the columns of dayz)
   # dividing by 30 gets us the adjustment when assuming a 30-day month
   # 
   # Note: I *think* we could do this on a daily scale by leveraging our met data 
   # by calculating day length from SWdown>0 per day and dividing by 12
   dayfact=NULL
   if(is.null(dayfact)){
     dayfact <- apply(dayz, 2, FUN=function(x){approx(0:50, x, min(lat, 50))$y})/30
   }
   # if(ncol(PE)==12){
   #   dayfact <- dayfact/30
   # }
   
   # Calculating adjusted PE; will do leap year adjustment in 
   # next step to keep matrices smaller
   # Will return PE in mm/mo
   PE <- t(apply(PE, 1, FUN=function(x){x * dayfact}))
   
   if(ncol(PE)==12){
     PE[yrs.leap,2] <- PE[yrs.leap,2] * 28/29
   }
   
   
   # # Convert PE from mm to in
   PE <- PE*mm2in
   # ------------------------------------------

   # ------------------------------------------
   # 7. Format & return output
   # ------------------------------------------
   row.names(PE) <- yrs
   return(PE)
   # ------------------------------------------
}
