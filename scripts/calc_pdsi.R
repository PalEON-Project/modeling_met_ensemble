##' Calculate PDSI
# ----------------------------------- 
# Description
# -----------------------------------
##' @title calc.pdsi
##' @family Meteorology
##' @author Christy Rollinson
##' @description Takes an input of monthly data formated in 1 file per year and calculates monthly PDSI
# ----------------------------------- 
# Parameters
# -----------------------------------
##' @param path.in File path to where the CF-standard .nc files should be housed
##' @param years.pdsi - which years to calculate PDSI for; if NULL (default), all available years will be used
##' @param years.calib - years to calibrate the PDSI against;
##' @param watcap - vector of length 2 indicating the water holding capacity of the upper and lower layers of the soil
##' @return list everything used to calculate PDSI
##'           1. Z  = Z-index (unitless); matrix; dim=c(nyr, 12)
##'           2. X  = PDSI value (unitless); matrix; dim=c(nyr, 12)
##'           3. XM = modified PDSI value (unitless); matrix; dim=c(nyr, 12)
##'           4. XH =  PDHI; defined as x4[t] = 0.897*x4[t-1] + Z/3
##'                    circumvents probability of switching between three 
##'                    alternative version of the index
##'           5. PE = potential evapotranspiration (in/mo); matrix; dim=c(nyr, 12)
##'           6. W  = calcualted average soil moisture (inches); matrix; dim=c(nyr, 12)
##'           7. R  = calculated monthly recharge (inches); matrix; dim=c(nyr, 12)
##'           8. RO = calculated monthly runoff (inches); matrix; dim=c(nyr, 12)
##'           9. S1 = effective preciptiation (inches); array; dim=c(nyr, 12, 10); max(0, p-f*pe)
##'          10. P  = precipitation input (inches); matrix; dim=c(nyr, 12)
##'          11. T  = temperature input (F); matrix; dim=c(nyr, 12)
##'          12. D  = daylength adjustment factor; matrix; dim=c(nyr, 12)
##' @export
# -----------------------------------

calc.pdsi <- function(path.in, years.pdsi, years.calib, watcap){
  
  # Get a list of all the files available to do the aggregation for
  files.met <- dir(path.in)
  yrs.files <- strsplit(files.met, "[.]")
  yrs.files <- matrix(unlist(yrs.files), ncol=length(yrs.files[[1]]), byrow=T)
  yrs.files <- as.numeric(yrs.files[,ncol(yrs.files)-1])
  
  # Subset to just the files we want to work with
  if(is.null(years.pdsi)){
    years.pdsi <- yrs.files
  }
  files.met <- files.met[which(yrs.files %in% years.pdsi)]
  years.pdsi <- yrs.files[which(yrs.files %in% years.pdsi)]
  
  
  # Loop through and extract our temperature, precipitation, & daylength (if available) values  
  mos <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  TEMP1 <- PRECIP1 <- dayfact <- matrix(nrow=length(years.pdsi), ncol=12)
  row.names(TEMP1) <- row.names(PRECIP1) <- row.names(dayfact) <- years.pdsi
  colnames(TEMP1) <- colnames(PRECIP1) <- colnames(dayfact) <- mos
  for(i in 1:length(files.met)){
    ncT <- ncdf4::nc_open(file.path(path.in, files.met[i]))
    
    if(i == 1) site.lat <- ncdf4::ncvar_get(ncT, "latitude")
    
    TEMP1[i,] <- ncdf4::ncvar_get(ncT, "air_temperature")
    PRECIP1[i,] <- ncdf4::ncvar_get(ncT, "precipitation_flux")
    if("daylength" %in% names(ncT$var)) dayfact[i,] <- ncdf4::ncvar_get(ncT, "daylength")
    
    ncdf4::nc_close(ncT)
  }
  
  # Make sure our matrices are ordered oldest to newest
  TEMP1   <- TEMP1[order(as.numeric(row.names(TEMP1))),]
  PRECIP1  <- PRECIP1[order(as.numeric(row.names(PRECIP1))),]
  dayfact <- dayfact[order(as.numeric(row.names(dayfact))),]
  
  # Convert Temp: Kelvin to Fahrenheit 
  C2F <- function(x){x*9/5 + 32}
  TEMP1 <- C2F(TEMP1-273.15)
  
  # Convert Precip: kg/m2/s = mm/mo 
  # Making a matrix to do the time conversion so we can factor in leap year
  sec2mo <- matrix(lubridate::days_in_month(1:12), nrow=nrow(PRECIP1), ncol=12, byrow=T)
  yrs.leap <- which(lubridate::leap_year(as.numeric(row.names(PRECIP1))))
  sec2mo[yrs.leap, 2] <- sec2mo[yrs.leap, 2]+1
  sec2mo <- sec2mo*60*60*24
  
  mm2in <- 1/25.4
  
  PRECIP1 <- (PRECIP1*sec2mo)/25.4 # mm to in
  
  # Calculate the proper daylength adjustment
  dpm <- lubridate::days_in_month(1:12)
  dayfact <- dayfact/12*dpm/30
  
  # Package Everythign for the function
  datmet <- list(Temp=TEMP1, Precip=PRECIP1)
  
  datother <- list()
  datother$pdsi.fun <- "."
  datother$metric <- F
  datother$lat <- site.lat
  datother$watcap <- list(awcs=watcap[1], awcu=watcap[2])
  datother$yrs.calib <- years.calib 
  datother$dayz      <- NULL
  datother$dayfact   <- dayfact
  
  # Run the actual PDSI calculation
  pdsi.out <- pdsi1(datmet, datother, metric=F, method.PE="Thornthwaite", snow=NULL, snowopts=NULL, penopts=NULL, datpen=NULL)
  pdsi.out$D <- dayfact
  
  return(pdsi.out)
}