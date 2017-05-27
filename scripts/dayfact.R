# Function to calculate the day factor adjustment for the Thornthwaite PE calculation
calc.dayfact <- function(timestep=c("daily", "monthly"), daylength=NULL, lat=NULL, dayz=NULL){
  # daylength = hours day length for a day or mean day length in hours for a month
  if(!is.null(daylength)){
    dayfact <- daylength/12
  } else if(timestep=="monthly" & is.null(daylength)){
    if(is.null(dayz)) stop("You have monthly data an no daylength.  Please provide Thornthwaite dayz table")
    
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
    
    # Note: I *think* we could do this on a daily scale by leveraging our met data 
    # by calculating day length from SWdown>0 per day and dividing by 12
    dayfact <- apply(dayz, 2, FUN=function(x){approx(0:50, x, min(lat, 50))$y})/30 
  } else {
    stop("You're working with non-monthly data and have not provided a daylength table... figure something out.")
  }
}


