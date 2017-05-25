PE.thorn <- function(T, yrc, lat, dayz){
   # Function to calculate potential evapotranspiration by Thornthwaite method
   # Translated from Ben Cook's matlab code
   # Note: Deleted Ben's hemisphere argument since this will 
   #       only work for N hemisphere righ tnow
   
   # -----------------------
   # Inputs:
   # -----------------------
   #  1. T = matrix of year (col1) and 12 monthly temperatures (mean values) 
   #         dim = nyear X 13
   #         units = degrees CELCIUS
   #  2. yrc = start & end years of 'calibration' period to compute heat index
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
   #      dim= nyear  x 13
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
   
   # Set up some constants
   daysmon=[31 28 31 30 31 30 31 31 30 31 30 31] # number of days in month
   anan=NaN
   
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

}
