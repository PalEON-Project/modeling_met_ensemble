# PDSI function
pdsi1 <- function(datmon, datother, snowinf, textin, kot, penopts, datpen){
  # --------------
  # Inputs:
  # --------------
  #  1. datmon: list, length=2
  #     1. P = precipf; data frame; units="in"; dim=c(nyr, 13)
  #     2. T = temperature; data frame; units="F"; dim=c(nyr, 13)
  #  2. datother: list, length=4
  #     1. lat = latitude; numeric, units="decimal degrees"; lenght=1
  #     2. watcap = water capacity; data frame; units="in"; dim=c(1, 2)
  #         watcap[,1] = awcs = awc surface layer (1")
  #         watcap[,2] = awcu = awc underlying layer (5")
  #     3. YRS = start & end years; data frame; units="years"; dim=c(2,2)
  #         YRS[1,] = desired PDSI window
  #         YRS[2,] = window for normals & CAFEC?? output
  #     4. dayz = lookup table for percentage of possible sunshine
  #  3. snowinf: list, length=3
  #     1. Tcrit = temperature to convert ppt to snow if kopt[[1]] = 1; default = 32F
  #     2. mons1 = months to convert precip to snow (typically oct to may)
  #     3. melttbl = lookup table to convert snow accumulation & snow melt
  #  4. textin: character string; site ID
  #  5. kopt: list; length=3 (general options)
  #     1. PE computation method
  #        1 = Thornthwaite
  #        2 = Penman-Monteith
  #        3 = Scaled Pan evaporation
  #     2. Mode of run
  #        1 = routine (no interaction)
  #        2 = exploratory/instructional (graphs, etc)
  #     3. snow option
  #        1 = all ppt as rain regardless of temp
  #        2 = convert rain to snow & re-distribute based on temperature & Tcrit (see snowinf)
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
  #  1. Z = Z index (unitless); data frame; dim=c(nyr, 13)
  #  2. X = PDSI value (unitless); data frame; dim=c(nyr, 13)
  #  3. XM = modified PDSI (unitless); data frame; dim=c(nyr, 13)
  #  4. W = avg soil moist (unit="in"); data frame; dim=c(nyr, 13)
  #  5. RO = monthly runoff (unit=??); dataframe; dim=c(nyr, 13)
  #  6. S1 = effective ppt (unit="in"); array; dim=c(nyr, 13, 10)
  #        = max(0, p - f*pe)
  #  7. P = precip input; data frame; dim=c(nyr, 13)
  #  8. T = temperature input; data frame; dim=c(nyr, 13)
  # --------------

  # --------------
  # Functions Called
  # --------------
  # 1. pdsix - calculates X values rom
  # 2. pethorn - calculate Thornthwaite PE
  # 3. snowppt (optional)
  # --------------
  
  # --------------
  # Workflow
  # --------------
  # 1. Check input, make pointers
  # 2. Compute PE based on desired option
  # 3. Redistribute monthly P if needed
  # 4. Get calib preiod PE, P, T, & calculate 'normals'
  # 5. Calculate soil moisture using AWC
  #     - 
  # 6. 
  # --------------
  
  return(datout)
}

