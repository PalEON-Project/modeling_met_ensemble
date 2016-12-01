PalEON Phase3 Met Ensemble README
Version: version 1
Date: 30 November, 2016

###################################
# PEOPLE
###################################
Dataset Author: Christine Rollinson 
                crollinson@gmail.com
                PalEON Postdoc 2014-2016 (Boston University)
                current affiliation: Morton Arboretum, Lisle, IL (as of 28 December 2016)

PI Contact:     Mike Dietze, dietze@bu.edu

###################################
# Dataset Description
###################################
General:
Hourly meteorological forcing data for 850-2015 A.D. for PalEON MIP Phase 3: Data Assimilation. Forcing includes the following variables: 
      Variable Name   Description                    Units          
   1) tair            air temperature                K
   2) precipf         precipitation rate             kg m-2 s-1     
   3) swdown          incident shortwave radiation   W m-2
   4) lwdown          incident longwave radiation    W m-2
   5) press           surface pressure               Pa
   6) qair            specific humidity              kg kg-1
   7) wind            wind speed                     m s-1

Variable dimensions: [lat, lon, time]
Dimension Units:
      Name    Description       Units
   1) lat     latitude          decimal degrees
   2) lon     longitude         decimal degrees
   3) time    simulation time   hours since 0850-01-01 00:00:00


More detailed Overview:
Hourly meteorological forcing was created through a two-step bias correction and temporal downscaling process by combining the following four datasets:
 1) NLDAS (training data, 1980-2015), 
 2) CRUNCEP (1901-1979), 
 3) CMIP5 historical simulations (1850-1900), 
 4) CMIP5 past millennium (p1000) simulations (850-1849).  

CMIP5 simulations include the following four GCMs had all required seven meteorological variables at daily resolution (except radiation) for the p1000 simulations and are included in our ensemble: bcc-csm1-1, CCSM4, MIROC-ESM, MPI-ESM-P.  Version 1 of the PalEON Phase 3 met ensemble uses output from the r1i1p1 scenario.  More information on CMIP5 can be found here: http://cmip-pcmdi.llnl.gov/cmip5/experiment_design.html and here: http://cmip-pcmdi.llnl.gov/cmip5/docs/cmip5_data_reference_syntax.pdf

The two steps for the creation of hourly meteorological forcing data are as follows:
 1) Spatial downscaling / climatological de-biasing & creation of continuous daily records — This step uses high-resolution LDAS (Land Data Assimilation Systems) data to statistically adjust the climatic mean of datasets to match the adjacent datasets. This climate adjustment takes into account the seasonal (day of year) cycle and corrects for both potential incorrect seasonal patterns in GCMs. The bias-adjustment is determined by comparing either the available overlapping period of datasets or closest available 30 years.  For example, CRUNCEP was downscaled by comparing the climatic data for 1980-2010 between CRUNCEP and NLDAS whereas the climatology from 1820-1849 of CMIP5 p1000 simulations was compared to 1850-1879 CMIP5 historical simulation that had already been bias-corrected. This step uses covariance among meteorological variables to estimate the daily variability of lwdown, press, and wind in the CRUNCEP dataset (all assumed constant prior to 1950) and daily variability of wind, short-, and longwave radiation in the CMIP5 datasets (only available at monthly resolution).
 2) Temporal Downscaling to hourly records — This step generates hourly meteorological values based on the observed hourly relationship between mean daily value, hour, and the last hour of the previous day seen in the NLDAS training data.  This step is applied in filter starting with the last day in the dataset (2015-12-31) and applied day by day towards the beginning of the dataset (850-01-01). This filter approach smooths over abrupt transitions at splice points and should reduce impossible jumps in meteorological forcing that might cause problems in ecosystem model simulations.

Uncertainty is propagated through time and hierarchical ensemble members are generated in both steps above.  See ‘Dataset Organization’ (below) for information on how this is indicated in the file structure.  

Code Location (Public Repository): 
https://github.com/PalEON-Project/modeling_met_ensemble


###################################
# Dataset Organization
###################################
Dataset Location: Cyverse (de.cyverse.org)
Cyverse File Path: /iplant/home/crollinson/paleon/phase3_met_drivers/1hr/
Current Version Directory: met_phase3_1hr_v1
Sites Available: HARVARD

The general file path structure is
[Cyverse File Path]/[Current Version]/[SITE]/[GCM]/[Ensemble ID]/[Ensemble ID]_[YEAR].nc

The naming scheme for Ensemble IDs is as follows:
[SITE]_[GCM]_[resolution]_[day ensemble ID]_[hour ensemble ID]

Within the directory met_phase3_1hr_[version], there is a directory for each site available.  Within that site directory is a directory for each GCM.  Each GCM then has directories for each ensemble member.  Within each ensemble member directory are files for each years with all variables for that year in a single netcdf file.

Each netcdf file has the following seven (7) variables at hourly resolution: air temperature, precipitation rate, incident shortwave radiation, incident longwave radiation, surface pressure, specific humidity, wind speed.  See ‘Dataset Description’ (above) for netcdf variable codes and units.  All variables should be stored with the following three dimensions: latitude, longitude, time.  (See ‘Dataset Description’ above for more information).  Note: The length of time the time dimension will have 24 more observations in leap years compared to non-leap years.

###################################
# Other Information: Met Variable Name Translations
###################################
Below is a table showing the naming crosswalk between the different met products that underly the PalEON MIP Phase 3 meteorological forcing data:

   PalEON                CF Standard                                CMIP5                 CRUNCEP        NLDAS
1) tair (K)              air_temperature                            tas, tasmax, tasmin   tair           N2-m_above_ground_Temperature
2) precipf (kg m-2 s-1)  precipitation_flux                         pr                    rain           Precipitation_hourly_total
3) swdown (W m-2)        surface_downwelling_shortwave_flux_in_air  rsds                  swdown         SW_radiation_flux_downwards_surface
4) lwdown (W m-2)        surface_downwelling_longwave_flux_in_air   rlds                  lwdown         LW_radiation_flux_downwards_surface
5) press (Pa)            air_pressure                               psl                   press          Pressure
6) qair (kg kg-1)        specific_humidity                          huss                  qair           N2-m_above_ground_Specific_humidity
7) wind (m s-1)          eastward_wind, northward_wind              ua, va                uwind, vwind   N10-m_above_ground_Zonal_wind_speed, N10-m_above_ground_Meridional_wind_speed
