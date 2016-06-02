#!/bin/bash
# ########################################################################### #
# Script to extract the paleon domain from CRUNCEP to make the downscaling for 
# individual sites a lot easier
#
# This code was based on Bjorn Brooks' extract_cruncep_gunzip.sh and extract_cruncep_gunzip.sh
# scripts (see Bjorn folder for originals)
#
# Author: Christy Rollinson, crollinson@gmail.com
# Date: 2 June 2016
# ########################################################################### #

module load cdo
module load nco

CDO="/project/earth/packages/cdo-1.6.3rc2/bin/cdo"

# CDO="/share/apps/cdo-1.5.1/bin/cdo"
dir_raw=/projectnb/dietzelab/paleon/met_ensemble/data/full_raw/cru_ncep/analysis
dir_out=/projectnb/dietzelab/paleon/met_ensemble/data/paleon_domain

# Make the CRUNCEP folders
if [ ! -d ${dir_out}/paleon1/cruncep/ ]; then
   mkdir ${dir_out}/paleon1/cruncep/
fi

if [ ! -d ${dir_out}/alaska/cruncep/ ]; then
   mkdir ${dir_out}/alaska/cruncep/
fi


# Variables we need for PalEON 
variables=(tair precipf swdown lwdown qair psurf wind)

# Some constants from a file that Bjorn made (but still holds true)
YR_BEG_CNCP=1901
YR_END_CNCP=2010

# Bounding box of the PalEON 1 Domain
MINLAT_P1=35
MAXLAT_P1=50
MINLON_P1=-100
MAXLON_P1=-60

# Bounding Box of Alaska (more or less)
MINLAT_AK=55
MAXLAT_AK=72
MINLON_AK=-165
MAXLON_AK=-140

# -----------------------------------------------------------
# Looping through each variable!
# -----------------------------------------------------------
for VAR in ${variables[@]}
do
echo ${VAR}

# -----------------------------------------
# Going through the PalEON 1 domain
# -----------------------------------------
echo 'Paleon 1'
# Make the output directory if necessary and move to it
if [ ! -d ${dir_out}/paleon1/cruncep/${VAR} ]; then
   mkdir ${dir_out}/paleon1/cruncep/${VAR}
fi
pushd ${dir_out}/paleon1/cruncep/${VAR}

# Going through each year
for (( yr=${YR_BEG_CNCP}; yr<=${YR_END_CNCP}; yr++ )); do

   # unzip the files; doing some renaming along the way
   #  Keep in min dthat 
   if [[ ${VAR} == "psurf" ]]; then
     cp ${dir_raw}/press/cruncep_press_${yr}.nc.gz .
     gunzip -f cruncep_press_${yr}.nc.gz
     mv cruncep_press_${yr}.nc cruncep_psurf_${yr}.nc
     ncrename -v press,psurf cruncep_${VAR}_${yr}.nc
   elif [[ ${VAR} == "precipf" ]]; then
     cp ${dir_raw}/rain/cruncep_rain_${yr}.nc.gz .
     gunzip -f cruncep_rain_${yr}.nc.gz
     mv cruncep_rain_${yr}.nc cruncep_precipf_${yr}.nc
     ncrename -v rain,precipf cruncep_${VAR}_${yr}.nc
   elif [[ ${VAR} == "wind" ]]; then
     cp ${dir_raw}/uwind/cruncep_uwind_${yr}.nc.gz .
     gunzip -f cruncep_uwind_${yr}.nc.gz
     cp ${dir_raw}/vwind/cruncep_vwind_${yr}.nc.gz .
     gunzip -f cruncep_vwind_${yr}.nc.gz     

     echo "Calculating wind speed from u,v components"
     # Copy uwind file here
     cp cruncep_uwind_${yr}.nc.gz cruncep_${variable}_${year}.nc
     # Add vwind to uwind file 
     ncks -A cruncep_vwind_${yr}.nc.gz cruncep_${variable}_${year}.nc

     # Calculate wind from u and v components
     ncap -O -s 'Wind=abs(vwind)/sin(atan(abs(vwind/uwind)))' \
     cruncep_${variable}_${yr}.nc cruncep_${variable}_${yr}.nc

     # Exclude undesired variables
     ncks -O -x -v uwind,vwind cruncep_${variable}_${yr}.nc \
     cruncep_${variable}_${yr}.nc

     # Update annotations
     ncatted -a units,Wind,o,c,'m s-1' \
     -a long_name,Wind,o,c,'Wind speed at the lowest atmosphere level' \
     -a standard_name,Wind,o,c,Wind \
     -a description,Wind,o,c,'Wind=abs(V)/sin(atan(abs(V/U)))' \
     cruncep_${variable}_${yr}.nc
   else
     cp ${dir_raw}/${VAR}/cruncep_${VAR}_${yr}.nc.gz .
     gunzip -f cruncep_${VAR}_${yr}.nc.gz
   fi

# 
#   # Edit NetCDF attributes to CF standard naming
#   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf"  || ${VAR} == "wind" ]]; then
#     ncrename -d latitude,lat -v latitude,lat -d longitude,lon \
#     -v longitude,lon cruncep_${VAR}_${yr}.nc
#   fi

   # Edit NetCDF attributes to CF standard naming
   ncatted -O -a units,lon,o,c,degrees_east \
   cruncep_${VAR}_${yr}.nc

   ncatted -O -a units,lat,o,c,degrees_north \
   cruncep_${VAR}_${yr}.nc


   # Parse desired geographic domain
#   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf" ]]; then
#     # lwdown and psurf are not on same grid as other variables so
#     #+ parse according to different grid
#     ${CDO} \
#     sellonlatbox,${MINLON},`echo "${MAXLON}-0.5" | bc`,${MINLAT},`echo "${MAXLAT}-0.5" | bc` \
#     cruncep_${VAR}_${yr}.nc cruncep_${VAR}_${yr}_tmp_1.nc
#   else
     ${CDO} sellonlatbox,${MINLON_P1},${MAXLON_P1},${MINLAT_P1},${MAXLAT_P1} cruncep_${VAR}_${yr}.nc \
     cruncep_${VAR}_${yr}_tmp_1.nc
#   fi

   # Flipping the lat & lon for some reason
   ${CDO} invertlat cruncep_${VAR}_${yr}_tmp_1.nc \
   cruncep_${VAR}_${yr}_tmp_2.nc

# Defining some of the dimensions and ranges
cat << EOF > 0_5grid
gridtype = lonlat
xsize = $[(MAXLON_P1-MINLON_P1)*2] # Number of longitudinal cells; 2 cells per degree
ysize = $[(MAXLAT_P1-MINLAT_P1)*2] # Number of latitudinal cells; 2 cells per degree
xfirst = `echo "${MINLON_P1}+0.25" | bc` # xmin (midpoint_
xinc = 0.5 # x resolution
yfirst = `echo "${MINLAT_P1}+0.25" | bc` # ymin (midpoint)
yinc = 0.5 # y resolution
EOF

   # Modifying the metadata for the grid with what we defined
   ${CDO} setgrid,0_5grid cruncep_${VAR}_${yr}_tmp_2.nc \
   cruncep_${VAR}_${yr}_tmp_3.nc

   # setting the time stamp
   ${CDO} -O settaxis,${yr}-01-01,00:00,6hour cruncep_${VAR}_${yr}_tmp_3.nc \
   cruncep_${VAR}_${yr}.nc

#   if [[ ${VAR} == "precipf" ]]; then
#     ncrename -v rain,precipf cruncep_${VAR}_${yr}.nc
#   elif [[ ${VAR} == "psurf" ]]; then
#     ncrename -v press,psurf cruncep_${VAR}_${yr}.nc
#   fi

   rm -f cruncep_${VAR}_${yr}_tmp_*.nc

   echo cruncep_${VAR}_${yr}.nc

rm -f 0_5grid

done

popd
# -----------------------------------------

# -----------------------------------------
# Going through the Alaska domain
# -----------------------------------------
echo 'Alaska'

# Make the output directory if necessary and move to it
if [ ! -d ${dir_out}/alaska/cruncep/${VAR} ]; then
   mkdir ${dir_out}/alaska/cruncep/${VAR}
fi
pushd ${dir_out}/alaska/cruncep/${VAR}

# Going through each year
for (( yr=${YR_BEG_CNCP}; yr<=${YR_END_CNCP}; yr++ )); do

   # unzip the files; doing some renaming along the way
   #  Keep in min dthat 
   if [[ ${VAR} == "psurf" ]]; then
     cp ${dir_raw}/press/cruncep_press_${yr}.nc.gz .
     gunzip -f cruncep_press_${yr}.nc.gz
     mv cruncep_press_${yr}.nc cruncep_psurf_${yr}.nc
     ncrename -v press,psurf cruncep_${VAR}_${yr}.nc
   elif [[ ${VAR} == "precipf" ]]; then
     cp ${dir_raw}/rain/cruncep_rain_${yr}.nc.gz .
     gunzip -f cruncep_rain_${yr}.nc.gz
     mv cruncep_rain_${yr}.nc cruncep_precipf_${yr}.nc
     ncrename -v rain,precipf cruncep_${VAR}_${yr}.nc
   elif [[ ${VAR} == "wind" ]]; then
     cp ${dir_raw}/uwind/cruncep_uwind_${yr}.nc.gz .
     gunzip -f cruncep_uwind_${yr}.nc.gz
     cp ${dir_raw}/vwind/cruncep_vwind_${yr}.nc.gz .
     gunzip -f cruncep_vwind_${yr}.nc.gz     

     echo "Calculating wind speed from u,v components"
     # Copy uwind file here
     cp cruncep_uwind_${yr}.nc.gz cruncep_${variable}_${year}.nc
     # Add vwind to uwind file 
     ncks -A cruncep_vwind_${yr}.nc.gz cruncep_${variable}_${year}.nc

     # Calculate wind from u and v components
     ncap -O -s 'Wind=abs(vwind)/sin(atan(abs(vwind/uwind)))' \
     cruncep_${variable}_${yr}.nc cruncep_${variable}_${yr}.nc

     # Exclude undesired variables
     ncks -O -x -v uwind,vwind cruncep_${variable}_${yr}.nc \
     cruncep_${variable}_${yr}.nc

     # Update annotations
     ncatted -a units,Wind,o,c,'m s-1' \
     -a long_name,Wind,o,c,'Wind speed at the lowest atmosphere level' \
     -a standard_name,Wind,o,c,Wind \
     -a description,Wind,o,c,'Wind=abs(V)/sin(atan(abs(V/U)))' \
     cruncep_${variable}_${yr}.nc

   else
     cp ${dir_raw}/${VAR}/cruncep_${VAR}_${yr}.nc.gz .
     gunzip -f cruncep_${VAR}_${yr}.nc.gz
   fi

# 
#   # Edit NetCDF attributes to CF standard naming
#   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf"  || ${VAR} == "wind" ]]; then
#     ncrename -d latitude,lat -v latitude,lat -d longitude,lon \
#     -v longitude,lon cruncep_${VAR}_${yr}.nc
#   fi

   # Edit NetCDF attributes to CF standard naming
   ncatted -O -a units,lon,o,c,degrees_east \
   cruncep_${VAR}_${yr}.nc

   ncatted -O -a units,lat,o,c,degrees_north \
   cruncep_${VAR}_${yr}.nc


   # Parse desired geographic domain
#   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf" ]]; then
#     # lwdown and psurf are not on same grid as other variables so
#     #+ parse according to different grid
#     ${CDO} \
#     sellonlatbox,${MINLON},`echo "${MAXLON}-0.5" | bc`,${MINLAT},`echo "${MAXLAT}-0.5" | bc` \
#     cruncep_${VAR}_${yr}.nc cruncep_${VAR}_${yr}_tmp_1.nc
#   else
     ${CDO} sellonlatbox,${MINLON_AK},${MAXLON_AK},${MINLAT_AK},${MAXLAT_AK} cruncep_${VAR}_${yr}.nc \
     cruncep_${VAR}_${yr}_tmp_1.nc
#   fi

   # Flipping the lat & lon for some reason
   ${CDO} invertlat cruncep_${VAR}_${yr}_tmp_1.nc \
   cruncep_${VAR}_${yr}_tmp_2.nc

# Defining some of the dimensions and ranges
cat << EOF > 0_5grid
gridtype = lonlat
xsize = $[(MAXLON_AK-MINLON_AK)*2] # Number of longitudinal cells; 2 cells per degree
ysize = $[(MAXLAT_AK-MINLAT_AK)*2] # Number of latitudinal cells; 2 cells per degree
xfirst = `echo "${MINLON_AK}+0.25" | bc` # xmin (midpoint_
xinc = 0.5 # x resolution
yfirst = `echo "${MINLAT_AK}+0.25" | bc` # ymin (midpoint)
yinc = 0.5 # y resolution
EOF

   # Modifying the metadata for the grid with what we defined
   ${CDO} setgrid,0_5grid cruncep_${VAR}_${yr}_tmp_2.nc \
   cruncep_${VAR}_${yr}_tmp_3.nc

   # setting the time stamp
   ${CDO} -O settaxis,${yr}-01-01,00:00,6hour cruncep_${VAR}_${yr}_tmp_3.nc \
   cruncep_${VAR}_${yr}.nc

#   if [[ ${VAR} == "precipf" ]]; then
#     ncrename -v rain,precipf cruncep_${VAR}_${yr}.nc
#   elif [[ ${VAR} == "psurf" ]]; then
#     ncrename -v press,psurf cruncep_${VAR}_${yr}.nc
#   fi

   rm -f cruncep_${VAR}_${yr}_tmp_*.nc

   echo cruncep_${VAR}_${yr}.nc

rm -f 0_5grid

done

popd
# -----------------------------------------


done
# -----------------------------------------------------------

 