#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @brief extract_cruncep_gunzip.sh is a script for inflating CRUNCEP          #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

#ulimit -v 10000000000 # set memory usage limit to 10 GB

VAR=$1
YR_BEG_CNCP=`grep '^YR_BEG_CNCP=' constants.txt | sed 's/.*=//'`
YR_END_CNCP=`grep '^YR_END_CNCP=' constants.txt | sed 's/.*=//'`
CDO="/share/apps/cdo-1.5.1/bin/cdo"
MINLAT=`grep '^MINLAT=' constants.txt | sed 's/.*=//'`
MAXLAT=`grep '^MAXLAT=' constants.txt | sed 's/.*=//'`
MINLON=`grep '^MINLON=' constants.txt | sed 's/.*=//'`
MAXLON=`grep '^MAXLON=' constants.txt | sed 's/.*=//'`

if [[ -z ${VAR} ]]; then
   echo "extract_cruncep_gunzip.sh requires one argument"
   echo
   echo "SYNOPSIS"
   echo "   ./extract_cruncep_gunzip.sh [variable]"
   echo
   echo "ACCEPTABLE VARIABLE ARGUMENTS"
   echo "lwdown precipf psurf qair swdown tair uwind vwind"
   echo
   echo "EXAMPLES"
   echo "   ./extract_cruncep_gunzip.sh lwdown"
   exit
fi

echo ${VAR}
if [ ! -d ../cruncep/${VAR} ]; then
   mkdir ../cruncep/${VAR}
fi

pushd ../cruncep/${VAR}

for (( yr=${YR_BEG_CNCP}; yr<=${YR_END_CNCP}; yr++ )); do

   if [[ ${VAR} == "psurf" ]]; then
      ln -s ../input/cruncep_press_${yr}.nc.gz .
      gunzip -f cruncep_press_${yr}.nc.gz
      mv cruncep_press_${yr}.nc cruncep_psurf_${yr}.nc
   elif [[ ${VAR} == "precipf" ]]; then
##      ln -s ../input/cruncep_rain_${yr}.nc.gz .
      cp ../input/cruncep_rain_${yr}.nc.gz .
      gunzip -f cruncep_rain_${yr}.nc.gz
      mv cruncep_rain_${yr}.nc cruncep_precipf_${yr}.nc
   else
      ln -s ../input/cruncep_${VAR}_${yr}.nc.gz .
      gunzip -f cruncep_${VAR}_${yr}.nc.gz
   fi

   # Edit NetCDF attributes to CF standard naming
   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf" ]]; then
      ncrename -d latitude,lat -v latitude,lat -d longitude,lon \
      -v longitude,lon cruncep_${VAR}_${yr}.nc
   fi

   # Edit NetCDF attributes to CF standard naming
   ncatted -O -a units,lon,o,c,degrees_east \
   cruncep_${VAR}_${yr}.nc

   ncatted -O -a units,lat,o,c,degrees_north \
   cruncep_${VAR}_${yr}.nc

   # Parse desired geographic domain
   if [[ ${VAR} == "lwdown" || ${VAR} == "psurf" ]]; then
      # lwdown and psurf are not on same grid as other variables so
      #+ parse according to different grid
      ${CDO} \
      sellonlatbox,${MINLON},`echo "${MAXLON}-0.5" | bc`,${MINLAT},`echo "${MAXLAT}-0.5" | bc` \
      cruncep_${VAR}_${yr}.nc cruncep_${VAR}_${yr}_tmp_1.nc
   else
      ${CDO} sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} cruncep_${VAR}_${yr}.nc \
      cruncep_${VAR}_${yr}_tmp_1.nc
   fi

   ${CDO} invertlat cruncep_${VAR}_${yr}_tmp_1.nc \
   cruncep_${VAR}_${yr}_tmp_2.nc

cat << EOF > 0_5grid
gridtype = lonlat
xsize = $[(MAXLON-MINLON)*2]
ysize = $[(MAXLAT-MINLAT)*2]
xfirst = `echo "${MINLON}+0.25" | bc`
xinc = 0.5
yfirst = `echo "${MINLAT}+0.25" | bc`
yinc = 0.5
EOF

   ${CDO} setgrid,0_5grid cruncep_${VAR}_${yr}_tmp_2.nc \
   cruncep_${VAR}_${yr}_tmp_3.nc

   ${CDO} -O settaxis,${yr}-01-01,00:00,6hour cruncep_${VAR}_${yr}_tmp_3.nc \
   cruncep_${VAR}_${yr}.nc

   if [[ ${VAR} == "precipf" ]]; then
      ncrename -v rain,precipf cruncep_${VAR}_${yr}.nc
   elif [[ ${VAR} == "psurf" ]]; then
      ncrename -v press,psurf cruncep_${VAR}_${yr}.nc
   fi

   rm -f cruncep_${VAR}_${yr}_tmp_*.nc

   echo cruncep_${VAR}_${yr}.nc

done

rm -f 0_5grid
popd
