#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @brief extract_cruncep_parse.sh updates the NetCDF annotations to CF        #
# standards and also calculates wind speed from the u and v components        #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

### qsub arguments
#$ -wd /home/scratch/dietze_lab/Getson_PalEON/ATM_Script_Run/Scripts/
#$ -j y
#$ -m eas
#$ -S /bin/bash
#$ -M bjorn@climatemodeling.org

YR_BEG_CNCP=`grep '^YR_BEG_CNCP=' constants.txt | sed 's/.*=//'`
YR_END_CNCP=`grep '^YR_END_CNCP=' constants.txt | sed 's/.*=//'`
CDO="/share/apps/cdo-1.5.1/bin/cdo"

if [ ! -s idx_latlon.txt ]
then

	ncks -d time,0,0 ../cruncep/lwdown/cruncep_lwdown_1901.nc | grep '^time\[0\]=' | \
	grep -v '1e+34' | grep -v ' day as ' | awk '{print $2, $3}' | \
	sed 's/lon\[//;s/lat\[//;s/\]=/ /g' | awk '{print $1, $3, $2, $4;}' | \
	grep -v '^hours' > idx_latlon.txt

fi

NE=`wc -l idx_latlon.txt | awk '{print $1;}'`
LATIDX=(`awk '{printf "%02i ", $1;}' idx_latlon.txt`)
LONIDX=(`awk '{printf "%02i ", $2;}' idx_latlon.txt`)

for variable in lwdown psurf precipf qair swdown tair wind
do

   pushd /home/scratch/dietze_lab/Getson_PalEON/ATM_Script_Run/cruncep/${variable}

   for (( year=${YR_BEG_CNCP};year<=${YR_END_CNCP};year++ )); do

      # Fix missing value according to CF standards
      if [[ ${variable} == "lwdown" || ${variable} == "psurf" ]]; then
         ncatted -O -a _FillValue,,o,f,'-1e+34' cruncep_${variable}_${year}.nc
      fi

      if [[ ${variable} == "wind" ]]; then

         echo "Calculating wind speed from u,v components"
         # Copy uwind file here
         cp /home/scratch/dietze_lab/Getson_PalEON/ATM_Script_Run/cruncep/uwind/cruncep_uwind_${year}.nc cruncep_${variable}_${year}.nc
         # Add vwind to uwind file 
         ncks -A /home/scratch/dietze_lab/Getson_PalEON/ATM_Script_Run/cruncep/vwind/cruncep_vwind_${year}.nc cruncep_${variable}_${year}.nc

         # Calculate wind from u and v components
         ncap -O -s 'Wind=abs(vwind)/sin(atan(abs(vwind/uwind)))' \
         cruncep_${variable}_${year}.nc cruncep_${variable}_${year}.nc

         # Exclude undesired variables
         ncks -O -x -v uwind,vwind cruncep_${variable}_${year}.nc \
         cruncep_${variable}_${year}.nc

         # Update annotations
         ncatted -a units,Wind,o,c,'m s-1' \
         -a long_name,Wind,o,c,'Wind speed at the lowest atmosphere level' \
         -a standard_name,Wind,o,c,Wind \
         -a description,Wind,o,c,'Wind=abs(V)/sin(atan(abs(V/U)))' \
         cruncep_${variable}_${year}.nc

      fi

   #   echo "Parsing input and target data for year: ${year}"
   #   if [ -d ${year} ]; then
   #      rm -rf ${year}
   #   fi
   #   mkdir ${year}
   #
   #   for (( latlon=0; latlon<${NE}; latlon++ )); do
   #
   #      ncks -O -d lat,${LATIDX[latlon]},${LATIDX[latlon]} \
   #      -d lon,${LONIDX[latlon]},${LONIDX[latlon]} \
   #      cruncep_${variable}_${year}.nc \
   #      ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.tmp
   #
   #      ncwa -O -a lat,lon \
   #      ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.tmp \
   #      ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.nc
   #
   #      rm ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.tmp
   #
   #      # Convert to ASCII text
   #      if [[ $[(year-1904)%4] == 0 ]]; then # Parse out leap day values
   #
   #         ${CDO} outputext ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.nc \
   #         | grep "^  " | awk 'NR<237 || NR>240{print $0}' > \
   #         ${year}/${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.csv
   #
   #      else # No leap days to exclude so parse normally
   #
   #         ${CDO} outputext ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.nc \
   #         | grep "^  " > \
   #         ${year}/${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.csv
   #
   #      fi
   #
   #      rm -f ${year}_${LONIDX[latlon]}lon_${LATIDX[latlon]}lat.nc
   #
   #   done

   done

done

popd
