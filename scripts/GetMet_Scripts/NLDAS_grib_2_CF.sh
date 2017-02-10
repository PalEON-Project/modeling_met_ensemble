#!/bin/sh

# Script to converting NLDAS from annoying (but compact) grib format to netcdf to make 
# it easier to query
#
# Alternative is to call cdo from R and have it pull individual points, but that's still 
# a pain; I think having things as netcdf will be a bit more flexible down the road.

# Steps needed (looping through each file):
# 1. Convert .grb to .nc
# 2. Assign variable names to make sure they make sense
# 3. Compress move .grb file to folder & compress if everything went okay


# Define the working directory & move there if we're not already
dir_NLDAS=/Volumes/Celtis/Meteorology/LDAS/NLDAS_FORA0125_H.002
dir_grib=GRIB_original
dir_ncdf=netcdf

yr_start=1979
yr_end=2017
yrs=(`seq ${yr_start} 1 ${yr_end}`)


cd $dir_NLDAS

# # ---------------------------------
# # Move grib files to folder just to keep everything straight
# # Note: Too many files to do this one at a time
# # ---------------------------------
# 
# echo ${yrs[@]}
# 
# for YEAR in ${yrs[@]}
# do
#  echo ${YEAR}
#  
#  mv *${YEAR}1*.grb GRIB_original/
#  mv *${YEAR}01*.grb GRIB_original/
#  mv *${YEAR}02*.grb GRIB_original/
#  mv *${YEAR}03*.grb GRIB_original/
#  mv *${YEAR}0*.grb GRIB_original/
# done
# # ---------------------------------


# ---------------------------------
# Converting files from .grb to netcdf to make them easier to use even though it comes
# at a cost of size & time
#
# ----------
# Workflow
# ----------
#  1. Start loop by year
#  2. Make folder for year in GRIB & NETCDF directories
#  3. Start loop by month
#  4. Make folder for month in GRIB & NETCDF directories
#  5. Identify days in the month (rather than prescribe a priori to make leap year easier)
#  6. Start loop by day of month
#  7. Make folder for day in GRIB directories
#  8. Identify all files for the day
#  9. Combine all files for day into 1 file & change to netcdf
# 10. move grib files to appropriate folder to clean things up a bit
# 11. change incomprehensible var names to CF
# 12. At end of month, compress grib directory
# ----------
# ---------------------------------

# NLDAS_FORA0125_H.A19790130.0300.002.grb
prefix=NLDAS_FORA0125_H.A

mkdir -p $dir_ncdf # make directory for our output

mos=(01 02 03 04 05 06 07 08 09 10 11 12)


#  1. Start loop by year
for YEAR in ${yrs[@]}
do
	#  2. Make folder for year in GRIB & NETCDF directories
	mkdir -p ${dir_grib}/${YEAR}
	mkdir -p ${dir_ncdf}/${YEAR}
	
	#  3. Start loop by month
	for MONTH in ${mos[@]}
	do
		#  4. Make folder for month in GRIB & NETCDF directories
		mkdir -p ${dir_grib}/${YEAR}/${MONTH}
		mkdir -p ${dir_ncdf}/${YEAR}/${MONTH}

		#  5. Identify days in the month (rather than prescribe a priori to make leap year easier)
		pushd ${dir_grib}
			files_mo=(`ls ${prefix}${YEAR}${MONTH}*`)
		popd
		
		# Really dumb way of doing this is to loop through all files and extract the day 
		# and then find the unique values
		days_mo=($(
			for ((i=0; i < ${#files_mo[@]}; i++))
			do
				today=(`echo ${files_mo[i]} | cut -c25-26`)
				echo $today
			done
			))
		
		# a hack way to make unique work on an array instead of a list from
		# http://stackoverflow.com/questions/13648410/how-can-i-get-unique-values-from-an-array-in-bash
		days_mo=($(echo "${days_mo[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

		# nday=${#days_mo[@]}
		
		# # -------------
		# # Snippet of code to determine if we have a leap year
		# # -------------
		# dpm=(31 28 31 30 31 30 31 31 30 31 30 31) #days per month
		# if  [[ $MONTH == "02" && ( $(($YEAR % 400)) -eq 0 || ( $(( $YEAR % 4 )) -eq 0 && $(( $YEAR % 100 )) -gt 0 ) ) ]]
		# then
		#     echo "${YEAR} is Leap Year"
		# else
		#     echo "${YEAR} is not a Leap Year"
		# fi
		# # -------------

		
		#  6. Start loop by day of month
		for DAY in ${days_mo[@]}
		do
			# Print where we are
			echo ${YEAR}-${MONTH}-${DAY}
			
			#  7. Make folder for day in GRIB directories
			mkdir -p ${dir_grib}/${YEAR}/${MONTH}/${DAY}
			
			#  8. Identify all files for the day
			pushd ${dir_grib}
				files_day=(`ls ${prefix}${YEAR}${MONTH}${DAY}*`)

				#  9. Combine all files for day into 1 file & change to netcdf
				cdo -f nc copy ${files_day[@]} ../${dir_ncdf}/${YEAR}/${MONTH}/TEMP.nc
			
				# 10. move grib files to appropriate folder to clean things up a bit
				mv ${prefix}${YEAR}${MONTH}${DAY}*.grb ${YEAR}/${MONTH}/${DAY}

			popd

			# 11. change incomprehensible var names to CF
			#     vars on separate lines just to make it easier to keep track of
			pushd ${dir_ncdf}/${YEAR}/${MONTH}
# 				cdo chname,var11,air_temperature,\
# 				var51,specific_humidity,\
# 				var1,air_pressure,\
# 				var33,eastward_wind,\
# 				var34,northward_wind,\
# 				var205,surface_downwelling_longwave_flux_in_air,\
# 				var153,convective_precipitation_fraction,\
# 				var157,specific_convective_available_potential_energy,\
# 				var228,water_potential_evaporation_amount,\
# 				var61,precipitation_amount,\
# 				var204,downwelling_shortwave_flux_in_air \
# 				TEMP.nc ${prefix}${YEAR}${MONTH}${DAY}.nc

				cdo chname,var11,air_temperature,var51,specific_humidity,var1,air_pressure,var33,eastward_wind,var34,northward_wind,var205,surface_downwelling_longwave_flux_in_air,var153,convective_precipitation_fraction,var157,specific_convective_available_potential_energy,var228,water_potential_evaporation_amount,var61,precipitation_amount,var204,downwelling_shortwave_flux_in_air TEMP.nc ${prefix}${YEAR}${MONTH}${DAY}.nc

			
				rm -f TEMP.nc # Clean up
			popd
			
		done # end day loop
		
		# 12. At end of month, compress grib directory
		pushd ${dir_grib}/${YEAR}
			tar -jcvf ${prefix}${YEAR}${MONTH}.tar.bz2 ${MONTH}
			rm -rf ${MONTH}
		popd
		
	done # end month loop
done # end year loop
# ---------------------------------
