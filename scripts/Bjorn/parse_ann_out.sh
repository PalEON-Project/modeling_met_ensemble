#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @brief parse_ann_out.sh regrids the ASCII output from ANN downscaling       #
# into arrays that can be easily written to NetCDF/HDF files                  #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

### qsub arguments
#$ -wd /home/scratch/dietze_lab/bjorn/ccsm4_sd/scripts/
#$ -j y
#$ -S /bin/bash
#$ -m eas
#$ -M bjorn@climatemodeling.org

SPY=1460			# No. 6 hour samples per year
NTEST=1000			# No. of testing patterns (years of data)
YR_BEG=`grep '^YR_BEG=' constants.txt | sed 's/.*=//'`	# Start yr
BASEDIR=`grep '^BASEDIR=' constants.txt | sed 's/.*=//'`
ANNDIR=${BASEDIR}/output/ann/
FGC=`wc -l ${BASEDIR}/output/cruncep/cncp_idx_latlon.txt | awk '{print $1;}'`

for variable in tair swdown qair lwdown psurf precipf wind
do

	d_dir=${ANNDIR}/${variable}/ # Destination of output
	pushd ${d_dir}/

	echo "Splitting each output_testing_${variable}_[grid_cell] into annual files"
	for (( i=1; i<=${FGC}; i++ ))
	do

		ii=`printf "%04i" ${i}`
		echo "Splitting grid cell ${ii}"

		split -d -a 4 -l ${SPY} ${d_dir}/output_testing_${variable}_${ii} \
		output_testing_${variable}_${ii}_

		rm -f ${d_dir}/output_testing_${variable}_${ii}

	done
	echo "Done splitting"

	# The routine below assumes that the file names identify grid cells
	#+ using type format %04i
	echo "Pasting annual output_testing_${variable}_* files into gridded arrays"
	if [[ -d ${d_dir}/netcdf/ ]]
	then
		rm -rf ${d_dir}/netcdf/
	fi
	mkdir ${d_dir}/netcdf/

	for (( yr=1; yr<=${NTEST}; yr++ ))
	do
		# Make ${year} variable with same type format as file names from split
		year=`printf "%04i" $[yr-1]`
		year_cal=`printf "%04i" $[yr+YR_BEG-1]` # calendar year
		echo "Pasting year ${year_cal}"

		# Clean-up files if any
		junk_files=$(ls output_testing_${variable}_${year}.tmp output_testing_${variable}_${year}.tmp2 output_testing_${variable}_${year}.txt 2> /dev/null | wc -l)
		if [[ ${junk_files} != 0 ]]; then

			rm -f output_testing_${variable}_${year}.tmp \
			output_testing_${variable}_${year}.tmp2 \
			output_testing_${variable}_${year}.txt

		fi

		if [[ $[FGC/1000] == 0 ]] # If < 1000 files the proceed normally
		then

			paste -d' ' output_testing_${variable}_*_${year} > \
			${d_dir}/netcdf/output_testing_${variable}_${year_cal}.txt

		else # If > 1000 files loop in chunks to avoid too many files error

			for (( i=0; i<=$[FGC/1000]; i++ ))
			do

				if [[ ${i} == 0 ]]
				then

					paste -d' ' \
					output_testing_${variable}_${i}*_${year} > \
					output_testing_${variable}_${year}.tmp

				else

					paste -d' ' \
					output_testing_${variable}_${year}.tmp \
					output_testing_${variable}_${i}*_${year} > \
					output_testing_${variable}_${year}.tmp2

					mv output_testing_${variable}_${year}.tmp2  \
					output_testing_${variable}_${year}.tmp

				fi

				rm -f output_testing_${variable}_${i}*_${year}

			done

			mv output_testing_${variable}_${year}.tmp \
			${d_dir}/netcdf/output_testing_${variable}_${year_cal}.txt

		fi

		rm -f output_testing_${variable}_[0-9]*_${year}

	done

	rm -f output_testing_${variable}_[0-9]*

	popd

done

echo "Done!"
