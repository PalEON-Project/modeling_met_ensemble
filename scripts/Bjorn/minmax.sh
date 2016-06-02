#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @breif minmax.sh extracts the min and max values for all variables          #
#                                                                             #
# Requires:                                                                   #
# standard linux utilities                                                    #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                            #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

### qsub arguments
#$ -wd /home/scratch/dietze_lab/bjorn/ccsm4_sd/scripts/
#$ -j y
#$ -m eas
#$ -S /bin/bash
#$ -M bjorn@climatemodeling.org

BASEDIR=`grep '^BASEDIR=' constants.txt | sed 's/.*=//'`/output/
MODEL=`grep '^MODEL=' constants.txt | sed 's/.*=//'`
ANN_DIR=${BASEDIR}/ann/
CNCP_INPUT=${BASEDIR}/cruncep/ # Path to parsed CRUNCEP data
COARSE_INPUT=${BASEDIR}/${MODEL}/ # Path to parsed COARSE data

if [ ! -d ${ANN_DIR} ]
then
	mkdir -p ${ANN_DIR}
fi

pushd ${ANN_DIR}

for variable in lwdown psurf precipf qair swdown tair wind
do

	cat ${COARSE_INPUT}/${variable}/[0-9]*/testing_${variable}.original > \
	testing_${variable}.minmax.tmp

	cat ${CNCP_INPUT}/${variable}/[0-9]*/input_${variable}.original > \
	input_${variable}.minmax.tmp

	cat \
	${CNCP_INPUT}/${variable}/[0-3]*/[0-9]*/desired_${variable}.original > \
	desired_${variable}.minmax.tmp1

	cat desired_${variable}.minmax.tmp1 \
	${CNCP_INPUT}/${variable}/[4-9]*/[0-9]*/desired_${variable}.original > \
	desired_${variable}.minmax.tmp2

	echo "Search 6-hourly ${variable} data over all divisions"
	awk 'NR==1{min=max=$1;next;} \
	NR>1{max=(max>$1)?max:$1; min=(min>$1)?$1:min;} \
	END{printf "%0.9f %0.9f\n", min, max;}' \
	desired_${variable}.minmax.tmp2 \
	> desired_${variable}.minmax

	echo "Search monthly ${variable} input data over all divisions"
	awk 'NR==1{min=max=$1;next;} \
	NR>1{max=(max>$1)?max:$1; min=(min>$1)?$1:min;} \
	END{printf "%0.9f %0.9f\n", min, max;}' \
	input_${variable}.minmax.tmp \
	> input_${variable}.minmax

	echo "Search monthly ${variable} testing data over all divisions"
	awk 'NR==1{min=max=$1;next;} \
	NR>1{max=(max>$1)?max:$1; min=(min>$1)?$1:min;} \
	END{printf "%0.9f %0.9f\n", min, max;}' \
	testing_${variable}.minmax.tmp \
	> testing_${variable}.minmax

	rm -f desired_${variable}.minmax.tmp[1,2] \
	input_${variable}.minmax.tmp testing_${variable}.minmax.tmp

	echo "Calculating output_testing_${variable}.minmax"
	MU_TST=`awk '{printf "%0.9f", ($2+$1)/2;}' testing_${variable}.minmax`
	MU_INP=`awk '{printf "%0.9f", ($2+$1)/2;}' input_${variable}.minmax`
	MU_DES=`awk '{printf "%0.9f", ($2+$1)/2;}' desired_${variable}.minmax`
	MU_BIAS=`echo "${MU_INP} ${MU_DES}" | awk '{printf "%0.9f", $1-$2;}'`
	D_INP=`awk '{printf "%0.9f", $2-$1;}' input_${variable}.minmax`
	D_DES=`awk '{printf "%0.9f", $2-$1;}' desired_${variable}.minmax`
	# Given the above means and ranges of input and desired min's/max's
	#+ calcualate the unknown min's/max's
	OT_MIN=`echo "${MU_TST} ${MU_BIAS} ${D_DES}" | \
	awk '{printf "%0.9f", ($1-$2) - ($3*0.5);}'`

	OT_MAX=`echo "${MU_TST} ${MU_BIAS} ${D_DES}" | \
	awk '{printf "%0.9f", ($1-$2) + ($3*0.5);}'`

	if [[ ${variable} == "swdown" || ${variable} == "precipf" ]]
	then

		echo "Readjusing range to zero baseline for ${variable}"
		OT_MIN="0.000000000"
		echo "${OT_MIN} ${OT_MAX}" > output_testing_${variable}.minmax

	# Check to make sure OT_MIN is a positive value. If not then constrain
	#+ range relative to desired data. This is an issue for variables with
	#+ min/max's very close to zero (e.g., qair)
	elif [[ ${OT_MIN} < 0 ]]
	then

		echo -n "Negative ${variable} min value detected. Defaulting to "
		echo "desired minmax values. refit.R will readjust to fit mean later"
		OT_MIN=`awk '{printf "%0.9f", $1;}' desired_${variable}.minmax`
		OT_MAX=`awk '{printf "%0.9f", $2;}' desired_${variable}.minmax`
		echo "${OT_MIN} ${OT_MAX}" > output_testing_${variable}.minmax

	else

		echo "${OT_MIN} ${OT_MAX}" > output_testing_${variable}.minmax

	fi

done

echo "Done"
