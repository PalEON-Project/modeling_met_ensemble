#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @brief parse_cruncep.sh is the final step in parsing the CRUNCEP NetCDF     #
# data to ASCII files that can be read by the ANN downscaling software        #
#                                                                             #
# Requires:                                                                   #
# standard linux utilities, resequence.R                                      #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

### qsub arguments
#$ -wd /home/scratch/dietze_lab/bjorn/mpi_sd2/scripts/
#$ -j y
#$ -m eas
#$ -S /bin/bash
#$ -M bjorn@climatemodeling.org

BASEDIR=`grep '^BASEDIR=' constants.txt | sed 's/.*=//'`
MODEL=`grep '^MODEL=' constants.txt | sed 's/.*=//'`
MODEL_CAP=`echo ${MODEL} | tr [a-z] [A-Z]`-ESM-P
OUTDIR=${BASEDIR}/output/cruncep/
CNCP_INPUT=${BASEDIR}/../cruncep/ # Path to CRUNCEP archive
CDO=/share/apps/cdo-1.5.1/bin/cdo     # Path to CDO on ebi-cluster
# No. of pieces to divide the domain into
CGC=`ls -d output/${MODEL}/lwdown/[0-9]*/ | wc -l`
SPY=1460			# No. 6 hour samples per year
MPY=12				# No. of months per year
MINLAT=`grep '^MINLAT=' constants.txt | sed 's/.*=//'`
MAXLAT=`grep '^MAXLAT=' constants.txt | sed 's/.*=//'`
MINLON=`grep '^MINLON=' constants.txt | sed 's/.*=//'`
MAXLON=`grep '^MAXLON=' constants.txt | sed 's/.*=//'`
YR_BEG_CNCP=`grep '^YR_BEG_CNCP=' constants.txt | sed 's/.*=//'` # Start yr
YR_END_CNCP=`grep '^YR_END_CNCP=' constants.txt | sed 's/.*=//'` # Final yr

if [ ! -d ${OUTDIR} ]
then
	mkdir -p ${OUTDIR}
fi

if [[ ! -s ${OUTDIR}/mask.nc || ! -s ${OUTDIR}/cncp_idx_latlon.txt ]]
then

	echo "Creating mask for oceans and lakes from lwdown file"
	ncks -v precipf,lat,lon,time -d time,0,0 \
	${CNCP_INPUT}/precipf/cruncep_precipf_1901.nc ${OUTDIR}/mask.tmp1

	ncks -x -v mask ${OUTDIR}/mask.tmp1 ${OUTDIR}/mask.tmp2
	${CDO} setmissval,'-99999' ${OUTDIR}/mask.tmp2 ${OUTDIR}/mask.tmp3
	${CDO} gtc,0 ${OUTDIR}/mask.tmp3 ${OUTDIR}/mask.nc
	rm -f ${OUTDIR}/mask.tmp*

	ncks -v Incoming_Long_Wave_Radiation,lat,lon,time \
	-d time,0,0 ${CNCP_INPUT}/lwdown/cruncep_lwdown_1901.nc \
	${OUTDIR}/mask.tmp1

	ncks -x -v mask ${OUTDIR}/mask.tmp1 ${OUTDIR}/mask.tmp2
	${CDO} setmissval,'-99999' ${OUTDIR}/mask.tmp2 ${OUTDIR}/mask.tmp3
	${CDO} gtc,0 ${OUTDIR}/mask.tmp3 ${OUTDIR}/mask2.nc

	${CDO} ifthen ${OUTDIR}/mask.nc ${OUTDIR}/mask2.nc ${OUTDIR}/mask3.nc
	mv ${OUTDIR}/mask3.nc ${OUTDIR}/mask.nc
	rm -f ${OUTDIR}/mask.tmp* ${OUTDIR}/mask2.nc

	echo "Creating cncp_idx_latlon.txt"
	ncks ${OUTDIR}/mask.nc | grep '^time\[' | grep lat | \
	grep -v -e '-99999' | awk '{print $2, $3;}' | \
	sed 's/ lon.*=/ /;s/lat.*=//' > ${OUTDIR}/cncp_idx_latlon.txt

fi

FGC=`wc -l ${OUTDIR}/cncp_idx_latlon.txt | awk '{print $1}'` # No. of fine cells
# Find lat,lon bounds for resequence.R
ncks -d time,0,0 \
../mpi_esm_p/past1000/rlds_Amon_MPI-ESM-P_past1000_r1i1p1_085001-184912.nc\
 foo1.nc


${CDO} sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} foo1.nc foo2.nc
NLAT=`ncdump foo2.nc | head | grep 'lat =' | sed 's/.* = //;s/ .*//'`
NLON=`ncdump foo2.nc | head | grep 'lon =' | sed 's/.* = //;s/ .*//'`
# No inflation
CLAT_PM=`echo "${MAXLAT} ${MINLAT} ${NLAT}" | awk '{print (($1-$2)/$3)/2;}'`
# No inflation
CLON_PM=`echo "${MAXLON} ${MINLON} ${NLON}" | awk '{print (($1-$2)/$3)/2;}'`
rm -f foo1.nc foo2.nc
# Update resequence.R lon,lat bounds
sed -i "s/clon_pm=.*/clon_pm=${CLON_PM}/;s/clat_pm=.*/clat_pm=${CLAT_PM}/;\
s/fgc=.*/fgc=${FGC}/;s/model=.*/model='${MODEL}'/;\
s/flon_min=.*/flon_min=${MINLON}/;s/flon_max=.*/flon_max=${MAXLON}/;\
s/flat_min=.*/flat_min=${MINLAT}/;s/flat_max=.*/flat_max=${MAXLAT}/" \
resequence.R

for variable in lwdown psurf precipf qair swdown tair wind
do

	echo -n "Step 1. Parse ${variable} CRUNCEP 6 hourly data from "
	echo "NetCDF files to ASCII monthly means"
	rm -rf ${OUTDIR}/${variable}
	mkdir -p ${OUTDIR}/${variable}

	# Copy parsing scripts to appropriate directory and edit to
	#+ process ${variable}
	sed "s/variable=.*/variable='${variable}'/" resequence.R > \
	${OUTDIR}/${variable}/resequence.R

	cp knn.R ${OUTDIR}/${variable}/

	pushd ${OUTDIR}/${variable}

	for ((yr=${YR_BEG_CNCP}; yr<=${YR_END_CNCP}; yr++))
	do

		# Reset missing value to fit in single precision
		${CDO} setmissval,'-99999' \
		${CNCP_INPUT}/${variable}/cruncep_${variable}_${yr}.nc \
		${yr}.tmp1

		# Remove mask variable
		ncks -O -x -v mask ${yr}.tmp1 ${yr}.tmp1

		# Mask out oceans from previously generated standard mask file
		${CDO} ifthen ../mask.nc ${yr}.tmp1 ${yr}.tmp2

		# Permute order to (lat,lon,time) so that the output
		#+ resequence from ncks will advance one grid cell at a
		#+ time across all time steps
		# Ignore warnings
		ncpdq ${yr}.tmp2 -a lat,lon,time -o ${yr}.tmp3

		# Output 6 hourly data for ANN to be used as 'desired' data
		if [[ $[(yr-1904)%4] == 0 ]] # Parse out leap day values
		then
			ncks -d time,0,235 -d time,240,1463 ${yr}.tmp3 | \
			grep '^lat\[' | grep time | grep -v -e '-99999' | \
			sed 's/.*=//;s/ .*//' > ${yr}.txt
		else
			ncks ${yr}.tmp3 | grep '^lat\[' | grep time | \
			grep -v -e '-99999' | sed 's/.*=//;s/ .*//' > \
			${yr}.txt
		fi

		# If variable is precipf then convert mm/6hr to mm/s (kg/m2/s)
		if [[ ${variable} == "precipf" ]]
		then

			echo "Converting precip values from mm/6h to mm/s"

			awk '{printf "%0.9f\n", $1/6/60/60;}' \
			${yr}.txt > ${yr}.tmp

			mv ${yr}.tmp ${yr}.txt

		fi

		echo -n "nrows in 6hourly: "
		echo "`wc -l ${yr}.txt | awk -v SPY=${SPY} '{print $1/SPY;}'`"

		# Output monthly means for step 2
		${CDO} monmean ${yr}.tmp2 ${yr}.tmp4

		# Permute order to (lat,lon,time) so that the output
		#+ sequence from  ncks will advance one grid cell at a time
		#+ across all time steps

		# Ignore warnings
		ncpdq ${yr}.tmp4 -a lat,lon,time -o ${yr}.tmp5

		ncks ${yr}.tmp5 | grep '^lat\[' | grep time | \
		grep -v -e '-99999' | sed 's/.*=//;s/ .*//' > \
		${yr}_monthly.txt

		# If variable is precipf then convert mm/6hr to mm/s (kg/m2/s)
		if [[ ${variable} == "precipf" ]]
		then

			echo "Converting precip values from mm/6h to mm/s"

			awk '{printf "%0.9f\n", $1/6/60/60;}' \
			${yr}_monthly.txt > \
			${yr}_monthly.tmp

			mv ${yr}_monthly.tmp ${yr}_monthly.txt

		fi

		echo "nrows in monthly data: \
		`wc -l ${yr}_monthly.txt|awk -v MPY=${MPY} '{print $1/MPY;}'`"
		rm -f ${yr}.tmp*

		# Transpose columns and rows in ${yr}.txt and
		#+ ${yr}_monthly.txt
R --no-save -q << EOF
fname='${yr}.txt'
data=read.table(fname,header=F)
nrow=nrow(data)/${SPY}
new_data=array(data[,1],c(${SPY},nrow))
ofname='${yr}_array.txt'
write.table(new_data, file = ofname, sep = ' ', 
col.names= FALSE, row.names = FALSE)
rm(list=ls())

# Repeat parsing for monthly data
fname='${yr}_monthly.txt'
data=read.table(fname,header=F)
nrow=nrow(data)/${MPY}
new_data=array(data[,1],c(${MPY},nrow))
ofname='${yr}_monthly_array.txt'
write.table(new_data, file = ofname, sep = ' ', 
col.names= FALSE, row.names = FALSE)
EOF
		rm -f ${yr}_monthly.txt ${yr}.txt

	done

	echo "Step 2. Generate mon. means corresponding to res. of coarse data"
	# Set up necessary directories. After resequence.R writes to
	#+ these directories then they will be removed and parsed to the
	#+ ANN processing directory.
	for ((i=1; i<=${CGC}; i++))
	do
		ii=`printf "%03i\n" ${i}`
		if [ -d ${ii} ]
		then
			rm -rf ${ii}/
		fi
		if [ ! -d ${ii} ]
		then
			mkdir -p ${ii}
		fi
	done

R --no-save -q << EOF
source("resequence.R")
EOF

	# Each column for each year output from resequence.R corresponds
	#+ to mean values for 12 months of data for one of the 62 divisions
	#+ of the total no. of CRUNCEP grid cells. These can be parsed
	#+ directly into 62 different directories as ANN input_training
	#+ & input_validation data. However the 6 hourly data has to
	#+ parsed using the ordering of the CRUNCEP cells saved in
	#+ cncp_recentered_idx.txt.
	rm -f [0-9][0-9][0-9][0-9]_monthly_array.txt \
	[0-9][0-9][0-9][0-9]_array.txt

	rm -f [0-9][0-9][0-9]/input_${variable}.original
	rm -f [0-9][0-9][0-9]/desired_${variable}.original
	rm -rf [0-9][0-9][0-9]/[0-9][0-9][0-9][0-9]/
	for ((yr=${YR_BEG_CNCP}; yr<=${YR_END_CNCP}; yr++))
	do

		echo "Parsing 6 hourly and monthly data for ANN year: ${yr}"

		for ((i=1; i<=${CGC}; i++))
		do

			test=`awk -v COL=${i} 'NR==6{print $COL;}' \
			${OUTDIR}/cncp_resequence_idx.txt`

			# If this coarse division has fine grid cells
			#+ in it then proceed
			if [[ ${test} > 0 ]]
			then

				ii=`printf "%03i\n" ${i}`

				# Parse the coarse means to the ANN
				#+ analysis directory
				awk -v GC=${i} 'NR>5{printf "%0.9f\n", $GC;}' \
				${yr}_coarse_monthly.out >> \
				${ii}/input_${variable}.original

				# Parse corresponding fine grid (6 hourly)
				#+ data to its own sub-dir
				gc=\
			$(ls ${ii}/[0-9][0-9][0-9][0-9]_desired_${yr}.out)
				ngc=`echo ${gc[*]} | wc -w`

				for ((j=0; j<${ngc}; j++))
				do

					jj=`echo ${gc[j]} | \
					sed "s/${ii}\///;s/_desired_.*//"`

					if [[ ! -d ${ii}/${jj}/ ]]
					then
						mkdir -p ${ii}/${jj}/
					fi

					# Parse corresponding fine, 6 hourly
					#+ data to the same place
					awk 'NR>5{printf "%0.9f\n", $1;}' \
					${ii}/${jj}_desired_${yr}.out >> \
					${ii}/${jj}/desired_${variable}.original

					rm -f ${ii}/${jj}_desired_${yr}.out

				done

			fi

		done
	done
	rm -f [0-9][0-9][0-9][0-4]_monthly.txt
	rm -f [0-9][0-9][0-9][5-9]_monthly.txt
	rm -f [0-9][0-9][0-9][0-4]_coarse_monthly.out
	rm -f [0-9][0-9][0-9][5-9]_coarse_monthly.out resequence.R knn.R
	popd
done

echo "Done!"
