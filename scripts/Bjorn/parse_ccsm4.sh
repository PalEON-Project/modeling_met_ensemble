#!/bin/bash
# ########################################################################### #
# /file                                                                       #
# @brief parse_ccsm4.sh parses the coarse model's history files for           #
# downscaling                                                                 #
#                                                                             #
# Requires:                                                                   #
# standard linux utilities, cdo, nco utilities (e.g., ncatted)                #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

### qsub arguments
#$ -S /bin/bash
#$ -wd /home/scratch/dietze_lab/bjorn/ccsm4_sd/scripts/
#$ -j y
#$ -m eas
#$ -M bjorn@climatemodeling.org

BASEDIR=`grep '^BASEDIR=' constants.txt | sed 's/.*=//'`
MODEL=`grep '^MODEL=' constants.txt | sed 's/.*=//'`
MODEL_CAP=`echo ${MODEL} | tr [a-z] [A-Z]`
COARSE_DIR=${BASEDIR}/../${MODEL}/past1000/ # Path to COARSE history files
OUTDIR=${BASEDIR}/output/${MODEL}/
MINLAT=`grep '^MINLAT=' constants.txt | sed 's/.*=//'`
MAXLAT=`grep '^MAXLAT=' constants.txt | sed 's/.*=//'`
MINLON=`grep '^MINLON=' constants.txt | sed 's/.*=//'`
MAXLON=`grep '^MAXLON=' constants.txt | sed 's/.*=//'`
CDO=/share/apps/cdo-1.5.1/bin/cdo # Path to CDO on ebi-cluster

if [[ -d ${OUTDIR} ]]
then
	echo "Removing ${OUTDIR}"
	rm -rf ${OUTDIR}
fi

#create land/sea mask
rm -f mask.nc mask.tmp*
ncks -d time,0,0 \
../${MODEL}/past1000/rlds_Amon_CCSM4_past1000_r1i1p1_085001-185012.nc mask.tmp1

${CDO} sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} mask.tmp1 mask.tmp2
${CDO} ltc,265 mask.tmp2 mask.tmp3
ncwa -a time,nb2 mask.tmp3 mask.nc
rm -f mask.tmp*

mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/lwdown
mkdir -p ${OUTDIR}/precipf
mkdir -p ${OUTDIR}/psurf
mkdir -p ${OUTDIR}/qair
mkdir -p ${OUTDIR}/swdown
mkdir -p ${OUTDIR}/tair
mkdir -p ${OUTDIR}/wind

pushd ${OUTDIR}/

ln -s ${BASEDIR}/mask.nc .
ncks mask.nc | grep '^lat\[' | grep 'lon\[' | grep -v '=0 W' | sed \
's/ rlds.*//;s/lon\[[0-9]\]=//;s/lon\[[0-9][0-9]\]=//;s/lat\[[0-9]\]=//;s/lat\[[0-9][0-9]\]=//' \
> ${MODEL}_idx_latlon.txt

echo "Parsing longwave radiation"
${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
${COARSE_DIR}/rlds_Amon_${MODEL_CAP}_past1000_r1i1p1_085001-185012.nc \
lwdown.tmp

${CDO} ifthen mask.nc lwdown.tmp lwdown.nc
rm -f lwdown.tmp
ncks -O -v lon,lat,time,rlds lwdown.nc lwdown.nc
ncrename -v rlds,lwdown lwdown.nc
${CDO} splityear lwdown.nc lwdown/
rm -f lwdown.nc

for year in `ls lwdown/*.nc | sed 's/lwdown\///;s/\.nc//'`
do
	${CDO} splitmon lwdown/${year}.nc lwdown/${year}_
	rm -f lwdown/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks lwdown/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > lwdown/${year}_${month}.txt

		rm -f lwdown/${year}_${month}.nc

	done
done

echo "Parsing precipitation flux"
${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
${COARSE_DIR}/pr_Amon_${MODEL_CAP}_past1000_r1i1p1_085001-185012.nc \
precipf.tmp

${CDO} ifthen mask.nc precipf.tmp precipf.nc
rm -f precipf.tmp
ncks -O -v lon,lat,time,pr precipf.nc precipf.nc
ncrename -v pr,precipf precipf.nc
${CDO} splityear precipf.nc precipf/
rm -f precipf.nc

for year in `ls precipf/*.nc | sed 's/precipf\///;s/\.nc//'`
do
	${CDO} splitmon precipf/${year}.nc precipf/${year}_
	rm -f precipf/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks precipf/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > precipf/${year}_${month}.txt

		rm -f precipf/${year}_${month}.nc

	done
done


echo "Parsing surface pressure"
${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
${COARSE_DIR}/ps_Amon_${MODEL_CAP}_past1000_r1i1p1_085001-185012.nc \
psurf.tmp

${CDO} ifthen mask.nc psurf.tmp psurf.nc
rm -f psurf.tmp
ncks -O -v lon,lat,time,ps psurf.nc psurf.nc
ncrename -v ps,psurf psurf.nc
${CDO} splityear psurf.nc psurf/
rm -f psurf.nc

for year in `ls psurf/*.nc | sed 's/psurf\///;s/\.nc//'`
do
	${CDO} splitmon psurf/${year}.nc psurf/${year}_
	rm -f psurf/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks psurf/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > psurf/${year}_${month}.txt

		rm -f psurf/${year}_${month}.nc

	done
done


echo "Parsing humidity"
j=1
for i in ${COARSE_DIR}/hus_Amon_${MODEL_CAP}_past1000_r1i1p1_*.nc
do

	ncks -O -d plev,2,2 ${i} qair.nc.${j}.tmp1

	${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
	qair.nc.${j}.tmp1 qair.nc.${j}.tmp2

	ncks -O -v lon,lat,time,hus qair.nc.${j}.tmp2 qair.nc.${j}

	rm -f qair.nc.${j}.tmp[1,2]

	j=$[j+1]

done

ncrcat qair.nc.[0-9]* qair.tmp
${CDO} ifthen mask.nc qair.tmp qair.nc
rm -f qair.tmp
ncrename -v hus,qair qair.nc
rm -f qair.nc.[0-9]*
${CDO} splityear qair.nc qair/
rm -f qair.nc

for year in `ls qair/*.nc | sed 's/qair\///;s/\.nc//'`
do
	${CDO} splitmon qair/${year}.nc qair/${year}_
	rm -f qair/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks qair/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > qair/${year}_${month}.txt

		rm -f qair/${year}_${month}.nc

	done
done


echo "Parsing shortwave radiation"
${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
${COARSE_DIR}/rsds_Amon_${MODEL_CAP}_past1000_r1i1p1_085001-185012.nc \
swdown.tmp

${CDO} ifthen mask.nc swdown.tmp swdown.nc
rm -f swdown.tmp
ncks -O -v lon,lat,time,rsds swdown.nc swdown.nc
ncrename -v rsds,swdown swdown.nc
${CDO} splityear swdown.nc swdown/
rm -f swdown.nc

for year in `ls swdown/*.nc | sed 's/swdown\///;s/\.nc//'`
do
	${CDO} splitmon swdown/${year}.nc swdown/${year}_
	rm -f swdown/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks swdown/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > swdown/${year}_${month}.txt

		rm -f swdown/${year}_${month}.nc

	done
done


echo "Parsing 2m air temp"
${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
${COARSE_DIR}/tas_Amon_${MODEL_CAP}_past1000_r1i1p1_085001-185012.nc \
tair.tmp

${CDO} ifthen mask.nc tair.tmp tair.nc
rm -f tair.tmp
ncks -O -v lon,lat,time,tas tair.nc tair.nc
ncrename -v tas,tair tair.nc
${CDO} splityear tair.nc tair/
rm -f tair.nc

for year in `ls tair/*.nc | sed 's/tair\///;s/\.nc//'`
do
	${CDO} splitmon tair/${year}.nc tair/${year}_
	rm -f tair/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks tair/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > tair/${year}_${month}.txt

		rm -f tair/${year}_${month}.nc

	done
done


echo "Parsing u & v wind components"
for vec in u v # Loop over u and v components
do

	for i in ${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_*.nc
	do

		ncks -v ${vec}a -d plev,0,0 ${i} ${i}.tmp1

		${CDO} \
		-O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
		${i}.tmp1 ${i}.tmp2

	done

	# Concatenate all time steps
	ncrcat \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_*.nc.tmp2 \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp1

	rm -f ${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_*.nc.tmp*

	# Remove orphaned plev dimension
	ncwa -a plev \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp1 \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp2

	# Extract geographic domain
	${CDO} -O -r sellonlatbox,${MINLON},${MAXLON},${MINLAT},${MAXLAT} \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp2 \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp3

	# Mask oceans
	${CDO} ifthen mask.nc \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp3 \
	${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.nc

	rm -f ${COARSE_DIR}/${vec}a_Amon_CCSM4_past1000_r1i1p1_all.tmp*

done

# Append ua variable to file with va variable
ncks -A ${COARSE_DIR}/ua_Amon_CCSM4_past1000_r1i1p1_all.nc \
${COARSE_DIR}/va_Amon_CCSM4_past1000_r1i1p1_all.nc

echo "Calculating wind speed"
ncap -O -s 'Wind=abs(va)/sin(atan(abs(va/ua)))' \
${COARSE_DIR}/va_Amon_CCSM4_past1000_r1i1p1_all.nc wind.tmp1
rm -f ${COARSE_DIR}/[u,v]a_Amon_CCSM4_past1000_r1i1p1_all.nc

ncks -O -v Wind wind.tmp1 wind.nc

ncatted -a units,Wind,o,c,'m s-1' \
-a long_name,Wind,o,c,'Wind speed at the lowest atmosphere level' \
-a standard_name,Wind,o,c,Wind \
-a description,Wind,o,c,'Wind=abs(V)/sin(atan(abs(V/U)))' wind.nc

${CDO} splityear wind.nc wind/
rm -f wind.nc wind.tmp1

for year in `ls wind/*.nc | sed 's/wind\///;s/\.nc//'`
do
	${CDO} splitmon wind/${year}.nc wind/${year}_
	rm -f wind/${year}.nc
	for (( mon=1; mon<=12; mon++ ))
	do

		month=`printf %02i ${mon}`
		ncks wind/${year}_${month}.nc | grep '^time\[' | \
		grep lat | grep -v '1e+20' | sed 's/.*=//;s/ .*//' | \
		awk '{printf "%0.9f\n", $1;}' > wind/${year}_${month}.txt

		rm -f wind/${year}_${month}.nc

	done

done


popd

CGC=`awk 'END{print NR;}' ${OUTDIR}/${MODEL}_idx_latlon.txt` # No. of pieces to divide the domain into
for variable in lwdown precipf psurf qair swdown tair wind
do

	pushd ${OUTDIR}/${variable}/

	for ((cgc=1; cgc<=${CGC}; cgc++))
	do

		ii=`printf "%03i\n" ${cgc}`
		if [ -d ${ii} ]
		then
			rm -rf ${ii}/
		fi
		if [ ! -d ${ii} ]
		then
			mkdir -p ${ii}
		fi

		# Concatenate all rows matching the desired grid cell
		cat [0-9][0-9][0-9][0-9]_[0,1][0-9].txt | \
		awk -v cgc=${cgc} -v CGC=${CGC} \
		'(NR-cgc)%CGC==0{print $0}' > ${ii}/testing_${variable}.original

	done

	rm -f [0,1][0-9][0-9][0-9]_[0,1][0-9].txt

	popd

done

echo "Done!"
