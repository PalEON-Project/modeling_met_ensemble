#!/bin/bash
# Simple loop for submitting parse_ann_in.sh to job scheduler

OUTDIR=`grep '^BASEDIR=' constants.txt | sed 's/.*=//'`/output/cruncep/
for variable in tair swdown qair lwdown psurf precipf wind
do

	ls -d ${OUTDIR}/${variable}/[0-9][0-9][0-9]/0* | \
	sed "s/.*${variable}\///;s/\// /" > gc.txt

	ls -d ${OUTDIR}/${variable}/[0-9][0-9][0-9]/1* | \
	sed "s/.*${variable}\///;s/\// /" >> gc.txt

	ngc=`awk 'END{print NR;}' gc.txt`

	for (( i=1; i<=${ngc}; i++ ))
	do

		DIVISION=`awk -v ROW=${i} 'NR==ROW{printf "%i", $1;}' gc.txt`
		FGCID=`awk -v ROW=${i} 'NR==ROW{printf "%i", $2;}' gc.txt`

		echo -n \
		"Submitting parse_ann_in.sh ${variable} ${DIVISION} ${FGCID}."

		qsub parse_ann_in.sh ${variable} ${DIVISION} ${FGCID}

	done

	rm -f gc.txt

done
