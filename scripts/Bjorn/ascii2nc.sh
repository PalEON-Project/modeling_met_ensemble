#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# @brief ascii2nc.sh submits ascii2nc.R to the job scheduler                  #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date September 11, 2012                                                    #
# ########################################################################### #

### qsub arguments
#$ -wd /home/scratch/dietze_lab/bjorn/ccsm4_sd/scripts/
#$ -j y
#$ -S /bin/bash
#$ -m eas
#$ -M bjorn@climatemodeling.org

echo "Creating netcdf files"
R --no-save -q << EOF
source("ascii2nc.R")
EOF

#for variable in tair swdown qair lwdown psurf precipf wind
for variable in psurf precipf wind
do
	for i in output/ann/${variable}/netcdf/*.nc
	do
		ncrename -a missing_value,_FillValue ${i}
	done
done

# ebi-cluster has an old version R and netcdf fill values will
#+ have to be specified manually:
#+ ncrename -a missing_value,_FillValue output/ann/VARNAME/netcdf/FNAME.nc
