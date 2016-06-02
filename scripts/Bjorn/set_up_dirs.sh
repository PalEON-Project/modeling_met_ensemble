#!/bin/bash
# ########################################################################### #
# \file                                                                       #
# diag.sh                                                                     #
# @brief Replicates directories and data necessary prior downscaling          #
#                                                                             #
# Contact:                                                                    #
# @author Bjorn Brooks, bjorn@climatemodeling.org                             #
#                                                                             #
# @date August 24, 2012                                                       #
# ########################################################################### #

MODEL=`grep '^MODEL=' constants.txt | sed 's/.*=//'`

pushd ../

ln -s /home/scratch/dietze_lab/bjorn/cruncep/ .
ln -s /home/scratch/dietze_lab/bjorn/cmip5/${MODEL}/ .
ln -s /home/scratch/dietze_lab/bjorn/mpi_sd/normalize/ .
ln -s /home/scratch/dietze_lab/bjorn/mpi_sd/neural_networks/ .

popd

echo "Directories created"
