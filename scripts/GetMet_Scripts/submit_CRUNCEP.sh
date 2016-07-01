#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/GetMet_Scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=2:00:00
#$ -N site_CRU
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/GetMet_Scripts/

R CMD BATCH extract.CRUNCEP.R
