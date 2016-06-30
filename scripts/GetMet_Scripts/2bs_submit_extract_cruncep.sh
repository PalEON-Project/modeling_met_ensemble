#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/GetMet_Scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N CRUNCEP
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/GetMet_Scripts/

sh 2b_extract_cruncep.sh
