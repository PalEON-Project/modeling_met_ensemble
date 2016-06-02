#!bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N Extract_CRU
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

sh 2b_extract_cruncep.sh
