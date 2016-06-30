#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -N NLDAS
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 3_extract_NLDAS.R
