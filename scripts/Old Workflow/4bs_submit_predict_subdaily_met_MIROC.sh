#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -pe omp 12
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -N diel_MIROC
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 4b_predict_subdaily_met_MIROC.R
