#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -pe omp 12
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -N predict_diel
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 4_predict_subdaily_met.R
