#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -N Debias_CCSM4
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 2a_bias_correct_CCSM4.R
