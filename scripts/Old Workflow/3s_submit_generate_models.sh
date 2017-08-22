#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=48:00:00
#$ -N subdaily_models
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 3_generate_subdaily_models.R
