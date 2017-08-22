#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -q "geo*"
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -N Debias_BCC
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

R CMD BATCH 2d_bias_correct_BCC.R
