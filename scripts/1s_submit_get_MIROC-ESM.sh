#!bin/bash
#$ -wd /projectnb/dietzelab/paleon/met_ensemble/scripts/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=24:00:00
#$ -N MIROC-ESM
#cd /projectnb/dietzelab/paleon/met_ensemble/scripts/

sh 1_get_data_MIROC-ESM_ceda.sh