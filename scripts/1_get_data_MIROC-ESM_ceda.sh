#!bin/bash
# --------------------------
# Script Information
# --------------------------
# Purpose: Download the data for a CMIP5 ensemble: MIROC-ESM
# Creator: Christy Rollinson 6 May 2016
# Contact: crollinson@gmail.com
# --------------------------

# --------------------------
# Description:
# --------------------------
# Download a CMIP5 ensemble member that has the variables we need (for PalEON) from the
# past 1000 experiment (850-1849) and historical runs (1850-2005)
#
# Note: This particular model and variables have been pre-screened and paths hard coded 
#       based on what was available at the following URL on 6 May 2016
#       http://browse.ceda.ac.uk/browse/badc/cmip5/data/cmip5/output1
#
#       I'm sure there are better ways to scrape the website and set up a dynamic system 
#       to get the highest temporal resolution possible for each variable, but I'm too lazy
#       to do that right now
# 
#       Right now I'm only pulling the basic r1i1p1 simulations because I think these were
#       the baseline runs (if I remember correclty)
#
# Note: ceda_user and cedapw are defined in my .bashrc script because these shouldn't be 
#       publicly available on github! (see lines 58, 85, 124)
# --------------------------

# --------------------------
# Workflow overview
# --------------------------
# 1. Get data from past 1000 experiment (850-1849)
#    a. Whatever's available at daily resolution
#    b. remaining variables at monthly resolution
# 2. Get data from historical experiment (1850-2005)
#    a. Daily resolution (sub-day data only started at 1950 & we can get better data for then)
# --------------------------

# get my ceda username and password from my .bashrc scrpt
source /usr2/postdoc/crolli/ceda_creds.sh

# Set the baseline working directory
out_base=/projectnb/dietzelab/paleon/met_ensemble/data/full_raw/MIROC-ESM/
# 
# # --------------------------
# # Past 1000 Experiment: 850-1849
# # --------------------------
# #    Make sure we have our base directory
# #    outdir=${out_base}p1000/
# #    if [ ! -d ${outdir} ]
# #    then
# #       mkdir -p ${outdir}
# #    fi
# # 
# #    -----------
# #    Getting daily data
# #    -----------
# #       echo "---------- Loop 1: p1000 daily ----------"
# #       Define file path & variables of interest
# #       path_p1000_day=/badc/cmip5/data/cmip5/output1/MIROC/MIROC-ESM/past1000/day/atmos/day/r1i1p1/latest
# #       var_p1000_day=(tas tasmax tasmin pr psl huss sfcwind)
# # 
# #       Make a directory for the daily data so we can keep things separate
# #       if [ ! -d ${outdir}day/ ]
# #       then
# #          mkdir -p ${outdir}day/
# #       fi
# # 
# #       Change to the daily p1000 working directory
# #       pushd ${outdir}day/
# # 
# #       loop through daily variables
# #       for var in ${var_p1000_day[@]}
# #       do
# #          echo $var
# #          wget -r -nH --cut-dirs=14 ftp://${ceda_user}:${cedapw}@ftp.ceda.ac.uk/${path_p1000_day}/${var}
# #       done
# #       popd # leave the daily directory
# #    -----------
# #    
# #    
# #    -----------
# #    loop through monthly
# #    -----------
# #       echo "---------- Loop 2: p1000 monthly ----------"
# #       Define the file path and variable list for monthly data
# #       path_p1000_mo=/badc/cmip5/data/cmip5/output1/MIROC/MIROC-ESM/past1000/mon/atmos/Amon/r1i1p1/latest
# #       var_p1000_mo=(ps rlds rsds)
# # 
# #       Make a directory for the monthly data so we can keep things separate
# #       if [ ! -d ${outdir}month/ ]
# #       then
# #          mkdir -p ${outdir}month/
# #       fi
# # 
# #       Change to the monthly p1000 working directory
# #       pushd ${outdir}month/
# # 
# #       loop through daily variables
# #       for var in ${var_p1000_mo[@]}
# #       do
# #          echo $var
# #          wget -r -nH --cut-dirs=14 ftp://${ceda_user}:${cedapw}@ftp.ceda.ac.uk/${path_p1000_mo}/${var}
# #       done
# #       popd # leave the monthly directory
# #    -----------
# # --------------------------
# # 
# # 
# --------------------------
# Historical experiment: 1850-2010
# NOTE: only have daily for the key 1850-1900 period!
# --------------------------
   Make sure we have our base directory
   outdir=${out_base}historical/
   if [ ! -d ${outdir} ]
   then
      mkdir -p ${outdir}
   fi

   # -----------
   # Getting daily data
   # -----------
      echo "---------- Loop 3: historical daily ----------"
      Define file path & variables of interest
      path_hist_day=/badc/cmip5/data/cmip5/output1/MIROC/MIROC-ESM/historical/day/atmos/day/r1i1p1/latest
      var_hist_day=(tas tasmin tasmax pr psl huss rlds rsds sfcWind)

      Make a directory for the daily data so we can keep things separate
      if [ ! -d ${outdir}day/ ]
      then
         mkdir -p ${outdir}day/
      fi

      Change to the daily p1000 working directory
      pushd ${outdir}day/

      loop through daily variables
      for var in ${var_hist_day[@]}
      do
         echo $var
         wget -r -nH --cut-dirs=14 ftp://${ceda_user}:${cedapw}@ftp.ceda.ac.uk/${path_hist_day}/${var}
      done
      popd # leave the daily directory
   # -----------

   # -----------
   # Getting sub-daily data
   # -----------
      echo "---------- Loop 4: historical sub-daily ----------"
      # Define file path & variables of interest
      path_hist_3h=/badc/cmip5/data/cmip5/output1/MIROC/MIROC-ESM/historical/3hr/atmos/3hr/r1i1p1/latest
      var_hist_3h=(tas pr ps huss rlds rsds)

      # Make a directory for the daily data so we can keep things separate
      if [ ! -d ${outdir}3h/ ]
      then
         mkdir -p ${outdir}3h/
      fi

      # Change to the daily p1000 working directory
      pushd ${outdir}3h/

      # loop through daily variables
      for var in ${var_hist_3h[@]}
      do
         echo $var
         wget -r -nH --cut-dirs=14 ftp://${ceda_user}:${cedapw}@ftp.ceda.ac.uk/${path_hist_3h}/${var}
      done
      popd # leave the daily directory
   # -----------

# --------------------------
