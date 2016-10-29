# -----------------------------------
# Script Information
# -----------------------------------
# Purpose: Create statistical models to predict subdaily meteorology from daily means
# Creator: Christy Rollinson, 28 October 2016
# Contact: crollinson@gmail.com
# -----------------------------------

# -----------------------------------
# Description
# -----------------------------------
# Apply the statistical models from step 3 to convert the daily, bias-corrected met 
# files from step 2 (daily means) and predict subdaily values.  This gets done by 
# filtering backwards in time starting with the present (where the trianing data is).
#
# There are ways to improve this and speed it up, but hopefully this works for now.
# We whould also probably think about applying this filter approach to the bias-
# correction step to avoid abrupt and unreasonable jumps in climate.
# -----------------------------------


# -----------------------------------
# Workflow
# -----------------------------------
# 0. Load libraries, set up file paths, etc
# ----- Loop through by ensemble member ----------
#    1. Load and format prediction data (1 file from step 2)
#       1.1 Load output file from bias correction (bias ensemble member)
#       1.2 select year we're working with (single file for output)
#    ----- Loop through by year ----------
#      2. Predict subdaily values for whole year, filtering backwards in time
#      3. Write annual output into .nc files 
#         - separate file for each year/ensemle member; 
#         - all met vars in one annual file (similar to pecan met structure)
#    ----- recycle steps 2 & 3 for all years in file ----------
# ----- recycle step 1 for all files for ensemble member ----------
# -----------------------------------


# -----------------------------------
# 0. Load libraries, set up file paths, etc
# -----------------------------------
# -----------------------------------


# -----------------------------------
# 1. Load bias-correction file
# -----------------------------------

# -----------------------------------
