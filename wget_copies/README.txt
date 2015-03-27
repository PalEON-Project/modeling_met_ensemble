README for folder: wget_copies
Christy Rollinson, crollinson@gmail.com
date created: 27 March 2014
date modified: 27 March 2014

This folder houses a copy of all wget files used to download CMIP5 met data.  Original wget files used to download data is stored in the appropriate folder under “data” which is ignored from github.

wget scripts were generated at http://pcmdi9.llnl.gov/esgf-web-fe/live with the following search constraints (one or more used as appropriate):
experiment: past1000, historical
realm: atmos
time frequency: day, 6 hrly, 3 hrly
model: MRI-CGCM3

Please adhere to the following naming convention for wget scripts:
wget_PROJECT_MODEL_EXPERIMENT_TIMEFREQ_REALM_VARS_ENSEMBLE_VERSIONDATE.sh
