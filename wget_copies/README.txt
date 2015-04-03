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
variables: tas, tasmin (day only), tasmax (day only), ps (or psl if necessary), huss, pr, ps, rlds, rsds, sfcWind (or uas & vas if necessary)

Please adhere to the following naming convention for wget scripts:
wget_PROJECT_MODEL_EXPERIMENT_TIMEFREQ_REALM_VARS_ENSEMBLE_VERSIONDATE.sh

To run the scripts, you will need to register with earth systems grid (you should have already done this to generate the wget script).  to execute:
$ sh [script name].sh -i

The -i is necessary to ignore some trusted certificate stuff that pops up


NOTE:  Sometimes random things happen and not all the files are downloaded.  After downloading through the wget script, scan through for files sizes of 0.  If there’s only a handful, you can go to the cmip5 website, select the variable you need (select “filter over text”) and then click the http link to download the file.  Note that you cannot use Safari for this download because of their built-in security settings (similar to how you have to run the wget script with -i).  Firefox works fine