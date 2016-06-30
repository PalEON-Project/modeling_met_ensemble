#!/bin/bash
# Get NLDAS data for individual sites and convert to .nc
# LDAS page: http://disc.sci.gsfc.nasa.gov/hydrology/overview
# NLDAS ftp: ftp://hydro1.sci.gsfc.nasa.gov

# Code to go from .grb to .nc found here:
# http://p-martineau.com/converting-grib-to-netcdf/
# 
# #!/bin/bash
# #the script must be placed outside of the folder where the data is
# #therefore script.sh folder_grib and folder_netcdf are in the same location
# folderin='folder_grib' #folder where the grib are
# folderout='folder_netcdf' #folder to store netcdf
# cd $folderin # cd to that folder
# #now loop through all the grib files in there
# for file in *.grib
# do
# echo 'test '$file #display the current selected file
# ncl_convert2nc $file -o "../"$folderout -e grb #convert it using ncl
# done

module load cdo
module load nco

CDO="/project/earth/packages/cdo-1.6.3rc2/bin/cdo"

# FTP URL
nldas_base=ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/
dir_raw=/projectnb/dietzelab/paleon/met_ensemble/data/full_raw/NLDAS
dir_nc=nc_files
nldas_out=/projectnb/dietzelab/paleon/met_ensemble/data/paleon_domain

# Make the output directories
mkdir -p ${nldas_out}/paleon1/NLDAS
mkdir -p ${nldas_out}/alaska/NLDAS

# Bounding box of the PalEON 1 Domain
MINLAT_P1=35
MAXLAT_P1=50
MINLON_P1=-100
MAXLON_P1=-60

# # Bounding Box of Alaska (more or less)
# # NOTE: ALASKA is not in the NLDAS, so we'll have to get the GLDAS for that one
# MINLAT_AK=55
# MAXLAT_AK=72
# MINLON_AK=-165
# MAXLON_AK=-140

# Min and max years of NLDAS dataset
minyr=1979
maxyr=2015

mkdir -p ${dir_raw}/$dir_nc
for (( YEAR=$minyr; YEAR<=$maxyr; YEAR++ )); do
   mkdir -p ${dir_raw}/${dir_nc}/${YEAR}
   
   paleon1=${nldas_out}/paleon1/NLDAS/${YEAR}
   alaska=${nldas_out}/alaska/NLDAS/${YEAR}
   
   mkdir -p $paleon1
   mkdir -p $alaska
   # Download all files for that year
   wget -r -nH --cut-dirs=4 ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_FORA0125_H.002/$YEAR
   
   pushd $YEAR
      days=(*)
   popd

   for DAY in ${days[@]}; do
      # Reconstruct the file structure in our nc folder
      path_now=${dir_raw}/${dir_nc}/${YEAR}/${DAY}
      mkdir -p $path_now
      
      pushd ${YEAR}/${DAY}
         mkdir -p ${paleon1}/${DAY}
         mkdir -p ${alaska}/${DAY}
         
         for FILE in *.grb; do
            # convert to .nc
            ncl_convert2nc $FILE -o  ${path_now} -e grb #convert it using ncl
            
            # crop for paleon 1
            ${CDO} sellonlatbox,${MINLON_P1},${MAXLON_P1},${MINLAT_P1},${MAXLAT_P1} ${path_now}/${FILE}.nc \
            ${paleon1}/${DAY}/${FILE}.nc
         done # end file
      popd 
   done # end day    
   
   # Assuming the year processed okay, get rid of the pain-in-the-butt .grib file
   # and compress the full .nc file
   rm -rf $YEAR
   
   pushd ${dir_raw}/${dir_nc}
      tar -jcvf ${YEAR}.tar.bz2 $YEAR
      rm -rf $YEAR
   popd
done # end year