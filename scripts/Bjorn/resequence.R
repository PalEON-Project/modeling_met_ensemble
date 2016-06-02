## \file 
# @file resequence.R
# @author Author: Bjorn Brooks (bjorn@climatemodeling.org)
# @brief Octave routine that resequences the lat,lon ordering of the
# finely gridded data (CRUNCEP) by finding all the coordinate points of
# the fine grid that fit within the corresponding rectangular coordinate
# box of the coarse data. Box may have a differing numbers of fine grid
# cells that fit within them due to lakes and coastlines.
# @date April 26, 2012

rm (list=ls())
source("knn.R")
model='ccsm4'
variable='lwdown'
# half of the longitudinal width of the average grid cell box
clon_pm=
# half of the latitudinal height of each grid cell box
clat_pm=
# No. of fine (CRUNCEP) grid cells in domain
fgc=1630
mpy=12 # months per year

flon_min=-100
flon_max=-60
flat_min=35
flat_max=50

# Load target grid lat lon indexes
fname='../cncp_idx_latlon.txt' # Lat,Lon coordinates from fine grid
fll=read.table(fname, header=FALSE) # Lat,Lon coordinates from fine grid

# Load input grid lat lon indexes
fname=paste('../../', model, '/', model, '_idx_latlon.txt', sep='')#coarse grid
cll=read.table(fname, header=FALSE)
cll[,2]=cll[,2]-360 # Adjust longitude so it corresponds to CRUNCEP values

ne_fll=nrow(fll) # No. of elements in fine grid
ne_cll=nrow(cll) # No. of elements in coarse grid
fll_row=1:ne_fll  # Row number in fll locating its lat,lon coordinates

flat=fll[,1] # Fine grid latitude values
flon=fll[,2] # Fine grid longitude values
seq_idx=array(0,c(ne_fll,ne_cll)) # Initialize array

# Step through each column (which represents a coarse grid cell) and
#+ assign all corresponding fine grid cells.  Corners and edges will
#+ include everything within the coarse grid cell and out to the 
#+ edge/corner so that no fine grid cells are skipped.
for (i in 1:ne_cll) {

	clat=cll[i,1]
	clon=cll[i,2]

	# If this coarse grid cell is in the SW corner then assign all
	#+ fine grid cells that are within its coordinate bounds ...
	if (clat == min(cll[,1]) && clon == min(cll[,2])) { # SW corner
		tmp=fll_row[flat >= flat_min & flat < clat+clat_pm & 
		flon >= flon_min & flon < clon+clon_pm]
	}
	else if (clat == max(cll[,1]) && clon == max(cll[,2])) { # NE corner
		tmp=fll_row[flat >= clat-clat_pm & flat <= flat_max & 
		flon >= clon-clon_pm & flon <= flon_max]
	}
	else if (clat == min(cll[,1]) && clon == max(cll[,2])) { # SE corner
		tmp=fll_row[flat >= flat_min & flat < clat+clat_pm & 
		flon >= clon-clon_pm & flon <= flon_max]
	}
	else if (clat == max(cll[,1]) && clon == min(cll[,2])) { # NW corner
		tmp=fll_row[flat >= clat-clat_pm & flat <= flat_max & 
		flon >= flon_min & flon < clon+clon_pm]
	}
	else if (clat == min(cll[,1])) { # Remaining S cells
		tmp=fll_row[flat >= flat_min & flat < clat+clat_pm & 
		flon >= clon-clon_pm & flon < clon+clon_pm]
	}
	else if (clat == max(cll[,1])) { # Remaining N cells
		tmp=fll_row[flat >= clat-clat_pm & flat <= flat_max & 
		flon >= clon-clon_pm & flon < clon+clon_pm]
	}
	else if (clon == min(cll[,2])) { # Remaining W cells
		tmp=fll_row[flat >= clat-clat_pm & flat < clat+clat_pm & 
		flon >= flon_min & flon < clon+clon_pm]
	}
	else if (clon == max(cll[,2])) { # Remaining E cells
		tmp=fll_row[flat >= clat-clat_pm & flat < clat+clat_pm & 
		flon >= clon-clon_pm & flon <= flon_max]
	}
	else { # All other central cells
		tmp=fll_row[flat >= clat-clat_pm & flat < clat+clat_pm & 
		flon >= clon-clon_pm & flon < clon+clon_pm]
	}
	# seq_idx is the grid cell index numbers for fine grid cells. Each
	#+ number represents a fine grid cell and is placed in a column
	#+ whos number identifies the coarse grid cell it belongs to
	if (length(tmp)>0) {
		seq_idx[1:length(tmp),i]=tmp
	}
	# Remove rows that have already been selected to avoid duplicates
	#fll_row=setdiff(fll_row,tmp)
	fll_row[fll_row=tmp]=0

}

# Collect up all the leftover fine (fll) grid cells not yet assigned
#+ and assign them to the nearest coarse cell (cll)
for (i in 1:length(fll_row)) { # Loop over remaining unassigned fine grid cells
	# Find the column index locating the coarse grid cell this fine
	#+ grid cell should be assigned to
        fll_i=matrix(fll[fll_row[i],],c(1,2)) # Current lat lon coord's
	# If this fll value is zero then it was already assigned a
	#+ corresponding ccll and it will be skipped
	if (! is.na(fll_i[1,1]>0)) {
		cll_col_idx=knn(cll, fll_i, 1)

		# Insert this leftover grid cell into seq_idx
		tarr=seq_idx[,cll_col_idx] # Temp array for easy subsetting
		tarr=append(tarr,fll_row[i]) # Add new idx
		tarr=tarr[tarr>0]            #+ after last non-zero element
		# Rewrite seq_idx with original values plus new value
		seq_idx[1:length(tarr),cll_col_idx]=tarr
	}
}
rm(tarr)

# Save the fll (CRUNCEP) lat,lon row indexes just in case
if ( variable == 'lwdown') {
	fname='../cncp_resequence_idx.txt'
	print('resequence.m: Saving row indexes.')
	print('              Note that n-th column lists the rows from')
	print('              cncp_idx_latlon.txt that correspond to n-th')
	print('              row (grid cell) of MODEL_idx_latlon.txt')
	write.table(seq_idx, file = fname, sep = ' ', 
	col.names= FALSE, row.names = FALSE)
}

for (year in 1901:2010) {

	# Find spatial mean of monthly means
	#+ Load the already-monthly-averaged 0.5-degree CRUNCEP data before
	#+ finding the spatial means corresponding to the coarsely gridded data
	#+ The data being loaded has the following array structure:
	#+ month_1  cell_1_value  cell_2_value  ... cell_fgc_value
	#+ month_2  cell_1_value  cell_2_value  ... cell_fgc_value
	fname=paste(year, '_monthly_array.txt', sep='')
	var_data=read.table(fname, header=FALSE)

	fll_out_mm=array(0,c(mpy,ne_cll))
	for (mo in 1:mpy) {

		for (i in 1:ne_cll) {

			# Calculte averages of fine grids corresponding to
			#+ location of each coarse cell. Note that fll grid
			#+ cells have been re-sequenced according to seq_idx,
			#+ which can have different numbers of fll grid cells
			#+ per cll grid box. The original sequence will have
			#+ to be restored afer ANN analysis using seq_idx.
			# Some seq_idx cols may have no indx row nums
			if (seq_idx[1,i]>0) {

				# index of non-null rows
				idx=seq_idx[seq_idx[,i]>0,i]
				var_tmp=var_data[mo,idx]

				# Calcuate mean across grid cells for
				#+ each time step and save corresponding
				#+ month and lat,lon coordinates
				# mean of selected fll gridcells for that month
				if (length(var_tmp) > 1) {
					fll_out_mm[mo,i]=apply(var_tmp,1,mean)
				} else {
					fll_out_mm[mo,i]=var_tmp
				}
			}

			else { # If no reference to any rows set mean to -99999

				fll_out_mm[mo,i]=-99999

			}

		}
	}

	fname=paste(year, '_coarse_monthly.out', sep='')
	print(year)

	write.table(fll_out_mm, file = fname, sep = ' ', 
	col.names= FALSE, row.names = FALSE)

	# Find spatial mean of 6 hourly values ('desired' data for ANN)
	#+ and write to file
	fname=paste(year, '_array.txt', sep='')
	var_data_sixhourly=read.table(fname, header=FALSE)
	for (i in 1:ne_cll) {

		# Only process cols of seq_idx that reference index row nums
		if (seq_idx[1,i] > 0) {

			# Loop over all rows in this col of seq_idx and
			#+ output each fll (6 hourly) grid cell to its own
			#+ file but in the corresponding cll direcctory
			for (j in 1:length(seq_idx[,i])) {

				if (seq_idx[j,i]>0) {

					fll_out=
					var_data_sixhourly[,seq_idx[j,i]]

					fname=paste(sprintf('%03i',i), '/', 
					sprintf('%04i',seq_idx[j,i]), 
					'_desired_', year, '.out', sep='')

					# Save each array of 6 hourly
					#+ data from a singe grid cell to
					#+ the directory corresponding to
					#+ its cll (coarse lat lon) box
					write.table(fll_out, file = fname, 
					sep = ' ', col.names= FALSE, 
					row.names = FALSE)

				}

			}

		}

	}
}
