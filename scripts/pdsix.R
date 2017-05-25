# Calculate PDSI from Palmer Z Scores
# Translated from code written by Ben Cook
#     from Table 12, Palmer, 1965
# Author: Christy Rollinson (crollinson@mortonarb.org); translated from B. Cook matlab code
#
# Inputs:
#   1. z = Palmer 'Z index"; vector with of monthly values; length evenly divisible by 12
# Outputs
#   1. x  = PDSI; vector with same length as Z
#   2. xm = modified PDSI
#   3. x4 = PDHI; defined as x4[t] = 0.897*x4[t-1] + Z/3
#           circumvents probability of switching between three 
#           alternative version of the index

# Workflow: 
# 1. Check inputs; setup vectors etc
# 2. Initialize values for first month
# 3. Loop through remaining time series
#	 - Identify wet/dry/normal periods
# 4. Format & return outputs
pdsix <- function(z) {

	# ------------------------------------------
	# 1. Check inputs; set up pieces
	# ------------------------------------------
	# Insert checks to make sure z is divisible by 12
	
	# Setting up place holder variables
	nz = length(z)
	Uw=rep(NA, length=nz) # Effective wetness
	Ud=rep(NA, length=nz) # Effective dryness
	V=rep(NA, length=nz) # numerator for probability of ending wet/dry
	Q=rep(NA, length=nz) # denominator for probability of ending wet/dry
	Ze=rep(NA, length=nz) # Z-value needed to end wet/dry
	Pe=rep(NA, length=nz) # Probability htat drought has ended
	x1=rep(0, length=nz)
	x2=rep(0, nz)
	x3=rep(0, nz)
	x4=rep(0, nz)
	
	# making a shortcut for z3
	z3=z/3
	
 	# Some Allocation vectors
 	LLd = rep(0, nz)
 	LLw = rep(0, nz)
 	LLn = vector(1, nz)
 	pullx1 = rep(1, nz)
 	pullx2 = rep(1, nz)
 	pullnorm = rep(1, nz)
	# ------------------------------------------


	# ------------------------------------------
	# 2. Initialize the First month of the time series
	# ------------------------------------------
	Uw[1] = NA # effective wetness; applies only with pre-existing drought
	Ud[1] = NA # effective dryness; applies only with pre-existing wet
	V[1] = 0 # numerator for probability (end of wet/dry)
	Ze[1] = NA # z-value needed to end wet/dry
	Q[1] = NA # denominator
	Pe[1] = 0 # probability that drought/wet has ended
	x1[1] = max(0, z3[1])
	x2[1] = min(0, z3[1])
	x3[1] = z3[1]
	x4[1] = z3[1]
	
	if(z3[1] <= -1.0 ){
		status="dry"
		LLd[1]=1 # Logical for dry period month
		newdry=1 # Started new dry period
		firstwet=0 # have not encountered wet yet
		nump=0 # Number of active previous months for probability
	} else if(z3[1] >= 1.0) {
		status="wet"
		LLw[1] = 1
		newwet=1
		firstdry=0
		nump=0
	} else {
		status="normal"
		LLn[1]=1
		x3[1]=0
		newnorm=1
		pullnorm[1]=1
	}
	# ------------------------------------------


	# ------------------------------------------
	# 3. Looping through the rest of the time series and doing the calculation
	# ------------------------------------------
	for(i in 2:nz){
		x4[i] = 0.897*x4[i-1] + z3[i]
		
		if(status == "normal") {
			x1[i] = max(0, 0.897*x1[i-1] + x3[i])
			x2[i] = min(0, 0.897*x2[i-1] + x3[i])
			x3[i] = 0
			
			# Figure out the new state
			if(x1[i] < 1.0 & x2[I] >-1.0 ) {  # No new drought 
				LLn[i] = 1
				pullnorm[i] = 1
			} else if(x1[i] >= 1.0 ) { # new wet period 
				status = "wet"; newwet=1; firstdry=0; nump=0
				LLw[i] = 1
				x3[i] = x1[i]
			} else if(x2[I] <= -1.0) {# new dry period 
				status = "dry"; newdry=1; firstwet=0; nump=0
				LLd[i] = 1
				x3[i] = x2[i]
			}
		} else if(status == "wet") { 
			x1[i] = max(0, 0.897*x1[i-1] + z3[i])
			x2[i] = min(0, 0.897*x2[i-1] + z3[i])
			x3[i] = 0.897*x3[i-1] + z3[i]
			Ud[i] = Z[i] - 0.15
			if(newwet==1) {
				newwet=0
				x1[i] = 0
			}
			if(firstdry==0){ # have not previously hit effective dry Ud
				if(Ud[i] >= 0){
					Ud[i] = NaN
					Pe[i] = 0
				} else {
					firstdry=1 # Flag that we have now hit the first effective dry Ud
					Ze[i] = -2.691 * x3[i-1] + 1.5 # Z-value needed to end wet period
					# Calculate probability that wet period has ended
					V[i] = Ud[i]
					Q[i] = Ze[i]
					Pe[i] = V[i] * 100/Q[i]
				}
			} else { # firstdry==1
				Ze[i] = -2.691 * x3[i-1] + 1.5 # Z-value needed to end wet
				
				# Calculate probability that wet period has ended
				V[i] = V[i-1] + Ud[i]
				Q[i] = V[i-1] + Ze[i]
				Pe[i] = V[i] * 100/Q[i]				
			} # End first dry case
			
			# Bound out probability of ending wet period to 0 - 100
			if(Pe[i] < 0) {
				Pe[i] = 0
			} else if (Pe[i] >= 100 ){
				Pe[i] = 100
			}
			
			if(Pe[i] == 0) { # zero chance that our wet period is over
				LLw[i] = 1; nump=0
				x1[i] = 0
				if(Pe[i-1]>0){ # If we fizzled out when trying to end the period, x1 & x2 get reset to 0
					x1[i] = 0
					x2[i] = 0
					firstdry=0; Ze[i] = NaN; V[i] = 0; Q[i] = NaN; nump=0
				}
			} else if(Pe[i] = 100 ) { # wet period ends!
				if(x2[i] <= -1.0 ) { # going straight into a drought
					status="dry"; LLd[i] = 1; newdry=1; firstwet=0
					x3[i] = x2[i]
					pullx2[(i-nump):i] <- 1
				} else if (x1[i] >= 1.0) { 
					stop("Something's Wrong: Ended wet period with start of a new one!")
				} else { # going to a "normal" stateof being
					status="normal"; LLn[i] = 1
					x3[i] = 0
					pullnorm[(i-nump):i] <- 1
				}
				nump=0
			}else { # We're getting closer to ending our wet period, but it's not over yet
				nump=nump+1
			}
		} else if (status == "dry") {
			x1[i] = max(0, 0.897*x1[i-1] + z3[i])
			x2[i] = min(0, 0.897*x2[i-1] + z3[i])
			x3[i] = 0.897*x3[i-1] + z3[i]
			Uw[i] = Z[i] + 0.15

			if(newdry==1) {
				newdry=0
				x2[i] = 0
			}
			if(firstwet==0){ # have not previously hit effective wet U2
				if(Uw[i] >= 0){
					Uw[i] = NaN
					Pe[i] = 0
				} else {
					firstwet=1 # Flag that we have now hit the first effective wet Ud
					Ze[i] = -2.691 * x3[i-1] - 1.5 # Z-value needed to end wet period
					# Calculate probability that wet period has ended
					V[i] = Uw[i]
					Q[i] = Ze[i]
					Pe[i] = V[i] * 100/Q[i]
				}
			} else { # firstwet==1
				Ze[i] = -2.691 * x3[i-1] - 1.5 # Z-value needed to end wet
				
				# Calculate probability that wet period has ended
				V[i] = V[i-1] + Uw[i]
				Q[i] = V[i-1] + Ze[i]
				Pe[i] = V[i] * 100/Q[i]				
			} # End first wet case
			
			# Bound out probability of ending wet period to 0 - 100
			if(Pe[i] < 0) {
				Pe[i] = 0
			} else if (Pe[i] >= 100 ){
				Pe[i] = 100
			}
			
			if(Pe[i] == 0) { # zero chance that our dry period is over
				LLd[i] = 1; nump=0 # Ben's code has LLw[i], but I htink we should be inverse of what we did for wet
				x2[i] = 0
				if(Pe[i-1]>0){ # If we fizzled out when trying to end the period, x1 & x2 get reset to 0
					x1[i] = 0
					x2[i] = 0
					firstwet=0; Ze[i] = NaN; V[i] = 0; Q[i] = NaN; nump=0
				}
			} else if(Pe[i] = 100 ) { # wet period ends!
				if(x1[i] >= 1.0 ) { # going straight into a wet period
					status="wet"; LLw[i] = 1; newwet=1; firstdry=0
					x3[i] = x1[i]
					pullx1[(i-nump):i] <- 1
				} else if (x2[i] <= -1.0) { 
					stop("Something's Wrong: Ended dry period with start of a new one!")
				} else { # going to a "normal" stateof being
					status="normal"; LLn[i] = 1
					x3[i] = 0
					pullnorm[(i-nump):i] <- 1
				}
				nump=0
			}else { # We're getting closer to ending our wet period, but it's not over yet
				nump=nump+1
			}			
		} # End selecting normal/wet/dry
	} # End time loop
	# ------------------------------------------

	# ------------------------------------------
	# 4. Format & return outputs
	# ------------------------------------------
	# From Ben Cook: 
	#    column vectors of monthly data x1,x2, and x3 have now been calculated.  
	#    Now must substitute values from series x1 and x2 for those in x3 when 
	#    appropriate
	
	
	x=x3 # Initialize storage vector of final pdsi
	# Putting in backgroudn wet & dry
	if(any(pullx1==1)) x[which(pullx1==1)] <- x1[which(pullx1==1)]
	if(any(pullx2==1)) x[which(pullx2==1)] <- x2[which(pullx2==1)]
	
	# Adding in the nomral periods; substitute into x whichever of x1 & x2 is stronger
	Lwet = which(pullnorm==1 & abs(x1) >= abs(x2))
	Ldry = which(pullnorm==1 & abs(x2) >  abs(x1))
	if(length(Lwet)>0) x[Lwet] <- x1[Lwet]
	if(length(Ldry)>0) x[Ldry] <- x2[Ldry]
	
	# From Ben Cook:
	# Modified Palmer Index.
	# The series xm is a modified index that may differ from x when the
	# probability of ending a drought or wet period is greater than zero but
	# less than 100%.  For such months, the modified index is a weighted average
	# of x3 and either x1 or x2, depending on whether the existing period as
	# indicated by x3 is a drought (x3<0) or a wet period (x3>0).  The weights are
	# the probabilities (Pe) of ending the drought or wet period. If this probability
	# is  prob=100 * Pe, the weighted value for dry period is:
	# 
	#  x = prob*x2 + (1-prob)* x3;
	#
	# For a wet period, it is:
	#
	#  x = prob*x1 + (1-prob)*x3

	prob=Pe/100 # The probability of ending a drought or wet period
	L7=prob>0 & prob<1 # True/False
	L7d = which(L7 & x3<0) # True/False
	L7w = which(L7 & x3>0) # True/False
	
	xm=x
	if(length(L7d)>0) xm[L7d] <- prob[L7d] * x1[L7d] + (1-prob[L7d])*x3[L7d]
	if(length(L7w)>0) xm[L7w] <- prob[L7w] * x2[L7w] + (1-prob[L7w])*x3[L7w]


	# Return outputs; x4 calculated in loops as we go & no formatting needed
	return(list(x=x, xm=xm, x4=x4))
	# ------------------------------------------


} # End function