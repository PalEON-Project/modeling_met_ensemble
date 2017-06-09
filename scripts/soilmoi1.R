# Calculate end of month soil moisture from soil moisture accounting
# Translated from code written by Ben Cook; which look slike it was translated by Dave Meko in 1997
# Author: Christy Rollinson (crollinson@mortonarb.org); translated from B. Cook matlab code
#
#
# Notes from Ben Cook:
# Source: Palmer, Wayne C., 1965.  Meteorological Drought; Research Paper No. 45.
# US Dept of Commersce, Washington, DC.  Coded from equations and method described
# on pages 6-11.
#
# soilmoi1.m carries out the actual monthly accounting of soil moisture and other
# variables.  The code was intended for inches & months, but since it's just accounting,
# units just need to be internally consistent with AWC matching the unit of P & PE
#
# Notes: I think units for p, pe, and awc need to match (mm/time), but actual units
#        don't matter (in other words in/mo should work the same)
# Inputs:
#   1. p    = monthly precipitation in (inches), jan - dec; length = 12*nyrs
#   2. pe   = potential evapotranspiration in (inches); length = 12*nyrs
#   3. awcs = available water capacity (inches) in surface (sfc) layer
#   4. awcu = available water capacity (inches) in underlying layer
#   5. ssgo = starting soil moisture (inches) in sfc layer
#   6. sugo = starting soil moisture (inches) in under layer
# Outputs (meat, class=list)
#   1. dels  = soil moisture change in sfc layer
#   2. delu  = soil moisture change in under layer
#   3. del   = total soil moisture change in both layers
#   4. ss1   = starting soil moisture, sfc
#   5. su1   = starting soil moisture, under
#   6. s1    = starting soil moisture, combined
#   7. ss2   = ending soil moisture, sfc
#   8. su2   = ending soil moisture, under
#   9. s2    = ending soil moisture, combined
#  10. smean = mean soil moisture, combined (s1 + s2)/2
#  11. r     = recharge, combined
#  12. pr    = potential recharge
#  13. ro    = runoff
#  14. pro   = potential runoff
#  15. loss  = loss
#  16. ploss = potential loss
#  17. et    = estimated actual evapotranspiration
#
# Workflow: 
# 1. Check Inputs, Create placeholders for outputs
# 2. Caluclate water balance
# 3. Initialize soil moisture
# 4. Loop through time to calculate soil moisture
#    4.1 Surface Layer Dynamics
#    4.2 Underlying Layer Dynamics
# 5. Calculate some additional variables for total soil
# 6. Calculate recharge, loss, and runoff
# 7. Format & return output
#

calc.soilmoist <- function(p, pe, awcs, awcu, ssgo, sugo) {
	# ------------------------------------------
	# 1. Check Inputs, Create placeholders for outputs
	# ------------------------------------------
	ntime = length(p)
	
	# Insert formatting checks when I'm not being lazy!
	
	# Creating place holders for variables
	# Soil Moisture
	ss1  <- rep(NaN, ntime)
	ss2  <- rep(NaN, ntime)
	su1  <- rep(NaN, ntime)
	su2  <- rep(NaN, ntime)

	# Recharge
	rs   <- rep(NaN, ntime)
	ru   <- rep(NaN, ntime)
	
	# Runoff
	ro   <- rep(NaN, ntime)

	# Net loss to Evapotranspiration
	es   <- rep(NaN, ntime)
	eu   <- rep(NaN, ntime)

	# Change in soil moisture
	dels <- rep(NaN, ntime)
	delu <- rep(NaN, ntime)	
	# ------------------------------------------

	# ------------------------------------------
	# 2. Caluclate water balance
	# ------------------------------------------
	d   <- pe - p # deficit = excess of potential evapotraspiration over precipitation
	awc <- awcu + awcs # Combined water capacity
	# ------------------------------------------

	# ------------------------------------------
	# 3. Initialize soil moisture
	# ------------------------------------------
	# Start things with the values provided
	ss1this = ssgo
	su1this = sugo
	# ------------------------------------------

	# ------------------------------------------
	# 4. Loop through time to calculate soil moisture
	#    4.1 Surface Layer Dynamics
	#    4.2 Underlying Layer Dynamics
	#    4.3 Calculate some additional variables for total soil
	#    4.4 Calculate recharge, loss, and runoff
	# ------------------------------------------
	for(i in 1:ntime){
		dthis = d[i] # pe-p for right now
		
		# -------------------------
		# 4.1 Surface Layer Dynamics
		# -------------------------
		sempty = awcs - ss1this # how much the sfc layer could take in
		
		if(dthis >= 0){ # if pe exceeds precip, we're going to lose soil moisture
		  dels[i] <- -dthis # tentatively set the soil moisture to pe-pe
			if(dthis > ss1this) { # if pe - p exceeds what we have in the sfc layer, get rid of what we have
				dels[i] <- -ss1this
			}
			
			rs[i] = 0 # No excess, so no recharge
			ro[i] = 0 # no excess, so no runoff
			
			# Net Loss from sfc layer
			if(dels[i] < 0){
				es[i] <- -dels[i]
			} else {
				es[i] <- 0
			}

			excess = 0
		} else { # ppt exceeds pe, so our soils will get wetter (or stay at capacity)
		  dels[i]   <- min(sempty, -dthis) # either all the precip, or as much as the soils can take in
			rs[i]     <- dels[i] # surface recharge
			excess <- -dthis - dels[i] #
			es[i]     <- 0			
		} # End surface balance ifelse
		
		ss1[i]  <- ss1this # save our starting point
		ss2[i]  <- ss1this + dels[i] #
		ss1this <- ss2[i] # Next starting point will be our current end
		# -------------------------

		# -------------------------
		# 4.2 Underlying Layer Dynamics
		# -------------------------
		uempty <- awcu - su1this # how much the under layer could take in
		
		if(excess<=0){ # no moisture input from above
		  eu[i] <- (dthis - es[i]) * (su1this/awc) # "loss" from the under layer
			eu[i] <- min(eu[i], su1this)
			
			if(eu[i] < 0) eu[i] = 0 # no negative values allowed
			ru[i] = 0 # no recharge
			ro[i] = 0 # no runoff
			delu[i] = -eu[i] # change in under soil moisture			
		} else { # There is some moisture input from above
		  eu[i] = 0 # no loss from underlying layer
			delu[i] = min(uempty, excess) # change is how much it could take or how much there is
			ru[i] = delu[i] # setting the recharge
			if(excess > uempty) { # We have more than the soil can take --> runoff!
				ro[i] <- excess - uempty
			} else { # no runoff because we can take it all
				ro[i] <- 0
			}
		}
		
		su1[i] <- su1this # Save our starting point
		su2[i] <- su1this + delu[i] # save our ending point
		su1this <- su2[i] # This ending point is the next time step's starting point
		# -------------------------
		if(is.na(ss1[i]) | is.na(su1[i])) stop("su1 is na")
	}
	# ------------------------------------------

	# -------------------------
	# 5. Calculate some additional variables for the soil
	# -------------------------
	del   <- delu + dels # total change in soil moisture, combining layers
	et    <- p - ro - del # evapotranspiration
	r     <- rs + ru # total recharge, combined layers
	loss  <- es + eu # total losses
	s1    <- ss1 + su1 # total starting soil moisture
	s2    <- ss2 + su2 # total ending soil moisture
	smean <- (s1 + s2)/2 # mean starting and ending soil moisture
	# -------------------------

	# -------------------------
	# 6. Calculate recharge, loss, and runoff
	# -------------------------
	pr <- awc - s1 # Potential recharge
	
	# Potential losses
	dope <- c(pe, ss1)
	plosss <- min(dope)
	plossu <- (pe - plosss) * (su1/awc)
	ploss  <- plosss + plossu
	
	pro <- awc - pr # potential runoff 
	# -------------------------

	# ------------------------------------------
	# 7. Format & return output
	# ------------------------------------------
	meat <- list()
	meat[["dels" ]] <- dels
	meat[["delu" ]] <- delu
	meat[["del"  ]] <- del
	
	meat[["ss1"  ]] <- ss1
	meat[["su1"  ]] <- su1
	meat[["s1"   ]] <- s1
	
	meat[["ss2"  ]] <- ss2
	meat[["su2"  ]] <- su2
	meat[["s2"   ]] <- s2
	meat[["smean"]] <- smean
	
	meat[["r"    ]] <- r
	meat[["pr"   ]] <- pr
	
	meat[["ro"   ]] <- ro
	meat[["pro"  ]] <- pro
	
	meat[["loss" ]] <- loss
	meat[["ploss"]] <- ploss
	
	meat[["et"   ]] <- et
	
	return(meat)
	# ------------------------------------------

}