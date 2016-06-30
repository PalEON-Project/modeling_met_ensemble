# Function to downscale data

# what we need:
# 1. Training Dataset -- the dataset used to for downscaling
# 2. Downscale Dataset -- the dataset that needs to be spatially & temporally downscaled

# 1. Determine if we have an overlapping period that should be used for the adjustment
# 1.a. If yes, use all overlapping days/years
# 1.b. If no, use the closest 30 years and assume they have similar means

# 2. Format the downscale data to have the temporal resolution of the training data
#    2.1. Center data points for the window they apply to
#    2.3. Add NAs around it?

# 3. Create the downscaling model
#    Variable notation
#           data=response dataset
#           y.down  = the observed value
#           y.train = mean observed value in training dataset for the spatial & temporal resolution we want in th eend
#           time = time stamp (hour) of training data
#           temp.down = timestep of the dataset to be downscaled (time-specific bias correction factor)
#    3.1. Basic bias-correction factor applied daily over the coarse of the year -- corrects seasonal cycles
#           y.down ~ s(doy) + y.train # Gets the basic bias-correction factor
#    3.2. temporal cycle from training data
#           (y.train - predict3.1) ~ s(time, by=day) + factor(day) # get the dirunal cycle for each day
#            ** predict.3.1 = predicted daily mean
#    3.2. Calculate timestep anomalies
#           y.down.mean = predict3.1 + predict3.2 # Daily bias correction + diurnal cycle
#           
#           y.anom[temp.res] = y.down.mean - y.down[temp.res]  # sub-daily temperatures now expressed as daily anomalies?
#           y.down(?) ~ s(time, by=temp.down) + y.train(?) + factor(temp.down)
#            ** everything on the right-hand side of the equation is from the trainign dataset
# 
#    
#    3.2. Calculate the timestep anomaly for the overlap period
# 
#    3.3. Predict the mean diurnal cycle over the full 
#    Basic idea: y.mean ~ s(time, by=day.of.year) + resid.day.year +