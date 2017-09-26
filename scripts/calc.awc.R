# Function to Calculate AWC from soil texture 
# Pulled from ED2 Code, ed_params.f90; lines 4521-4553
#    According to documentation: this is from Cosby et al 1984, Table 4, equation 1
#    NML (whoever that is) said it should be saturated moisture potential over moisture potential
#
# Note: in a couple test runs, these give us similar, but not identical values to Web Soil Survey
#.      looks like it under-shoots awc in highly sandy soils

calc.awc <- function(sand, clay){
	# sand = faction sand (0-1)
	# clay = fraction clay (0-1)
	# fc = field capacity (m^3/m^3)
	# pwp = permanent wilting point at -1.5 MPa (m^3/m^3)
	# awc = available water capacity = field capacity - permanent wilting point 
  #       units: m3/m3
	#       https://en.wikipedia.org/wiki/Available_water_capacity
	
	# Calculating percentage silt
	silt = 1 - sand - clay

	# Constants; mix of from consts_cos.F90 & ed_params.f90
	fieldcp_k = 0.1 # hydraulic conductivity at field capacity (mm/day)
	wdns = 1.000e3 # liquid water density (kg/m3)
	day_sec = 86400 # num of seconds in day (s/day)
	hr_sec = 3600 # num seconds in an hour (s/hr)
	soilwp_MPa = -1.5 # matric potential at wilting point (MPa)
	grav = 9.80665 # Gravity acceleration (m/s3)

	# Calculated things we need to know
	slbs = 3.10 + 15.7*clay - 0.3*sand # B exponent; unitless
	slpots = -1. * (10^(2.17 - 0.63*clay - 1.58*sand)) * 0.01 # Soil moisture potential at saturation (m)
	slmsts = (50.5 - 14.2*sand - 3.7*clay) / 100. # Soil moisture at saturation (m^3/m^3)
	slcons = 10^(-0.60 + 1.26*sand - 0.64*clay) * 0.0254/hr_sec # Hydraulic conductivity at saturation (m/2)
	
	# Calculating Field Capacity
	fc = slmsts * (fieldcp_k/wdns/day_sec)/slcons^(1/(2*slbs + 3))
	pwp = slmsts * (slpots / (soilwp_MPa * wdns / grav))^(1/slbs)

	# Cacluate awc
	awc <- abs(fc - pwp)
	
	return(awc)
}
