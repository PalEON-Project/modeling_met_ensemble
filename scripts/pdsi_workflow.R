# PDSI workflow following Ben Cook

# 1. Set UP
#    A. Set up File paths
#    B. Read in Variables
#       1. Lat, Lon
#       2. Extract Soil Info: Depth, layer/moisture (soilm_layer)
#       3. Radiation: LWdown, LWup, SWdown, SWup
#            - Rnet = (LWdown + SWdown) - (LWup + SWup)
#       4. Vapor Pressure: press, humidity (qair)
#           - Vapor Pressure (e_all) = press*qair / 0.6213 # In Pascals
#       5. Temp & Precip: tair, precipf
#    C. Set Standardization interval; NOAA calibrates 1931-1990
# 2. Unit Conversions
#    A. Tair = C (from K)
#    B. precip (mm/day)
#    C. Penman Monteith: ea (VP, Pa), press (Pa), Rnet (W/m2)
# 3. Define some variables
#     k1 = days per month, non-leap year
#     k2 = days per month, leap year
#     P  = precip
#     L  = leapyear (yes/no); use to conver P to in/mo (form mm/day)
#     F  = repate days per month
#     P  = convert precip to in per mo = 0.3937*F*P
#     More leap year conversions
#     T  = convert temperature to Farenheit  = 32+9/5*T
# 4. Run PDSI1

# datmon (List)
#  1. Precip (in/mo)
#  2. Temperature (F)

# datother (List)
#  1. Latitude
#  2. Soil Moisture Capacity Depths: 1" & 5" (Used for bucket model for soil moist); length=2
#     - awcs 
#     - awcu
#  3. year range for PDSI & standard: yrgo yrsp, stdp1 stdp2; dim=c(2,2)
#      - yrgo  = 
#      - yrsp  = 
#      - stdp1 = 1931
#      - stdp2 = 1990
#  4. Daylength factor (dayz); from black-box matlab file
#      - percentage of possible sunshine

# datout (List from pdsi1)
#  1. Z = Z index
#  2. X = PDSI (Thornthwaite)
#  3. XM = Modified PDSI
#  4. W = Monthly soil moisture (in)
#  5. RO = model runoff
#  6. S1 = effective precip; max(0, p-r*pe)
#      precip minus PE
#  7. Precipitation
#  8. Temperature
#  9.
# 10.
# 11. Thorthtwaite PE

# Options
#  - snowinf
#  - kopt
#  - penopts
#  - datpen