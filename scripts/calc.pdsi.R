# % BEN'S CODE
# %  Pull out model diagnostics data for each grid cell, one at a 
# %  time.
# %       tsurf should be in units of deg C
# %		prec should be in units of mm/day
# % Make these unit corrections, if necessary.
# What we need
# 1. tair  - deg C (k-273)
# 2. precipf - mm/day
# 3. Rnet = (lwdown + swdown) - (lwup + swup) ## THIS IS A PROBLEM -- CAN WE GET SOME HAND-WAVY NUMBER FOR UP??

temp <- 

tmp_pdsi=tas_all(:,:,n_lat,n_lon)-273;
pre_pdsi=pr_all(:,:,n_lat,n_lon);

% Other variables needed for Penman-Monteith ET 
%    ea and press should be in pascals, Rnet in W/m2
ea_pdsi=e_all(:,:,n_lat,n_lon);
press_pdsi=press_all(:,:,n_lat,n_lon);      
Rnet_pdsi=Rnet_all(:,:,n_lat,n_lon);      

% Check to see if the current gridcell contains any land. If it 
% does, then perform the PDSI calculation. 
if isnan(landfrac(n_lat,n_lon))==0;

% This is a table with daylength factors. The PDSI calculation
% needs access to this hard-coded table, so do not change.
pf3 = 'daylennh.mat';

% OTHER HARD CODED FILE NAMES
pf4       = 'palmer.xls';         % an excel spreadsheet to hold the output
lat_curr  = abs(lat(n_lat));      % latitude in decimal degrees
site_name = 'lon_1_11_lat_33_37'; % name used for output labeling only

% Get daylength factor
dayz = load(pf3);
dayz=dayz.dayz;

% BECAUSE INPUT PRECIPITATION IS NOT MONTHLY TOTAL, WILL NEED NUMBER OF
% DAYS IN YEAR TO CONVERT FROM MM/DAY TO INCHES/HR
k1 =[31 28 31 30 31 30 31 31 30 31 30 31];   % days in month, nonleap
k2 = [31 29 31 30 31 30 31 31 30 31 30 31];  % days in month, leap

% P is going to be the ultimate precipitation input matrix. It
% will ultimately be in the form of a 13-column MATRIX (P) with
% the year vector as Column 1 and inches of precipitation for
% Jan, Feb, ... Dec as the other 12 columns.
P=pre_pdsi;

% My original year vector (yrs) is column oriented; change it
% to a single column
yr = yrs';

% Sets L=1 if leap year; this is needed to identify leap years 
% for conversion of P in mm/day to P in inches/month
L = leapyr(yr);  

% Check to make sure I have 12 months of precipitation
[mP,nP]=size(P);
if nP~=12;
error('Not 12 mos in P');
end

% Dupe days per month vector for nonleap months
F = repmat(k1,mP,1); 

% Convert precipitation to inches per month.
P = 0.03937.*(F .* P);

% Make corrections for a Leap year, if necessary
Pfeb = P(:,2);
Pleap = Pfeb(L) * 29/28;
Pold =P;
P(L,2)=Pleap;
P=[yr P];                
Pold=[yr Pold];

% Store precipitation (inches/month) in input cell variable "datmon"
datmon{1} = P;

% Temperature input data matrix will have same form as
% precipitation, with year in first column and each month of
% temperature after
T=tmp_pdsi;

% Check to make sure I have 12 months of temperature
[mT,nT]=size(T); 
if nT~=12;
error('Not 12 mos in T');
end

% Convert to Temperature from Celsius to Fahrenheit
T=32 + (9/5)*T;

% Add on the Year Vector
T=[yr T];

% Store temperature (F) in "datmon", the input cell variable
datmon{2}=T;

% Check to make sure T and P have identical years
yrP = P(:,1);
yrT = T(:,1);
if ~all(yrP==yrT);
error('P and T not identical years');
end

% Now, set the beginning and end years to calculate PDSI over.
% Here, I'm just setting it to all years.
yrgo = yrP(1);
yrsp = yrP(end);

% Now, load other input data into a different cell array,
% "datother"

% Latitude
datother{1}=lat_curr;

% Soil moisture capacities (inches) for upper and lower levels in the
% PDSI bucket model. 1 inch and 5 inches are the standard.
datother{2}=[1 5];

% Vector containing year range to calculate PDSI and year range
% for the standardization step.
datother{3}=[yrgo yrsp; stdp1 stdp2];

% Daylength factor
datother{4}=dayz;

% Other options. These are the default, because alternative
% options have not been coded or implemented yet.
snowinf=[];
kopt=[1 1 1];
penopts=[];
datpen=[];

% Site Label.
textin{1}=site_name; % load site label into cell     

%% Thornthwaite PDSI

% This uses the default version of Meko's PDSI code (pdsi1). I
% did modify it to include the Thornthwaite PET in the output.

% Calculate the PDSI, and pull out the Thornthwaite PET and
% PDSI from the output cell arrays.
datout=pdsi1(datmon,datother,snowinf,textin,kopt,penopts,datpen);
pe_thorn=datout{11}; % TH PE
pdsi_th=datout{2};

% Reformat the PDSI output and save to my main output variable.
PDSI_th(1:length(yrs),1:12,n_lat,n_lon) = pdsi_th(:,2:13);

%% Penman-Monteith PDSI

% This uses my own modifications and functions to calculate the
% Penman-Monteith version of the PDSI.

% First, calculate PET (mm/day) using my own Penman-Monteith
% functin, based on the FAO formulation:            
%       [PET,VPD,RH] = penmont_vpd_model(Ta,Rnet,ea,press,u)
%
% For these applications, I generally just set windspeed to 1
% m/s, to simplify the data requirements.            
[PET_pm,VPD,RH] = penmont_vpd_model(tmp_pdsi,Rnet_pdsi,ea_pdsi,press_pdsi,1);

% Convert PET to inches/month from mm/day
PE_new_pm = 0.03937.*(F .* PET_pm); % Convert PET to inches

% Format Penman-Monteith for the PDSI Calculation
PE_new_pm=[yr PE_new_pm];

% Calculate Penman-Monteith PDSI. This uses "pdsi1_netrad", a
% modification of pdsi1 I made that allows the user to define
% their own PET. So really, you could sub in anything you want.
datout=pdsi1_netrad(datmon,datother,snowinf,textin,kopt,penopts,datpen,PE_new_pm);
pdsi_pm=datout{2};

% Reformat the PDSI output and save to my main output variable.
PDSI_pm(1:length(yrs),1:12,n_lat,n_lon) = pdsi_pm(:,2:13);

else

% For Cells with NO land, make their PDSI values "NaN"
PDSI_pm(1:length(yrs),1:12,n_lat,n_lon) = nan;
PDSI_th(1:length(yrs),1:12,n_lat,n_lon) = nan;

%soilm_total(1:length(yrs),1:12,n_lat,n_lon) = nan;
%soilm_layer(:,1:length(yrs),1:12,n_lat,n_lon) = nan;

end

