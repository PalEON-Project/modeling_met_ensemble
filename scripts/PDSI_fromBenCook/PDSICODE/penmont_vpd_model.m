function [PET,VPD,RH] = penmont_vpd_model(Ta,Rnet,ea,press,u)

% [PET,VPD,RH] = penmont_vpd_model(Ta,Rnet,ea,press,u)
%
% Modified to take data from a model
%
% This function will use the Penman-Monteith method to calculate reference
% (potential) evapotranspiration.
% 
% Inputs:
%       Ta      = temperature, degrees C
%       Rnet    = surface net radiation, W/m2
%       ea      = actual vapor pressure, Pascals
%       press   = surface pressure, Pascals
%       u       = wind speed, m/s
%
% Outputs:
%       PET  = potential evapotranspiration (mm/day)       
%       VPD  = vapor presssure deficit (kPa)       
%       RH   = relative humidity (fraction)
%
% Written by Benjamin I. Cook
% 
% Based on:
%       Xu, C-Y., and V. P. Singh. "Cross comparison of empirical equations 
%               for calculating potential evapotranspiration with data 
%               from Switzerland." Water Resources Management,
%               16.3 (2002): 197-219.
%
%           FAO Document:
%           http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents
%
% For Tetens, both above and below zero:
%       http://cires.colorado.edu/~voemel/vp.html

% ground heat flux, set to zero (works fine on >monthly timescales, and is 
% accurate if one calculates Rnet as SH+LH)
gflux=0;    

% Calculate the latent heat of vaporization (MJ kg-1)
lambda_lv=2.501-(2.361e-3).*Ta;

% Convert Pressure to kPa
press=press./1000;

% Convert actual vapor pressure to kPa
ea=ea./1000;

% Calculate Saturation vapor pressure (kPa)
es=0.611.*exp((17.27.*Ta)./(Ta+237.3));   

% Slope of the vapor pressure curve (kPa C-1)
delta_vpc=(4098.*es)./((Ta+237.3).^2);

% Psychometric constant (kPa C-1)
psych_const=0.00163.*(press./lambda_lv);

% Net radiation
Rnet=Rnet./11.6;  % convert W/m2 to MJ/m2/d

% Calculate VPD
VPD=es-ea;

% Calculate Relative Humidity
RH = ea./es;

% Potential Evapotranspiration (mm/day)
PET=(0.408.*delta_vpc.*(Rnet-gflux)+psych_const.*(900./(Ta+273)).*u.*(es-ea))./(delta_vpc+psych_const.*(1+0.34.*u));









