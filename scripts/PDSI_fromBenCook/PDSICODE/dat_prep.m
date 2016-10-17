% dat_prep
%
% Example of a script that runs the functions to compute PDSI and related
% variables.   This example is for gridded monthly precipitation and
% temperature in Tunisia, at 1 11E, 33 37N
%
% This input ascii T and P happens to be in mm/day.  The
% script converts T to degrees F and P
% to inches/month, which are the units required by the functions.  Modify
% the sections of this script to deal with your data. 
%
%
% This script puts the output into an xls file, with sheets corresponding
% to PDSI, PDMI, PDHI.  The xls is pre-formatted so that the data show 2
% digits to right of decimal point, which is reasonable for these
% variables.   For convenience, if using new P and T data and writing
% another xls file, just blank out the data in the existing file (e.g.,
% palmer.xls) and re-save the file with the desired xls filename.  That
% wayt the sheet names and formatting is preserved when your newly
% generattred pdsi is written to the xls. 
%
% Comments preced each section 


%-- CLEAR THE WORKSPACE
clear all;
close all;
clc;

pf3 = 'daylennh.mat';  % a table with the daylength factor;  you need access to this hard-coded table. Do not Change.

%--- HARD CODE SOME FILE NAMES 
pf4 = 'palmer.xls'; % an excel spreadsheet to hold the output
pf1 ='icru2_pre_eh_1-11E_33-37N_n.dat'; % monthly P in mm/day
pf2='icru2_tmp_eh_1-11E_33-37N_n.dat'; % Monthly T deg C
lat = 33+37/60; % latitude in decimal degrees
site_name='lon_1_11_lat_33_37'; % name used for output labeling only
% Set input required for labeling only



%--- Get daylength factor
dayz = load(pf3);
dayz=dayz.dayz;

%--- BECAUSE INPUT PRECIPITATION IS NOT MONTHLY TOTAL, WILL NEED NUMBER OF
% DAYS IN YEAR TO CONVERT FROM MM/DAY TO INCHES/HR
k1 =[31 28 31 30 31 30 31 31 30 31 30 31]; % days in month, nonleap
k2 = [31 29 31 30 31 30 31 31 30 31 30 31]; % days in month, leap


%--- LOAD P AND CONVERT TO INCHES PER MONTH;  WANT 13-COL MATRIX P WITH
%YEAR AS  COL 1 AND INCHES OF PRECIP FOR JAN, FEB, ... DEC AS THE OTHER 12
%COLS

eval(['P = load(''' pf1 ''',''-ascii'');']); % year and 12 mos precip, as mm/day
yr = P(:,1);
L = leapyr(yr);  % 1 if leap year; % need to identify leap years for conversion of P in mm/day to P in inches/mnonth

P(:,1)=[]; % strip year off P
[mP,nP]=size(P);
if nP~=12;
    error('Not 12 mos in P');
end
F = repmat(k1,mP,1); % dupe the vector of day in nonleap mos
P = 25.4*(F .* P); % precip in inches
Pfeb = P(:,2);
Pleap = Pfeb(L) * 29/28;
Pold =P;
P(L,2)=Pleap;
P=[yr P];
Pold=[yr Pold];


%--- BEGIN STORING INPUT DATA IN CELL VARIABLE datmon

datmon{1}=P; % store precip, inches per month



%--- LOAD T AND CONVERT TO DEGREES F.  WANT 13-COL MATRIX WITH YR AS COL 1,
% JAN-DEC MEAN MONTHLY T IN DEG F AS COLS 2-13

eval(['T = load(''' pf2 ''',''-ascii'');']); % year and 12 mos T, deg C
yr = T(:,1);

T(:,1)=[]; % strip year off T
[mT,nT]=size(T);
if nT~=12;
    error('Not 12 mos in T');
end
T=32 + (9/5)*T;
T=[yr T];
datmon{2}=T; % store T, deg F


%--- Check that identical years in P and T
yrP = P(:,1);
yrT = T(:,1);
if ~all(yrP==yrT);
    error('P and T not identical years');
end

yrgo = yrP(1);
yrsp = yrP(end);




%--- LOAD DATA INTO CELL datother

datother{1}=lat;
datother{2}=[1 5];  % soil moisture capacity (in) in upper and lower levels
datother{3}=[yrgo yrsp; yrgo yrsp];
datother{4}=dayz;


%--- HANDLING OF SNOWMELT NOT YET CODED (THIS WOULD BE AN EXTENSION OF THE
% PDSI);  SET INPUT DATA THAT DEALS WITH SNOWMELT AND OTHER FACTORS NOT YET
% CODED
snowinf=[];
kopt=[1 1 1];
penopts=[];
datpen=[];



%--- CALL FUNCTION pdsi1 TO COMPUTE THE PDSI AND RELATED OUTPUTS

textin{1}=site_name; % load site label into cell 

datout=pdsi1(datmon,datother,snowinf,textin,kopt,penopts,datpen);


%--- WRITE DESIRED OUTPUT TO AND EXCEL SPREADSHEET

xlswrite(pf4,datout{2},'Palmer_Index','A1');
xlswrite(pf4,datout{3},'Modified Palmer Index','A1');
xlswrite(pf4,datout{10},'hydrologic index','A1'); %---modified 20080312

