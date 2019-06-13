function [TUV]=convert_AmSeas_NC_to_TUVstruct(file,dtime,conf)
%% function written on April 14, 2015 by Hugh Roarty to convert a HYCOM
%% netCDF file to TUV structured array

%% read in the nc variables
tot.t=ncread(file, 'time');
tot.lat = ncread(file, 'lat');
tot.lon = ncread(file, 'lon');
I = find(tot.lon >= conf.Totals.grid_lims(1) & tot.lon <= conf.Totals.grid_lims(2)); % lon in degrees east
J = find(tot.lat >= conf.Totals.grid_lims(3) & tot.lat <= conf.Totals.grid_lims(4)); % lat
lons = tot.lon(I); 
lons = lons-360;
lats = tot.lat(J);
[LON,LAT] = meshgrid(lons, lats);
LON=LON';
LAT=LAT';

%% get the time units from the nc variable
% for example:
% time.units: 'hours since 2010-01-25T22:00:00Z'
% time.str: '2010-01-25T22:00:00'
% time.matlab: 7.3416e+05
 
time.units=ncreadatt(file,'time','units');
ind=strfind(time.units,'since');
time.str=time.units(ind+6:ind+28);
time.matlab=datenum(time.str,'yyyy-mm-dd HH:MM:SS');

%% calculate the time index so you can pull the correct data 
%% convert the hours in the nc file to matlab time, first divide by 24 to 
%% convert hours to days and then add the reference time
t0=tot.t/24+time.matlab;

%% calculate the difference between the matlab time of the nc files (t0)
%% and the matlab time of the analysis (dtime)
tDiff = abs(t0 - dtime(1));
tDiff2= abs(t0 - dtime(end));

%% find the minimum of the difference and that's the time_index that you need
[~,tlen] = min(tDiff);
time_index=tlen;%example of last day;cd
[~,tlen2] = min(tDiff2);
time_index2=tlen2;%example of last day;cd

ti=time_index2-time_index+1;

%% pull out the u and v from the nc files
tot.u=ncread(file, 'water_u', [I(1) J(1) 1 time_index], [I(end)-I(1)+1 J(end)-J(1)+1 1 ti]);
tot.v=ncread(file, 'water_v', [I(1) J(1) 1 time_index], [I(end)-I(1)+1 J(end)-J(1)+1 1 ti]);

%% convert the m/s of the u and v to cm/s
tot.u = double(tot.u)*100;
tot.v = double(tot.v)*100;


%% create the empty TUV stuctured array
TUV = TUVstruct( [length(LON(:)) 1], 1 );

%% get the size of the returned data so you can reshape it into the TUV structure
d=size(tot.u);

%% Populate the TUV struct with the data from the .nc file

TUV.TimeStamp=t0(time_index:time_index2);
TUV.LonLat(:,1)=LON(:);
TUV.LonLat(:,2)=LAT(:);

if numel(d)>2
    TUV.U=reshape(tot.u,[d(1)*d(2) d(4)]);
    TUV.V=reshape(tot.v,[d(1)*d(2) d(4)]);
else
    TUV.U=reshape(tot.u,[d(1)*d(2) 1]);
    TUV.V=reshape(tot.v,[d(1)*d(2) 1]);
end

TUV.LON=LON;
TUV.LAT=LAT;
TUV.Ug=tot.u;
TUV.Vg=tot.v;

% TUV.ErrorEstimates.Type='OIuncert';
% TUV.ErrorEstimates.Uerr=tot.u_err(:);
% TUV.ErrorEstimates.Verr=tot.v_err(:);
% TUV.ErrorEstimates.UerrUnits='cm2/s2';
% TUV.ErrorEstimates.VerrUnits='cm2/s2';
end
