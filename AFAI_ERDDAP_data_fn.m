function [ERD] = AFAI_ERDDAP_data_fn(f,dtime,lims)
%AFAI_ERDDAP_data_fn extracts AFAI data from NOAA erddap server
%   Usage [ERD] = AFAI_ERDDAP_data_fn(f,dtime,lims)
%   f       the url of the erddap server
%   dtime   matlab time for when you want to extract the single timestamp
%   lims    [lon_min lon_max lat_min lat_max]
%
%   Output
%   structured array with fields
%   LON     matrix of longitudes
%   LAT     matrix of latitudes
%   AFAI    matrix of AFAI values

% read inthe variables
AFAI_time=ncread(f,'time');
lon=double(ncread(f,'longitude'));
lat=double(ncread(f,'latitude'));

% to convert the time to matlab time
AFAI_matlab_time = datenum(1970,1,1,0,0,AFAI_time);

%% calculate the difference between the matlab time of the nc files (t0)
%% and the matlab time of the analysis (dtime)
tDiff = abs(AFAI_matlab_time - dtime(1));

%% find the minimum of the difference and that's the time_index that you need
[~,tlen] = min(tDiff);
ti=tlen+1;%example of last day;cd

I = find(lon >= lims(1) & lon <= lims(2)); % lon in degrees east
J = find(lat >= lims(3) & lat <= lims(4)); % lat

start=[I(1) J(1) ti];
count=[I(end)-I(1)+1 J(end)-J(1)+1 1];
% stride
ERD.AFAI=ncread(f,'AFAI',start,count);

lon2=lon(I);
lat2=lat(J);
[LON,LAT] = meshgrid(lon2, lat2);
ERD.LON=LON';
ERD.LAT=LAT';

ERD.dtime=AFAI_matlab_time(ti);


end

