function [LON,LAT,Z] = read_in_etopo_bathy(f)
%read_in_etopo_bathy Summary of this function goes here
%   Detailed explanation goes here

% ncdisp(f);

% read in the lon and lat info from the nc file
lon=ncread(f,'lon');
lat=ncread(f,'lat');

[LON,LAT]=meshgrid(lon,lat);

LON=LON';
LAT=LAT';

Z=ncread(f,'Band1');
end

