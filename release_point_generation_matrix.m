function [wp3]=release_point_generation_matrix(wp,resolution,range)

%  Usage: [wp3]=release_point_generation_matrix(wp,resolution,range)
%
% This function creates a matrix of latitude and longitude points for 
% input into a trajectory model.
%
% Inputs
% ------
% wp = array of [latitude longitude] that is the origin of the matrix, it is
%       also the south west corner of the matrix, e.g. wp1=[17+45/60 -68-00/60];
% resolutiuon = distance between each release point (km)
% range = how far from the release point to generate drifters (km)
%
% Outputs
% -------
% wp3 = matrix of release points for input into the drifter model

wp2=wp;

count=round(range/resolution);

%% course2ll is a function john kerfoot wrote to give a new lat and lon 
%% based on existing lat and lon, course in degrees and a distance in km
%% this loop creates the drifter points of constant latitude
for ii=2:count
    wp(ii,:)=course2ll([wp(ii-1,1) wp(ii-1,2) 90 resolution]);
end

%%  this loop creates the drifetrs with constant longitude

for ii=2:count
    wp2(ii,:)=course2ll([wp2(ii-1,1) wp2(ii-1,2) 0 resolution]);
end

[X,Y]=meshgrid(wp(:,2),wp2(:,1));

wp3=[Y(:) X(:)];

%[X,Y]=meshgrid(unique(wp3(:,2)),unique(wp3(:,1)));

end

