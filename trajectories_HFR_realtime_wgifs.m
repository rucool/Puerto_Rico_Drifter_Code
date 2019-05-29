%% Script to View particle trajectories, as measured by HF radar, of Sargassum seaweed near Puerto Rico.
% Written by Joseph Anarumo on 3/19/2019

%% Clear workspace

close all
clear all
clc

compType=computer;

if ~isempty(strmatch('PCWIN64',compType))
    root='L:';
else
    root='/home';
end

%% Add directories to MATLAB path for following functions.

addpath([root '/codaradm/HFR_Progs-2_1_3beta/matlab/general']);
%addpath([root '/codaradm/operational_scripts/totals']);
add_subdirectories_to_path([root '/codaradm/HFR_Progs-2_1_3beta/matlab/'],{'CVS','private','@'});


%% Load HFR data
f = 'http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/PRVI/6km/hourly/RTV/HFRADAR_Puerto_Rico_and_the_US_Virgin_Islands_6km_Resolution_Hourly_RTV_best.ncd';

%% Set date & time variables.

dnow= ((round(now*24))/24);
% dtime = dnow-12:1/24:dnow-10;
dtime = datenum(2019,05,23):1/24:datenum(2019,05,25);
%% Set Latitude & Longitude for figure window limits.

conf.HourPlot.axisLims=[-68 -65 17 18.25];
conf.HourPlot.DomainName='Puerto Rico';
conf.HourPlot.Print=false;
conf.HourPlot.grid=1;

%% Create a Directory for all Processed Imagery 
%AFAI process imagery directories are created within the function
%AFAI_img2coordhighres_v1.m

%Directory for all k-means clustered algorithem imagery
HFR_dir = [pwd '\HFR_animations\'];
if ~exist(HFR_dir, 'dir')
    mkdir(HFR_dir);
end

%PNG images to be turned into animations
HFR_imgs = [HFR_dir  'HFR_imgs\' ];
if ~exist(HFR_imgs, 'dir')
    mkdir(HFR_imgs);
end

%High Res Particle Tracking folder
gifs = [HFR_dir  'HFR_gifs\' ];
if ~exist(gifs, 'dir')
    mkdir(gifs);
end


%conf.Totals.DomainName

conf.Totals.grid_lims=[conf.HourPlot.axisLims(1) conf.HourPlot.axisLims(2) ...
    conf.HourPlot.axisLims(3) conf.HourPlot.axisLims(4)];

[TUVcat] = convert_NN_NC_to_TUVstruct(f,dtime,conf);

%% Grid the total data onto a rectangular grid
[TUVcat,dim]=gridTotals(TUVcat,0,0);

X=TUVcat.LonLat(:,1);
Y=TUVcat.LonLat(:,2);
U=TUVcat.U;
V=TUVcat.V;
tt=TUVcat.TimeStamp;
tspan=TUVcat.TimeStamp(1):1/24:TUVcat.TimeStamp(end);


wp1=[17+30/60+0/3600 -67-30/60-00/3600];
%wp1=[35.25 -75];
resolution=10;
range=200;
 
[wp]=release_point_generation_matrix(wp1,resolution,range);

drifter=[wp(:,2) wp(:,1)];

[X1,Y1]=meshgrid(unique(X),unique(Y));

[r,c]=size(X1);

U1=reshape(U,r,c,length(TUVcat.TimeStamp));
V1=reshape(V,r,c,length(TUVcat.TimeStamp));

%% generate the particle trajectories
[x,y,ts]=particle_track_ode_grid_LonLat(X1,Y1,U1,V1,tt,tspan,drifter);

[r1,c1]=size(x);

%% determine which position estimates are Nans to use in the coloring of the 
%% particles
color_flag=~isnan(x);

%% When the drifter leaves the domain the position turns into nans
%% replace the nans with the last known position
xx=x;
for ii=1:c1
    tf=find(isnan(x(:,ii)),1,'first');
    x(tf:end,ii)=x(tf-1,ii);
    y(tf:end,ii)=y(tf-1,ii);
end

%% plot the results

%WHERE THE PNG FILES ARE BEING SAVED TO
%conf.Plot.BaseDir='C:\Users\Joe Anarumo\Documents\Puerto_Rico\Sargassum\Test_trajectories\Realtime_HFR\';

dir_sav=datestr(dnow,'yyyymmdd');
conf.Plot.BaseDir=HFR_imgs;

f1=[root '/jpa104/caricoos/etopo1_Puerto_Rico.nc'];
 
[LON,LAT,Z] = read_in_etopo_bathy(f1);
% bathy=load ('/Users/roarty/Documents/GitHub/HJR_Scripts/data/bathymetry/puerto_rico/puerto_rico_6second_grid.mat');
% ind2= bathy.depthi==99999;
% bathy.depthi(ind2)=NaN;
bathylines=[ -50 -100 -500 -1000 -2000 -3000 -4000 -5000];

close all

for ii=1:r1
hold on
m_proj('albers equal-area','lat',conf.HourPlot.axisLims(3:4),'long',conf.HourPlot.axisLims(1:2),'rect','on');
m_gshhs_f('patch',[240,230,140]./255);
m_grid('box','fancy','tickdir','in','xaxisloc','bottom','yaxisloc','left');

%% plot bathymetry
[cs, h1] = m_contour(LON,LAT, Z,bathylines);
clabel(cs,h1,'fontsize',8,'Color',[0.8 0.8 0.8]);
set(h1,'LineColor',[0.8 0.8 0.8])


%% plot entire track in gray
if ii>=2 
    m_plot(x(1:ii,:),y(1:ii,:),'-','Color',[0.5 0.5 0.5]);
end 

%% plot last 6 hours in black 
if ii>6
    m_plot(x(ii-6:ii,:),y(ii-6:ii,:),'-','Color','k','LineWidth',2);
elseif ii>1
    m_plot(x(1:ii,:),y(1:ii,:),'-','Color','k','LineWidth',2);
end

%% plot last location of drifter in red
m_plot(x(ii,:),y(ii,:),'ro','MarkerFaceColor','r');

if ~isempty(find(isnan(xx(ii,:)), 1)) 
    m_plot(x(ii,find(isnan(xx(ii,:)))),y(ii,find(isnan(xx(ii,:)))),'bo','MarkerFaceColor','b');
end 

%% add zeros before the number string
N=append_zero(ii);

%%-------------------------------------------------
%% Add title string

conf.HourPlot.TitleString = [' Particle Trajectories: ', ...
                            datestr(tspan(ii),'mm/dd/yyyy HH:MM'),' ',TUVcat.TimeZone(1:3)];

hdls.title = title( conf.HourPlot.TitleString, 'fontsize', 12,'color',[0 0 0] );

timestamp(1,'trajectories_HFR.m')

if ~exist(conf.Plot.BaseDir, 'dir')
        mkdir(conf.Plot.BaseDir)
    end

print(1,'-dpng','-r100',[ conf.Plot.BaseDir conf.HourPlot.DomainName '_' datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'])

close all
clear conf.HourPlot.TitleString

end

%% TURN DAILY FILES INTO ANIMATED GIFS

for ii=1:length(tspan)
    if ii==1
        file_NameTemp = [ 'Puerto Rico' '_'  datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'];
file_name = file_NameTemp;
    else
        file_NameTemp = [ 'Puerto Rico' '_'  datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'];
        file_name = [file_name {file_NameTemp}];
    end
end
%WHERE THE PNG FILES CAN BE FOUND
file_path = HFR_imgs;
file_name=sort(file_name);

file_name2 = [ 'Puerto Rico' '_'  datestr(tspan(1),'yyyy_mm_dd') '.gif'];
%WHERE ANIMATED GIFS ARE BEING SAVED TO
file_path2 = gifs;

loops=65535;    %Forever loop
delay=0;        %delay between images

h = waitbar(0,['0% done'],'name','Progress') ;
for i=1:length(file_name)
    if strcmpi('gif',file_name{i}(end-2:end))
        [M  c_map]=imread([file_path,file_name{i}]);
    else
        a=imread([file_path,file_name{i}]);
        [M  c_map]= rgb2ind(a,256);
    end
    if i==1
        imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay)
    elseif i==length(file_name)
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
    else
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
    end
    waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;
end
close(h);
