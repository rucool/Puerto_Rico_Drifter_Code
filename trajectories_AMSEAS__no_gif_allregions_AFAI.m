

%% Script to generate particle trajectories of Sargassum seaweed near Puerto Rico.
% Written by Joseph Anarumo on 3/19/2019

%% Clear workspace

close all
clear all
clc

tic

% Determine which computer you are working on

compType=computer;

if ~isempty(strmatch('PCWIN64',compType))
    root='/Volumes';
else
    root='/home';
end

%% Add directories to MATLAB path for following functions.

% addpath([root '/codaradm/HFR_Progs-2_1_3beta/matlab/external_matlab_packages/m_map/'])
% addpath([root '/codaradm/HFR_Progs-2_1_3beta/matlab/general']);
% addpath([root '/codaradm/HFR_Progs-2_1_3beta/matlab/trajectories']);
% addpath([root '/codaradm/HFR_Progs-2_1_3beta/matlab/totals']);
%add_subdirectories_to_path([root '/codaradm/HFR_Progs-2_1_3beta/matlab/'],{'CVS','private','@'});


%% Load NOAA data
f = 'https://ecowatch.ncddc.noaa.gov/thredds/dodsC/ncom_amseas_agg/AmSeas_Apr_05_2013_to_Current_best.ncd';

%% Set date & time variables.

dnow= ((round(now*24))/24);%nearest hour
dtime = dnow:1/24:dnow+1;
% dtime = datenum(2019,6,1):1/24:datenum(2019,6,2);

%% Create a Directory for all Processed Imagery

%Directory for all AMSEAS imagery
% AMSEAS_dir = '/Users/roarty/COOL/01_CODAR/CARICOOS/AMSEAS_animations/';
AMSEAS_dir = '/Volumes/hroarty/public_html/caricoos/animations/';
% if ~exist(AMSEAS_dir, 'dir')
%     mkdir(AMSEAS_dir);
% end

%PNG images to be turned into a

nimations
% AMSEAS_imgs = [AMSEAS_dir  'AMSEAS_imgs/' ];
% if ~exist(AMSEAS_imgs, 'dir')
%     mkdir(AMSEAS_imgs);
% end

% load the AFAI data
f2='https://cwcgom.aoml.noaa.gov/erddap/griddap/noaa_aoml_atlantic_oceanwatch_AFAI_3D';
lims=[-68 -65 17 19];
dtime2=datestr(dnow-1);

[ERD]=AFAI_ERDDAP_data_fn(f2,dtime2,lims);




%% Set Latitude & Longitude for figure window limits for each study region.

sites={'PR_North','PR_East','PR_South','PR_LaPar','PR_West','VI_North','VI_South','PR_NorthEast'};
% sites={'PR_NorthEast'};
 
for gg=1:length(sites)
 
n=sites{gg};

conf.Plot.BaseDir=[AMSEAS_dir n '/'];
 
switch n
    case 'PR_North'
        conf.HourPlot.axisLims=[-66-20/60 -65-54/60 18+25/60 18+45/60];
        conf.HourPlot.DomainName='Puerto_Rico_North';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_North.mat';
        conf.HourPlot.Title='Puerto Rico North ';
    case 'PR_East'
        conf.HourPlot.axisLims=[-66 -65-34/60 17+50/60 18+10/60];
        conf.HourPlot.DomainName='Puerto_Rico_East';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_East.mat';
        conf.HourPlot.Title='Puerto Rico East ';
    case 'PR_South'
        conf.HourPlot.axisLims=[-66-50/60 -66-24/60 17+40/60 18];
        conf.HourPlot.DomainName='Puerto_Rico_South';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_South.mat';
        conf.HourPlot.Title='Puerto Rico South ';
    case 'PR_LaPar'
        conf.HourPlot.axisLims=[-67-13/60 -66-47/60 17+40/60 18];
        conf.HourPlot.DomainName='Puerto_Rico_LaParguera';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_LaParguera.mat';
        conf.HourPlot.Title='Puerto Rico La Parguera ';
    case 'PR_West'
        conf.HourPlot.axisLims=[-67-36/60 -67-10/60 17+50/60 18+10/60];
        conf.HourPlot.DomainName='Puerto_Rico_West';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_West.mat';
        conf.HourPlot.Title='Puerto Rico West ';
    case 'VI_North'
        conf.HourPlot.axisLims=[-65-14/60 -64-48/60 18+20/60 18+40/60];
        conf.HourPlot.DomainName='Virgin_Islands_North';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Virgin_Islands_North.mat';
        conf.HourPlot.Title='Puerto Rico North ';
    case 'VI_South'
        conf.HourPlot.axisLims=[-65-14/60 -64-48/60 18 18+20/60];
        conf.HourPlot.DomainName='Virgin_Islands_South';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Virgin_Islands_South.mat';
        conf.HourPlot.Title='Puerto Rico South ';
    case 'PR_NorthEast'
        conf.HourPlot.axisLims=[-65-40/60 -65-18/60 18+8/60 18+26/60];
        conf.HourPlot.DomainName='Puerto_Rico_NorthEast';
        conf.HourPlot.coast='/Users/roarty/data/coastline/Puerto_Rico_NorthEast.mat.mat';
        conf.HourPlot.Title='Puerto Rico NorthEast ';
end
 


% conf.HourPlot.axisLims=[-66 -65-36/60 17+50/60 18+10/60];

conf.HourPlot.Print=false;
conf.HourPlot.grid=1;

%% Create a Directory for all Processed Imagery in region and day folers

%Sort by region
% regionfolder = [AMSEAS_imgs n '/'];
% if ~exist(regionfolder, 'dir')
%     mkdir(regionfolder);
% end

% %Sort by date
% dir_save = datestr(dtime(1), 'yyyymmdd');
% dayfolder = [regionfolder dir_save '/'];
% if ~exist(dayfolder, 'dir')
%     mkdir(dayfolder);
% end

%High Res Particle Tracking folder
% gifs = [AMSEAS_dir  'AMSEAS_gifs/' ];
% if ~exist(gifs, 'dir')
%     mkdir(gifs);
% end
% 
% % Also Sort gifs by region
% regionfolder2 = [gifs n '/'];
% if ~exist(regionfolder2, 'dir')
%     mkdir(regionfolder2);
% end

%conf.Totals.DomainName

conf.Totals.grid_lims=[(conf.HourPlot.axisLims(1)+360) (conf.HourPlot.axisLims(2)+360) ...
    conf.HourPlot.axisLims(3) conf.HourPlot.axisLims(4)];

[TUVcat] = convert_AmSeas_NC_to_TUVstruct(f,dtime,conf);


%% Grid the total data onto a rectangular grid
[TUVcat,dim]=gridTotals(TUVcat,0,0);

X=TUVcat.LonLat(:,1);
Y=TUVcat.LonLat(:,2);
U=TUVcat.U;
V=TUVcat.V;
tt=TUVcat.TimeStamp(1:end);
tspan=TUVcat.TimeStamp(1):1/24:TUVcat.TimeStamp(end);

%% Define the points where we will release drifters
%wp1=[17+52/60+0/3600 -65-55/60-0/3600];
wp1=[conf.HourPlot.axisLims(3)+2/60 conf.HourPlot.axisLims(1)+5/60];%lat lon
%wp1=[35.25 -75];
resolution=5;
range=25;
 
[wp]=release_point_generation_matrix(wp1,resolution,range);

drifter=[wp(:,2) wp(:,1)];%[lon lat]

% color code the drifters by AFAI concentration
% dc=interp2(ERD.LON,ERD.LAT,ERD.AFAI,wp(:,2),wp(:,1));
dc=griddata(ERD.LON,ERD.LAT,ERD.AFAI,wp(:,2),wp(:,1));

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
% conf.Plot.BaseDir='C:\Users\Joe Anarumo\Documents\Puerto_Rico\Sargassum\Test_trajectories\Realtime_AMSEAS\';




% f1=[root '/jpa104/caricoos/etopo1_Puerto_Rico.nc'];
f1='/Users/roarty/data/bathymetry/etopo1_Puerto_Rico.nc';
 
[LON,LAT,Z] = read_in_etopo_bathy(f1);
% bathy=load ('/Users/roarty/Documents/GitHub/HJR_Scripts/data/bathymetry/puerto_rico/puerto_rico_6second_grid.mat');
% ind2= bathy.depthi==99999;
% bathy.depthi(ind2)=NaN;
bathylines=[ -50 -100 -500 -1000 -2000 -3000 -4000 -5000];

close all

for ii=1:r1
hold on
% m_proj('albers equal-area','lat',conf.HourPlot.axisLims(3:4),'long',conf.HourPlot.axisLims(1:2),'rect','on');
% m_gshhs_f('patch',[240,230,140]./255);
% m_grid('box','fancy','tickdir','in','xaxisloc','bottom','yaxisloc','left');

plotBasemap(conf.HourPlot.axisLims(1:2),conf.HourPlot.axisLims(3:4),conf.HourPlot.coast,'Mercator','patch',[240,230,140]./255,'edgecolor','k');

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
% m_plot(x(ii,:),y(ii,:),'ro','MarkerFaceColor','r');
m_scatter(x(ii,:),y(ii,:),50,dc,'filled');

if ~isempty(find(isnan(xx(ii,:)), 1)) 
    m_plot(x(ii,find(isnan(xx(ii,:)))),y(ii,find(isnan(xx(ii,:)))),'ko','MarkerFaceColor','k');
end 

%% add zeros before the number string
N=append_zero(ii);

colormap(jet)
caxis([-0.004 0.006])
colorbar

%%-------------------------------------------------
%% Add title string

conf.HourPlot.TitleString = [conf.HourPlot.Title,' Particle Trajectories: ', ...
                            datestr(tspan(ii),'mm/dd/yyyy HH:MM'),' ',TUVcat.TimeZone(1:3)];

hdls.title = title( conf.HourPlot.TitleString, 'fontsize', 12,'color',[0 0 0] );

timestamp(1,'trajectories_AMSEAS.m')

% if ~exist(conf.Plot.BaseDir, 'dir')
%     mkdir(conf.Plot.BaseDir)
% end

% print(1,'-dpng','-r100',[ conf.Plot.BaseDir conf.HourPlot.DomainName '_' datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'])
print(1,'-dpng','-r100',[ conf.Plot.BaseDir conf.HourPlot.DomainName '_' sprintf('%02d',ii) '.png'])

close all
clear conf.HourPlot.TitleString

end


%% TURN DAILY FILES INTO ANIMATED GIFS

% for ii=1:length(tspan)
%     if ii==1
%         file_NameTemp = [ conf.HourPlot.DomainName '_' datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'];
% file_name = file_NameTemp;
%     else
%         file_NameTemp = [ conf.HourPlot.DomainName '_' datestr(tspan(ii),'yyyy_mm_dd_HHMM') '.png'];
%         file_name = [file_name {file_NameTemp}];
%     end
% end
% %WHERE THE PNG FILES CAN BE FOUND
% file_path = dayfolder;
% file_name=sort(file_name);
% 
% file_name2 = [ 'Puerto Rico' '_'  datestr(tspan(1),'yyyy_mm_dd') '.gif'];
% %WHERE ANIMATED GIFS ARE BEING SAVED TO
% file_path2 = regionfolder2;
% 
% loops=65535;    %Forever loop
% delay=0;        %delay between images
% 
% h = waitbar(0,['0% done'],'name','Progress') ;
% for i=1:length(file_name)
%     if strcmpi('gif',file_name{i}(end-2:end))
%         [M  c_map]=imread([file_path,file_name{i}]);
%     else
%         a=imread([file_path,file_name{i}]);
%         [M  c_map]= rgb2ind(a,256);
%     end
%     if i==1
%         imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay)
%     elseif i==length(file_name)
%         imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
%     else
%         imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
%     end
%     waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;
% end
% close(h);
end

toc

exit
