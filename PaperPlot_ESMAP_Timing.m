clc;clear all;close all;

%the goal of this figure is to compliment the Screening figure in Fig. 2
%showing the timing of ESMAP and a time series of a point, Little Washita,
%showing when SMAP obs are vs. when Esoil is estimates, and then when
%precip screening occurs.

%get precip from older ESMAP product that is in mm/interval:
%load in ESMAP precipitation
ESMAP_precip=ncread('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip_screen_precip.nc','precip');
%load in ESMAP slength:
slength=ncread('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip_screen_precip.nc','slength');

%load in ESMAP Esoil
% esoil=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','esoil_screened');
esoil=ncread('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip_screen_precip.nc','esoil_screened');
%load in ESMAP time
time=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','time');
%load in ESMAP grid:
lat=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','lat');
lon=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','lon');
[LAT,LON]=meshgrid(lat,lon);
latvec=LAT(:);
lonvec=LON(:);
%define test point, Little Washita: (34.88,-98.09)
lat=34.88;lon=-98.09;
IDX=knnsearch([latvec,lonvec],[lat,lon]);
[r,c]=ind2sub(size(LAT),IDX);
site_lat=latvec(IDX); site_lon=lonvec(IDX);
%define slength and esoil from point:
ntimes=length(time);
start_date=datenum([2015 3 31]);
dates=start_date+time;

%initialize data:
site_esoil=[];
site_slength=[];
site_precip=[];
for t=1:ntimes
    site_esoil=[site_esoil;esoil(r,c,t)];
    site_slength=[site_slength;slength(r,c,t)];
    site_precip=[site_precip;ESMAP_precip(r,c,t)];
end
site_slength_original=site_slength;

idx_nan=isnan(site_esoil);
site_esoil(idx_nan)=[];
site_slength(idx_nan)=[];
valid_dates=dates;
valid_dates(idx_nan)=[];


ndata=length(site_esoil);
store_obs_dates=[];
for i=1:ndata
    valid_date=valid_dates(i);
    current_slength=site_slength(i);
    if rem(current_slength,2)==0
        offset=6/24;
    else
        offset=18/24;
    end
    valid_date=valid_date+offset;
    SMAP_obs1=valid_date-current_slength/2;
    SMAP_obs2=valid_date+current_slength/2;
    store_obs_dates=[store_obs_dates;SMAP_obs1,SMAP_obs2];
end

%store obs dates with screened precip intervals:
ndata=length(site_slength_original);
store_screened_dates=[];
for i=1:ndata
    date=dates(i);
    current_slength=site_slength_original(i);
    if rem(current_slength,2)==0
        offset=6/24;
    else
        offset=18/24;
    end
    date=date+offset;
    SMAP_obs1=date-current_slength/2;
    SMAP_obs2=date+current_slength/2;
    
% %     if i==1 || site_precip(i)>=2
    if site_precip(i)>=2
% %         if i==1
% %             store_screened_dates=[store_screened_dates;datenum([2015 4 1 6 0 0]),datenum([2015 4 3 6 0 0])];
% %         else
% %             store_screened_dates=[store_screened_dates;SMAP_obs1,SMAP_obs2];
% %         end
        store_screened_dates=[store_screened_dates;SMAP_obs1,SMAP_obs2];
    end
end


%Load in original obs timeseries
SMAP_obs_raw_QC_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP/%.15g/%.15g/SMAP_raw_QC.xls',site_lat,site_lon+360);
SMAP_raw_QC=xlsread(SMAP_obs_raw_QC_filename);
idx=isnan(SMAP_raw_QC(:,4));
SMAP_raw_QC(idx,:)=[];
SMAP_raw_QC_dates=datenum(SMAP_raw_QC(:,1:3))+6/24;
%trim dates to month period (April 1-May 1 2015)
start_date=datenum([2015 9 1]);
end_date=datenum([2015 10 1]);
idx=find(SMAP_raw_QC_dates>=start_date & SMAP_raw_QC_dates<end_date);
SMAP_raw_QC_dates=SMAP_raw_QC_dates(idx);
SMAP_raw_QC=SMAP_raw_QC(idx,:);


SMAP_obs_raw_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP/%.15g/%.15g/SMAP_raw.xls',site_lat,site_lon+360);
SMAP_raw=xlsread(SMAP_obs_raw_filename);
idx=isnan(SMAP_raw(:,4));
SMAP_raw(idx,:)=[];
SMAP_raw_dates=datenum(SMAP_raw(:,1:3))+6/24;
idx=find(SMAP_raw_dates>=start_date & SMAP_raw_dates<end_date);
SMAP_raw_dates=SMAP_raw_dates(idx);
SMAP_raw=SMAP_raw(idx,:);

%Load in forcing (precip) time series
Forcing_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP/%.15g/%.15g/NLDAS_Forcing_raw',site_lat,site_lon+360);
Forcing=importdata(Forcing_filename);
precip=Forcing(:,8);
Forcing_dates=datenum([Forcing(:,1:4),zeros(length(Forcing(:,1)),2)]);
%include local time offset:
Forcing_dates=Forcing_dates+7/24;
%trime to plot month:
idx=find(Forcing_dates>=start_date & Forcing_dates<end_date);
Forcing_dates=Forcing_dates(idx);
precip=precip(idx);

%calculate ESMAP exact date:
ESMAP_dates=mean(store_obs_dates');
idx=find(ESMAP_dates>=start_date & ESMAP_dates<end_date);
ESMAP_dates=ESMAP_dates(idx);
site_slength=site_slength(idx);

%Create plot:
f=figure;
yyaxis left
hold on
%first plot vertical lines with original SMAP Obs:
for i=1:length(SMAP_raw_dates)
    SMAP_raw=xline(SMAP_raw_dates(i),'linewidth',4,'color','r');  
end

%second plot vertical lines with SMAP Obs valid after QC flag:
for i=1:length(SMAP_raw_QC_dates)
    SMAP_QC=xline(SMAP_raw_QC_dates(i),'linewidth',4,'color','k');
end

%third plot vertical dashed lines with ESMAP recording:
for i=1:length(ESMAP_dates)
    ESMAP=plot(linspace(ESMAP_dates(i),ESMAP_dates(i),10),linspace(0,20,10),'--r','linewidth',1);
end
set(gca,'YTick',[])
ylim([0 1])
set(gca,'fontsize',30)
%highlight in green the intervals screened for precip
for i=1:length(store_screened_dates)
    x_Left=ones(10,1)*store_screened_dates(i,1);
    x_Right=ones(10,1)*store_screened_dates(i,2);
    Vertical=linspace(0,1,10);
    X=[x_Left;x_Right];
    inBetween=[Vertical,fliplr(Vertical)];
    h0=fill(X,inBetween,[0.85 0.85 0.85]);
end


%fourth plot precip atop:
yyaxis right
right_color= [0 0 1]; %blue for precip
b=bar(Forcing_dates,precip,'EdgeColor','b','FaceColor','b','LineWidth',0.01,'BarWidth',0.01);
axR=gca;
axR.YColor=right_color;
ylabel({'Precipitation'; 'mm hr^{-1}'},'Fontsize',36)
hold off
set(gca,'ydir','reverse')
ylim([0 10])
xlim([start_date end_date])
set(gca,'fontsize',30)
xticks([start_date:3:end_date])
datetick('x','mm/dd/yy','keeplimits','keepticks')
xtickangle(45)
% legend([SMAP_QC,ESMAP],{'   SMAP retreival screened based on quality flag','   E-SMAP E_{soil}'},'fontsize',24)
title('Little Washita','fontsize',30)
grid on
grid minor
box on


% % %make sure this aligns with original SMAP data:
% % %get grid:
% % lat=ncread('/Volumes/LabShare/Global/SMAP/SMAP.2015.04.30.nc','lat');
% % lon=ncread('/Volumes/LabShare/Global/SMAP/SMAP.2015.04.30.nc','lon');
% % [LAT,LON]=meshgrid(lat,lon);
% % latvec=LAT(:);
% % lonvec=LON(:);
% % idx=knnsearch([latvec,lonvec],[site_lat,site_lon]);
% % [r,c]=ind2sub(size(LAT),idx);
% % indir='/Volumes/LabShare/Global/SMAP/';
% % store_qual=[];
% % store_SMAP=[];
% % for year=2015:2015
% %     for month=4:4
% %         for day=1:30
% %             filename=sprintf('SMAP.%04d.%02d.%02d.nc',year,month,day);
% %             data=ncread([indir,filename],'SMAP',[r,c],[1,1]);
% %             qual_flag=ncread([indir,filename],'qual_flag',[r,c],[1,1]);
% %             store_SMAP=[store_SMAP;data];
% %             store_qual=[store_qual;qual_flag];
% %         end
% %     end
% % end
           