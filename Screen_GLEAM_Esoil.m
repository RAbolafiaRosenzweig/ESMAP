clc;clear all;close all;


%define Noah in directory:
indir='/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil/';
outdir='/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/';
%Load in GLEAM files:
GLEAM_Esoil_2015=ncread([indir,'Eb_9km_2015_GLEAM_v3.3a.nc'],'Eb');
GLEAM_2015_time=ncread([indir,'Eb_9km_2015_GLEAM_v3.3a.nc'],'time');  
lat=ncread([indir,'Eb_9km_2015_GLEAM_v3.3a.nc'],'lat');
lon=ncread([indir,'Eb_9km_2015_GLEAM_v3.3a.nc'],'lon');
    
GLEAM_Esoil_2016=ncread([indir,'Eb_9km_2016_GLEAM_v3.3a.nc'],'Eb');
GLEAM_2016_time=ncread([indir,'Eb_9km_2016_GLEAM_v3.3a.nc'],'time'); 

GLEAM_Esoil_2017=ncread([indir,'Eb_9km_2017_GLEAM_v3.3a.nc'],'Eb');
GLEAM_2017_time=ncread([indir,'Eb_9km_2017_GLEAM_v3.3a.nc'],'time'); 

GLEAM_Esoil_2018=ncread([indir,'Eb_9km_2018_GLEAM_v3.3a.nc'],'Eb');
GLEAM_2018_time=ncread([indir,'Eb_9km_2018_GLEAM_v3.3a.nc'],'time'); 

%Load in ESMAP product that uses quality flags:
ESMAP_filename='/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip.nc';
ESMAP_esoil_QC=ncread(ESMAP_filename,'esoil_screened');
ESMAP_time_QC=ncread(ESMAP_filename,'time');
ESMAP_slength_QC=ncread(ESMAP_filename,'slength');
ESMAP_lat_QC=ncread(ESMAP_filename,'lat');
ESMAP_lon_QC=ncread(ESMAP_filename,'lon');

[LAT,LON]=meshgrid(ESMAP_lat_QC,ESMAP_lon_QC);
[nrow ncol]=size(LAT);

lat_vec=reshape(LAT,nrow*ncol,1);
lon_vec=reshape(LON,nrow*ncol,1);

ESMAP_domain_points=[lat_vec,lon_vec];
ESMAP_domain_points=round(ESMAP_domain_points,5);

%define dates:
start_date=datenum([2015 3 31]);
dates=start_date+ESMAP_time_QC;
%trim at end of 2018, end of GLEAM record:
idx=find(dates>=datenum([2019 1 1]));
dates(idx)=[];
dates=double(dates);
date_vec=datevec(dates);

%Load in point list for ESMAP product post QC flags:
Points=importdata('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
% Points=importdata('/projects/roab9675/SMAP/9km_points/SMAP_9km_Points_iter999999');
Points=round(Points,5);

%find intersection between entire domain and vlaid points:
[C,ia,ib]=intersect(ESMAP_domain_points,Points,'rows');
[ia_row,ia_col]=ind2sub(size(LAT),ia);
%the point list shrinks by exclusion of points outside of domain:
%lat[25-50], lon[-125 - -67]

npoints=length(C);


for s=1:npoints
    s
    site_Esoil=ESMAP_esoil_QC(ia_row(s),ia_col(s),:);
    site_slength=ESMAP_slength_QC(ia_row(s),ia_col(s),:);
    site_lat=LAT(ia_row(s),ia_col(s));
    site_lon=LON(ia_row(s),ia_col(s));
    store_GLEAM_Esoil=[];
    for d=1:length(date_vec)
         %check if site has data at this time:
         current_Esoil=site_Esoil(1,1,d);
         idx=isnan(current_Esoil);   
        %if so calculate GLEAM ET during this overpass interval:
        if idx==0
            year=date_vec(d,1);
            month=date_vec(d,2);
            day=date_vec(d,3);
            date=datenum([year month day]);
            current_length=site_slength(1,1,d);
            start_date=double(date-current_length/2);
            end_date=double(date+current_length/2);
            b=datevec(end_date);
            %if the hour of b=12 --> truncated at hour 18
            %if the hour of b=0  --> truncated at hour 6
            if b(4)==12
                date=date+18/24;
                start_date=double(date-current_length/2);
                end_date=double(date+current_length/2);
            elseif b(4)==0
                date=date+6/24;
                start_date=double(date-current_length/2);
                end_date=double(date+current_length/2);
            else
                disp('issue with hour')
                ERRORNOW
            end

            %end one day before 'end date'. This makes the GLEAM interval
            %equal length to the SMAP interval:
            end_date=end_date-30/24;
            start_date=start_date-6/24;
            timesteps=start_date:end_date;
            timesteps_datevec=datevec([timesteps]');
            
            idx=find(timesteps_datevec(:,1)>=2019);
            timesteps(idx)=[];
            timesteps_datevec(idx,:)=[];
            
            GLEAM_Esoil_total=[];
            for ts=1:length(timesteps)
                GLEAM_year=timesteps_datevec(ts,1);
                GLEAM_month=timesteps_datevec(ts,2);
                GLEAM_day=timesteps_datevec(ts,3);
                if GLEAM_year==2015
                    GLEAM_Esoil=GLEAM_Esoil_2015;
                    GLEAM_time=GLEAM_2015_time;
                elseif GLEAM_year==2016
                    GLEAM_Esoil=GLEAM_Esoil_2016;
                    GLEAM_time=GLEAM_2016_time;
                elseif GLEAM_year==2017
                    GLEAM_Esoil=GLEAM_Esoil_2017;
                    GLEAM_time=GLEAM_2017_time;
                elseif GLEAM_year==2018
                    GLEAM_Esoil=GLEAM_Esoil_2018;
                    GLEAM_time=GLEAM_2018_time;
                end
                GLEAM_dates=datenum([1970 1 1])+GLEAM_time;
                idx_date=find(GLEAM_dates==datenum([GLEAM_year,GLEAM_month,GLEAM_day]));
                GLEAM_Esoil=GLEAM_Esoil(:,:,idx_date);
                GLEAM_Esoil_site=GLEAM_Esoil(ia_col(s),ia_row(s)); %mm/day
                GLEAM_Esoil_total=[GLEAM_Esoil_total;GLEAM_Esoil_site];
            end
            store_GLEAM_Esoil=[store_GLEAM_Esoil;year,month,day,nanmean(GLEAM_Esoil_total)];
            
        end
    end
    %export:
    if length(store_GLEAM_Esoil)>0
        outfilename=[outdir,sprintf('%.15g/%.15g/GLEAM_Esoil.csv',site_lat,site_lon+360)];
        disp(outfilename)
        dlmwrite(outfilename,store_GLEAM_Esoil,'delimiter',',','precision',8);
    end
end
            
            