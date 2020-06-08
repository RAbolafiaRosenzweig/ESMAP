clc;clear all;close all;

%define Mosaic in directory:
indir='/scratch/summit/roab9675/SMAP/Validation_Data/NLDAS2/Mosaic/';
outdir='/scratch/summit/roab9675/SMAP/Gridded_ESMAP_blank/';

LH_vap= 2.501e6; % J kg^-1

%Load in ESMAP product that uses quality flags:
fileloc=getenv('SLURM_SCRATCH');
ESMAP_filename=[fileloc,'/ESMAP_QC.nc'];

ESMAP_Esoil_QC=ncread(ESMAP_filename,'esoil_screened');
ESMAP_time_QC=ncread(ESMAP_filename,'time');
ESMAP_slength_QC=ncread(ESMAP_filename,'slength');
ESMAP_lat_QC=ncread(ESMAP_filename,'lat');
ESMAP_lon_QC=ncread(ESMAP_filename,'lon');

[LAT,LON]=meshgrid(ESMAP_lat_QC,ESMAP_lon_QC);
[nrow ncol]=size(LAT);

lat_vec=reshape(LAT,nrow*ncol,1);
lon_vec=reshape(LON,nrow*ncol,1);

ESMAP_domain_points=[lat_vec,lon_vec];
ESMAP_domain_points(:,2)=ESMAP_domain_points(:,2)+360;
ESMAP_domain_points=round(ESMAP_domain_points,5);

%define dates:
start_date=datenum([2015 3 31]);
dates=start_date+ESMAP_time_QC;
dates=double(dates);
date_vec=datevec(dates);



%Load in point list for ESMAP product post QC flags:
% Points=importdata('/scratch/summit/roab9675/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
Points=importdata('/projects/roab9675/SMAP/9km_points/SMAP_9km_Points_iter999999');
Points=round(Points,5);



%find intersection between entire domain and vlaid points:
[C,ia,ib]=intersect(ESMAP_domain_points,Points,'rows');
%the point list shrinks by exclusion of points outside of domain:
%lat[25-50], lon[-125 - -67]
[ia_row,ia_col]=ind2sub(size(LAT),ia);

npoints=length(C);
%define UTC time adjustment column (the number required to ADD to local
%time o get UTC time)
UTC_offset=zeros(npoints,1);

idx_5=find(C(:,1)<=-82.5);
idx_6=find(C(:,1)>-82.5 & C(:,1)<=-97.5);
idx_7=find(C(:,1)>-97.5 & C(:,1)<=-112.5);
idx_8=find(C(:,1)>-112.5);

UTC_offset(idx_5)=5;
UTC_offset(idx_6)=6;
UTC_offset(idx_7)=7;
UTC_offset(idx_8)=8;

%loop through sites:
for s=1:npoints
    site_esoil=ESMAP_Esoil_QC(ia_row(s),ia_col(s),:);
    site_slength=ESMAP_slength_QC(ia_row(s),ia_col(s),:);
    site_UTC_offset=UTC_offset(s);
    site_lat=LAT(ia_row(s),ia_col(s));
    site_lon=LON(ia_row(s),ia_col(s));
    store_Mosaic_Esoil=[];
    for d=1:length(date_vec)
        %check if site has data at this time:
        current_esoil=site_esoil(1,1,d);
        idx=isnan(current_esoil);
        %if so calculate Mosaic esoil during this overpass interval:
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
            
            %Convert to UTC for NLDAS Time:
            start_date=start_date;
            end_date=end_date;
            
            start_datevec=datevec(start_date);
            end_datevec=datevec(end_date);
            assert(start_datevec(1,4)==6,'start overpass is not 6 am');
            assert(end_datevec(1,4)==6,'end overpass is not 6 am');
            
            start_datevec=start_datevec(:,1:3);
            start_date=datenum(start_datevec);
            start_date=start_date+1;
            end_datevec=end_datevec(:,1:3);
            end_date=datenum(end_datevec);
            
            %loop through start and end date accujulating Mosaic Esoil:
            timesteps=start_date:end_date;
            timesteps_datevec=datevec([timesteps]');
            Mosaic_Esoil_total=[];
            for ts=1:length(timesteps)
                Mosaic_year=timesteps_datevec(ts,1);
                Mosaic_month=timesteps_datevec(ts,2);
                Mosaic_day=timesteps_datevec(ts,3);
                validation_dir=sprintf('/gpfs/summit/scratch/roab9675/SMAP/Validation_Data/NLDAS2/MOSAIC_%dhr_Offset/',site_UTC_offset);
                Mosaic_filename=[validation_dir,sprintf('%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.nc',Mosaic_year,Mosaic_month,Mosaic_year,Mosaic_month,Mosaic_day)];
                Mosaic_Esoil_site=ncread(Mosaic_filename,'var199_EVBS',[ia_row(s),ia_col(s),1],[1,1,1]);
                Mosaic_Esoil_site=Mosaic_Esoil_site/LH_vap*86400; %W/m^2 --> mm/day
                Mosaic_Esoil_total=[Mosaic_Esoil_total;Mosaic_Esoil_site];
            end
            store_Mosaic_Esoil=[store_Mosaic_Esoil;year,month,day,nanmean(Mosaic_Esoil_total)];
        end
    end
    %export:
    if length(store_Mosaic_Esoil)>0
        outfilename=[outdir,sprintf('%.15g/%.15g/Mosaic_Esoil.csv',site_lat,site_lon+360)];
        disp(outfilename)
        dlmwrite(outfilename,store_Mosaic_Esoil,'delimiter',',','precision',8);
    end
end
