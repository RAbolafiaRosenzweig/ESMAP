clc;clear all;close all;

%This scripts checks to see if considering Esoil/ET only during 'valid'
%ESMAP intervals is consisten with Esoil/ET if not screening were to have
%occurred in the GLEAM product.


%%====================================================================================
% % %Using the non-screened data, create a spatial map of Esoil/ET over CONUS:
% %
% % %load in monthly GLEAM ET for the ESMAP 04/2015-03/31/2019
indir_ET='/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/ET/Monthly/';
indir_Esoil='/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil/Monthly/';

%initialize variable to plot:
store_monthly_ratio=[];
store_monthly_esoil=[];
for year=2015:2018
    if year==2015
        startmonth=4;
    else
        startmonth=1;
    end
    
    for month=startmonth:12
        %load ET data
        filename_ET=[indir_ET,sprintf('E_9km_%04d%02d_GLEAM_v3.3a.nc',year, month)];
        GLEAM_ET=ncread(filename_ET,'E');
        GLEAM_ET=GLEAM_ET';
        %load esoil data
        filename_Esoil=[indir_Esoil,sprintf('Eb_9km_%04d%02d_GLEAM_v3.3a.nc',year, month)];
        GLEAM_Esoil=ncread(filename_Esoil,'Eb');
        GLEAM_Esoil=GLEAM_Esoil';
        %calculate Esoil/ET ratio
        GLEAM_Esoil_ET_ratio=GLEAM_Esoil./GLEAM_ET;
        
        %store the ratio:
        store_monthly_ratio=cat(3,store_monthly_ratio,GLEAM_Esoil_ET_ratio);
        store_monthly_esoil=cat(3,store_monthly_esoil,GLEAM_Esoil);
    end
end
%%====================================================================================

%Now store the ratio for the screened version of GLEAM, mathcing ESMAP
%overpasses

Points=importdata('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
Points(:,2)=Points(:,2)+360;
Points=round(Points,5);
npoints=length(Points);
%create grid to store outputs on:
lat=ncread(filename_ET,'lat');
lon=ncread(filename_ET,'lon');
[LAT,LON]=meshgrid(lat,lon);
[nrow ncol]=size(LAT);
LAT_vec=round(reshape(LAT,nrow*ncol,1),5);
LON_vec=round(reshape(LON,nrow*ncol,1),5);
[ia,ib]=ismember([LAT_vec,LON_vec],Points,'rows');
IDX=find(ia==1);
assert(length(IDX)==length(Points),'points missing');
empty_grid=nan(nrow,ncol);
empty_grid_esoil=nan(nrow,ncol);

for i=1:npoints
    i
    lat=Points(i,1);
    lon=Points(i,2);
    filename_ET=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/GLEAM_ET.csv',lat,lon);
    filename_Esoil=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/GLEAM_Esoil.csv',lat,lon);
    if exist(filename_ET,'file')>0 && i~=23459 && i~=28977 && i~=32918 && i~=34232 && i~=40511 && i~=62119
        GLEAM_ET=csvread(filename_ET);
        GLEAM_Esoil=csvread(filename_Esoil);
        Esoil_ET_ratio=GLEAM_Esoil(:,4)./GLEAM_ET(:,4);
        %concat to monthly
        [u,~,j]=unique(GLEAM_ET(:,1:2),'rows','stable');
        Monthly_screened_ratio=[accumarray(j,Esoil_ET_ratio,[],@nanmean)];
        Monthly_screened_Esoil=[accumarray(j,GLEAM_Esoil(:,4),[],@nanmean)];
        %I checked to make sure the first site had all months:
        if i==1
            store_months=u;
        end
        %store the mean monthly ratio of Esoil/ET
        empty_grid(IDX(i))=mean(Monthly_screened_ratio);
        empty_grid_esoil(IDX(i))=mean(Monthly_screened_Esoil);
        %remove months of the the unscreened GLEAM data that does not have
        %corresponding ESMAP data:
        
        if length(Monthly_screened_ratio)~=45
            [r,c]=ind2sub([nrow,ncol],IDX(i));
            Unscreened=store_monthly_ratio(r,c,:);
            for m=1:45
                current_month=store_months(m,:);
                [ia,ib]=ismember(current_month,u,'rows');
                if ia==0
                    store_monthly_ratio(r,c,m)=NaN;
                    store_monthly_esoil(r,c,m)=NaN;
                end
            end
        end
    end
    
end
%write out variable for future plotting:
if exist('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/','dir')==0
    CMD='mkdir -p /Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/';
    system(CMD);
end

dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_Screened.csv',empty_grid,'delimiter',',','precision',15);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_Screened.csv',empty_grid_esoil,'delimiter',',','precision',15);

%write out variable for future plotting:
mean_monthly_ratio=nanmean(store_monthly_ratio,3);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_noScreening.csv',mean_monthly_ratio,'delimiter',',','precision',15);
mean_monthly_esoil=nanmean(store_monthly_esoil,3);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_noScreening.csv',mean_monthly_esoil,'delimiter',',','precision',15);

