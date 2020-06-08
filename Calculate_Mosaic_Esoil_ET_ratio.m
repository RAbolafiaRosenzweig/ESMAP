clc;clear all;close all;

%This scripts checks to see if considering Esoil/ET only during 'valid'
%ESMAP intervals is consisten with Esoil/ET if not screening were to have
%occurred in the Mosaic product.


%%====================================================================================
%This script corrects the conversion of Esoiland ET from Noah and Mosaic.
%In the original converion, Lv was assumed to be constant whereas here, I
%use equation 4 from Henderson-Sellers, 1986 where Lv is a function of temp.
Lv=@(T) ( 2500.82 - 2.358*(T-273) ) * 10^3; %J/kg (T in K)

% % %load in monthly Mosaic ET for the ESMAP 04/2015-03/31/2019
indir_monthly='/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Monthly/';

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
        filename=[indir_monthly,sprintf('NLDAS_MOSAIC0125_H.A%04d%02d.nc',year, month)];
        Mosaic_ET=ncread(filename,'var57_EVP');
        Mosaic_ET=Mosaic_ET*24; %mm/hr --> mm/day
        %load esoil data
        Mosaic_Esoil=ncread(filename,'var199_EVBS');
        forcing_filename=sprintf('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/FORCING_Monthly/NLDAS_FORB0125_H.A%04d%02d.nc',year,month);
        Temp=ncread(forcing_filename,'var11_TMP');
        Lvp=Lv(Temp);
        Mosaic_Esoil=Mosaic_Esoil./Lvp*86400*-1; %W/m^2 --> mm/day
        %calculate Esoil/ET ratio
        Mosaic_Esoil_ET_ratio=Mosaic_Esoil./Mosaic_ET;
        
        %store the ratio:
        store_monthly_ratio=cat(3,store_monthly_ratio,Mosaic_Esoil_ET_ratio);
        store_monthly_esoil=cat(3,store_monthly_esoil,Mosaic_Esoil);
        
    end
end
%%====================================================================================

%Now store the ratio for the screened version of Mosaic, mathcing ESMAP
%overpasses

Points=importdata('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
Points=round(Points,5);

npoints=length(Points);
%create grid to store outputs on:
lat=ncread(filename,'lat');
lon=ncread(filename,'lon');
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
    lon=Points(i,2)+360;
    filename_ET=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/Mosaic_ET.csv',lat,lon);
    filename_Esoil=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/Mosaic_Esoil.csv',lat,lon);
    if exist(filename_ET,'file')>0 && i~=23459 && i~=28977
        Mosaic_ET=csvread(filename_ET);
        Mosaic_Esoil=csvread(filename_Esoil);
        Esoil_ET_ratio=Mosaic_Esoil(:,4)./Mosaic_ET(:,4);
        %concat to monthly
        [u,~,j]=unique(Mosaic_ET(:,1:2),'rows','stable');LAT
        Monthly_screened_ratio=[accumarray(j,Esoil_ET_ratio,[],@nanmean)];
        Monthly_screened_Esoil=[accumarray(j,Mosaic_Esoil(:,4),[],@nanmean)];
        %I checked to make sure the first site had all months:
        if i==1
            store_months=u;
        end
        %store the mean monthly ratio of Esoil/ET
        empty_grid(IDX(i))=mean(Monthly_screened_ratio);
        empty_grid_esoil(IDX(i))=mean(Monthly_screened_Esoil);
        %remove months of the the unscreened Mosaic data that does not have
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
if exist('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio','dir')==0
    CMD='mkdir -p /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio';
    system(CMD);
end

dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_Screened.csv',empty_grid,'delimiter',',','precision',15);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_Screened.csv',empty_grid_esoil,'delimiter',',','precision',15);

%write out variable for future plotting:
mean_monthly_ratio=nanmean(store_monthly_ratio,3);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_noScreening.csv',mean_monthly_ratio,'delimiter',',','precision',15);
mean_monthly_esoil=nanmean(store_monthly_esoil,3);
dlmwrite('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_noScreening.csv',mean_monthly_esoil,'delimiter',',','precision',15);

