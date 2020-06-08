clc;clear all;close all;


%This script corrects the conversion of Esoiland ET from Noah and Mosaic.
%In the original converion, Lv was assumed to be constant whereas here, I
%use equation 4 from Henderson-Sellers, 1986 where Lv is a function of temp.
Lv=@(T) ( 2500.82 - 2.358*(T-273) ) * 10^3; %J/kg (T in K)
%The old Lv was assumed to be: 2.501e6
Lv_original=2.501e06;

%load in ESMAP points
Points=importdata('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
Points(:,2)=Points(:,2)+360;
Points=round(Points,5);
npoints=length(Points);

for j=1:2
    if j==1
        model='Noah';
        multiplier=1;
    elseif j==2
        model='Mosaic';
        multiplier=-1;
    end
    
    for i=1:npoints
        i
        lat=Points(i,1);
        lon=Points(i,2);
        %Load in NLDAS2 forcing
        forcing_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP/%.15g/%.15g/NLDAS_Forcing_raw',lat,lon);
        forcing_data=importdata(forcing_filename);
        %get temp time series for points
        temp=forcing_data(:,6);
        %Convert from hourly to daily temp:
        [u,~,j]=unique(forcing_data(:,1:3),'rows','stable');
        Daily_Temp=[u,accumarray(j,temp,[],@nanmean)];
        %Load in Noah time series match to ESMAP:
        Esoil_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/%s_Esoil.csv',lat,lon,model);
        ET_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/%s_ET.csv',lat,lon,model);
        Ec_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/%s_Ec.csv',lat,lon,model);
        Trans_filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/%s_Trans.csv',lat,lon,model);
        if exist(ET_filename,'file')>0
            Esoil=csvread(Esoil_filename);
            ET=csvread(ET_filename);
            Ec=csvread(Ec_filename);
            Trans=csvread(Trans_filename);
            %temporally match forcing eith Noah time series:
            [ia,ib]=ismember(u,Esoil(:,1:3),'rows');
            Daily_Temp=Daily_Temp(ia,:);
            %calc new Lv
            Lv_new=Lv(Daily_Temp(:,4));
            Esoil(:,4)=Esoil(:,4)*Lv_original./Lv_new *multiplier;
            Ec(:,4)=Ec(:,4)*Lv_original./Lv_new *multiplier;
            Trans(:,4)=Trans(:,4)*Lv_original./Lv_new *multiplier;
            
%             ET(:,4)=ET(:,4)*Lv_original./Lv_new;
%             ET(:,4)=ET(:,4)*multiplier;
            dlmwrite(Esoil_filename,Esoil,'delimiter',',','precision',8);
            dlmwrite(Ec_filename,Ec,'delimiter',',','precision',8);
            dlmwrite(Trans_filename,Trans,'delimiter',',','precision',8);
%             dlmwrite(ET_filename,ET,'delimiter',',','precision',8);
        end
        
    end
end
