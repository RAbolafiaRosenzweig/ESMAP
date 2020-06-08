clc;clear all;close all;

%define bounding domain box
latlim=[24.9,50.1];
lonlim=[-125.1,-66.9];

%%=======================================================================================================================
%Mosaic
Mosaic_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_noScreening.csv');
Screened_Mosaic=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_Screened.csv');
%calculate the Esoil/EST ratio for screened vs. unscreened:
Screened_unScreened_ratio=Screened_Mosaic./Mosaic_og;
Mosaic_Scale_Factor=1./Screened_unScreened_ratio;
idx=isinf(Mosaic_Scale_Factor);
Mosaic_Scale_Factor(idx)=NaN;

%Noah
Noah_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_noScreening.csv');
Screened_Noah=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_Screened.csv');
%calculate the Esoil/EST ratio for screened vs. unscreened:
Screened_unScreened_ratio=Screened_Noah./Noah_og;
Noah_Scale_Factor=1./Screened_unScreened_ratio;

%GLEAM
%Plot data:
GLEAM_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_noScreening.csv');
Screened_GLEAM=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_Screened.csv');
%calculate the Esoil/EST ratio for screened vs. unscreened:
Screened_unScreened_ratio=Screened_GLEAM./GLEAM_og;
GLEAM_Scale_Factor=1./Screened_unScreened_ratio;
idx=isinf(GLEAM_Scale_Factor);
GLEAM_Scale_Factor(idx)=NaN;

%calculate ensemble mean scale factor:
SCALE_FACTOR_Ensemble_Mean=(Mosaic_Scale_Factor+Noah_Scale_Factor+GLEAM_Scale_Factor)./3;

%get grid to export to:
filename_ET='/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Monthly/NLDAS_MOSAIC0125_H.A201708.nc';
lat=ncread(filename_ET,'lat');
lon=ncread(filename_ET,'lon');
%%=======================================================================================================
%Export scale factors to netCDF data repository:

outfile='/Volumes/REESEN/SMAP/Data_Repository/ESMAP_SCREENING_SCALE_FACTOR.nc';
if exist(outfile,'file')>0
    CMD_rm=['rm ',outfile];
    system(CMD_rm);
end

%create latitude:
nccreate(outfile,'lat','Dimensions',{'lat' length(lat)});
ncwriteatt(outfile, 'lat', 'standard_name', 'latitude');
ncwriteatt(outfile, 'lat', 'long_name', 'latitude');
ncwriteatt(outfile, 'lat', 'units', 'degrees north');
ncwriteatt(outfile, 'lat', '_CoordinateAxisType', 'Lat');
%create longitude:
nccreate(outfile,'lon','Dimensions',{'lon' length(lon)});
ncwriteatt(outfile, 'lon', 'standard_name', 'longitude');
ncwriteatt(outfile, 'lon', 'long_name', 'longitude');
ncwriteatt(outfile, 'lon', 'units', 'degrees north');
ncwriteatt(outfile, 'lon', '_CoordinateAxisType', 'Lon');

%create Mosaic data
nccreate(outfile,'Mosaic_Scale_Factor','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
ncwriteatt(outfile, 'Mosaic_Scale_Factor', 'long_name', 'Continuous/Screened Mosaic Esoil');

%create Noah data
nccreate(outfile,'Noah_Scale_Factor','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
ncwriteatt(outfile, 'Noah_Scale_Factor', 'long_name', 'Continuous/Screened Noah Esoil');

%create GLEAM data
nccreate(outfile,'GLEAM_Scale_Factor','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
ncwriteatt(outfile, 'GLEAM_Scale_Factor', 'long_name', 'Continuous/Screened GLEAM Esoil');

%create Ensemnle Mean data
nccreate(outfile,'Scale_Factor','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
ncwriteatt(outfile, 'Scale_Factor', 'long_name', 'Mean scale factor from Noah, Mosaic and GLEAM derivations');

%write variables to file
ncwrite(outfile,'lat',lat);
ncwrite(outfile,'lon',lon);
ncwrite(outfile,'Mosaic_Scale_Factor',Mosaic_Scale_Factor);
ncwrite(outfile,'Noah_Scale_Factor',Noah_Scale_Factor);
ncwrite(outfile,'GLEAM_Scale_Factor',GLEAM_Scale_Factor);
ncwrite(outfile,'Scale_Factor',SCALE_FACTOR_Ensemble_Mean);
