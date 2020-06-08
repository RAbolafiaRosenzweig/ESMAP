clc;clear all;close all;

start_date=[2015 4 1];
end_date=[2019 3 31];

datelist=datenum(start_date):datenum(end_date);

%Define path to soil moisture, lat and lon (SMAP)
DATAFIELD_NAME= 'Soil_Moisture_Retrieval_Data_AM/soil_moisture';
DATAFIELD_NAME_qual= 'Soil_Moisture_Retrieval_Data_AM/retrieval_qual_flag';
DATAFIELD_NAME_surf= 'Soil_Moisture_Retrieval_Data_AM/surface_flag';
Lat_NAME= 'Soil_Moisture_Retrieval_Data_AM/latitude';
Lon_NAME = 'Soil_Moisture_Retrieval_Data_AM/longitude';

for i=1:length(datelist)
    %define the current date:
    current_date=datevec(datelist(i))
    year=current_date(1);
    month=current_date(2);
    day=current_date(3);
    %define the directory of the SMAP file in hdf format:
    FILE_PATH_SMAP = sprintf('/Volumes/LabShare/SMAP/data_raw/n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP_E.002/%04d.%02d.%02d',year,month,day);
    SMAP_FILE=dir([FILE_PATH_SMAP '/*.h5']);
    if length(SMAP_FILE)>0
    FILE_NAME_SMAP=SMAP_FILE.name;
    %define file ID
    file_id = H5F.open ([FILE_PATH_SMAP '/' FILE_NAME_SMAP], 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

    % Open the dataset.
    data_id = H5D.open (file_id, DATAFIELD_NAME);
    qual_id = H5D.open (file_id, DATAFIELD_NAME_qual);
    surf_id = H5D.open (file_id, DATAFIELD_NAME_surf);
    lat_id=H5D.open(file_id, Lat_NAME);
    lon_id=H5D.open(file_id, Lon_NAME);
    
    % Read the dataset.
    data=H5D.read (data_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    latitude=H5D.read(lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    %turn latitude into a row vector instead of grid:
    lat=max(latitude); %use max to avoid -999 fill values

    longitude=H5D.read(lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    %turn longitude into a row vector instead of grid:
    lon=max(longitude');
    
    qual_flag=H5D.read(qual_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    surf_flag=H5D.read(surf_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    %define time as days since 3/31/2015:
    time=datenum(current_date)-datenum(start_date)-1;
    %write the file to netcdf:
    outfile=sprintf('/Volumes/LabShare/Global/SMAP/SMAP.%04d.%02d.%02d.nc',year,month,day);
    %if file already exists remove the old version:
    if exist(outfile,'file')==2
        CMD=['rm ',outfile];
        system(CMD);
    end
    %create variables and file
    
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
    %create time:
    nccreate(outfile,'time','datatype','single','Dimensions',{'time' 1});
    ncwriteatt(outfile, 'time', 'long_name', 'Time variable');
    ncwriteatt(outfile, 'time', 'units', 'days since 2015-3-31');
    ncwriteatt(outfile, 'time', '_CoordinateAxisType', 'Time');
    %create SMAP data
    nccreate(outfile,'SMAP','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
    ncwriteatt(outfile, 'SMAP', 'units', 'cm3/cm3');
    ncwriteatt(outfile, 'SMAP', 'missing_value', -9999.0);
    %create quality and surface flags
    nccreate(outfile,'surf_flag','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
    nccreate(outfile,'qual_flag','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});

    %write variables to file
    ncwrite(outfile,'lat',lat);
    ncwrite(outfile,'lon',lon);
    ncwrite(outfile,'SMAP',data);
    ncwrite(outfile,'surf_flag',surf_flag);
    ncwrite(outfile,'qual_flag',qual_flag);
    
%     ncdisp(outfile);
    end
end
    
    