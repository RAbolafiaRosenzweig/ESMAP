clc;clear all;close all;

SMAP_dir='/Volumes/LabShare/Global/SMAP/';
SMAP_QC_dir='/Volumes/LabShare/Global/SMAP/QC/';

days=[31,28,31,30,31,30,31,31,30,31,30,31];
for year=2015:2019
    if year==2015
        startmonth=4;
    else
        startmonth=1;
    end
    
    if year==2019
        endmonth=3;
    else
        endmonth=12;
    end
    
    for month=startmonth:endmonth
        daysinmonth=days(month);
        if month==2 && year==2016
            daysinmonth=29;
        end 
        
        for day=1:daysinmonth
            [year month day]
            
            filename=sprintf('SMAP.%04d.%02d.%02d.nc',year,month,day);
            if exist([SMAP_dir,filename],'file')>0
            qual_flag=ncread([SMAP_dir,filename],'qual_flag');
            SMAP=ncread([SMAP_dir,filename],'SMAP');
            lat=ncread([SMAP_dir,filename],'lat');
            lon=ncread([SMAP_dir,filename],'lon');
            surf_flag=ncread([SMAP_dir,filename],'surf_flag');
            time=ncread([SMAP_dir,filename],'time');
            
            %only keep qual flag=0 or 8 based on email from Andrew
            %Badger(02/25/2020)
            idx_bad=find(qual_flag~=0 & qual_flag~=8);
            SMAP(idx_bad)=-9999;
            
            
            
            %write to new netcdf file
            outfile=[SMAP_QC_dir,filename];
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
            ncwrite(outfile,'SMAP',SMAP);
            ncwrite(outfile,'surf_flag',surf_flag);
            ncwrite(outfile,'qual_flag',qual_flag);
            end
        end
    end
end
