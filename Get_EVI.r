rm(list=ls())

library(MASS)
library(ncdf4)
library(fields)
library(maps)
library(hydroGOF)
#library(plotrix)
library(rworldmap)
library(RColorBrewer)
library(humidity)

esmap.dir = "/Volumes/LabShare/SMAP/Gridded_ESMAP/"
evi.dir = "/Volumes/LabShare/Global/MOD13A2/EVI_dir/"

smap.datelist = seq.Date(as.Date("2015/04/01"),as.Date("2019/3/31"),"days") #Change when extending SMAP
months = as.numeric(format(smap.datelist,'%m'))
years = as.numeric(format(smap.datelist,'%Y'))
days = as.numeric(format(smap.datelist,'%d'))

#Open file to get dimensions
tmp.file = paste(evi.dir,"EVI.2015.04.07.nc", sep='')
p.nc = nc_open(tmp.file)
evi.lon = ncvar_get(p.nc,'lon')
evi.lat = ncvar_get(p.nc,'lat')
nc_close(p.nc)

south = 25
north = 50
east = -60
west = -130

chosen.lat = which(evi.lat >= south & evi.lat <= north)
chosen.lon = which(evi.lon >= west & evi.lon <= east)

evi.lon = evi.lon[chosen.lon]
evi.lat = evi.lat[chosen.lat]

nlon = length(evi.lon)
nlat = length(evi.lat)
ntime = length(smap.datelist)

evi.data = array(NA, c(nlon,nlat,ntime))

for(i in 1:ntime){
	tmp.file = paste(evi.dir,"EVI.",years[i],".",sprintf("%02d", months[i]),".",sprintf("%02d", days[i]),".nc", sep='')
	if(file.exists(tmp.file) == TRUE){
	
		p.nc = nc_open(tmp.file)
		tmp.data = ncvar_get(p.nc,'Band1')
		nc_close(p.nc)
		
		tmp.data = tmp.data[chosen.lon, chosen.lat]
		
		evi.data[,,i] = tmp.data
		print(tmp.file)
		}
	}

evi.lon = evi.lon + 360
smap.lat = list.files(esmap.dir)
smap.lat = as.numeric(smap.lat)
nlat = length(smap.lat)
for(i in 1:nlat){
	smap.lon = list.files(paste(esmap.dir,smap.lat[i],"/",sep=''))
	lat.loc = which.min(abs(evi.lat - smap.lat[i]))
	nlon = length(smap.lon)
	smap.lon = as.numeric(smap.lon)
	for(j in 1:nlon){
		lon.loc = which.min(abs(evi.lon - smap.lon[j]))
		tmp.evi = evi.data[lon.loc,lat.loc,]
		if(min(tmp.evi,na.rm=TRUE)<0){
			tmp.evi = evi.data[(lon.loc+1),lat.loc,]
			}
		if(min(tmp.evi,na.rm=TRUE)<0){
			tmp.evi = evi.data[(lon.loc-1),lat.loc,]
			}	
		to.write = cbind(years, months, days, tmp.evi)
		new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/EVI_raw",sep='')
		write.table(to.write,new.file,sep="\t", col.names = F, row.names = F)
		}

	}
	