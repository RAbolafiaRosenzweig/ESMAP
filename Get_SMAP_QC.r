rm(list=ls())

library(MASS)
library(ncdf4)
#library(fields)
#library(maps)
#library(hydroGOF)
#library(plotrix)
#library(rworldmap)
#library(RColorBrewer)
#library(humidity)

smap.dir = "/scratch/summit/roab9675/SMAP/QC/"
smap.datelist = seq.Date(as.Date("2015/04/01"),as.Date("2019/3/31"),"days")
months = as.numeric(format(smap.datelist,'%m'))
years = as.numeric(format(smap.datelist,'%Y'))
days = as.numeric(format(smap.datelist,'%d'))

ndays = length(days)
print(ndays)

tmp.file = paste(smap.dir,"SMAP.2015.04.01.nc",sep='')
p.nc = nc_open(tmp.file)
lon = ncvar_get(p.nc,'lon')
lat = ncvar_get(p.nc,'lat')
nc_close(p.nc)

nlon = length(lon)
nlat = length(lat)

sm.data = array(NA, c(nlon, nlat, ndays))


for(i in 1:ndays){
	tmp.file = paste(smap.dir,"SMAP.",years[i],".",sprintf("%02d", months[i]),".",sprintf("%02d", days[i]),".nc",sep='')
	#only perform if the tmp.file existis:
	if (file.exists(tmp.file)){
	p.nc = nc_open(tmp.file)
	tmp.data = ncvar_get(p.nc,'SMAP')
	nc_close(p.nc)
	
	sm.data[,,i] = tmp.data
	print(tmp.file)
	}
	}

south = 25
north = 50
east = -60
west = -130

#nldas 25 to 53N, -67 to -125W

south.loc = which.min(abs(lat - south))
north.loc = which.min(abs(lat - north))
east.loc = which.min(abs(lon - east))
west.loc = which.min(abs(lon - west))

sm.data = sm.data[west.loc:east.loc,south.loc:north.loc,]
lon = lon[west.loc:east.loc] + 360
lat = lat[south.loc:north.loc]
nlon = length(lon)
nlat = length(lat)

esmap.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP/"
for(i in 1:nlat){
	for(j in 1:nlon){
		tmp.sm = sm.data[j,i,]
		nvalid = length(which(is.na(tmp.sm)==FALSE))
		if(nvalid >= 5){
				print(lat[i])
				print(lon[j])
				new.dir = paste(esmap.dir,round(lat[i],5),"/",round(lon[j],5),"/",sep='')
				dir.create(new.dir,recursive=TRUE)
				new.file = paste(new.dir,"SMAP_raw_QC",sep='')
				to.write = cbind(years, months, days, tmp.sm)
				write.table(to.write,new.file,sep="\t", col.names = F, row.names = F)
			}
		}
	}

