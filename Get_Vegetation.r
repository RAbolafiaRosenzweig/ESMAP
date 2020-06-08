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
veg.dir = "/Volumes/LabShare/SMAP/ESMAP_for_Ronnie/"

#Open file to get dimensions
tmp.file = paste(veg.dir,"AVHRR_8km_LANDCOVER_1981_1994.GLOBAL.nc", sep='')
p.nc = nc_open(tmp.file)
veg.lon = ncvar_get(p.nc,'lon')
veg.lat = ncvar_get(p.nc,'lat')
veg.data = ncvar_get(p.nc,'Band1')
nc_close(p.nc)

veg.lon = veg.lon + 360

smap.lat = list.files(esmap.dir)
smap.lat = as.numeric(smap.lat)
#remove NAs: (other files in directory that are not coordinates)
idx<-which(is.na(smap.lat))
if (length(idx)>0) {
  smap.lat<-smap.lat[-idx]
}
nlat = length(smap.lat)

for(i in 1:nlat){
	smap.lon = list.files(paste(esmap.dir,smap.lat[i],"/",sep=''))
	lat.loc = which.min(abs(veg.lat - smap.lat[i]))
	nlon = length(smap.lon)
	smap.lon = as.numeric(smap.lon)
	for(j in 1:nlon){
		lon.loc = which.min(abs(veg.lon - smap.lon[j]))
		tmp.veg = veg.data[lon.loc,lat.loc]
		if(tmp.veg == 0 | tmp.veg == 13){
			tmp.data = veg.data[,lat.loc]
			tmp.data[which(tmp.data == 0)] = NA
			tmp.data[which(tmp.data == 13)] = NA
			lon.loc = which(is.na(tmp.data) == FALSE)[which.min(abs(lon.loc - which(is.na(tmp.data) == FALSE)))]
			tmp.veg = veg.data[lon.loc,lat.loc]
			}	
		new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/Vegetation_raw",sep='')
		write.table(tmp.veg,new.file,sep="\t", col.names = F, row.names = F)
		
		to.screen = c(smap.lat[i], smap.lon[j], tmp.veg)
		print(to.screen)
		}

	}
