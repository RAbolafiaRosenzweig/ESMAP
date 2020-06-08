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
soil.dir = "/Volumes/LabShare/SMAP/ESMAP_for_Ronnie/"

#Open file to get dimensions
soil.file = paste(soil.dir,"NLDAS_Soil_Class.txt", sep='')
soil.raw = read.table(soil.file, fill=TRUE)
soil.lon = soil.raw[,3] + 360
soil.lat = soil.raw[,4]
soil.class = soil.raw[,5:20]

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
	nlon = length(smap.lon)
	smap.lon = as.numeric(smap.lon)
	for(j in 1:nlon){
		
		lat.final = which(abs(soil.lat - smap.lat[i]) == min(abs(soil.lat - smap.lat[i])))
		lon.final = which(abs(soil.lon - smap.lon[j]) == min(abs(soil.lon - smap.lon[j])))
		tmp.soil = soil.class[intersect(lat.final,lon.final),]
		if(dim(tmp.soil)[1] >= 2) tmp.soil = colMeans(tmp.soil)
		final.soil = which.max(tmp.soil)
		if(final.soil== 14){
			tmp.soil[14] = NA
			final.soil = which.max(tmp.soil)
			}	
		new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/Soil_raw",sep='')
		write.table(final.soil,new.file,sep="\t", col.names = F, row.names = F)
		
		#all.soil = which(tmp.soil >= 1)
		#new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/Soil_All_raw",sep='')
		#write.table(all.soil,new.file,sep="\t", col.names = F, row.names = F)
		
		to.screen = c(smap.lat[i], smap.lon[j], final.soil)
		print(to.screen)
		}

	}