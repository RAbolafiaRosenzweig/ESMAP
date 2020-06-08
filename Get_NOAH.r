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
nldas.dir = "/Volumes/LabShare/NLDAS2/NOAH/"

smap.datelist = seq.Date(as.Date("2015/04/01"),as.Date("2019/4/1"),"days") #Change when extending SMAP, do + 1 day for time conversion
months = as.numeric(format(smap.datelist,'%m'))
years = as.numeric(format(smap.datelist,'%Y'))
days = as.numeric(format(smap.datelist,'%d'))
nsmap = length(days)

#Open file to get dimensions
tmp.file = paste(nldas.dir,"2015/04/NLDAS_NOAH0125_H..A20150401.0000.nc", sep='')
p.nc = nc_open(tmp.file)
nldas.lon = ncvar_get(p.nc,'lon')
nldas.lat = ncvar_get(p.nc,'lat')
nc_close(p.nc)

nnlon = length(nldas.lon)
nnlat = length(nldas.lat)

nldas.lon = nldas.lon + 360

ntime = length(smap.datelist)*24

smap.lat = list.files(esmap.dir)
smap.lat = as.numeric(smap.lat)
nlat = length(smap.lat)
			
for(i in 1:nsmap){

	sw.data = array(NA, c(nnlon,nnlat,24))
	lw.data = array(NA, c(nnlon,nnlat,24))
	years.data = array(years[i], c(24))
	months.data = array(months[i], c(24))
	days.data = array(days[i], c(24))
	hours.data = array(NA, c(24))
	
	for(j in 0:23){
	
		hours.data[(j+1)] = j
		
		tmp.file = paste(nldas.dir,years[i],"/",sprintf("%02d", months[i]),"/NLDAS_NOAH0125_H..A",years[i],sprintf("%02d", months[i]),sprintf("%02d", days[i]),".",sprintf("%02d", j),"00.nc", sep='')
		if(file.exists(tmp.file) == TRUE){
	
			p.nc = nc_open(tmp.file)
			tmp.sw = ncvar_get(p.nc,'var111_NSWRS')
			tmp.lw = ncvar_get(p.nc,'var112_NLWRS')
			nc_close(p.nc)
			
			sw.data[,,(j+1)] = tmp.sw
			lw.data[,,(j+1)] = tmp.lw
			}
		}
			
	for(k in 1:nlat){
		smap.lon = list.files(paste(esmap.dir,smap.lat[k],"/",sep=''))
		lat.loc = which.min(abs(nldas.lat - smap.lat[k]))
		nlon = length(smap.lon)
		smap.lon = as.numeric(smap.lon)
				
		for(l in 1:nlon){
			lon.loc = which.min(abs(nldas.lon - smap.lon[l]))
				
			to.write.sw = sw.data[lon.loc,lat.loc,]
			to.write.lw = lw.data[lon.loc,lat.loc,]
					
			if(is.na(to.write.lw[1]) == TRUE){
				lon.loc = which(is.na(tmp.lw[,lat.loc]) == FALSE)[which.min(abs(lon.loc - which(is.na(tmp.lw[,lat.loc]) == FALSE)))]
				to.write.sw = sw.data[lon.loc,lat.loc,]
				to.write.lw = lw.data[lon.loc,lat.loc,]
				}
							
			new.file = paste(esmap.dir,smap.lat[k],"/",smap.lon[l],"/NLDAS_NOAH_raw",sep='')
			if(i == 1) file.remove(new.file)
			to.write = cbind(years.data, months.data, days.data, hours.data, to.write.sw, to.write.lw)
			write.table(to.write,new.file,append=TRUE,sep="\t", col.names = F, row.names = F)
					
			to.screen = c(i, nsmap, k, nlat, l, nlon)
			print(to.screen)
			}

		}	
	print(tmp.file)
	}
		
	
