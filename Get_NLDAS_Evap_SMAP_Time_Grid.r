rm(list=ls())

library(MASS)
library(ncdf)
library(fields)
library(maps)
library(hydroGOF)
library(plotrix)
library(rworldmap)
library(RColorBrewer)
library(humidity)

nldas.dir = "/Volumes/LabShare/NLDAS2/NOAH_DAILY/"
esmap.dir = "/Volumes/LabShare/SMAP/Gridded_ESMAP/"
smap.points.file = "/Users/andrewbadger/Documents/SMAP/ESMAP_Calc/SMAP_36km_Points"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]
npoints = length(smap.lon)

smap.file = paste(esmap.dir,smap.lat[1],"/",smap.lon[1],"/SMAP_raw",sep='')
smap.raw = read.table(smap.file)
smap.raw = as.matrix(smap.raw)
smap.year = smap.raw[,1]
smap.month = smap.raw[,2]
smap.day = smap.raw[,3]
ntime = length(smap.year)

tmp.file = "/Volumes/LabShare/NLDAS2/NOAH_DAILY/2015/05/NLDAS_NOAH0125_H.A20150516.nc"
p.nc = open.ncdf(tmp.file)
nldas.lon = get.var.ncdf(p.nc,'lon') + 360
nldas.lat = get.var.ncdf(p.nc,'lat')
close.ncdf(p.nc)

nlon = length(nldas.lon)
nlat = length(nldas.lat)

nldas.data = array(NA, c(nlon,nlat,ntime))

for(i in 1:ntime){
	tmp.file = paste(nldas.dir,smap.year[i],"/",sprintf("%02d", smap.month[i]),"/NLDAS_NOAH0125_H.A",smap.year[i],sprintf("%02d", smap.month[i]),sprintf("%02d",smap.day[i]),".nc", sep='')
	if(file.exists(tmp.file) == TRUE){
		p.nc = open.ncdf(tmp.file)
		tmp.e = get.var.ncdf(p.nc,'var57_EVP')
		close.ncdf(p.nc)
			
		nldas.data[,,i] = 24 * tmp.e/10
		}
	}

tmp.file = paste("/Users/andrewbadger/Documents/SMAP/ESMAP_Calc/SMAP_36km/SMAP_36km.nc",sep='')
p.nc = open.ncdf(tmp.file)
new.lon = get.var.ncdf(p.nc,'lon')
new.lat = get.var.ncdf(p.nc,'lat')
close.ncdf(p.nc)

south = 25
north = 50
east = -67
west = -125

south.loc = which.min(abs(new.lat - south))
north.loc = which.min(abs(new.lat - north))
east.loc = which.min(abs(new.lon - east))
west.loc = which.min(abs(new.lon - west))
new.lon = new.lon[west.loc:east.loc] + 360
new.lat = new.lat[south.loc:north.loc]

new.nlon = length(new.lon)
new.nlat = length(new.lat)

new.evap = array(NA, c(new.nlon,new.nlat,ntime))
for(p in 1:npoints){
	smap.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/SMAP_raw",sep='')
	
	smap.raw = read.table(smap.file)
	smap.raw = as.matrix(smap.raw)
	smap.year = smap.raw[,1]
	smap.month = smap.raw[,2]
	smap.day = smap.raw[,3]
	smap.val = smap.raw[,4]
	smap.valid = which(is.na(smap.val)==FALSE)
	nint = length(smap.valid) - 1
	calc.time = array(NA, c(nint))
	calc.esoil = array(NA, c(nint))
	
	lat.loc = which.min(abs(nldas.lat - smap.lat[p]))
	lon.loc = which.min(abs(nldas.lon - smap.lon[p]))
	nldas.e = nldas.data[lon.loc,lat.loc,]
	for(i in 1:nint){
		calc.time[i] = floor((smap.valid[i] + smap.valid[i+1])/2)
		calc.esoil[i] = sum(nldas.e[smap.valid[i]:(smap.valid[i+1]-1)])
		}
	lat.loc = which.min(abs(new.lat - smap.lat[p]))
	lon.loc = which.min(abs(new.lon - smap.lon[p]))
	new.evap[lon.loc,lat.loc,calc.time] = calc.esoil
	}

filename = paste(esmap.dir,"/NLDAS_ESOIL.nc",sep='')
var = paste("ESOIL",sep='')
units = paste("cm",sep='')
time = seq(1:ntime)
nc_filename = paste(filename,sep='')
nc_var = var
nc_var_units = units
nc_dim1 = "time"
nc_dim2 = "lat"
nc_dim3 = "lon"
nc_dim1_units = "time"
nc_dim2_units = "degrees_north"
nc_dim3_units = "degrees_east"
nc_vartype = "double"
nc_missing = -9999
# Define some straightforward dimensions
y = dim.def.ncdf( nc_dim3, nc_dim3_units, new.lon)
x = dim.def.ncdf( nc_dim2, nc_dim2_units, new.lat)
t = dim.def.ncdf( nc_dim1, nc_dim1_units, time, unlim=TRUE)
# Make a variable with those dimensions.
variable = var.def.ncdf(nc_var, nc_var_units, list(y,x,t), nc_missing, prec = nc_vartype)
# Create a netCDF file with this variable
ncnew <- create.ncdf( nc_filename, variable)
# Write some values to this variable on disk.
put.var.ncdf(ncnew, variable, new.evap,)
close.ncdf(ncnew)

	
		