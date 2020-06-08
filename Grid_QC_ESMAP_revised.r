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

esmap.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP/"
smap.qc.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP_blank/"
smap.points.file = "/projects/roab9675/SMAP/Point_Check/ESMAP_Calc/Esoil_QC_has_data"
xtl.dir = "/scratch/summit/roab9675/SMAP/ESMAP_Calc/"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]
npoints = length(smap.lon)

smap.file = paste(smap.qc.dir,smap.lat[1],"/",smap.lon[1],"/SMAP_raw_QC",sep='')
smap.raw = read.table(smap.file)
smap.raw = as.matrix(smap.raw)
smap.year = smap.raw[,1]
smap.month = smap.raw[,2]
smap.day = smap.raw[,3]
ntime = length(smap.year)

tmp.file = paste("/scratch/summit/roab9675/SMAP/ESMAP_Calc/SMAP_9km/SMAP_9km.nc",sep='')
p.nc = nc_open(tmp.file)
new.lon = ncvar_get(p.nc,'lon')
new.lat = ncvar_get(p.nc,'lat')
nc_close(p.nc)

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
xres = mean(new.lon[2:new.nlon] - new.lon[1:(new.nlon-1)])
yres = mean(new.lat[2:new.nlat] - new.lat[1:(new.nlat-1)])

smap.time = array(NA, c(new.nlon,new.nlat,ntime))
smap.length = array(NA, c(new.nlon,new.nlat,ntime))
smap.evi = array(NA, c(new.nlon,new.nlat,ntime))
smap.pre = array(NA, c(new.nlon,new.nlat,ntime))
smap.del = array(NA, c(new.nlon,new.nlat,ntime))
smap.trans = array(NA, c(new.nlon,new.nlat,ntime))
smap.qbot = array(NA, c(new.nlon,new.nlat,ntime))
smap.esoil = array(NA, c(new.nlon,new.nlat,ntime))

for(p in 1:npoints){
	hydrus.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/ESMAP_QC_DATA_new",sep='')
	
	hydrus.raw = read.table(hydrus.file)
	hydrus.raw = as.matrix(hydrus.raw)
	calc.time = hydrus.raw[,1]
	
	hydrus.time = hydrus.raw[,1]
	hydrus.length = hydrus.raw[,2]
	hydrus.evi = hydrus.raw[,3]
	hydrus.pre = hydrus.raw[,4] #cm
	hydrus.pre = hydrus.pre*10 #mm
	hydrus.del = hydrus.raw[,5]*10 #mm
	hydrus.trans = hydrus.raw[,6]*10 #mm
	hydrus.qbot = hydrus.raw[,7]*10 #mm
	hydrus.esoil = hydrus.raw[,8]*10 #mm
	old.hydrus.esoil = hydrus.esoil
	hydrus.esoil = 0 - hydrus.del - hydrus.trans + hydrus.qbot
	hydrus.esoil = hydrus.esoil/hydrus.length #mm/day
	hydrus.del = hydrus.del/hydrus.length #mm/day
	hydrus.trans = hydrus.trans/hydrus.length #mm/day
	hydrus.qbot = hydrus.qbot/hydrus.length #mm/day
	
	lat.loc = which.min(abs(new.lat - smap.lat[p]))
	lon.loc = which.min(abs(new.lon - smap.lon[p]))
	
	smap.time[lon.loc,lat.loc,calc.time] = hydrus.time
	smap.length[lon.loc,lat.loc,calc.time] = hydrus.length
	smap.evi[lon.loc,lat.loc,calc.time] = hydrus.evi
	smap.pre[lon.loc,lat.loc,calc.time] = hydrus.pre
	smap.del[lon.loc,lat.loc,calc.time] = hydrus.del
	smap.trans[lon.loc,lat.loc,calc.time] = hydrus.trans
	smap.qbot[lon.loc,lat.loc,calc.time] = hydrus.qbot
	smap.esoil[lon.loc,lat.loc,calc.time] = hydrus.esoil
	print(paste(p," of ",npoints," have been gridded!",sep=''))
	}

smap.esoil.pos = smap.esoil
smap.esoil.screened = smap.esoil
tmp = which(smap.esoil < 0)
smap.esoil.pos[tmp] = NA
smap.esoil.screened[tmp] = NA
tmp = which(smap.pre > 0.2)	
smap.esoil.screened[tmp] = NA

for(i in 1:10){

	if(i == 1) filename = paste(esmap.dir,"ESMAP_QC_TIME_revised.nc",sep='')
	if(i == 2) filename = paste(esmap.dir,"ESMAP_QC_LENGTH_revised.nc",sep='')
	if(i == 3) filename = paste(esmap.dir,"ESMAP_QC_EVI_revised.nc",sep='')
	if(i == 4) filename = paste(esmap.dir,"ESMAP_QC_PRE_revised.nc",sep='')
	if(i == 5) filename = paste(esmap.dir,"ESMAP_QC_DELSMAP_revised.nc",sep='')
	if(i == 6) filename = paste(esmap.dir,"ESMAP_QC_TRANS_revised.nc",sep='')
	if(i == 7) filename = paste(esmap.dir,"ESMAP_QC_QBOT_revised.nc",sep='')
	if(i == 8) filename = paste(esmap.dir,"ESMAP_QC_ESOIL_RAW_revised.nc",sep='')
	if(i == 9) filename = paste(esmap.dir,"ESMAP_QC_ESOIL_noNEG_revised.nc",sep='')
	if(i == 10) filename = paste(esmap.dir,"ESMAP_QC_ESOIL_SCREENED_revised.nc",sep='')
	
	if(i == 1) xfilename = paste(xtl.dir,"ESMAP_QC_TIME_revised.xtl",sep='')
	if(i == 2) xfilename = paste(xtl.dir,"ESMAP_QC_LENGTH_revised.xtl",sep='')
	if(i == 3) xfilename = paste(xtl.dir,"ESMAP_QC_EVI_revised.xtl",sep='')
	if(i == 4) xfilename = paste(xtl.dir,"ESMAP_QC_PRE_revised.xtl",sep='')
	if(i == 5) xfilename = paste(xtl.dir,"ESMAP_QC_DELSMAP_revised.xtl",sep='')
	if(i == 6) xfilename = paste(xtl.dir,"ESMAP_QC_TRANS_revised.xtl",sep='')
	if(i == 7) xfilename = paste(xtl.dir,"ESMAP_QC_QBOT_revised.xtl",sep='')
	if(i == 8) xfilename = paste(xtl.dir,"ESMAP_QC_ESOIL_RAW_revised.xtl",sep='')
	if(i == 9) xfilename = paste(xtl.dir,"ESMAP_QC_ESOIL_noNEG_revised.xtl",sep='')
	if(i == 10) xfilename = paste(xtl.dir,"ESMAP_QC_ESOIL_SCREENED_revised.xtl",sep='')

	if(i == 1) var = paste("stime",sep='')
	if(i == 2) var = paste("slength",sep='')
	if(i == 3) var = paste("evi",sep='')
	if(i == 4) var = paste("pre",sep='')
	if(i == 5) var = paste("delsmap",sep='')
	if(i == 6) var = paste("trans",sep='')
	if(i == 7) var = paste("qbot",sep='')
	if(i == 8) var = paste("esoil",sep='')
	if(i == 9) var = paste("esoil",sep='')
	if(i == 10) var = paste("esoil",sep='')
	
	if(i == 1) units = paste("days",sep='')
	if(i == 2) units = paste("days",sep='')
	if(i == 3) units = paste("unitless",sep='')
	if(i == 4) units = paste("mm/day",sep='')
	if(i == 5) units = paste("mm/day",sep='')
	if(i == 6) units = paste("mm/day",sep='')
	if(i == 7) units = paste("mm/day",sep='')
	if(i == 8) units = paste("mm/day",sep='')
	if(i == 9) units = paste("mm/day",sep='')
	if(i == 10) units = paste("mm/day",sep='')
	
	if(i == 1) nc.chosen = smap.time
	if(i == 2) nc.chosen = smap.length
	if(i == 3) nc.chosen = smap.evi
	if(i == 4) nc.chosen = smap.pre
	if(i == 5) nc.chosen = smap.del
	if(i == 6) nc.chosen = smap.trans
	if(i == 7) nc.chosen = smap.qbot
	if(i == 8) nc.chosen = smap.esoil
	if(i == 9) nc.chosen = smap.esoil.pos
	if(i == 10) nc.chosen = smap.esoil.screened
	
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
	y = ncdim_def( nc_dim3, nc_dim3_units, new.lon)
	x = ncdim_def( nc_dim2, nc_dim2_units, new.lat)
	t = ncdim_def( nc_dim1, nc_dim1_units, time, unlim=TRUE)
	# Make a variable with those dimensions.
	variable = ncvar_def(nc_var, nc_var_units, list(y,x,t), nc_missing, prec = nc_vartype)
	# Create a netCDF file with this variable
	ncnew <- nc_create( nc_filename, variable)
	# Write some values to this variable on disk.
	ncvar_put(ncnew, variable, nc.chosen,)
	nc_close(ncnew)
	
	xtl.write = array(NA, c(10,1))
	xtl.write[1,1] = paste("DSET ",nc_filename,sep='')
	xtl.write[2,1] = "DTYPE netcdf"
	xtl.write[3,1] = "UNDEF -9999 missing_value"
	xtl.write[4,1] = paste("XDEF ",new.nlon," LINEAR ",new.lon[1]+360," ",xres,sep='')
	xtl.write[5,1] = paste("YDEF ",new.nlat," LINEAR ",new.lat[1]," ",yres,sep='')
	xtl.write[6,1] = paste("TDEF ",ntime,"  LINEAR 00z1apr2015 1dy",sep='')
	xtl.write[7,1] = "ZDEF 1 LEVELS 1"
	xtl.write[8,1] = "vars 1"
	xtl.write[9,1] = paste(nc_var,"=>",nc_var," 1 t,y,x ",nc_var_units,sep='')
	xtl.write[10,1] = "endvars"
	xtl.name = xfilename
	write(xtl.write, xtl.name)
	}


