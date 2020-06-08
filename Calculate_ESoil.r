rm(list=ls())

library(MASS)
#library(ncdf4)
#library(fields)
#library(maps)
#library(hydroGOF)
#library(plotrix)
#library(rworldmap)
#library(RColorBrewer)
#library(humidity)

esmap.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP/"
esmap.qc.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP_blank/"
smap.points.file = "/projects/roab9675/SMAP/Point_Check/ESMAP_Calc/Esoil_has_data"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]
npoints = length(smap.lon)

for(p in 1:npoints){
	smap.file = paste(esmap.qc.dir,smap.lat[p],"/",smap.lon[p],"/SMAP_raw_QC",sep='')
	q.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/OBS_NODE.OUT",sep='')
	for.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/ATMOSPH.IN",sep='')
	trans.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Transpiration_Estimates",sep='')
	hydrus.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/A_LEVEL.OUT",sep='')
	t.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/T_LEVEL.OUT",sep='')
	ks.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/PROFILE.OUT",sep='')
	
	smap.raw = read.table(smap.file)
	smap.raw = as.matrix(smap.raw)
	smap.year = smap.raw[,1]
	smap.month = smap.raw[,2]
	smap.day = smap.raw[,3]
	smap.val = smap.raw[,4]
	
	trans.raw = read.table(trans.file)
	trans.raw = as.matrix(trans.raw)
	trans.val = trans.raw[,7]
	
	lhead = grep('h        theta    Temp',readLines(q.file))
	lend = grep('end',readLines(q.file))
	q.raw = readLines(q.file)
	q.raw = q.raw[(lhead+1):(lend-1)]
	nq = length(q.raw)
	nq.sel = seq(1,nq,by=1)
	nq = length(nq.sel)
	q.time = array(NA, c(nq))
	q.time.start = array(NA, c(nq))
	q.time.start[1] = 0
	q.bot.raw = array(NA, c(nq))
	for(i in 1:nq){
		tmp = q.raw[nq.sel[i]]
		q.time[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[1]])
		q.time.start[(i+1)] = q.time[i]
		q.bot.raw[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[19]])	
		}
	time.diff = q.time[1:nq] - q.time.start[1:nq]
	q.bot.tot = time.diff*q.bot.raw
		
	lhead = grep('tAtm',readLines(for.file))
	lend = grep('end',readLines(for.file))
	for.raw = readLines(for.file)
	for.raw = for.raw[(lhead+1):(lend-1)]
	nfor = length(for.raw)
	for.time = array(NA, c(nfor))
	for.pre = array(NA, c(nfor))
	for.evi = array(NA, c(nfor))
	for(i in 1:nfor){
		tmp = for.raw[i]
		for.time[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[1]])
		for.pre[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[2]])
		for.evi[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[4]])
		}
	
	lhead = grep('Time',readLines(hydrus.file))
	lend = grep('end',readLines(hydrus.file))
	hydrus.raw = readLines(hydrus.file)
	hydrus.raw = hydrus.raw[(lhead+3):(lend-1)]
	nhydrus = length(hydrus.raw)
	hydrus.time = array(NA, c(nhydrus))
	hydrus.evap = array(NA, c(nhydrus))
	for(i in 1:nfor){
		tmp = hydrus.raw[i]
		hydrus.time[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[1]])
		hydrus.evap[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[4]])
		}
		
	lhead = max(grep('Time',readLines(t.file)))
	lend = grep('end',readLines(t.file))
	t.raw = readLines(t.file)
	t.raw = t.raw[(lhead+3):(lend-1)]
	nt = length(t.raw)
	t.time = array(NA, c(nt))
	t.vol = array(NA, c(nt))
	for(i in 1:nt){
		tmp = t.raw[i]
		t.time[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[1]])
		t.vol[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[17]])
		}
		
	lhead = max(grep('depth',readLines(ks.file)))
	ks.raw = readLines(ks.file)
	ks.raw = ks.raw[(lhead+3)]
	ks.val = as.numeric(strsplit(ks.raw," ")[[1]][which(strsplit(ks.raw," ")[[1]] != "")[6]])
	
	smap.valid = which(is.na(smap.val)==FALSE)
	smap.times = which(is.na(smap.val)==FALSE)*24-18
	nint = length(smap.times) - 1
	calc.time = array(NA, c(nint))
	calc.length = array(NA, c(nint))
	calc.qbot = array(NA, c(nint))
	calc.pre = array(NA, c(nint))
	calc.delsmap = array(NA, c(nint))
	calc.trans = array(NA, c(nint))
	calc.evi = array(NA, c(nint))
	calc.esoil = array(NA, c(nint))
	calc.hydruse = array(NA, c(nint))
	for(i in 1:nint){
		calc.time[i] = floor((smap.valid[i] + smap.valid[i+1])/2)
		calc.length[i] = smap.valid[i+1] - smap.valid[i]
		calc.pre[i] = sum(for.pre[which(for.time >= smap.times[i] & for.time < smap.times[i+1])]) 
		calc.delsmap[i] = (smap.val[smap.valid[i+1]] - smap.val[smap.valid[i]])*5
		calc.trans[i] = sum(trans.val[smap.valid[i]:(smap.valid[i+1]-1)])/10
		calc.evi[i] = mean(for.evi[which(for.time >= smap.times[i] & for.time < smap.times[i+1])]) 
		calc.hydruse[i] = sum(hydrus.evap[which(hydrus.time >= smap.times[i] & hydrus.time < smap.times[i+1])])
		calc.qbot[i] = sum(q.bot.tot[which(q.time >= smap.times[i] & q.time < smap.times[i+1])])
		calc.esoil[i] = calc.pre[i] - calc.qbot[i] - calc.trans[i] - calc.delsmap[i]
		}
		
		tmp = which(calc.qbot == 0)
		print(length(tmp))
		if(length(tmp) >= 1) calc.qbot[tmp] = NA
		if(length(tmp) >= 1) calc.esoil[tmp] = NA
		
		to.write = cbind(calc.time, calc.length, calc.evi, calc.pre, calc.delsmap, calc.trans, calc.qbot, calc.esoil)
		new.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/ESMAP_QC_DATA_new",sep='')
		write.table(to.write,new.file,sep="\t", col.names = F, row.names = F)
		print(new.file)
	}
