rm(list=ls())

library(MASS)
#library(fields)
#library(maps)
#library(hydroGOF)
#library(plotrix)
#library(rworldmap)
#library(RColorBrewer)
#library(humidity)

esmap.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP/"
smap.points.file = "/projects/roab9675/SMAP/Point_Check/Good_Points"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]
npoints = length(smap.lon)

bad.file = "/projects/roab9675/SMAP/Point_Check/TMP_BAD_Points2"
if(file.exists(bad.file)==TRUE) file.remove(bad.file)
good.file = "/projects/roab9675/SMAP/Point_Check/TMP_GOOD_Points2"
if(file.exists(good.file)==TRUE) file.remove(good.file)
	#npoints	
for(p in 1:npoints){
	old.dir = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_SpinUp/",sep='')
	new.dir = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Hydrus_FINAL/",sep='')
    if(dir.exists(new.dir)==TRUE) unlink(new.dir,recursive=TRUE)
	dir.create(new.dir,recursive=TRUE)
	print(new.dir)
	#Open NOD_INF and get new SM
	nod.file = paste(old.dir,"NOD_INF.OUT", sep='')
	opt1 = grep('Time:     35064',readLines(nod.file))
	opt2 = grep('Time:     8802.7600',readLines(nod.file))
	opt3 = grep('Time:     8727.7300',readLines(nod.file))
	starts = grep('1    -0.0000',readLines(nod.file))
	ends = grep('101  -100.0000',readLines(nod.file))
	
	opt = c(opt1,opt2,opt3)[which(c(opt1,opt2,opt3) >= 1000)[1]]
	if(is.na(opt) == TRUE){
		write(paste(smap.lat[p]," ",smap.lon[p],sep=''),bad.file,append=TRUE)
		
		#Make Selector.in file
		soil.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Soil_raw",sep='')
		soil.raw = read.table(soil.file)
		
		smap.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/SMAP_raw",sep='')
		smap.raw = read.table(smap.file)
		ndays = dim(smap.raw)[1]
		nhours = ndays*24
		sel.file = paste(old.dir,"SELECTOR.IN",sep='')
		if(file.exists(sel.file)==TRUE) file.remove(sel.file)
		write("Pcp_File_Version=4",sel.file,append=TRUE)
		write("*** BLOCK A: BASIC INFORMATION *****************************************",sel.file,append=TRUE)
		write("Heading",sel.file,append=TRUE)
		write("dry to wet",sel.file,append=TRUE)
		write("LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)",sel.file,append=TRUE)
		write("cm",sel.file,append=TRUE)
		write("hours",sel.file,append=TRUE)
		write("mmol",sel.file,append=TRUE)
		write("lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse",sel.file,append=TRUE)
		write(" t     f     f      t     f     f      f     t       t         t         f",sel.file,append=TRUE)
		write("lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy",sel.file,append=TRUE)
 		write(" f       f       f       f       f       t       f       f       f       f",sel.file,append=TRUE)
		write("NMat    NLay  CosAlpha",sel.file,append=TRUE)
  		write("  1       1       1",sel.file,append=TRUE)
		write("*** BLOCK B: WATER FLOW INFORMATION ************************************",sel.file,append=TRUE)
		write("MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)",sel.file,append=TRUE)
 	 	write("  300    0.001      1",sel.file,append=TRUE)
		write("TopInf WLayer KodTop InitCond",sel.file,append=TRUE)
 		write(" t     t      -1       t",sel.file,append=TRUE)
		write("BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep",sel.file,append=TRUE)
 		write(" f     f     t     f     -1      f      0",sel.file,append=TRUE)
    	write("    hTab1   hTabN",sel.file,append=TRUE)
        write("        0       0",sel.file,append=TRUE)
    	write("    Model   Hysteresis",sel.file,append=TRUE)
      	write("      0          0",sel.file,append=TRUE)
   		write("   thr     ths    Alfa      n         Ks       l",sel.file,append=TRUE)
		if(soil.raw == 1){ # 	Sand - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000107 * 360000
			l = 0.5
			}
		if(soil.raw == 2){ # 	Loamy sand - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.0000141 * 360000
			l = 0.5
			}
		if(soil.raw == 3){ # 	Sandy loam - medium
			thr = 0.010
			ths = 0.439
			alfa = 0.0314
			n = 1.1804
			Ks = 0.00000523 * 360000
			l = 0.5
			}
		if(soil.raw == 4){ # 	Silt loam - medium
			thr = 0.010
			ths = 0.439
			alfa = 0.0314
			n = 1.1804
			Ks = 0.00000281 * 360000
			l = 0.5
			}
		if(soil.raw == 5){ # 	Silt - medium fine
			thr = 0.010
			ths = 0.430
			alfa = 0.0083
			n = 1.2539
			Ks = 0.00000281 * 360000
			l = 0.5
			}
		if(soil.raw == 6){ # 	Loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000338 * 360000
			l = 0.5
			}
		if(soil.raw == 7){ # 	Sandy clay loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000445 * 360000
			l = 0.5
			}
		if(soil.raw == 8){ # 	Silty clay loam - medium fine
			thr = 0.010
			ths = 0.430
			alfa = 0.0083
			n = 1.2539
			Ks = 0.00000204 * 360000
			l = 0.5
			}
		if(soil.raw == 9){ # 	Clay loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000245 * 360000
			l = 0.5
			}
		if(soil.raw == 10){ # 	Sandy clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000722 * 360000
			l = 0.5
			}
		if(soil.raw == 11){ # 	Silty clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000134 * 360000
			l = 0.5
			}
		if(soil.raw == 12){ # 	Clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.000000974 * 360000
			l = 0.5
			}
		if(soil.raw == 13){ # 	Organic materials
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000338 * 360000
			l = 0.5
			}
		if(soil.raw == 15){ # 	Bedrock
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.0000000974 * 360000
			l = 0.5
			}
		if(soil.raw == 16){ # 	Other
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000134 * 360000
			l = 0.5
			}
		write(paste("  ",round(thr,digits=5),"    ",round(ths,digits=5),"    ",round(alfa,digits=5),"    ",round(n,digits=5),"    ",round(Ks,digits=5),"    ",round(l,digits=5)),sel.file,append=TRUE)
		write("*** BLOCK C: TIME INFORMATION ******************************************",sel.file,append=TRUE)
		write("        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL",sel.file,append=TRUE)
		nouts = 468
		write(paste("      0.024     0.00024        0.96     1.3     0.7     3     7   ",round(nouts,digits=0),sep=''),sel.file,append=TRUE)
		write("      tInit        tMax",sel.file,append=TRUE)
		write(paste("          0       ",nhours,sep=''),sel.file,append=TRUE)
		write("  lPrintD  nPrintSteps tPrintInterval lEnter",sel.file,append=TRUE)
		write("     t           1             1       t",sel.file,append=TRUE)
		write("TPrint(1),TPrint(2),...,TPrint(MPL)",sel.file,append=TRUE)
		time = seq(24,nhours,length.out=nouts)
		for(k in 1:(nouts/6)){
			write(paste("     ",round(time[(k*6-5)],digits=2),"     ",round(time[(k*6-4)],digits=2),"     ",round(time[(k*6-3)],digits=2),"     ",round(time[(k*6-2)],digits=2),"     ",round(time[(k*6-1)],digits=2),"     ",round(time[(k*6-0)],digits=2),sep=''),sel.file,append=TRUE)
			}
		write("*** BLOCK G: ROOT WATER UPTAKE INFORMATION *****************************",sel.file,append=TRUE)
		write("     Model  (0 - Feddes, 1 - S shape)  cRootMax    OmegaC",sel.file,append=TRUE)
		write("        0                                   1",sel.file,append=TRUE)
		write("       P0       P2H       P2L       P3          r2H        r2L",sel.file,append=TRUE)
		write("       -1     -8000    -12000    -24000   0.0208333  0.00416667",sel.file,append=TRUE)
		write("POptm(1),POptm(2),...,POptm(NMat)",sel.file,append=TRUE)
		write("     -25 ",sel.file,append=TRUE)
		write("*** END OF INPUT FILE 'SELECTOR.IN' ************************************",sel.file,append=TRUE)
		write("    ",sel.file,append=TRUE)
		}
		
	if(is.na(opt) == FALSE){
		write(paste(smap.lat[p]," ",smap.lon[p],sep=''),good.file,append=TRUE)
				
		lstart = starts[which(starts > opt)[1]]
		lend = ends[which(ends > opt)[1]]

		nod.raw = readLines(nod.file)
		nod.chosen = nod.raw[lstart:lend]
	
		newsm = array(NA, c(101))
		for(i in 1:101){
			tmp = nod.chosen[i]
			newsm[i] = as.numeric(strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[4]])
			
			}
	
		#Create new PROFILE.DAT file
	
		pro.file = paste(old.dir,"PROFILE.DAT", sep='')
		pro.raw = readLines(pro.file)
		for(i in 1:101){
			tmp = pro.raw[(i+5)]
			value = strsplit(tmp," ")[[1]][which(strsplit(tmp," ")[[1]] != "")[3]]
			new = gsub(value, format(newsm[i], scientific=TRUE), tmp)
			pro.raw[(i+5)] = new
			}
		new.file = paste(new.dir,"PROFILE.DAT", sep='')
		write(pro.raw, new.file)
		
		#Make Selector.in file
		soil.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/Soil_raw",sep='')
		soil.raw = read.table(soil.file)
		
		smap.file = paste(esmap.dir,smap.lat[p],"/",smap.lon[p],"/SMAP_raw",sep='')
		smap.raw = read.table(smap.file)
		ndays = dim(smap.raw)[1]
		nhours = ndays*24
		sel.file = paste(new.dir,"SELECTOR.IN",sep='')
		if(file.exists(sel.file)==TRUE) file.remove(sel.file)
		write("Pcp_File_Version=4",sel.file,append=TRUE)
		write("*** BLOCK A: BASIC INFORMATION *****************************************",sel.file,append=TRUE)
		write("Heading",sel.file,append=TRUE)
		write("dry to wet",sel.file,append=TRUE)
		write("LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)",sel.file,append=TRUE)
		write("cm",sel.file,append=TRUE)
		write("hours",sel.file,append=TRUE)
		write("mmol",sel.file,append=TRUE)
		write("lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse",sel.file,append=TRUE)
		write(" t     f     f      t     f     f      f     t       t         t         f",sel.file,append=TRUE)
		write("lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy",sel.file,append=TRUE)
 		write(" f       f       f       f       f       t       f       f       f       f",sel.file,append=TRUE)
		write("NMat    NLay  CosAlpha",sel.file,append=TRUE)
  		write("  1       1       1",sel.file,append=TRUE)
		write("*** BLOCK B: WATER FLOW INFORMATION ************************************",sel.file,append=TRUE)
		write("MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)",sel.file,append=TRUE)
 	 	write("  10    0.001      1",sel.file,append=TRUE)
		write("TopInf WLayer KodTop InitCond",sel.file,append=TRUE)
 		write(" t     t      -1       t",sel.file,append=TRUE)
		write("BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep",sel.file,append=TRUE)
 		write(" f     f     t     f     -1      f      0",sel.file,append=TRUE)
    	write("    hTab1   hTabN",sel.file,append=TRUE)
        write("        0       0",sel.file,append=TRUE)
    	write("    Model   Hysteresis",sel.file,append=TRUE)
      	write("      0          0",sel.file,append=TRUE)
   		write("   thr     ths    Alfa      n         Ks       l",sel.file,append=TRUE)
		if(soil.raw == 1){ # 	Sand - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000107 * 360000
			l = 0.5
			}
		if(soil.raw == 2){ # 	Loamy sand - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.0000141 * 360000
			l = 0.5
			}
		if(soil.raw == 3){ # 	Sandy loam - medium
			thr = 0.010
			ths = 0.439
			alfa = 0.0314
			n = 1.1804
			Ks = 0.00000523 * 360000
			l = 0.5
			}
		if(soil.raw == 4){ # 	Silt loam - medium
			thr = 0.010
			ths = 0.439
			alfa = 0.0314
			n = 1.1804
			Ks = 0.00000281 * 360000
			l = 0.5
			}
		if(soil.raw == 5){ # 	Silt - medium fine
			thr = 0.010
			ths = 0.430
			alfa = 0.0083
			n = 1.2539
			Ks = 0.00000281 * 360000
			l = 0.5
			}
		if(soil.raw == 6){ # 	Loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000338 * 360000
			l = 0.5
			}
		if(soil.raw == 7){ # 	Sandy clay loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000445 * 360000
			l = 0.5
			}
		if(soil.raw == 8){ # 	Silty clay loam - medium fine
			thr = 0.010
			ths = 0.430
			alfa = 0.0083
			n = 1.2539
			Ks = 0.00000204 * 360000
			l = 0.5
			}
		if(soil.raw == 9){ # 	Clay loam - coarse
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000245 * 360000
			l = 0.5
			}
		if(soil.raw == 10){ # 	Sandy clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000722 * 360000
			l = 0.5
			}
		if(soil.raw == 11){ # 	Silty clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000134 * 360000
			l = 0.5
			}
		if(soil.raw == 12){ # 	Clay - fine
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.000000974 * 360000
			l = 0.5
			}
		if(soil.raw == 13){ # 	Organic materials
			thr = 0.025
			ths = 0.403
			alfa = 0.0383
			n = 1.3774
			Ks = 0.00000338 * 360000
			l = 0.5
			}
		if(soil.raw == 15){ # 	Bedrock
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.0000000974 * 360000
			l = 0.5
			}
		if(soil.raw == 16){ # 	Other
			thr = 0.010
			ths = 0.614
			alfa = 0.0265
			n = 1.1033
			Ks = 0.00000134 * 360000
			l = 0.5
			}
		write(paste("  ",round(thr,digits=5),"    ",round(ths,digits=5),"    ",round(alfa,digits=5),"    ",round(n,digits=5),"    ",round(Ks,digits=5),"    ",round(l,digits=5)),sel.file,append=TRUE)
		write("*** BLOCK C: TIME INFORMATION ******************************************",sel.file,append=TRUE)
		write("        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL",sel.file,append=TRUE)
		nouts = 468
		write(paste("      0.024     0.00024        0.96     1.3     0.7     3     7   ",round(nouts,digits=0),sep=''),sel.file,append=TRUE)
		write("      tInit        tMax",sel.file,append=TRUE)
		write(paste("          0       ",nhours,sep=''),sel.file,append=TRUE)
		write("  lPrintD  nPrintSteps tPrintInterval lEnter",sel.file,append=TRUE)
		write("     t           1             1       t",sel.file,append=TRUE)
		write("TPrint(1),TPrint(2),...,TPrint(MPL)",sel.file,append=TRUE)
		time = seq(24,nhours,length.out=nouts)
		for(k in 1:(nouts/6)){
			write(paste("     ",round(time[(k*6-5)],digits=2),"     ",round(time[(k*6-4)],digits=2),"     ",round(time[(k*6-3)],digits=2),"     ",round(time[(k*6-2)],digits=2),"     ",round(time[(k*6-1)],digits=2),"     ",round(time[(k*6-0)],digits=2),sep=''),sel.file,append=TRUE)
			}
		write("*** BLOCK G: ROOT WATER UPTAKE INFORMATION *****************************",sel.file,append=TRUE)
		write("     Model  (0 - Feddes, 1 - S shape)  cRootMax    OmegaC",sel.file,append=TRUE)
		write("        0                                   1",sel.file,append=TRUE)
		write("       P0       P2H       P2L       P3          r2H        r2L",sel.file,append=TRUE)
		write("       -1     -8000    -12000    -24000   0.0208333  0.00416667",sel.file,append=TRUE)
		write("POptm(1),POptm(2),...,POptm(NMat)",sel.file,append=TRUE)
		write("     -25 ",sel.file,append=TRUE)
		write("*** END OF INPUT FILE 'SELECTOR.IN' ************************************",sel.file,append=TRUE)
		write("    ",sel.file,append=TRUE)
		
		#Copy Descript and Hydrus1D.Dat
		file.copy(paste(old.dir,"DESCRIPT.TXT",sep=''),paste(new.dir,"DESCRIPT.TXT",sep=''))
		file.copy(paste(old.dir,"HYDRUS1D.DAT",sep=''),paste(new.dir,"HYDRUS1D.DAT",sep=''))
		file.copy(paste(old.dir,"ATMOSPH.IN",sep=''),paste(new.dir,"ATMOSPH.IN",sep=''))
		}
	print(new.dir)
	}
