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
hydrus.dir = "/Volumes/LabShare/SMAP/ESMAP_for_Ronnie/Hudrus_Files/"
smap.points.file = "/Volumes/LabShare/SMAP/ESMAP_for_Ronnie/SMAP_9km_Points"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]
npoints = length(smap.lon)
for(i in 1:npoints){

		working.dir = paste(esmap.dir,smap.lat[i],"/",smap.lon[i],"/",sep='')
		#check if directory exists - if it does not, then it does not have at least 5 valid SMAP measurements during this time
		if (dir.exists(working.dir)){
		new.dir = paste(working.dir,"Hydrus_SpinUp/",sep='')
		if (dir.exists(new.dir)) {
		  unlink(new.dir,recursive = TRUE)
		}
		dir.create(new.dir,recursive=TRUE)
		
		soil.file = paste(working.dir,"Soil_raw",sep='')
		soil.raw = read.table(soil.file)
		
		smap.file = paste(working.dir,"SMAP_raw",sep='')
		smap.raw = read.table(smap.file)
		smap.ic = smap.raw[which(is.na(smap.raw[,4])==FALSE)[1],4]
		ndays = dim(smap.raw)[1]
		nhours = ndays*24
		#Copy Descript and Hydrus1D.Dat
		file.copy(paste(hydrus.dir,"DESCRIPT.TXT",sep=''),paste(new.dir,"DESCRIPT.TXT",sep=''))
		file.copy(paste(hydrus.dir,"HYDRUS1D.DAT",sep=''),paste(new.dir,"HYDRUS1D.DAT",sep=''))

		#Make Selector.in file
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

		
		#Write Profile.Dat file
		pro.file = paste(new.dir,"PROFILE.DAT",sep='')
		if(file.exists(pro.file)==TRUE) file.remove(pro.file)
		write("Pcp_File_Version=4",pro.file,append=TRUE)
		write("    2",pro.file,append=TRUE)
		write("    1  0.000000e+000  1.000000e+000  1.000000e+000",pro.file,append=TRUE)
		write("    2 -1.000000e+002  1.000000e+000  1.000000e+000",pro.file,append=TRUE)
		write("  101    0    0    1 x         h      Mat  Lay      Beta           Axz            Bxz            Dxz          Temp          Conc ",pro.file,append=TRUE)
		write(paste("    1 -0.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE)         
		write(paste("    2 -1.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    3 -2.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    4 -3.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    5 -4.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    6 -5.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    7 -6.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    8 -7.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("    9 -8.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   10 -9.000000e+000  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   11 -1.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   12 -1.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   13 -1.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   14 -1.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   15 -1.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   16 -1.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   17 -1.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   18 -1.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   19 -1.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   20 -1.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   21 -2.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   22 -2.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   23 -2.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   24 -2.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   25 -2.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   26 -2.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   27 -2.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   28 -2.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   29 -2.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   30 -2.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   31 -3.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   32 -3.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   33 -3.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   34 -3.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   35 -3.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   36 -3.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   37 -3.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   38 -3.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   39 -3.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   40 -3.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   41 -4.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   42 -4.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   43 -4.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   44 -4.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   45 -4.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   46 -4.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   47 -4.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   48 -4.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   49 -4.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   50 -4.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   51 -5.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   52 -5.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   53 -5.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   54 -5.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   55 -5.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   56 -5.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   57 -5.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   58 -5.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   59 -5.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   60 -5.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   61 -6.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   62 -6.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   63 -6.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   64 -6.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   65 -6.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   66 -6.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   67 -6.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   68 -6.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   69 -6.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   70 -6.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   71 -7.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   72 -7.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   73 -7.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   74 -7.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   75 -7.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   76 -7.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   77 -7.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   78 -7.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   79 -7.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   80 -7.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   81 -8.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   82 -8.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   83 -8.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   84 -8.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   85 -8.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   86 -8.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   87 -8.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   88 -8.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   89 -8.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   90 -8.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   91 -9.000000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   92 -9.100000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   93 -9.200000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   94 -9.300000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   95 -9.400000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   96 -9.500000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   97 -9.600000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   98 -9.700000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("   99 -9.800000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("  100 -9.900000e+001  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write(paste("  101 -1.000000e+002  ",sprintf("%13.6e",smap.ic),"    1    1  1.000000e+000  1.000000e+000  1.000000e+000  1.000000e+000              ",sep=''),pro.file,append=TRUE) 
		write("   20",pro.file,append=TRUE)
		write("    1    2    3    4    5    6    7    8    9   10   11   21   31   41   51   61   71   81   91   101",pro.file,append=TRUE)
		write("   ",pro.file,append=TRUE)
		
		#Write Atmosp.IN file
		trans.file = paste(working.dir,"Transpiration_Estimates",sep='')
		trans.raw = read.table(trans.file)
		pet = trans.raw[,4]
		pet[which(pet <= 0)] = 0
		
		evi.file = paste(working.dir,"EVI_raw",sep='')
		evi.raw = read.table(evi.file)
		d.evi = evi.raw[,4]
		tmp.evi = na.approx(d.evi)
		d.evi[which(is.na(d.evi) == FALSE)[1]:tail(which(is.na(d.evi) == FALSE),1)] = tmp.evi
		beg.delta = d.evi[which(is.na(d.evi) == FALSE)[2]] - d.evi[which(is.na(d.evi) == FALSE)[1]]
		beg.missing = which(is.na(d.evi) == TRUE)[which(which(is.na(d.evi) == TRUE) < which(is.na(d.evi) == FALSE)[1])]
		beg.evi = d.evi[which(is.na(d.evi) == FALSE)[1]] - (which(is.na(d.evi) == FALSE)[1] - beg.missing)*beg.delta
		d.evi[beg.missing] = beg.evi
		end.delta = d.evi[tail(which(is.na(d.evi) == FALSE),1)] - d.evi[tail(which(is.na(d.evi) == FALSE),2)[1]]
		end.missing = which(is.na(d.evi) == TRUE)
		end.evi = d.evi[tail(which(is.na(d.evi) == FALSE),1)] + (end.missing - tail(which(is.na(d.evi) == FALSE),1))*end.delta
		d.evi[end.missing] = end.evi
		
		nldas.file = paste(working.dir,"NLDAS_Forcing_raw",sep='')
		nldas.raw = read.table(nldas.file)
		
		if((smap.lon[i]-360) <= -67.5 & (smap.lon[i]-360) > -82.5)	offset = 5
		if((smap.lon[i]-360) <= -82.5 & (smap.lon[i]-360) > -97.5)	offset = 6
		if((smap.lon[i]-360) <= -97.5 & (smap.lon[i]-360) > -112.5) offset = 7
		if((smap.lon[i]-360) <= -112.5 & (smap.lon[i]-360) > -127.5) offset = 8
		
		write.pet = rep((pet/24),each=24)/10
		write.evi = rep((d.evi/1000),each=24)
		write.pre = nldas.raw[offset:(nhours+offset-1),8]/10
		write.pre[which(is.na(write.pre)==TRUE)] = 0
		
		atm.file = paste(new.dir,"ATMOSPH.IN",sep='')
		if(file.exists(atm.file)==TRUE) file.remove(atm.file)
		
		write("Pcp_File_Version=4",atm.file,append=TRUE)
		write("*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************",atm.file,append=TRUE)
		write("   MaxAL                    (MaxAL = number of atmospheric data-records)",atm.file,append=TRUE)
		write(paste("  ",round(nhours,digits=0),sep=''),atm.file,append=TRUE)
		write(" DailyVar  SinusVar  lLay  lBCCycles lInterc lDummy  lDummy  lDummy  lDummy  lDummy",atm.file,append=TRUE)
		write("       f       f       t       f       f       f       f       f       f       f",atm.file,append=TRUE)
		write(" Extinction",atm.file,append=TRUE)
		write("    0.39",atm.file,append=TRUE)
		write(" hCritS                 (max. allowed pressure head at the soil surface)",atm.file,append=TRUE)
		write("      1",atm.file,append=TRUE)
		write("       tAtm        Prec       rSoil       rRoot      hCritA          rB          hB          ht    RootDepth",atm.file,append=TRUE)
		hlist = seq(1:nhours)
		to.write = cbind(round(hlist,digits=0),signif(write.pre, digits = 6),signif(write.pet, digits = 10),signif(write.evi, digits = 8),100000,0,0,0)
		write.table(to.write,atm.file,sep=" ", col.names = F, row.names = F, quote=FALSE,append=TRUE)
		write("end*** END OF INPUT FILE 'ATMOSPH.IN' **********************************",atm.file,append=TRUE)
		write(" ",atm.file,append=TRUE)
		
		print(working.dir)
		} else {
		  store_bad_coord=cbind(store_bad_coord,smap.points[i,])
		}
		
}



