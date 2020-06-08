rm(list=ls())

source("/Volumes/LabShare/SMAP/ESMAP_for_Ronnie/Calc_Roots.r")

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
print(esmap.dir)
smap.lat = list.files(esmap.dir)
smap.lat = as.numeric(smap.lat)
nlat = length(smap.lat)
for(i in 1:nlat){

	smap.lon = list.files(paste(esmap.dir,smap.lat[i],"/",sep=''))
	nlon = length(smap.lon)
	smap.lon = as.numeric(smap.lon)
	for(j in 1:nlon){
		working.dir = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/",sep='')
		evi.file = paste(working.dir,"EVI_raw",sep='')
		for.file = paste(working.dir,"NLDAS_Forcing_raw",sep='')
		noah.file = paste(working.dir,"NLDAS_NOAH_raw",sep='')
		smap.file = paste(working.dir,"SMAP_raw",sep='')	
		veg.file = paste(working.dir,"Vegetation_raw",sep='')
		soil.file = paste(working.dir,"Soil_raw",sep='')
		
		evi.raw = read.table(evi.file)
		for.raw = read.table(for.file)
		noah.raw = read.table(noah.file)
		smap.raw = read.table(smap.file)
		veg.raw = read.table(veg.file)
		soil.raw = read.table(soil.file)
		
		d.years = smap.raw[,1]
		d.months = smap.raw[,2]
		d.days = smap.raw[,3]
		d.smap = smap.raw[,4]
		n.d = length(d.days)
		
		h.years = noah.raw[,1]
		h.months = noah.raw[,2]
		h.days = noah.raw[,3]
		h.hours = noah.raw[,4]
		h.sw = noah.raw[,5]
		h.lw = noah.raw[,6]
		n.h = length(h.days)
		
		h.pres = for.raw[,5]
		h.temp = for.raw[,6]
		h.sh = for.raw[,7]
		
		tmp.smap = na.approx(d.smap)
		d.smap[which(is.na(d.smap) == FALSE)[1]:tail(which(is.na(d.smap) == FALSE),1)] = tmp.smap
		beg.delta = d.smap[which(is.na(d.smap) == FALSE)[2]] - d.smap[which(is.na(d.smap) == FALSE)[1]]
		beg.missing = which(is.na(d.smap) == TRUE)[which(which(is.na(d.smap) == TRUE) < which(is.na(d.smap) == FALSE)[1])]
		beg.smap = d.smap[which(is.na(d.smap) == FALSE)[1]] - (which(is.na(d.smap) == FALSE)[1] - beg.missing)*beg.delta
		d.smap[beg.missing] = beg.smap
		end.delta = d.smap[tail(which(is.na(d.smap) == FALSE),1)] - d.smap[tail(which(is.na(d.smap) == FALSE),2)[1]]
		end.missing = which(is.na(d.smap) == TRUE)
		end.smap = d.smap[tail(which(is.na(d.smap) == FALSE),1)] + (end.missing - tail(which(is.na(d.smap) == FALSE),1))*end.delta
		d.smap[end.missing] = end.smap
		
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
		
		if((smap.lon[j]-360) <= -67.5 & (smap.lon[j]-360) > -82.5)	offset = 5
		if((smap.lon[j]-360) <= -82.5 & (smap.lon[j]-360) > -97.5)	offset = 6
		if((smap.lon[j]-360) <= -97.5 & (smap.lon[j]-360) > -112.5) offset = 7
		if((smap.lon[j]-360) <= -112.5 & (smap.lon[j]-360) > -127.5) offset = 8
		
		t.min = array(NA, c(n.d))
		t.mean = array(NA, c(n.d))
		sh.mean = array(NA, c(n.d))
		p.mean = array(NA, c(n.d))
		evi.val = array(NA, c(n.d))
		vpd.val = array(NA, c(n.d))
		rh.val = array(NA, c(n.d))
		ra.val = array(NA, c(n.d))
		m.t = array(NA, c(n.d))
		m.vpd = array(NA, c(n.d))
		svp.list = array(NA, c(n.d))
		e.list = array(NA, c(n.d))
		pe.list = array(NA, c(n.d))
		d.sw = array(NA, c(n.d))
		d.lw = array(NA, c(n.d))
		
		if(veg.raw == 1){ # Evergreen Needleleaf Forest	
			t.open = 8.31 + 273.15
			t.close = -8 + 273.15
			v.close = 3000
			v.open = 650
			cl = 0.0032
			roots = Calc_Roots(7,2,5)
			}

		if(veg.raw == 2){ # Evergreen Broadleaf Forest
			t.open = 9.09 + 273.15
			t.close = -8 + 273.15
			v.close = 4000
			v.open = 1000
			cl = 0.0025
			roots = Calc_Roots(7,1,5)
			}
	
		if(veg.raw == 3){ # 	Deciduous Needleleaf Forest	
			t.open = 10.44 + 273.15
			t.close = -8 + 273.15
			v.close = 3500
			v.open = 650
			cl = 0.0032
			roots = Calc_Roots(7,2,5)
			}

		if(veg.raw == 4){ #	Deciduous Broadleaf Forest	
			t.open = 9.94 + 273.15
			t.close = -6 + 273.15
			v.close = 2900
			v.open = 650
			cl = 0.0028
			roots = Calc_Roots(6,2,5)
			}

		if(veg.raw == 5){ #	Mixed Forest
			t.open = 9.5 + 273.15
			t.close = -7 + 273.15
			v.close = 2900
			v.open = 650
			cl = 0.0025
			roots = Calc_Roots(7,1.5,5)
			}
	
		if(veg.raw == 6){ # 	Woodland
			t.open = 11.39 + 273.15
			t.close = -8 + 273.15
			v.close = 3500
			v.open = 650
			cl = 0.0065
			roots = Calc_Roots(7,1.5,5)
			}
	
		if(veg.raw == 7){ # 	Wooded Grassland
			t.open = 11.39 + 273.15
			t.close = -8 + 273.15
			v.close = 3500
			v.open = 650
			cl = 0.0065
			roots = Calc_Roots(7,1.5,5)
			}
	
		if(veg.raw == 8){ # 	Closed Shrubland
			t.open = 8.61 + 273.15
			t.close = -8 + 273.15
			v.close = 4300
			v.open = 650
			cl = 0.0065
			roots = Calc_Roots(7,1.5,5)
			}
			
		if(veg.raw == 9){ # 	Open Shrubland
			t.open = 8.8 + 273.15
			t.close = -8 + 273.15
			v.close = 4400
			v.open = 650
			cl = 0.0065
			roots = Calc_Roots(7,1.5,5)
			}
	
		if(veg.raw == 10){ # 	Grassland
			t.open = 12.02 + 273.15
			t.close = -8 + 273.15
			v.close = 4200
			v.open = 650
			cl = 0.007
			roots = Calc_Roots(11,2,5)
			}
	
		if(veg.raw == 11){ # 	Cropland	
			t.open = 12.02 + 273.15
			t.close = -8 + 273.15
			v.close = 4500
			v.open = 650
			cl = 0.007
			roots = Calc_Roots(6,3,5)
			}

		if(veg.raw == 12){ # 	Bare Ground	
			t.open = 12.02 + 273.15
			t.close = -8 + 273.15
			v.close = 4200
			v.open = 650
			cl = 0.007
			roots = Calc_Roots(1,1,5)
			}
		
		for(k in 1:n.d){
			chosen = which(h.years == d.years[k] & h.months == d.months[k] & h.days == d.days[k])
			
			t.min[k] = min(h.temp[chosen + offset],na.rm=TRUE)
			t.mean[k] = mean(h.temp[chosen + offset],na.rm=TRUE)
			sh.mean[k] = mean(h.sh[chosen + offset],na.rm=TRUE)
			p.mean[k] = mean(h.pres[chosen + offset],na.rm=TRUE)
			d.sw[k] = sum(h.sw[chosen + offset],na.rm=TRUE)
			d.lw[k] = sum(h.lw[chosen + offset],na.rm=TRUE)
			
			tmpk = t.mean[k]
			tmpmin = t.min[k]
			
			svp = SVP.ClaCla((tmpk))*100
			svp.list[k] = svp
			
			tmpsh = sh.mean[k]
			tmpp = p.mean[k]
			
			e = tmpp*tmpsh/0.622
			e.list[k] = e
			rh = 100*e/svp
			
			vpd.val[k] = (1 - rh/100)*svp
			rh.val[k] = rh
			
			tmpvpd = vpd.val[k]
			
			if(tmpmin >= t.open) m.t[k] = 1
			if(tmpmin < t.open & tmpmin > t.close) m.t[k] = (tmpmin - t.close)/(t.open - t.close)
			if(tmpmin < t.close) m.t[k] = 0.1
	
			if(tmpvpd <= v.open) m.vpd[k] = 1
			if(tmpvpd > v.open & tmpvpd < v.close) m.vpd[k] = (v.close - tmpvpd)/(v.close - v.open)
			if(tmpvpd >= v.close) m.vpd[k] = 0.1
			
			}
		
		lai.val = d.evi/1000
		resistance = 1/(cl * m.t * m.vpd * lai.val)
		
		if(soil.raw == 1){ # 	Sand
			REFSMC = 0.236
			WLTSMC = 0.010
			}
		if(soil.raw == 2){ # 	Loamy sand
			REFSMC = 0.383
			WLTSMC = 0.028
			}
		if(soil.raw == 3){ # 	Sandy loam
			REFSMC = 0.383
			WLTSMC = 0.047
			}
		if(soil.raw == 4){ # 	Silt loam
			REFSMC = 0.360
			WLTSMC = 0.084
			}
		if(soil.raw == 5){ # 	Silt
			REFSMC = 0.383
			WLTSMC = 0.084
			}
		if(soil.raw == 6){ # 	Loam
			REFSMC = 0.329
			WLTSMC = 0.066
			}
		if(soil.raw == 7){ # 	Sandy clay loam
			REFSMC = 0.314
			WLTSMC = 0.067
			}
		if(soil.raw == 8){ # 	Silty clay loam
			REFSMC = 0.387
			WLTSMC = 0.120
			}
		if(soil.raw == 9){ # 	Clay loam
			REFSMC = 0.382
			WLTSMC = 0.103
			}
		if(soil.raw == 10){ # 	Sandy clay
			REFSMC = 0.338
			WLTSMC = 0.100
			}
		if(soil.raw == 11){ # 	Silty clay
			REFSMC = 0.404
			WLTSMC = 0.126
			}
		if(soil.raw == 12){ # 	Clay
			REFSMC = 0.412
			WLTSMC = 0.138
			}
		if(soil.raw == 13){ # 	Organic materials
			REFSMC = 0.329
			WLTSMC = 0.066
			}
		if(soil.raw == 15){ # 	Bedrock
			REFSMC = 0.17
			WLTSMC = 0.006
			}
		if(soil.raw == 16){ # 	Other
			REFSMC = 0.283
			WLTSMC = 0.028
			}
			
		A = d.sw + d.lw
		
		# To calculate Penman Monteith
		# ET = (s * A + rho * cp * (esat - e) / ra) / (s + gamma * (1 + rs / ra))
		# s = d(esat)/dT, the slope of the curve relating saturated water vapor pressure (esat) to temperature (Pa/K) - CALCULATED
		#	s = 4098 * 0.6108 * exp((17.27 + T)/(T + 237.3)) / ((T + 237.3)^2)
		#	T = Temperature (C) - FROM NLDAS FORCING
		# A = net radiation (W/m2) - SOLAR NET AND LONGWAVE NET IN NLDAS NOAH
		# rho = air density (kg/m3) - 1.225
		# cp = specific heat capacity of air (J/kg/K) - 1005
		# ra = aerodynamic resistance (s/m) - CALCULATE	
		#	ra = (rc * rr) / (rc + rr)
		#	rc = resistance to convective heat transfer (s/m) - CALUCLATE
		#		rc = 107 / (((273.15 + T)/393.15)^1.75 * 101300/P)
		#		T = degrees (C) - FROM NLDAS FORCING
		#		P = Pressure (Pa) - FROM NLDAS FORCING
		#	rr = resistance to radiative heat transfer  (s/m)
		#		rr	= (rho * cp) / (4 * sigma * T^3) 
		#		sigma = Stefan Boltzmann -  5.670 * 10^-8
		#		T = Temperature (K) -  - FROM NLDAS FORCING (CONVERT TO K)
		# gamma = psychometric constant (Pa/K)
		# 	gamma = (Ma / Mw) * (cp * P / lambda)
		# 	Ma = molecular masses of dry air (kg/mol) - 0.02897
		# 	Mw = molecular masses of wet air (kg/mol) - 0.01802
		#	cp = specific heat capacity of air (J/kg/K) - 1005
		# 	P = atmospheric pressure (Pa) - FROM NLDAS FORCING
		# 	lambda = latent heat of evaporation (J/kg) - 2501000
		# rs = surface resistance (s/m) - CALCULLATED VIA MODIS
		# e = water vapor pressure (Pa) - CALCULATED
		# esat = staturation vapor pressure (Pa) - CALCULATED
			
		rho = 1.225
		cp = 1005
		Ma = 0.02897
		Mw = 0.01802
		lambda = 2501000
		sigma = 5.670 * 10^-8
		
		T.c = t.mean - 273.15	
		s = 4098 * 0.6108 * exp((17.27 + T.c)/(T.c + 237.3)) / ((T.c + 237.3)^2)	
		rc = 107 / (((273.15 + T.c)/393.15)^1.75 * 101300/p.mean)	
		rr	= (rho * cp) / (4 * sigma * t.mean^3)
		ra = (rc * rr) / (rc + rr)
		gamma = (Ma / Mw) * (cp * p.mean / lambda)
		lai.min = min(lai.val)
		lai.max = max(lai.val)
		Fc = (lai.val - lai.min)/(lai.max - lai.min)
		RH = rh.val/100
		tmp = which(RH < 0.7)
		if(length(tmp) >= 1) RH[tmp] = 0
		Fwet = RH^4
		Ac = A*Fc
		PM.ET = (s * A + rho * cp * (svp.list - e.list) / ra) / (s + gamma * (1 + resistance / ra))
		PM.ET.22 = (s * Ac * Fc + rho * cp * (svp.list - e.list) * Fc / ra)*(1 - Fwet) / (s + gamma * (1 + resistance / ra))
		PM.ET.s = PM.ET / lambda
		PM.ET.day = PM.ET.s * 60 * 60 * 24
		PM.ET.22.s = PM.ET.22 / lambda
		PM.ET.22.day = PM.ET.22.s * 60 * 60 * 24
		PM.ET.22.day.roots = PM.ET.22.day * roots
		gx = (d.smap - WLTSMC) / (REFSMC - WLTSMC)
		tmp = which(gx <= 0)
		gx[tmp] = 0
		tmp = which(gx >= 1)
		gx[tmp] = 1
	
		PM.ET.22.day.roots.sm = PM.ET.22.day.roots * gx
		
		new.file = paste(working.dir,"Transpiration_Estimates",sep='')
		to.write = cbind(d.years, d.months, d.days, PM.ET.day, PM.ET.22.day, PM.ET.22.day.roots, PM.ET.22.day.roots.sm)
		write.table(to.write,new.file,sep="\t", col.names = F, row.names = F)
		}
	}