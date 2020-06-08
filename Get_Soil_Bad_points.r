rm(list=ls())

#library(MASS)
#library(ncdf4)
#library(fields)
#library(maps)
#library(hydroGOF)
#library(plotrix)
#library(rworldmap)
#library(RColorBrewer)
#library(humidity)

no_secondary.file = paste("/projects/roab9675/SMAP/Point_Check/No_Secondary_Altered_Points")
if(file.exists(no_secondary.file)==TRUE) file.remove(no_secondary.file)
secondary.file = paste("/projects/roab9675/SMAP/Point_Check/Secondary_Altered_Points")
if(file.exists(secondary.file)==TRUE) file.remove(secondary.file)

esmap.dir = "/scratch/summit/roab9675/SMAP/Gridded_ESMAP/"
soil.dir = "/projects/roab9675/SMAP/Data/"

#Open file to get dimensions
soil.file = paste(soil.dir,"NLDAS_Soil_Class.txt", sep='')
soil.raw = read.table(soil.file, fill=TRUE)
soil.lon = soil.raw[,3] + 360
soil.lat = soil.raw[,4]
soil.class = soil.raw[,5:20]

smap.points.file = "/projects/roab9675/SMAP/Point_Check/Altered_Points"
smap.points = read.table(smap.points.file)
smap.points = as.matrix(smap.points)
smap.lat = smap.points[,1]
smap.lon = smap.points[,2]


#remove NAs: (other files in directory that are not coordinates)
idx<-which(is.na(smap.lat))
if (length(idx)>0) {
smap.lat<-smap.lat[-idx]
}
nlat = length(smap.lat)
no_secondary=c()
secondary=c()
for(i in 1:nlat){
    
		lat.final = which(abs(soil.lat - smap.lat[i]) == min(abs(soil.lat - smap.lat[i])))
		lon.final = which(abs(soil.lon - smap.lon[i]) == min(abs(soil.lon - smap.lon[i])))
		tmp.soil = soil.class[intersect(lat.final,lon.final),]
		if(dim(tmp.soil)[1] >= 2) tmp.soil = colMeans(tmp.soil)
		final.soil = which.max(tmp.soil)
		if(final.soil== 14){
		  tmp.soil[14] = NA
		  final.soil = which.max(tmp.soil)
		}	
		##================================================================
		#                        Added section                           #
		##================================================================
		final.soil_first=final.soil
		#call NA on the first choice and move to second choice:
		tmp.soil[final.soil] = NA
		final.soil = which.max(tmp.soil)
		#if the secondary soil type is 0%, choose the next most similar soil type:
		if (tmp.soil[final.soil]==0){
		  if (final.soil_first==1 | final.soil_first==2 | final.soil_first==6 | final.soil_first==7 | final.soil_first==9){
		    final.soil=16
		  }
		  if (final.soil_first==5 | final.soil_first==8){
		    final.soil=16
		  }
		  if (final.soil_first==3 | final.soil_first==4){
		    final.soil=16
		  }
		  if (final.soil_first==13){
		    final.soil=16
		  }
		  if (final.soil_first==10 | final.soil_first==11 | final.soil_first==12 | final.soil_first==15 | final.soil_first==16){
		    final.soil=5
		  }
		  #store file to note you are using the next most similar soil type, rather than secondary:
          write(paste(smap.lat[i]," ",smap.lon[i],sep=' '),no_secondary.file,append=TRUE)

		} else {
          write(paste(smap.lat[i]," ",smap.lon[i],sep=' '),secondary.file,append=TRUE)
		}
		##================================================================
		#                End of Added section                           #
		##================================================================
		#write soil file:
		new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[i],"/Soil_updated",sep='')
		write.table(final.soil,new.file,sep="\t", col.names = F, row.names = F)
		
		#record points with secondary soil type:
		
		
		
		#record points with no secondary soil type:
		
		
		
		#record points with secondary soil type:
		
		#all.soil = which(tmp.soil >= 1)
		#new.file = paste(esmap.dir,smap.lat[i],"/",smap.lon[j],"/Soil_All_raw",sep='')
		#write.table(all.soil,new.file,sep="\t", col.names = F, row.names = F)
		
		to.screen = c(smap.lat[i], smap.lon[i], final.soil)
		print(to.screen)
		}

