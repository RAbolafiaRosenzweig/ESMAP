#install.packages("exactRankTests")
library(exactRankTests)

#Load in Screened and Unscreened Esoil data
#Mosaic
Mosaic_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_noScreening.csv',header=FALSE,sep=",");
Mosaic_og=as.matrix(Mosaic_og)
Mosaic_og=as.vector(Mosaic_og)
Screened_Mosaic=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_Screened.csv',header=FALSE,sep=",");
Screened_Mosaic=as.matrix(Screened_Mosaic)
Screened_Mosaic=as.vector(Screened_Mosaic)

Screened_Mosaic<-Screened_Mosaic[!is.na(Mosaic_og)]
Mosaic_og<-Mosaic_og[!is.na(Mosaic_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_Mosaic,y=Mosaic_og,alternative="greater",paired=TRUE)

#Noah
Noah_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_noScreening.csv',header=FALSE,sep=",");
Noah_og=as.matrix(Noah_og)
Noah_og=as.vector(Noah_og)
Screened_Noah=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_Screened.csv',header=FALSE,sep=",");
Screened_Noah=as.matrix(Screened_Noah)
Screened_Noah=as.vector(Screened_Noah)

Screened_Noah<-Screened_Noah[!is.na(Noah_og)]
Noah_og<-Noah_og[!is.na(Noah_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_Noah,y=Noah_og,alternative="greater",paired=TRUE)

#GLEAM
GLEAM_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_noScreening.csv',header=FALSE,sep=",");
GLEAM_og=as.matrix(GLEAM_og)
GLEAM_og=as.vector(GLEAM_og)
Screened_GLEAM=read.table('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_Screened.csv',header=FALSE,sep=",");
Screened_GLEAM=as.matrix(Screened_GLEAM)
Screened_GLEAM=as.vector(Screened_GLEAM)

Screened_GLEAM<-Screened_GLEAM[!is.na(GLEAM_og)]
GLEAM_og<-GLEAM_og[!is.na(GLEAM_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_GLEAM,y=GLEAM_og,alternative="greater",paired=TRUE)


##===============================================================================================
##Now test for Esoil/ET ratio
##===============================================================================================

#Load in Screened and Unscreened Esoil data
#Mosaic
Mosaic_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_noScreening.csv',header=FALSE,sep=",");
Mosaic_og=as.matrix(Mosaic_og)
Mosaic_og=as.vector(Mosaic_og)
Screened_Mosaic=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_Screened.csv',header=FALSE,sep=",");
Screened_Mosaic=as.matrix(Screened_Mosaic)
Screened_Mosaic=as.vector(Screened_Mosaic)

Screened_Mosaic<-Screened_Mosaic[!is.na(Mosaic_og)]
Mosaic_og<-Mosaic_og[!is.na(Mosaic_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_Mosaic,y=Mosaic_og,alternative="greater",paired=TRUE)

#Noah
Noah_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_ratio_noScreening.csv',header=FALSE,sep=",");
Noah_og=as.matrix(Noah_og)
Noah_og=as.vector(Noah_og)
Screened_Noah=read.table('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_ratio_Screened.csv',header=FALSE,sep=",");
Screened_Noah=as.matrix(Screened_Noah)
Screened_Noah=as.vector(Screened_Noah)

Screened_Noah<-Screened_Noah[!is.na(Noah_og)]
Noah_og<-Noah_og[!is.na(Noah_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_Noah,y=Noah_og,alternative="greater",paired=TRUE)

#GLEAM
GLEAM_og=read.table('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_noScreening.csv',header=FALSE,sep=",");
GLEAM_og=as.matrix(GLEAM_og)
GLEAM_og=as.vector(GLEAM_og)
Screened_GLEAM=read.table('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_Screened.csv',header=FALSE,sep=",");
Screened_GLEAM=as.matrix(Screened_GLEAM)
Screened_GLEAM=as.vector(Screened_GLEAM)

Screened_GLEAM<-Screened_GLEAM[!is.na(GLEAM_og)]
GLEAM_og<-GLEAM_og[!is.na(GLEAM_og)]

#determine if the mean is significantly different with wilcox rank sum test
wilcox.exact(x=Screened_GLEAM,y=GLEAM_og,alternative="greater",paired=TRUE)
