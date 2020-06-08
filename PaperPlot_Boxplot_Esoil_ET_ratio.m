clc;clear all;close all;

%define data to plot:

%Mosaic
Mosaic_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_noScreening.csv');
Screened_Mosaic=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_ratio_Screened.csv');
Mosaic_ratio=Screened_Mosaic./Mosaic_og;
Mosaic_ratio=Mosaic_ratio(:);

%Noah
Noah_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_ratio_noScreening.csv');
Screened_Noah=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_ratio_Screened.csv');
Noah_ratio=Screened_Noah./Noah_og;
Noah_ratio=Noah_ratio(:);

%GLEAM
GLEAM_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_noScreening.csv');
Screened_GLEAM=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_ratio_Screened.csv');
GLEAM_ratio=Screened_GLEAM./GLEAM_og;
GLEAM_ratio=GLEAM_ratio(:);

%trim negatives:
idx=find(Mosaic_ratio<0);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

idx=find(Noah_ratio<0);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

idx=find(GLEAM_ratio<0);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

%trim inf:
idx=isinf(Mosaic_ratio);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

idx=isinf(Noah_ratio);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

idx=isinf(GLEAM_ratio);
Mosaic_ratio(idx)=[];
Noah_ratio(idx)=[];
GLEAM_ratio(idx)=[];

%create boxplot of the data:
figure
hold on
b=boxplot([GLEAM_ratio,Mosaic_ratio,Noah_ratio],'Labels',{'GLEAM','Mosaic','Noah'},'color','k','Whisker',1,'Notch','on');
set(b,'linewidth',3)
ylabel('Valid/Continuous E_{soil}/ET Ratio','fontsize',30)
set(gca,'fontsize',34)
lines=findobj(b,'type','line','Tag','Median');
set(lines,'Color','b')
ylim([0 3])
xlim([0 4])
%plot 1.0 ref line
plot(linspace(0,4,10),linspace(1,1,10),'--k','linewidth',3)

sprintf('Esoil/ET increases by %.2f, %.2f, and %.2f percent for GLEAM, Mosaic and Noah, respectively',nanmedian(GLEAM_ratio),nanmedian(Mosaic_ratio),nanmedian(Noah_ratio))


%%================================================================================================================================



%Mosaic
Mosaic_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_noScreening.csv');
Screened_Mosaic=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Esoil_ET_ratio/Mosaic_esoil_Screened.csv');
Mosaic_esoil=Screened_Mosaic./Mosaic_og;
Mosaic_esoil=Mosaic_esoil(:);

%Noah
Noah_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_noScreening.csv');
Screened_Noah=csvread('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/Esoil_ET_ratio/Noah_esoil_Screened.csv');
Noah_esoil=Screened_Noah./Noah_og;
Noah_esoil=Noah_esoil(:);

%GLEAM
GLEAM_og=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_noScreening.csv');
Screened_GLEAM=csvread('/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil_ET_ratio/GLEAM_esoil_Screened.csv');
GLEAM_esoil=Screened_GLEAM./GLEAM_og;
GLEAM_esoil=GLEAM_esoil(:);

%trim negatives:
idx=find(Mosaic_esoil<0);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

idx=find(Noah_esoil<0);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

idx=find(GLEAM_esoil<0);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

%trim inf:
idx=isinf(Mosaic_esoil);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

idx=isinf(Noah_esoil);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

idx=isinf(GLEAM_esoil);
Mosaic_esoil(idx)=[];
Noah_esoil(idx)=[];
GLEAM_esoil(idx)=[];

%create boxplot of the data:
figure
hold on
b=boxplot([GLEAM_esoil,Mosaic_esoil,Noah_esoil],'Labels',{'GLEAM','Mosaic','Noah'},'color','k','Whisker',1,'Notch','on');
set(b,'linewidth',3)
ylabel('Valid/Continuous E_{soil} Ratio','fontsize',25)
set(gca,'fontsize',34)
lines=findobj(b,'type','line','Tag','Median');
set(lines,'Color','b')
ylim([0.5 1.5])
xlim([0 4])
%plot 1.0 ref line
plot(linspace(0,4,10),linspace(1,1,10),'--k','linewidth',3)

sprintf('Esoil increases by %.2f, %.2f, and %.2f percent for GLEAM, Mosaic and Noah, respectively',nanmedian(GLEAM_esoil),nanmedian(Mosaic_esoil),nanmedian(Noah_esoil))


%report stats of the IQR and whiskers:
prctile25_GLEAM=prctile(GLEAM_esoil,25)
prctile75_GLEAM=prctile(GLEAM_esoil,75)
lower_whisker_gleam=prctile25_GLEAM-(prctile75_GLEAM-prctile25_GLEAM)

prctile25_Mosaic=prctile(Mosaic_esoil,25)
prctile75_Mosaic=prctile(Mosaic_esoil,75)
lower_whisker_mosaic=prctile25_Mosaic-(prctile75_Mosaic-prctile25_Mosaic)

prctile25_Noah=prctile(Noah_esoil,25)
prctile75_Noah=prctile(Noah_esoil,75)
lower_whisker_noah=prctile25_Noah-(prctile75_Noah-prctile25_Noah)
