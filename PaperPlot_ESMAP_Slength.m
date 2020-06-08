clc;clear all;close all;

%load in ESMAP slength
slength=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','slength');
slength_vec=slength(:);
idx_nan=isnan(slength_vec);
slength_vec(idx_nan)=[];
idx1=find(slength_vec==1);
idx2=find(slength_vec==2);
idx3=find(slength_vec==3);

frequency_1=length(idx1)/length(slength_vec)*100;
frequency_2=length(idx2)/length(slength_vec)*100;
frequency_3=length(idx3)/length(slength_vec)*100;

figure
b1=bar([1,2,3],[frequency_1,frequency_2,frequency_3]);
xlim([0.5 3.5])
xticks([1 2 3])
xlabel({'E-SMAP Interval'; 'Length (days)'},'fontsize',25)
ylabel('Frequency (%)','fontsize',25)
set(gca,'fontsize',25)
b1.BarWidth=1;
b1.FaceColor=[0.0 0.75 0.75];
b1.LineWidth=3;