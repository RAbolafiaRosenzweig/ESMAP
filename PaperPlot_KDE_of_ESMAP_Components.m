clc;clear all;close all;

%The goal of this is to show that ESMAP esoil is dominated by dS:
delsmap=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','delsmap');
esoil=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','esoil_screened');
Et=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','transp');
precip=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','infilt');
qbot=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','qbot');
slength=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','slength');
time=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','time');


idx=isnan(esoil);
delsmap(idx)=NaN;
Et(idx)=NaN;
precip(idx)=NaN;
qbot(idx)=NaN;

mean_delsmap=nanmean(delsmap,3);
mean_esoil=nanmean(esoil,3);
ratio=-delsmap./esoil;
dS_Esoil_ratio=nanmedian(ratio,3);
% dS_Esoil_ratio=-mean_delsmap./mean_esoil;
%%%============================================================================================================%%%============================================================================================================
%Plot dS/Esoil fraction
%%%============================================================================================================%%%============================================================================================================R2=csvread('/Volumes/REESEN/SMAP/Validation_Data/ESMAP_Comparisson/GLEAM_Noah_R2.csv');
lat=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','lat');
lon=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','lon');
[LAT,LON]=meshgrid(lat,lon);

%calculate the Esoil/EST ratio for screened vs. unscreened:
PLOT_DATA=dS_Esoil_ratio;

%define bounding domain box
latlim=[24.9,50.1];
lonlim=[-125.1,-66.9];

valid_min=0.4;
valid_max=1;
% % Create a set of level ranges to be used in converting the data to a
% % geolocated image that has a color assigned to each range.
levels = linspace(valid_min, valid_max, 7);

%polarmap goes from blue to red (flip this)
levels = linspace(valid_min, valid_max, 7);

%polarmap goes from blue to red (flip this)
cmap1 = flipud(copper(12)); 
cmap2 = flipud(summer(12)); 

cmap = [cmap2([12,8,3,1],:);cmap1([3,7,12],:)];

% % % % % Create a color map.
% Convert the data to an geolcated image by setting a color for each level
% range.
Z = PLOT_DATA;
% 
% Clamp the min and max values to the level index.
Z(Z > levels(end)) = length(levels);
Z(Z < levels(1)) = 1;


% Assign Z as an indexed image with the index value corresponding to the
% level range.
for k = 1:length(levels) - 1
    Z(PLOT_DATA >= levels(k) & PLOT_DATA < levels(k+1)) = double(k) ;
end

% Plot the data.
figure

    xlim(lonlim);
    ylim(latlim);
    hold on
% %tightmap
colormap(cmap)
p1=geoshow(LAT,LON,Z,'displaytype','surface');

%Display outline of land area
landareas = shaperead('landareas.shp', 'UseGeo', true);
coast.lat = [landareas.Lat];
coast.long = [landareas.Lon];
geoshow(coast.lat, coast.long, 'Color', 'k','LineWidth',3)
%overlay states:
states=shaperead('usastatehi','UseGeoCoords',true,'BoundingBox',[lonlim' latlim']);
p2=geoshow(states,'FaceColor','none','Linewidth',3);
for i=1:49
    p2.Children(i).ZData=10*ones(length(p2.Children(i).YData),1);
end

caxis auto
levels=round(levels,4);
clevels =  cellstr(num2str(levels'));
clevels = [clevels]';
h = lcolorbar(clevels, 'Location', 'horizontal','fontsize',24,'title','SMAP Drying / E_{soil}');
idx=3:2:length(h.XTickLabel);
set(gca,'fontsize',28)
h.Position(2)=h.Position(2)-0.05;

%Regional bounds are from Xia et al., 2014
%define latlims for regions
SE_latlim=[25, 37.2];
NE_latlim=[37.2, 50];
GP_latlim=[25 , 37.2];
MW_latlim=[37.2,50];
NW_latlim=[44 50];
SW_latlim=[30 44];
%define lon lims for regions
SE_lonlim=[-94.1, -75];
NE_lonlim=[-88,-67];
GP_lonlim=[-108.2,-94.1];
MW_lonlim=[-108.2,-88];
NW_lonlim=[-125 -108.2];
SW_lonlim=[-125 -108.2];


%box for Southeast:
bottom_line_y=linspace(SE_latlim(1),SE_latlim(1),20);
bottom_line_x=linspace(SE_lonlim(1),SE_lonlim(2),20);
top_line_y=linspace(SE_latlim(2),SE_latlim(2),20);
left_line_x=linspace(SE_lonlim(1),SE_lonlim(1),20);
right_line_x=linspace(SE_lonlim(2),SE_lonlim(2),20);
left_line_y=linspace(SE_latlim(1),SE_latlim(2),20);

p3=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p4=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p5=plot(left_line_x,left_line_y,'-k','linewidth',5);
p6=plot(right_line_x,left_line_y,'-k','linewidth',5);

p3.ZData=10*ones(length(p3.YData),1);
p4.ZData=10*ones(length(p4.YData),1);
p5.ZData=10*ones(length(p5.YData),1);
p6.ZData=10*ones(length(p6.YData),1);

%box for Northeast:
bottom_line_y=linspace(NE_latlim(1),NE_latlim(1),20);
bottom_line_x=linspace(NE_lonlim(1),NE_lonlim(2),20);
top_line_y=linspace(NE_latlim(2),NE_latlim(2),20);
left_line_x=linspace(NE_lonlim(1),NE_lonlim(1),20);
right_line_x=linspace(NE_lonlim(2),NE_lonlim(2),20);
left_line_y=linspace(NE_latlim(1),NE_latlim(2),20);

p7=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p8=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p9=plot(left_line_x,left_line_y,'-k','linewidth',5);
p10=plot(right_line_x,left_line_y,'-k','linewidth',5);

p7.ZData=10*ones(length(p7.YData),1);
p8.ZData=10*ones(length(p8.YData),1);
p9.ZData=10*ones(length(p9.YData),1);
p10.ZData=10*ones(length(p10.YData),1);

%box for Greatplains:
bottom_line_y=linspace(GP_latlim(1),GP_latlim(1),20);
bottom_line_x=linspace(GP_lonlim(1),GP_lonlim(2),20);
top_line_y=linspace(GP_latlim(2),GP_latlim(2),20);
left_line_x=linspace(GP_lonlim(1),GP_lonlim(1),20);
right_line_x=linspace(GP_lonlim(2),GP_lonlim(2),20);
left_line_y=linspace(GP_latlim(1),GP_latlim(2),20);

p11=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p12=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p13=plot(left_line_x,left_line_y,'-k','linewidth',5);
p14=plot(right_line_x,left_line_y,'-k','linewidth',5);

p11.ZData=10*ones(length(p11.YData),1);
p12.ZData=10*ones(length(p12.YData),1);
p13.ZData=10*ones(length(p13.YData),1);
p14.ZData=10*ones(length(p14.YData),1);

%box for Midwest:
bottom_line_y=linspace(MW_latlim(1),MW_latlim(1),20);
bottom_line_x=linspace(MW_lonlim(1),MW_lonlim(2),20);
top_line_y=linspace(MW_latlim(2),MW_latlim(2),20);
left_line_x=linspace(MW_lonlim(1),MW_lonlim(1),20);
right_line_x=linspace(MW_lonlim(2),MW_lonlim(2),20);
left_line_y=linspace(MW_latlim(1),MW_latlim(2),20);

p15=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p16=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p17=plot(left_line_x,left_line_y,'-k','linewidth',5);
p18=plot(right_line_x,left_line_y,'-k','linewidth',5);

p15.ZData=10*ones(length(p15.YData),1);
p16.ZData=10*ones(length(p16.YData),1);
p17.ZData=10*ones(length(p17.YData),1);
p18.ZData=10*ones(length(p18.YData),1);

%box for Northwest:
bottom_line_y=linspace(NW_latlim(1),NW_latlim(1),20);
bottom_line_x=linspace(NW_lonlim(1),NW_lonlim(2),20);
top_line_y=linspace(NW_latlim(2),NW_latlim(2),20);
left_line_x=linspace(NW_lonlim(1),NW_lonlim(1),20);
right_line_x=linspace(NW_lonlim(2),NW_lonlim(2),20);
left_line_y=linspace(NW_latlim(1),NW_latlim(2),20);

p19=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p20=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p21=plot(left_line_x,left_line_y,'-k','linewidth',5);
p22=plot(right_line_x,left_line_y,'-k','linewidth',5);

p19.ZData=10*ones(length(p19.YData),1);
p20.ZData=10*ones(length(p20.YData),1);
p21.ZData=10*ones(length(p21.YData),1);
p22.ZData=10*ones(length(p22.YData),1);

%box for Southwest:
bottom_line_y=linspace(SW_latlim(1),SW_latlim(1),20);
bottom_line_x=linspace(SW_lonlim(1),SW_lonlim(2),20);
top_line_y=linspace(SW_latlim(2),SW_latlim(2),20);
left_line_x=linspace(SW_lonlim(1),SW_lonlim(1),20);
right_line_x=linspace(SW_lonlim(2),SW_lonlim(2),20);
left_line_y=linspace(SW_latlim(1),SW_latlim(2),20);

p23=plot(bottom_line_x,bottom_line_y,'-k','linewidth',5);
p24=plot(bottom_line_x,top_line_y,'-k','linewidth',5);
p25=plot(left_line_x,left_line_y,'-k','linewidth',5);
p26=plot(right_line_x,left_line_y,'-k','linewidth',5);

p23.ZData=10*ones(length(p23.YData),1);
p24.ZData=10*ones(length(p24.YData),1);
p25.ZData=10*ones(length(p25.YData),1);
p26.ZData=10*ones(length(p26.YData),1);

%format axis labels:
AC=gca;
xticks=AC.XTickLabel;

%x-axis (lon)
for i=1:length(xticks)
    current_tick=str2num(xticks{i});
    current_tick=abs(current_tick);
    current_label=[sprintf('%d',current_tick),'\circW'];
    xticks{i}=current_label;
end
xticklabels(xticks)

%y-axis (lat)
yticks=AC.YTickLabel;

%x-axis (lon)
for i=1:length(yticks)
    current_tick=str2num(yticks{i});
    current_tick=abs(current_tick);
    current_label=[sprintf('%d',current_tick),'\circN'];
    yticks{i}=current_label;
end
yticklabels(yticks)

% % % %explore low ratio around 30-degrees lat, for example
% % % %(lon=-106.9,lat=29.54) - precip is large
% % % %(lon=-103.6,lat=29.78) - precip is large
% % % %(lon=-107.1,lat=30.68) - 
% % % latvec=LAT(:);
% % % lonvec=LON(:);
% % % IDX=knnsearch([latvec,lonvec],[30.86,-107.1]);
% % % [r,c]=ind2sub(size(LAT),IDX);
% % % 
% % % delsmap_site=[];
% % % esoil_site=[];
% % % Et_site=[];
% % % precip_site=[];
% % % qbot_site=[];
% % % for i=1:length(time)
% % %     delsmap_site=[delsmap_site;delsmap(r,c,i)];
% % %     esoil_site=[esoil_site;esoil(r,c,i)];
% % %     Et_site=[Et_site;Et(r,c,i)];
% % %     precip_site=[precip_site;precip(r,c,i)];
% % %     qbot_site=[qbot_site;qbot(r,c,i)];
% % % end
% % % 
% % % idx=isnan(esoil_site);
% % % delsmap_site(idx)=[];
% % % esoil_site(idx)=[];
% % % Et_site(idx)=[];
% % % precip_site(idx)=[];
% % % qbot_site(idx)=[];
% % % site_time=time;
% % % site_time(idx)=[];
% % % site_time=site_time+datenum([2015 3 31]);
% % % site_datevec=datevec(double(site_time));
% % % 
% % % figure
% % % subplot(2,2,1)
% % % hold on
% % % p1=plot(site_time,esoil_site,'r-');
% % % p2=plot(site_time,Et_site,'g-');
% % % p3=plot(site_time,precip_site,'b-');
% % % p4=plot(site_time,qbot_site,'m-');
% % % p5=plot(site_time,delsmap_site,'k-');
% % % legend([p1 p2 p3 p4 p5],{'Esoil','Ets','I','qbot','dS'})
% % % xlim([site_time(1) site_time(end)])
% % % datetick('x','mm//dd/yyyy','keepticks','keeplimits')
% % % 
% % % ratio=-delsmap_site./esoil_site;
% % % subplot(2,2,2)
% % % plot(site_time,ratio)
% % % xlim([site_time(1) site_time(end)])
% % % datetick('x','mm//dd/yyyy','keepticks','keeplimits')
% % % ylim([-1 1])
% % % 
% % % %nearby site that doesnt have screened effect:
% % % %(lon=-106.6,lat=30.19)
% % % latvec=LAT(:);
% % % lonvec=LON(:);
% % % IDX=knnsearch([latvec,lonvec],[30.19,-106.6]);
% % % [r,c]=ind2sub(size(LAT),IDX);
% % % 
% % % delsmap_site_full=[];
% % % esoil_site_full=[];
% % % Et_site_full=[];
% % % precip_site_full=[];
% % % qbot_site_full=[];
% % % for i=1:length(time)
% % %     delsmap_site_full=[delsmap_site_full;delsmap(r,c,i)];
% % %     esoil_site_full=[esoil_site_full;esoil(r,c,i)];
% % %     Et_site_full=[Et_site_full;Et(r,c,i)];
% % %     precip_site_full=[precip_site_full;precip(r,c,i)];
% % %     qbot_site_full=[qbot_site_full;qbot(r,c,i)];
% % % end
% % % 
% % % idx=isnan(esoil_site_full);
% % % delsmap_site_full(idx)=[];
% % % esoil_site_full(idx)=[];
% % % Et_site_full(idx)=[];
% % % precip_site_full(idx)=[];
% % % qbot_site_full(idx)=[];
% % % site_full_time=time;
% % % site_full_time(idx)=[];
% % % site_full_time=site_full_time+datenum([2015 3 31]);
% % % site_full_datevec=datevec(double(site_full_time));
% % % idx_screened=find(site_full_datevec(:,2)==1 |site_full_datevec(:,2)==2|site_full_datevec(:,2)==3|site_full_datevec(:,2)==4|site_full_datevec(:,2)==11|site_full_datevec(:,2)==12);
% % % 
% % % subplot(2,2,3)
% % % hold on
% % % p1=plot(site_full_time,esoil_site_full,'r-');
% % % p2=plot(site_full_time,Et_site_full,'g-');
% % % p3=plot(site_full_time,precip_site_full,'b-');
% % % p4=plot(site_full_time,qbot_site_full,'m-');
% % % p5=plot(site_full_time,delsmap_site_full,'k-');
% % % legend([p1 p2 p3 p4 p5],{'Esoil','Ets','I','qbot','dS'})
% % % xlim([site_full_time(1) site_full_time(end)])
% % % datetick('x','mm//dd/yyyy','keepticks','keeplimits')
% % % 
% % % ratio_full=-delsmap_site_full./esoil_site_full;
% % % a=ratio_full;
% % % a(idx_screened,:)=[];
% % % ratio_full_screened=ratio_full(idx_screened);
% % % subplot(2,2,4)
% % % plot(site_full_time,ratio_full)
% % % xlim([site_full_time(1) site_full_time(end)])
% % % datetick('x','mm//dd/yyyy','keepticks','keeplimits')
% % % ylim([-1 1])


%Plot histogram of ESMAP components for each region:
%don't include northeast and northwest because limited data in these regions

%southwest:
idx_SW=find(LAT>= SW_latlim(1) & LAT <= SW_latlim(2) & LON>=SW_lonlim(1) & LON<=SW_lonlim(2));
%midwest:
idx_MW=find(LAT>= MW_latlim(1) & LAT <= MW_latlim(2) & LON>=MW_lonlim(1) & LON<=MW_lonlim(2));
%Great Plains:
idx_GP=find(LAT>= GP_latlim(1) & LAT <= GP_latlim(2) & LON>=GP_lonlim(1) & LON<=GP_lonlim(2));
%Southeast:
idx_SE=find(LAT>= SE_latlim(1) & LAT <= SE_latlim(2) & LON>=SE_lonlim(1) & LON<=SE_lonlim(2));
%Northwest:
idx_NW=find(LAT>= NW_latlim(1) & LAT <= NW_latlim(2) & LON>=NW_lonlim(1) & LON<=NW_lonlim(2));
%Northwest:
idx_NE=find(LAT>= NE_latlim(1) & LAT <= NE_latlim(2) & LON>=NE_lonlim(1) & LON<=NE_lonlim(2));

disp('NW')
a1=nanmedian(dS_Esoil_ratio(idx_NW))
disp('NE')
b1=nanmedian(dS_Esoil_ratio(idx_NE))
disp('MW')
c1=nanmedian(dS_Esoil_ratio(idx_MW))
disp('SW')
d1=nanmedian(dS_Esoil_ratio(idx_SW))
disp('GP')
e1=nanmedian(dS_Esoil_ratio(idx_GP))
disp('SE')
f1=nanmedian(dS_Esoil_ratio(idx_SE))


%initialize outputs:
esoil_SW=[];
esoil_MW=[];
esoil_GP=[];
esoil_SE=[];
esoil_NE=[];
esoil_NW=[];

Et_SW=[];
Et_MW=[];
Et_GP=[];
Et_SE=[];
Et_NE=[];
Et_NW=[];

qbot_SW=[];
qbot_MW=[];
qbot_GP=[];
qbot_SE=[];
qbot_NE=[];
qbot_NW=[];

precip_SW=[];
precip_MW=[];
precip_GP=[];
precip_SE=[];
precip_NE=[];
precip_NW=[];

delsmap_SW=[];
delsmap_MW=[];
delsmap_GP=[];
delsmap_SE=[];
delsmap_NE=[];
delsmap_NW=[];

for t=1:length(time)
    current_esoil=esoil(:,:,t);
    current_esoil=current_esoil(:);
    
    current_delsmap=delsmap(:,:,t);
    current_delsmap=current_delsmap(:);
    
    current_precip=precip(:,:,t);
    current_precip=current_precip(:);
    
    current_Et=Et(:,:,t);
    current_Et=current_Et(:);
    
    current_qbot=qbot(:,:,t);
    current_qbot=current_qbot(:);
    
    esoil_SW=[esoil_SW,current_esoil(idx_SW)];
    esoil_MW=[esoil_MW,current_esoil(idx_MW)];
    esoil_GP=[esoil_GP,current_esoil(idx_GP)];
    esoil_SE=[esoil_SE,current_esoil(idx_SE)];
    esoil_NE=[esoil_NE,current_esoil(idx_NE)];
    esoil_NW=[esoil_NW,current_esoil(idx_NW)];
    
    delsmap_SW=[delsmap_SW,current_delsmap(idx_SW)];
    delsmap_MW=[delsmap_MW,current_delsmap(idx_MW)];
    delsmap_GP=[delsmap_GP,current_delsmap(idx_GP)];
    delsmap_SE=[delsmap_SE,current_delsmap(idx_SE)];
    delsmap_NE=[delsmap_NE,current_delsmap(idx_NE)];
    delsmap_NW=[delsmap_NW,current_delsmap(idx_NW)];
    
    Et_SW=[Et_SW,current_Et(idx_SW)];
    Et_MW=[Et_MW,current_Et(idx_MW)];
    Et_GP=[Et_GP,current_Et(idx_GP)];
    Et_SE=[Et_SE,current_Et(idx_SE)];
    Et_NE=[Et_NE,current_Et(idx_NE)];
    Et_NW=[Et_NW,current_Et(idx_NW)];
    
    precip_SW=[precip_SW,current_precip(idx_SW)];
    precip_MW=[precip_MW,current_precip(idx_MW)];
    precip_GP=[precip_GP,current_precip(idx_GP)];
    precip_SE=[precip_SE,current_precip(idx_SE)];
    precip_NE=[precip_NE,current_precip(idx_NE)];
    precip_NW=[precip_NW,current_precip(idx_NW)];
    
    qbot_SW=[qbot_SW,current_qbot(idx_SW)];
    qbot_MW=[qbot_MW,current_qbot(idx_MW)];
    qbot_GP=[qbot_GP,current_qbot(idx_GP)];
    qbot_SE=[qbot_SE,current_qbot(idx_SE)];
    qbot_NE=[qbot_NE,current_qbot(idx_NE)];
    qbot_NW=[qbot_NW,current_qbot(idx_NW)];
end
esoil_SW=nanmean(esoil_SW');
esoil_MW=nanmean(esoil_MW');
esoil_GP=nanmean(esoil_GP');
esoil_SE=nanmean(esoil_SE');
esoil_NE=nanmean(esoil_NE');
esoil_NW=nanmean(esoil_NW');

delsmap_SW=nanmean(delsmap_SW');
delsmap_MW=nanmean(delsmap_MW');
delsmap_GP=nanmean(delsmap_GP');
delsmap_SE=nanmean(delsmap_SE');
delsmap_NE=nanmean(delsmap_NE');
delsmap_NW=nanmean(delsmap_NW');

Et_SW=nanmean(Et_SW');
Et_MW=nanmean(Et_MW');
Et_GP=nanmean(Et_GP');
Et_SE=nanmean(Et_SE');
Et_NE=nanmean(Et_NE');
Et_NW=nanmean(Et_NW');

precip_SW=nanmean(precip_SW');
precip_MW=nanmean(precip_MW');
precip_GP=nanmean(precip_GP');
precip_SE=nanmean(precip_SE');
precip_NE=nanmean(precip_NE');
precip_NW=nanmean(precip_NW');

qbot_SW=nanmean(qbot_SW');
qbot_MW=nanmean(qbot_MW');
qbot_GP=nanmean(qbot_GP');
qbot_SE=nanmean(qbot_SE');
qbot_NE=nanmean(qbot_NE');
qbot_NW=nanmean(qbot_NW');

%Plot for southwest region:
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_SW);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_SW);
[f_precip,xi_precip]=ksdensity(precip_SW);
[f_Et,xi_Et]=ksdensity(Et_SW);
[f_qbot,xi_qbot]=ksdensity(qbot_SW);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
% legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Southwest','fontsize',25)
xlim([-1.5 1.5])

%Midwest
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_MW);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_MW);
[f_precip,xi_precip]=ksdensity(precip_MW);
[f_Et,xi_Et]=ksdensity(Et_MW);
[f_qbot,xi_qbot]=ksdensity(qbot_MW);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
% legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Midwest','fontsize',25)
xlim([-1.5 1.5])

%Great Plains
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_GP);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_GP);
[f_precip,xi_precip]=ksdensity(precip_GP);
[f_Et,xi_Et]=ksdensity(Et_GP);
[f_qbot,xi_qbot]=ksdensity(qbot_GP);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
% legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Great Plains','fontsize',25)
xlim([-1.5 1.5])

%Southeast
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_SE);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_SE);
[f_precip,xi_precip]=ksdensity(precip_SE);
[f_Et,xi_Et]=ksdensity(Et_SE);
[f_qbot,xi_qbot]=ksdensity(qbot_SE);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
% legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Southeast','fontsize',25)
xlim([-1.5 1.5])


%Northeast
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_NE);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_NE);
[f_precip,xi_precip]=ksdensity(precip_NE);
[f_Et,xi_Et]=ksdensity(Et_NE);
[f_qbot,xi_qbot]=ksdensity(qbot_NE);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
% legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Northeast','fontsize',25)
xlim([-1.5 1.5])

%Northwest
figure
hold on
[f_esoil,xi_esoil]=ksdensity(esoil_NW);
[f_delsmap,xi_delsmap]=ksdensity(delsmap_NW);
[f_precip,xi_precip]=ksdensity(precip_NW);
[f_Et,xi_Et]=ksdensity(Et_NW);
[f_qbot,xi_qbot]=ksdensity(qbot_NW);

p1=plot(xi_esoil,f_esoil,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_delsmap,f_delsmap,'-','color',[0.75 0 0.75],'linewidth',3);
p3=plot(xi_precip,f_precip,'-r','linewidth',3);
p4=plot(xi_Et,f_Et,'-','color',[0 0.5 0],'linewidth',3);
p5=plot(xi_qbot,f_qbot,'-','color',[0 0.75 0.75],'linewidth',3);

a1=area(xi_esoil,f_esoil,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_delsmap,f_delsmap,'FaceColor',[0.75 0 0.75],'FaceAlpha',0.2);
a3=area(xi_precip,f_precip,'FaceColor','r','FaceAlpha',0.2);
a4=area(xi_Et,f_Et,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a5=area(xi_qbot,f_qbot,'FaceColor',[0 0.75 0.75],'FaceAlpha',0.2);

xlabel('mm day^{-1}','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
legend([p1 p2 p3 p4 p5],{'E_{soil}','SMAP Drying','I','E_{T}_{s}','q_{bot}'},'fontsize',20)
title('Northwest','fontsize',25)
xlim([-1.5 1.5])

sprintf('ratio of dS/Esoil southwest is: %.2f',nanmedian(delsmap_SW./esoil_SW))
sprintf('ratio of dS/Esoil midwest is: %.2f',nanmedian(delsmap_MW./esoil_MW))
sprintf('ratio of dS/Esoil northwest is: %.2f',nanmedian(delsmap_NW./esoil_NW))
sprintf('ratio of dS/Esoil great plains is: %.2f',nanmedian(delsmap_GP./esoil_GP))
sprintf('ratio of dS/Esoil southeast is: %.2f',nanmedian(delsmap_SE./esoil_SE))
sprintf('ratio of dS/Esoil northeast is: %.2f',nanmedian(delsmap_NE./esoil_NE))

sprintf('ratio of qbot+I/Et southwest is: %.2f',nanmedian( ((-qbot_SW)+(precip_SW) )./Et_SW))
sprintf('ratio of qbot+I/Et midwest is: %.2f',nanmedian(((-qbot_MW)+(precip_MW) )./Et_MW))
sprintf('ratio of qbot+I/Et northwest is: %.2f',nanmedian(((-qbot_NW)+(precip_NW) )./Et_NW))
sprintf('ratio of qbot+I/Et great plains is: %.2f',nanmedian(((-qbot_GP)+(precip_GP) )./Et_GP))
sprintf('ratio of qbot+I/Et southeast is: %.2f',nanmedian(((-qbot_SE)+(precip_SE) )./Et_SE))
sprintf('ratio of qbot+I/Et northeast is: %.2f',nanmedian(((-qbot_NE)+(precip_NE) )./Et_NE))

disp('NW')
a1=nanmedian(delsmap_NW);
a2=nanmedian(delsmap_NW./esoil_NW);
disp('NE')
b1=nanmedian(delsmap_NE);
b2=nanmedian(delsmap_NE./esoil_NE);
disp('MW')
c1=nanmedian(delsmap_MW);
c2=nanmedian(delsmap_MW./esoil_MW);
disp('SW')
d1=nanmedian(delsmap_SW);
d2=nanmedian(delsmap_SW./esoil_SW);
disp('GP')
e1=nanmedian(delsmap_GP);
e2=nanmedian(delsmap_GP./esoil_GP);
disp('SE')
f1=nanmedian(delsmap_SE);
f2=nanmedian(delsmap_SE./esoil_SE);

ds=[a1;b1;c1;d1;e1;f1];
ratio=[a2;b2;c2;d2;e2;f2];

disp('correlation between drying rate and the ratio of ds/esoil:')
corr(ds,ratio,'Type','Pearson')

