clc;clear all;close all;

%%====================================================================================

% % %Load in ESMAP grid:
% % lat=ncread('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip.nc','lat');
% % lon=ncread('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip.nc','lon');
% % 
% % %load in valid points
% % Points=importdata('/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_Points');
% % Points=round(Points,5);
% % 
% % npoints=length(Points);
% % %create grid to store outputs on:
% % [LAT,LON]=meshgrid(lat,lon);
% % [nrow ncol]=size(LAT);
% % LAT_vec=round(reshape(LAT,nrow*ncol,1),5);
% % LON_vec=round(reshape(LON,nrow*ncol,1),5);
% % [ia,ib]=ismember([LAT_vec,LON_vec],Points,'rows');
% % IDX=find(ia==1);
% % assert(length(IDX)==length(Points),'points missing');
% % 
% % validation_sets={'Mosaic','Noah','GLEAM'};
% % %initialize output:
% % Mean_Esoil.Mosaic=nan(nrow,ncol);
% % Mean_Esoil.Noah=nan(nrow,ncol);
% % Mean_Esoil.GLEAM=nan(nrow,ncol);
% % 
% % for j=1:length(validation_sets)
% %     validation_set=validation_sets{j};
% %     for i=1:npoints
% %         [j,i]
% %         lat=Points(i,1);
% %         lon=Points(i,2)+360;
% %         filename_Esoil=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP_blank/%.15g/%.15g/%s_Esoil.csv',lat,lon,validation_set);
% %         if exist(filename_Esoil,'file')>0 && i~=23459 && i~=28977 && i~=32918 && i~=34232 && i~=40511 && i~=62119
% %             Mosaic_Esoil=csvread(filename_Esoil);
% %             mean_Esoil=nanmean(Mosaic_Esoil(:,4));
% %             %store the mean Esoil
% %             if j==1
% %                 Mean_Esoil.Mosaic(IDX(i))=mean_Esoil;
% %             elseif j==2
% %                 Mean_Esoil.Noah(IDX(i))=mean_Esoil;
% %             elseif j==3
% %                 Mean_Esoil.GLEAM(IDX(i))=mean_Esoil;
% %             end
% %         end
% %     end
% % end
% % save('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mean_Gridded_Validation_Esoil.mat','Mean_Esoil');

%Load in ESMAP data:
ESMAP_esoil=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','esoil_screened');
ESMAP_esoil=double(ESMAP_esoil);
%calculate mean ESMAP esoil
ESMAP_esoil_mean_grid=nanmean(ESMAP_esoil,3);
%Load in mean gridded valdiation data:
Validation_data=load('/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mean_Gridded_Validation_Esoil.mat');
Mosaic_Esoil=Validation_data.Mean_Esoil.Mosaic;
Mosaic_Esoil=Mosaic_Esoil;
Noah_Esoil=Validation_data.Mean_Esoil.Noah;
GLEAM_Esoil=Validation_data.Mean_Esoil.GLEAM;

%define color bar titles:
TITLES={'E-SMAP - Mosaic E_{soil} mm day^{-1}','E-SMAP - Noah E_{soil} mm day^{-1}','E-SMAP - GLEAM E_{soil} mm day^{-1}'};
%Plot bias figure between each:
for i=1:3
    
    if i==1
        PLOT_DATA=ESMAP_esoil_mean_grid-Mosaic_Esoil;
    elseif i==2
        PLOT_DATA=ESMAP_esoil_mean_grid-Noah_Esoil;
    elseif i==3
        PLOT_DATA=ESMAP_esoil_mean_grid-GLEAM_Esoil;
    end
    Title=TITLES{i};
    %define bounding domain box
    latlim=[24.9,50.1];
    lonlim=[-125.1,-66.9];
    
    %%=======================================================================================================================
    %create grid to plot on:
    filename_ET='/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/Monthly/NLDAS_MOSAIC0125_H.A201708.nc';
    lat=ncread(filename_ET,'lat');
    lon=ncread(filename_ET,'lon');
    [LAT,LON]=meshgrid(lat,lon);
    
        valid_min=-1.3;
        valid_max=1.3;
        
    % % Create a set of level ranges to be used in converting the data to a
    % % geolocated image that has a color assigned to each range.
    levels = linspace(valid_min, valid_max, 7);
    cmap = redbluecmap(7);
    cmap=flipud(cmap);
    cmap(4,:)=[0.8947    0.9474    0.4000];
    
    % % % % % Create a color map.
    % Convert the data to an geolcated image by setting a color for each level
    % range.
    PLOT_DATA=double(PLOT_DATA);
    Z = PLOT_DATA;
    
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
    geoshow(LAT,LON,Z,'DisplayType','surface')
    %Display outline of land area
    landareas = shaperead('landareas.shp', 'UseGeo', true);
    coast.lat = [landareas.Lat];
    coast.long = [landareas.Lon];
    geoshow(coast.lat, coast.long, 'Color', 'k','linewidth',3)
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
    h = lcolorbar(clevels, 'Location', 'horizontal','fontsize',40,'title',Title);
    idx=3:2:length(h.XTickLabel);
    set(gca,'fontsize',40)
    h.Position(2)=h.Position(2)-0.06;
    %defin ocrrelation with longitude:
    a=LON(:);
    b=PLOT_DATA(:);
    idx=isnan(a);
    a(idx)=[];
    b(idx)=[];
    idx=isnan(b);
    b(idx)=[];
    a(idx)=[];
    [R,P]=corrcoef(a,b);
    sprintf('R2 between difference and lon is %.5f which has a p-value of %.5f  for %s',R(2)^2,P(2),Title)
    
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
end

%Include a KDE comparisson:
figure
hold on
ESMAP_esoil_mean_grid=ESMAP_esoil_mean_grid(:);
[f_esmap,xi_esmap]=ksdensity(ESMAP_esoil_mean_grid);
Mosaic_Esoil=Mosaic_Esoil(:);
[f_mosaic,xi_mosaic]=ksdensity(Mosaic_Esoil);
Noah_Esoil=Noah_Esoil(:);
[f_noah,xi_noah]=ksdensity(Noah_Esoil);
GLEAM_Esoil=GLEAM_Esoil(:);
[f_gleam,xi_gleam]=ksdensity(GLEAM_Esoil);

p1=plot(xi_esmap,f_esmap,'-','color',[0.25 0.25 0.25],'linewidth',3);
p2=plot(xi_mosaic,f_mosaic,'-','color',[0.75 0.75 0.0],'linewidth',3);
p3=plot(xi_noah,f_noah,'-','color',[0 0.5 0],'linewidth',3);
p4=plot(xi_gleam,f_gleam,'-r','linewidth',3);

a1=area(xi_esmap,f_esmap,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.2);
a2=area(xi_mosaic,f_mosaic,'FaceColor',[0.75 0.75 0.0],'FaceAlpha',0.2);
a3=area(xi_noah,f_noah,'FaceColor',[0 0.5 0],'FaceAlpha',0.2);
a4=area(xi_gleam,f_gleam,'FaceColor','r','FaceAlpha',0.2);

%include vertical lines for means:
xline(nanmedian(ESMAP_esoil_mean_grid),'--','linewidth',4,'color',[0.25 0.25 0.25]);
xline(nanmedian(Mosaic_Esoil),'--','linewidth',4,'color',[0.75 0.75 0.0]);
xline(nanmedian(Noah_Esoil),'--','linewidth',4,'color',[0 0.5 0]);
xline(nanmedian(GLEAM_Esoil),'--','linewidth',4,'color','r');

xlabel('E_{soil} mm/day','fontsize',24)
ylabel('Probability','fontsize',24)
set(gca,'fontsize',24)
legend([p1 p2 p3 p4],{'E-SMAP','Mosaic','Noah','GLEAM'},'fontsize',20)
xlim([0 2.5])


%Load in ESMAP data:
ESMAP_esoil=ncread('/Volumes/REESEN/SMAP/Data_Repository/Evaporation_SMAP.nc','esoil_screened');
ESMAP_esoil=double(ESMAP_esoil);
%calculate mean ESMAP esoil
ESMAP_esoil_mean_grid=nanmean(ESMAP_esoil,3);

PLOT_DATA=ESMAP_esoil_mean_grid;

%define bounding domain box
latlim=[24.9,50.1];
lonlim=[-125.1,-66.9];

valid_min=0.1;
valid_max=1.3;
% % Create a set of level ranges to be used in converting the data to a
% % geolocated image that has a color assigned to each range.
levels = linspace(valid_min, valid_max, 7);
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
h = lcolorbar(clevels, 'Location', 'horizontal','fontsize',34,'title','E-SMAP E_{soil} mm day^{-1}');
idx=3:2:length(h.XTickLabel);
set(gca,'fontsize',34)
h.Position(2)=h.Position(2)-0.06;

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
        