clc;clear all;close all;

%Load in all of the points that required the final run of hydrus to alter
%the original soil params:
Final_Altered_Points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Final_GOOD_Points_2');

%Load in file of points Altered in the spin up simulation that did the
%secondary soil type soil matching:
NN_Altered_Points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Secondary_Altered_Points');

%Load in file of points Altered in the spin up simulation that did not have a
%secondary soil type  - so matched with the next closest soil texture classification:
no_NN_Altered_Points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/No_Secondary_Altered_Points');

%define the list of final altered points that match with NN:
[ia,ib]=ismember(Final_Altered_Points,NN_Altered_Points,'rows');
Final_Altered_Points_NN=Final_Altered_Points(ia,:);

%define the list of final altered points that match with closest classification:
[ia,ib]=ismember(Final_Altered_Points,no_NN_Altered_Points,'rows');
Final_Altered_Points_no_NN=Final_Altered_Points(ia,:);

%export these files:
outdir='/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Soil_Param_Flags/';

dlmwrite([outdir,'Secodnary_params'],Final_Altered_Points_NN,'delimiter',' ','precision','%.15g');
dlmwrite([outdir,'no_Secondary_params'],Final_Altered_Points_no_NN,'delimiter',' ','precision','%.15g');


%read in the points that sill have not converged:
Bad_Points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/TMP_BAD_Points2');

store=[];
for i=1:length(Bad_Points)
    i
    lat=Bad_Points(i,1);
    lon=Bad_Points(i,2);
    filename=sprintf('/Volumes/REESEN/SMAP/Gridded_ESMAP/%.15g/%.15g/Soil_updated',lat,lon);
    st=importdata(filename);
    store=[store;st];
end