clc;clear all;close all;

%begin by looping through bad points:
% bad_points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/ESMAP_Calc/Esoil_bad_points');
% other_bad_points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Final_Run_Check/Final_BAD_Points_2');
% Points_save=[bad_points;other_bad_points];
% 
% dlmwrite('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Altered_Points',Points_save,'delimiter',' ','precision','%.15g');

Points=importdata('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/from_summit/Point_Check/Final_BAD_Points_2');

npoints=length(Points(:,1));
outdir='/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/for_summit/9km_points_Final/';
CMD=['rm /Volumes/REESEN/SMAP/ESMAP_for_Ronnie/for_summit/9km_points_Final/*iter*'];
system(CMD);

nruns=24;
start_idx=1;
for prun=1:nruns
    end_idx=start_idx+npoints/nruns;
    end_idx=round(end_idx);

    if end_idx>length(Points(:,1))
        end_idx=length(Points(:,1));
    end
    current_points=Points(start_idx:end_idx,:);

    outfilename=sprintf('SMAP_9km_Points_iter%03d',prun);

    dlmwrite([outdir,outfilename],current_points,'delimiter',' ','precision','%.15g');

    %write points:
    start_idx=end_idx+1;
    
    if end_idx==length(Points(:,1))
        break
    end
end

%confirm all points are accounted for:
store_data=[];
for prun=1:nruns
    current_filename=sprintf('/Volumes/REESEN/SMAP/ESMAP_for_Ronnie/for_summit/9km_points_Final/SMAP_9km_Points_iter%03d',prun);
    if (exist(current_filename,'file')>0)
    data=importdata(current_filename);
    store_data=[store_data;data];
    end
end

assert(length(store_data)==length(Points),'points were left out!!')
