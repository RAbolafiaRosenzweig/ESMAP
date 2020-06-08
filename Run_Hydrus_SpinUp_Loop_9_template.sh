#!/bin/bash

now=$(date +"%T")
echo "Current time : $now"

cd /scratch/summit/roab9675/SMAP/ESMAP_Calc/

#99999 is fill value for the point file - to be filled by MATLAB script
FILE=/projects/roab9675/SMAP/9km_points/SMAP_9km_Points_iter99999

#read -p "Press enter to continue"

count=0
while read line <&3 ; do
	
	stringpoint=( $line )
	stringlat=${stringpoint[0]}
	stringlon=${stringpoint[1]}
#check if the rounding should be 3 or 4 decimal places:

	((count++))
#   if [ $count -gt 50 ]; then
#       ((tmpcount=count-50))
#       rm /scratch/summit/roab9675/SMAP/ESMAP_Calc/hydrus_${stringlat}_${stringlon}.sh
#fi
Directory=/scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
if [ -d $Directory ]; then
rm -r /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
fi

    mkdir /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
	cp /scratch/summit/roab9675/SMAP/Gridded_ESMAP/${stringlat}/${stringlon}/Hydrus_SpinUp/*.IN /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
	cp /scratch/summit/roab9675/SMAP/Gridded_ESMAP/${stringlat}/${stringlon}/Hydrus_SpinUp/*.DAT /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
	cp /scratch/summit/roab9675/SMAP/Gridded_ESMAP/${stringlat}/${stringlon}/Hydrus_SpinUp/*.TXT /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
	cp /projects/roab9675/hydrus/hydrus_original/build/h1d_calc /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
#changes for case sensitivity:
mv /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/PROFILE.DAT /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/Profile.dat
mv /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/SELECTOR.IN /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/Selector.in
#create shell script:
	echo "/scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}" > /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/LEVEL_01.DIR
	echo "\#\!/bin/bash" >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh
	echo "cd /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}"  >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh
	echo "yes | ./h1d_calc"  >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh
	echo "cd /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/"  >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh
	echo "cp /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/*.OUT /scratch/summit/roab9675/SMAP/Gridded_ESMAP/${stringlat}/${stringlon}/Hydrus_SpinUp/"  >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh

#echo "mv /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/nohup_${stringlat}_${stringlon} /scratch/summit/roab9675/SMAP/Gridded_ESMAP/${stringlat}/${stringlon}/Hydrus_SpinUp/nohup_${stringlat}_${stringlon}"  >> /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/hydrus_${stringlat}_${stringlon}.sh

    cd /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
	chmod +x hydrus_${stringlat}_${stringlon}.sh
    ./hydrus_${stringlat}_${stringlon}.sh 2> /dev/null
    cd /scratch/summit/roab9675/SMAP/ESMAP_Calc/
    rm -r /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}
#submit=1
#	while [ $submit == 1 ] ; do
#		RUNS=$(ps | grep "hydrus_" | wc -l)
#		echo $RUNS
#		if [ "$RUNS" -lt "9" ]; then
#			nohup ./hydrus_${stringlat}_${stringlon}.sh > /scratch/summit/roab9675/SMAP/ESMAP_Calc/tmp_${stringlat}_${stringlon}/nohup_${stringlat}_${stringlon} &
#			sleep 2
#           ./hydrus_${stringlat}_${stringlon}.sh
#			((submit++))
#		else
#			echo "Too many sims running! - SLEEP 10 SECONDS"
#			sleep 10
#		fi
#	done	
	#RUNS='ps | grep "hydrus_" | wc -l'
	#if [ "$RUNS" -lt "7" ]

done 3< "$FILE"

