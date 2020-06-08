#!/usr/bin/perl
#

#use strict;
use warnings;

$startyear=2018;
$endyear=2019;
#this will match data to PGF 0.5-degree grid:
$mask_file="/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC_withPrecip.nc";

@days=(0,31,28,31,30,31,30,31,31,30,31,30,31);

for($year=$startyear; $year<=$endyear; $year++) {
    if ($year==2018) {
            $startmonth=1;
	  } else {
	    $startmonth=1;
	  }
    if ($year==2019){
      $endmonth=3;
    } else{
      $endmonth=12;
    }
  for ($month=$startmonth; $month<=$endmonth; $month++){
    $daysinmonth=$days[$month];
    
        print "$daysinmonth\n";
        
        if ($month==2 && $year==2016) {
            $daysinmonth=29;
	  }
	    #create output directory:
	    $outdir=sprintf("/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Noah/%04d/%02d/",$year,$month);
	 unless(-e $outdir) {
            $cmd="mkdir -p $outdir";
            print "$cmd\n";
            system($cmd);
	  }

    
	    for ($day=1;$day<=$daysinmonth;$day++) {
	      
	      for ($hour=0;$hour<=23;$hour++) {
		
	print "$year $month $day $hour\n";
		
  $IN_FILE=sprintf("/Volumes/LabShare/NLDAS2/NOAH/%04d/%02d/NLDAS_NOAH0125_H.A%04d%02d%02d.%02d00.nc",$year,$month,$year,$month,$day,$hour);



  #copy over file to new directory and remove un-needed variables          
  $Out_File=sprintf("NLDAS_NOAH0125_H..A%04d%02d%02d.%02d00.nc",$year,$month,$day,$hour);
  $cmd="cp $IN_FILE $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);
	
#perform remap:
        $Out_File_Weights=sprintf("Weights.nc");
      $cmd="cdo genbil,$mask_file $IN_FILE $outdir/$Out_File_Weights";
      print "$cmd\n";
      system($cmd);
            
      $cmd="cdo remap,$mask_file,$outdir/$Out_File_Weights $IN_FILE $outdir/$Out_File";
      print "$cmd\n";
	system($cmd);


	#remove unwanted variables:
  $cmd="ncks -C -O -x -v depth $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_2 $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_2_bnds $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_3 $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_3_bnds $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_4 $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_4_bnds $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v depth_bnds $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var111_NSWRS $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var112_NLWRS $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var121_LHTFL $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var122_SHTFK $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var145_PEVPR $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var148_AVSFT $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var151_LSOIL $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var155_GFLUX $outdir/$Out_File $outdir/$Out_File";
  system($cmd);

  $cmd="ncks -C -O -x -v var161_ASNOW $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var162_ARAIN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var179_ACOND $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var181_CCOND $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var182_LAI $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var198_SBSNO $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var203_RSMIN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var207_MSTAV $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var223_CNWAT $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var229_SNOHF $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var234_BGRUN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var235_SSRUN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var238_SNOWC $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var246_RCS $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var247_RCT $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var248_RCQ $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var249_RCSOL $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var250_UNKNOWN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

  $cmd="ncks -C -O -x -v var255_UNKNOWN $outdir/$Out_File $outdir/$Out_File";
  print "$cmd\n";
  system($cmd);

	$IN_FILE=$Out_File;
	$Out_File=sprintf("NLDAS_NOAH0125_H..A%04d%02d%02d.%02d00.nc",$year,$month,$day,$hour);
        }
     }
   }  
}
      


