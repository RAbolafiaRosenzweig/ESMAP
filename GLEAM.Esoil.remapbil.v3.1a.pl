#!/usr/bin/perl
#

#use strict;
use warnings;

$startyear=2015;
$endyear=2018;
#this will match data to PGF 0.5-degree grid:
$mask_file="/Volumes/REESEN/SMAP/Gridded_ncdf_Products/Final_Data/ESMAP_QC.nc";

for($year=$startyear; $year<=$endyear; $year++) {
    
  $IN_FILE=sprintf("/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Raw/%04d/Eb_%04d_GLEAM_v3.3a.nc",$year,$year);
  print "$IN_FILE\n";
            $outdir="/Volumes/REESEN/SMAP/Validation_Data/GLEAM_9km/Esoil";
   print "$outdir\n";


      unless(-e $outdir) {
                $cmd="mkdir -p $outdir";
                print "$cmd\n";
                system($cmd);
            }
      $Out_File_Weights=sprintf("Weights.nc");
            $cmd="cdo genbil,$mask_file $IN_FILE $outdir/$Out_File_Weights";
      print "$cmd\n";
      system($cmd);
            
      $Out_File=sprintf("Eb_9km_%04d_GLEAM_v3.3a.nc",$year);
      $cmd="cdo remap,$mask_file,$outdir/$Out_File_Weights $IN_FILE $outdir/$Out_File";
      print "$cmd\n";
  system($cmd);

  
      }
      


