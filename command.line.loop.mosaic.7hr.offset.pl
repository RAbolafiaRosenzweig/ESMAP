#!/usr/bin/perl
#
# To run, type:
#which perl
#echo $PATH
#perl


# <full path>/command.line.loop.pl
#

#$series1=shift; #PIXEL OUTPUT FIL
#$series2=shift; #PIXEL OUTPUT FIL
#$series3=shift; #PIXEL OUTPUT FIL

$startyear=2018;
$endyear=2018;
$startmonth=1;

@days=(0,31,28,31,30,31,30,31,31,30,31,30,31);

for($year=$startyear;$year<=$endyear;$year++) {

    for ($month=1;$month<=1;$month++) {
        
        $daysinmonth=$days[$month];
        print "$daysinmonth\n";
        
        if ($month==2 && $year==2016) {
            $daysinmonth=29;
        }
        
        $outdir=sprintf("/Volumes/REESEN/SMAP/Validation_Data/NLDAS2/MOSAIC_7hr_Offset/%04d/%02d/",$year,$month);
        
        unless(-e $outdir) {
            $cmd="mkdir -p $outdir";
            print "$cmd\n";
            system($cmd);
        }

        for ($day=1;$day<=2;$day++) {
            @files=();
            
            #if it is atleast the thrid day of the month
            if ($day>1) {
                $filename=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.0?00.nc",$year,$month,$year,$month,$day);
                $filename_01=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[012]00.nc",$year,$month,$year,$month,$day);
                $filename_02=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[3456789]00.nc",$year,$month,$year,$month,$day-1);
                $filename_03=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.2?00.nc",$year,$month,$year,$month,$day-1);
                
                @files=($filename,$filename_01,$filename_02,$filename_03);
            }
            #if it's the second day of the month
            if ($day==1&& $month>1) {
                $month_before=$month-1;
                $day_before=$days[$month_before];
                $filename=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.0?00.nc",$year,$month,$year,$month,$day);
                $filename_01=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[012]00.nc",$year,$month,$year,$month,$day);
                $filename_02=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[3456789]00.nc",$year,$month_before,$year,$month_before,$day_before);
                $filename_03=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.2?00.nc",$year,$month_before,$year,$month_before,$day_before);
                
                @files=($filename,$filename_01,$filename_02,$filename_03);
            }
            
            #if it's the second day of the month
            if ($day==1 && $month==1) {
                $month_before=12;
                $day_before=31;
                $year_before=$year-1;
                $filename=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.0?00.nc",$year,$month,$year,$month,$day);
                $filename_01=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[012]00.nc",$year,$month,$year,$month,$day);
                $filename_02=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.1[3456789]00.nc",$year_before,$month_before,$year_before,$month_before,$day_before);
                $filename_03=sprintf(" /Volumes/REESEN/SMAP/Validation_Data/NLDAS2/Mosaic/%04d/%02d/NLDAS_MOSAIC0125_H..A%04d%02d%02d.2?00.nc",$year_before,$month_before,$year_before,$month_before,$day_before);
                @files=($filename,$filename_01,$filename_02,$filename_03);
            }
        
            $cmd=sprintf("ncrcat -O @files $outdir/NLDAS_MOSAIC0125_H..A%04d%02d%02d.daily.nc",$year,$month,$day);
            print "$cmd\n";
            system($cmd);
            
            $cmd=sprintf("ncra -O $outdir/NLDAS_MOSAIC0125_H..A%04d%02d%02d.daily.nc $outdir/NLDAS_MOSAIC0125_H..A%04d%02d%02d.nc",$year,$month,$day,$year,$month,$day);
            print "$cmd\n";
            system($cmd);
        }
    }
}


