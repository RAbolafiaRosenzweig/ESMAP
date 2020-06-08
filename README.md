# ESMAP
A soil evaporation data set based on SMAP observed drying rates (ESMAP; Evaporation–Soil Moisture Active Passive). Journal articles: https://www.mdpi.com/2072-4292/10/12/1945


The below scripts were used to create the ESMAP soil evaporation product

Screen_SMAP_QC_Flag.m - Turns all SMAP retrievals that have an uncertain quality flag (flags other than 0 or 8) and turns the retrievals to the fill value (-9999). Writes new netcdf file, identical to originals.

Get_SMAP_QC.r - Extracts SMAP observations at each grid point and creates a directory structure for gridded ESMAPSMAP file columns. These SMAP obs have  been screened by the quality flag (Screen_SMAP_QC_Flag.m).
1 - Year 
2 - Month 
3 - Day 
4 - SMAP Observation

Get_EVI.r - Extracts EVI observations at each SMAP grid point 
SMAP file columns
1 - Year 
2 - Month 
3 - Day 
4 - EVI Observation
** Note that EVI needs to be divided by 1000 to get an LAI value

Get_Forcing.r - Extracts forcing file information at each SMAP grid pont
Forcing file columns
1 - Year 
2 - Month 
3 - Day 
4 - Hour 
5 - Pressure, Pa
6 - Temperature, K
7 - Specific Humidity, kg/kg
8 - Precipitation, kg/m^2 
** Note that NLDAS time is in UTC, while SMAP is at 6am local time
** Initial file took a while to process.   Suggest changing the dates to be for new time period wanted and deleting the file.remove statement

Get_NOAH.r - Extracts NOAH LSM information at each SMAP grid pont
Forcing file columns
1 - Year 
2 - Month 
3 - Day 
4 - Hour 
5 - Net Shortwave, W/m^2
6 - Net Longwave, W/m^2
** Note that NLDAS time is in UTC, while SMAP is at 6am local time
** Initial file took a while to process.   Suggest changing the dates to be for new time period wanted and deleting the file.remove statement

Get_Vegetation.r - Extracts UMD vegetation classification at each SMAP grid point
No columns, just one number
0	Water 
1	Evergreen Needleleaf Forest	
2	Evergreen Broadleaf Forest	
3	Deciduous Needleleaf Forest	
4	Deciduous Broadleaf Forest	
5	Mixed Forest	
6	Woodland	
7	Wooded Grassland	
8	Closed Shrubland	
9	Open Shrubland	
10	Grassland	
11	Cropland	
12	Bare Ground	
13	Urban and Built 
** Note that Water and Urban SHOULD NOT APPEAR IN FILES GENERATED
** Source: http://glcf.umd.edu/data/landcover/

Get_Soil.r - Extracts predominant soil classification type at each SMAP grid point
Reads NLDAS_Soil_Class.txt
Column 1	X Coordinate Index
Column 2	Y Coordinate Index
Column 3	Longitude (center of 1/8th-degree grid boxes)
Column 4	Latitude (center of 1/8th-degree grid boxes)
Columns 5-21	Number of Occurrences of Soil Classes 1-16 in Each 1/8th-Grid Box
Soil Types:
1	Sand
2	Loamy sand
3	Sandy loam
4	Silt loam
5	Silt
6	Loam
7	Sandy clay loam
8	Silty clay loam
9	Clay loam
10	Sandy clay
11	Silty clay
12	Clay
13	Organic materials
14	Water
15	Bedrock
16	Other
Sources: https://ldas.gsfc.nasa.gov/nldas/NLDASsoils.php
https://ldas.gsfc.nasa.gov/nldas/asc/soils/STEX_TAB.01
https://ral.ucar.edu/sites/default/files/public/product-tool/noah-multiparameterization-land-surface-model-noah-mp-lsm/SOILPARM.TBL_.txt
** Note that we are excluding water as a valid option
** Commented out option in code to produce a file will all soil classes present in gridbox

Calc_Roots.r - R Function that calculates root distribution for transpiration purposes, called by Calc_Mu_Transpiration.r
Uses parameters from CLM4 Technical Note
Cites Zeng (2001) 
http://www.cesm.ucar.edu/models/ccsm4.0/clm/CLM4_Tech_Note.pdf
Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530

Calc_Mu_Transpiration.r - Calculates the Mu transpiration quantities at each of the SMAP points
Uses the output from: 
	Get_SMAP.r
	Get_EVI.r
	Get_Forcing.r
	Get_NOAH.r
	Get_Vegetation.r
	Get_Soil.r
Calls Calc_Roots.r to get root density of surface layer
Due to NLDAS and SMAP on different timings (NLDAS is UTC and SMAP is 6am local time), use an offset to get NLDAS on SMAP time
	Longitude UTC Conversion
	Long_1	Long_2	Offset
	-67.5	-82.5	5
	-82.5	-97.5	6
	-97.5	-112.5	7
	-112.5	-127.5	8	
Outputs a file of transpiration estimates, columns are:
1 - Year 
2 - Month 
3 - Day 
4 - PM (Mu eqn. 1), mm
5 - PM (Mu eqn. 22), mm
6 - Root restriction (applied to Mu eqn. 22, column 5), mm
7 - SM restriction (applied to the root restriction ET, column 6), mm

Soil_Constants_Hydrus.txt
Source for Van Genuchten parameters:
	https://www.ars.usda.gov/ARSUserFiles/80420525/EnvironmentalTransport/CalcPTFFiles/PTF_Manual.version_3.0.pdf
Source for soil texture:
	soil texture triangle - https://upload.wikimedia.org/wikipedia/commons/6/65/SoilTextureTriangle.jpg
Ks comes from (in units of m/s, converted to cm/hr):
	Chen, F., & Dudhia, J. (2001). Coupling an advanced land surface–hydrology model with the Penn State–NCAR MM5 modeling system. Part I: Model implementation and sensitivity. Monthly Weather Review, 129(4), 569-585.
Organic Materials was given same properties as loam - they have same Ks
Bedrock was given same properties as clay - they had most similar Ks
Other was given same properties as Silty Clay - they have same Ks
Por connectivity (l) is set to 0.5 for all soils - HYDRUS documentation says it is 0.5 for many soils and I found no source to dispute or verify this

Make_Hydrus_Spinup_Files.r
Creates hydrus spinup file for the 9km grid

Create_Final_Run_Files.r
Views output from Hydrus spinup runs, determines whether the Hydrus simulation ran to completion at this point (good_point) or not (bad_point). Then this script creates either new final run files or new spinup files (larger maxIT) depending on if the point was identified as having a completed/fully-converged Hydrus simulation or not.

Get_Soil_Bad_points.r
Creates a new soil file (analogous to that created in Get_Soil.r) with an updated soil type for the 10,161 points that did not converge in the original attempt. The first attempt to alter the soil parameters for points that did not converge looks at the second primary soil type at the grid cell (Secondary_Altered_Points). If the grid cell only has 1 soil type defined, i.e. all other soil types are noted as 0% of the grid cell, the next closest soil type with a different soil parameter set is chosen, other is the default, but in some cases other's parameter set is the same as the original parameter set. In these cases, Silt is the chosen soil classification because it has the most similar soil parameters to 'other'.
1	Sand --> other
2	Loamy sand --> other
3	Sandy loam --> other
4	Silt loam --> other
5	Silt --> other -->other
6	Loam --> other --other
7	Sandy clay loam --> other
8	Silty clay loam --> other
9	Clay loam --> other
10	Sandy clay --> Silt
11	Silty clay --> Silt
12	Clay --> Silt
13	Organic materials --> other
14	Water 
15	Bedrock --> Silt
16	Other --> Silt

Create_Selector_Spinup_SoilMod_BadPoints.r
Create new Hydrus Selector file for the points that did not originally converge, but now have update soil parameters.

Run_Hydrus_SpinUp_Loop_9_template.sh
Template for scripts that run Hydrus spinup simulations. The template is iteratively altered by Create_Summit_Run_Hydrus_SpinUp_Loop_9.m, to create scripts that can run in parallel for different sub-regions in the ESMAP domain.

Run_Hydrus_Final_Loop_9_GoodPoints_Template.sh
Template for scripts that run Hydrus final simulations. The template is iteratively altered by Create_Summit_Final_Run_Loop_9.m, to create scripts that can run in parallel for different sub-regions in the ESMAP domain.

create_lat_lon_list_9km_final_summit_parallel_portions.m
Creates sub-regions from a large lat-lon coordinate list to run multiple smaller jobs in parallel on the supercomputer (Summit)

Get_NLDAS_Evap_SMAP_Time_Grid.r 
Creates gridded NLDAS Esoil on SMAP intervals at each of the 9km grid points, units are cm

Calculate_ESoil.r
Calculates water balance terms and then ESOIL
Outputs a file with the following information:
1 - Time, Truncated day of mid-point between observations
2 - Length of SMAP overpass interval [days]
3 - EVI, mean of EVI used in Hydrus forcing for SMAP overpass interval
4 - Precipitation [cm] summed for SMAP overpass interval, want less than 0.2 cm
5 - Delta SMAP[cm], change in SMAP for SMAP overpass interval
6 - Transpiration [cm] or SMAP overpass interval, calculate in Calc_Mu_Transpiration.r with the SM and root restriction
7 - Qbot [cm] bottom flux calculated in Hydrus
8 - Esoil [cm], calculated Esoil from water balance calculation

Define_Soil_Param_Flag.m
Outputs the coordinate list of the points that converged using the 'Secondary_params' and those that used the next closest parameter set to the original parameter set 'no_Secondary_Params'. Because the 'Silt' parameter set is most similar to the 'Other' parameter set, in these cases where the original parameter set was identical to other's the 'Silt' texture was chosen.

Grid_QC_ESMAP_revised.r
Creates gridded netcdf files for ESMAP outputs after screening for the SMAP qual_flag  

convert_hdf_to_netcdf.m
Converts SMAP native hdf files to netCDF files.
___________________________________________________________________________________________________________________________________________________________________________________________________________________________________

The below scripts use CDO and NCO to interpolate and aggregate data in netCDF format.

NLDAS2.Mosaic.remapbil.pl
Uses CDO bilinear interpolation to bring NLDAS2 Mosaic Data to the ESMAP 9km EASE Grid

NLDAS2.Noah.remapbil.pl
Uses CDO bilinear interpolation to bring NLDAS2 Noah Data to the ESMAP 9km EASE Grid

GLEAM.Esoil.remapbil.v3.1a.pl
Uses CDO bilinear interpolation to bring GLEAM Esoil Data to the ESMAP 9km EASE Grid

GLEAM.ET.remapbil.v3.1a.pl
Uses CDO bilinear interpolation to bring GLEAM ET Data to the ESMAP 9km EASE Grid

command.line.loop.mosaic.[5-8]hr.offset.pl
Uses NCO aggregate hourly NLDAS-2 Mosaic data to daily data accounting for a 5-8 hour UTC to local time difference.

command.line.loop.noah.[5-8]hr.offset.pl
Uses NCO aggregate hourly NLDAS-2 Noah data to daily data accounting for a 5-8 hour UTC to local time difference.

___________________________________________________________________________________________________________________________________________________________________________________________________________________________________

The below scripts were used as a part of the analysis and validation of the ESMAP data product

Screen_Mosaic_Esoil_Daily.m
Screens NLDAS-2 Mosaic's temporally continuous Esoil outputs to match the ESMAP Esoil gridded product's 'valid' intervals

Screen_Noah_Esoil_Daily.m
Screens NLDAS-2 Noah's temporally continuous Esoil outputs to match the ESMAP Esoil gridded product's 'valid' intervals

Screen_GLEAM_Esoil.m
Screens GLEAM's temporally continuous Esoil estimates to match the ESMAP Esoil gridded product's 'valid' intervals

Convert_Lv_of_T.m 
Converts Esoil, Ec and Trans NLDAS outputs from W/m2 to mm/day using a latent heat of vaporization formula from equation 4 from Henderson-Sellers, 1986 where Lv is a function of temp

PaperPlot_Screening.m 
Creates Figure 2a in the ESMAP data product's data descriptor manuscript

PaperPlot_ESMAP_Slength.m
Creates Figure 2b in ESMAP data product's data descriptor manuscript

PaperPlot_ESMAP_Timing.m
Creates Figure 2c in ESMAP data product's data descriptor manuscript

PaperPlot_KDE_of_ESMAP_Components.m
Plots all panels in Figure 3 in the ESMAP data product's data descriptor manuscript

Calculate_GLEAM_Esoil_ET_ratio.m 
Calculates the ratio of the mean Esoil/ET for screened/unscreened time series from GLEAM. Where the screened time series temporally matches the ESMAP time series and the unscreened is the continuous GLEAM time series. PaperPlot_GLEAM_Esoil_ET_ratio.m uses these calculations to creates Figure 4a in the ESMAP data product's data descriptor manuscript

Calculate_Mosaic_Esoil_ET_ratio.m 
Calculates the ratio of the mean Esoil/ET for screened/unscreened time series from NLDAS-2 Mosaic. Where the screened time series temporally matches the ESMAP time series and the unscreened is the continuous Mosaic time series. PaperPlot_Mosaic_Esoil_ET_ratio.m uses these calculations to creates Figure 4b in the ESMAP data product's data descriptor manuscript	

Calculate_Noah_Esoil_ET_ratio.m 
Calculates the ratio of the mean Esoil/ET for screened/unscreened time series from NLDAS-2 Noah. Where the screened time series temporally matches the ESMAP time series and the unscreened is the continuous Noah time series. PaperPlot_Noah_Esoil_ET_ratio.m uses these calculations to creates Figure 4c in the ESMAP data product's data descriptor manuscript	

PaperPlot_Boxplot_Esoil_ET_ratio.m
Creates Figure 4d in the ESMAP data product's data descriptor manuscript

Fig4_Sginificance_Testing.R
Employs a Paired Wilcoxon Rank Sum test, to test the significance of the difference in medians between Screened and unscreened evaluation data sets. Uses wilcox.exact package in R from exactRankTests

PaperPlot_ESMAP_Difference_Maps.m	
Creates all panels in Figure 5 in the ESMAP data product's data descriptor manuscript

Create_SCALE_FACTOR_files.m
Creates the corresponding scale factor file in the data repository, based on Fig. 4 of the data descriptor (Abolafia-Rosenzweig et al., 2020, in review)
