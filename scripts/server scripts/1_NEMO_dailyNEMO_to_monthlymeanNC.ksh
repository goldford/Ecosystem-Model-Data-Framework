#!/bin/ksh
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;

# by G Oldford 2021-2024
# processes NEMO NC to monthly means in monthly and in annual NC files
# the monthly files are converted by another script to ASC
# the annual files are used in calculations of PAR and then converted to ASC (annual files faster than month)
# command prompt example: 'ksh 1_NEMO_to_monthlymean_asannualNC_NEMO.ksh'
# WARNING: script doesn't purge temp files

# note - monmin monmax monstd monvar should also work
# https://code.mpimet.mpg.de/projects/cdo/embedded/index.html

# daily files - converts to mean monthly data in annual files
modelrun="216"
yearsindirs="F"
mos={"01","02","03","04","03","04","05","06","07","08","09","10","11","12"} 
yrs={1979..2018}

pathin1="/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN"
echo $pathin ;

pathtemp="../DATA/TEMP/"
pathout1="../DATA/SS1500-RUN"
pathout2="/NEMO_annual_NC/"
pathout3="/NEMO_monthly_NC/"
pathout_mo=$pathout1""$modelrun""$pathout3
pathout_yr=$pathout1""$modelrun""$pathout2

echo $pathout_mo ;
dir="${pathout_mo}" ;
if [[ ! -e $dir ]]; then
  mkdir $dir
fi

echo $pathout_yr ;
dir="${pathout_yr}" ;
if [[ ! -e $dir ]]; then
  mkdir $dir
fi


#T  centre of cube, "V" middle side "north" (y vel), "U" middle side "east-west" (x vel), "W" centre of upper  / lower centre
declare -a grids2D=("T")
declare -a grids3D=("T" "V" "U" "W")
declare -a grids3D_to_2D=("V" "U" "T")
declare -a grids3D_to_2D_W=("W") # different levels for upwelling

# /////////////////////////////////////////////////////////////////////
# 2D files
# calc monthly mean of already daily averaged fields stored in monthly NC files
echo 2D files to monmean ;
for y in $yrs ; do
  echo y: $y ;
  if [[ $yearsindirs = "T" ]];
  then
	  pathin2="/CDF/"$y"/"
  else
 	  pathin2="/CDF/"
  fi
  pathin=$pathin1""$modelrun""$pathin2
  for m in $mos ; do 
    echo m: $m ;
    for grid in "${grids2D[@]}" ; do
      echo ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_${grid}_2D_y${y}m${m}.nc ;
      cdo monmean ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_${grid}_2D_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_2D_y${y}m${m}.nc;
    done  
  done
  for grid in "${grids2D[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_2D_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_2D_y${y}.nc ;
  done
done

# ///////////////// MONTH MEAN ALL LEVS 3D ///////////////////////////
# 3D files, monthly means, monhtly and annual NC files, all fields, avg'd
echo 3D fields to mon mean ;
for y in $yrs ; do
  echo y: $y ;
  if [[ $yearsindirs = "T" ]];
  then
	  pathin2="/CDF/"$y"/"
  else
 	  pathin2="/CDF/"
  fi
  
  pathin=$pathin1""$modelrun""$pathin2
  
  for m in $mos ; do 
    for grid in "${grids3D[@]}" ; do
      echo ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_${grid}_y${y}m${m}.nc
      cdo monmean ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_${grid}_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ;
    done
  done
  for grid in "${grids3D[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}.nc ;
  done
done

# this section not done, just currently a copy of 10 m vertmean
# ///////////////// SELECTED LAYERS / DEPTHS //////////////////////
# mean 3D files to single lev
#echo extracting single lev from avgd monthly ;
#for y in $yrs ; do
#  echo y: $y ;
#   not complete - left if here 2022-08-26
#  for m in $mos ; do 
#    for grid in "${grids3D_to_2D[@]}" ; do
    
      # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
      # if 'levels' aren't integers then I have to type out the level values exactly, as below 
      # https://code.mpimet.mpg.de/boards/1/topics/786
#      echo ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
#      cdo sellevel,0.500000298,1.5000031,2.50001144,3.50003052,4.50007057,5.50015068,6.50031042,7.50062323,8.50123596,9.50243282,10.5047655 ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ; 
      
#      cdo vertmean ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}m${m}.nc ;
      
#    done
#  done
#  for grid in "${grids3D_to_2D[@]}" ; do
#    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}.nc ;
#  done                                                                                                                    
#done


# ///////////////// TOP 10 layers (~metres) //////////////////////
# mean 3D files to mean 2D
echo 3D to mean 2D grid top 10 lyrs;
for y in $yrs ; do
  echo y: $y ;
  for m in $mos ; do 
    for grid in "${grids3D_to_2D[@]}" ; do
    
      # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
      # if 'levels' aren't integers then I have to type out the level values exactly, as below 
      # https://code.mpimet.mpg.de/boards/1/topics/786
      echo ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
      cdo sellevel,0.500000298,1.5000031,2.50001144,3.50003052,4.50007057,5.50015068,6.50031042,7.50062323,8.50123596,9.50243282,10.5047655 ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ; 
      
      cdo vertmean ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}m${m}.nc ;
      
    done
  done
  for grid in "${grids3D_to_2D[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}.nc ;
  done                                                                                                                    
done

# ///////////////// full water col vertmean //////////////////////
# mean 3D files to mean 2D
echo 3D to mean 2D grid vertmean all depths;
for y in $yrs ; do
  echo y: $y ;
  for m in $mos ; do 
    for grid in "${grids3D_to_2D[@]}" ; do
    

      echo ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
      
      cdo vertmean ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_vertmean_y${y}m${m}.nc ;
      
    done
  done
  for grid in "${grids3D_to_2D[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_vertmean_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_vertmean_y${y}.nc ;
  done                                                                                                                    
done

# /////////////////////////////////////////////////////////////////////
# mean 3D files to mean 2D W (upwelling W grid has different levels)
echo 3D to mean 2D W grid - top 10 lyrs; 
for y in $yrs ; do
  echo y: $y ;
  for m in $mos ; do 
    for grid in "${grids3D_to_2D_W[@]}" ; do
    
      # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
      # if 'levels' aren't integers then I have to type out the level values exactly, as below 
      # https://code.mpimet.mpg.de/boards/1/topics/786
      echo ${pathout_mo}SalishSea1500-${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
      cdo sellevel,0,1.000001,2.000006,3.000019,4.000047,5.000104,6.000217,7.000441,8.000879,9.001736,10.00341 ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ; 
      
      cdo vertmean ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top10m_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}m${m}.nc ;
      
    done
  done
  for grid in "${grids3D_to_2D_W[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop10m_y${y}.nc ;
  done
done

# ///////////////// TOP 4 layers (~metres) //////////////////////
# mean 3D files to mean 2D 
echo 3D to mean 2D grid top 4 lyrs;
for y in $yrs ; do
  echo y: $y ;
  for m in $mos ; do 
    for grid in "${grids3D_to_2D[@]}" ; do
    
      # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
      # if 'levels' aren't integers then I have to type out the level values exactly, as below 
      # https://code.mpimet.mpg.de/boards/1/topics/786
      echo ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
      cdo sellevel,0.500000298,1.5000031,2.50001144,3.50003052 ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top4m_y${y}m${m}.nc ; 
      
      cdo vertmean ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top4m_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}m${m}.nc ;
      
    done
  done
  for grid in "${grids3D_to_2D[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}.nc ;
  done                                                                                                                    
done

# /////////////////////////////////////////////////////////////////////
# mean 3D files to mean 2D W (different levels)
echo 3D to mean 2D W grid top 4 lyrs; 
for y in $yrs ; do
  echo y: $y ;
  for m in $mos ; do 
    for grid in "${grids3D_to_2D_W[@]}" ; do
    
      # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
      # if 'levels' aren't integers then I have to type out the level values exactly, as below 
      # https://code.mpimet.mpg.de/boards/1/topics/786
      echo ${pathout_mo}SalishSea1500-${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc
      cdo sellevel,0,1.000001,2.000006,3.000019 ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_y${y}m${m}.nc ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top4m_y${y}m${m}.nc ; 
      
      cdo vertmean ${pathtemp}RUN${modelrun}_MonthlyMean_grid_${grid}_top4m_y${y}m${m}.nc ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}m${m}.nc ;
      
    done
  done
  for grid in "${grids3D_to_2D_W[@]}" ; do
    ncrcat ${pathout_mo}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}* ${pathout_yr}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_${grid}_avgtop4m_y${y}.nc ;
  done
done





# OLD

#for y in {1979..1979} ; do
#  for file in ${pathin}/SalishSea1500-${modelrun}_1d_grid_${grid}_y${y}*.nc ; do 
#
#    echo found file $file ;
#
#    # this should be the monthandday
#    #echo "${file:113:4}" ;
#    #date1="${file:113:4}" ;
#    
#    # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
#    # if 'levels' aren't integers then I have to type out the level values exactly, as below 
#    # https://code.mpimet.mpg.de/boards/1/topics/786
#    cdo sellevel,0.500000298,1.5000031,2.50001144,3.50003052,4.50007057,5.50015068,6.50031042,7.50062323,8.50123596,9.50243282,10.5047655 ${file} ${pathtemp}/${modelrun}_1d_grid_${grid}_selectedlevs_${y}${date1}.nc ; 
#    
#    cdo vertmean ${pathtemp}${modelrun}_1d_grid_${grid}_selectedlevs_${y}${date1}.nc ${pathtemp}${modelrun}_1d_grid_${grid}_vertmean_${y}${date1}.nc ;
#    #cdo vertmean ${file} ${pathout}/{$NEMOrun}_1d_grid_T_${date1}.nc ;
#    #cdo vertmean ${file} ${pathout}${y}/{$NEMOrun}_1d_grid_T_${date1}.nc ; # change if files moved to annual folders
#  done
#  
#  # concatenate a year of daily or hourly files into annual files
#  #echo concatenating y: $y ; 
#  #ncrcat ${pathtemp}${modelrun}_1d_grid_${grid}_vertmean_${y}* ${pathtemp}${modelrun}_1d_grid_${grid}_vertmean_${y}.nc ;
#  
#  echo calculating monmean y: $y ;
#  cdo monmean ${pathtemp}${modelrun}_1d_grid_${grid}_vertmean_${y}.nc ${pathout}SalishSea1500-${modelrun}_MonthlyVertMean_grid_${grid}_selectedlevs_${y}.nc ;
#
#  cd ${pathtemp} ; 
#  rm *.nc
#  cd - ; 
#  
#done
#
#
## calc monthly mean of daily averaged fields stored in monthly NC files
#for y in {1979..1979} ; do
#  echo y: $y ;
#  for m in {"01","02","03","04","03","04","05","06","07","08","09","10","11","12"} ; do 
#    echo m: $m ; 
#    echo ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_T_y${y}m${m}.nc ;
#    cdo monmean ${pathin}SalishSea1500-RUN${modelrun}_1d_grid_T_y${y}m${m}.nc ${pathout}SalishSea1500-RUN${modelrun}_MonthlyMean_grid_T_y${y}m${m}.nc;
#done 
#done


# if issues with nan or nodata or weird values look at cdo monavg

# old code bleow
#for y in {1979..1983} ; do
#  echo y: $y ; 
#  echo N1: $N1 ;
#  if [[ $n -eq 1 ]] ; then
#    N2=$((N1 + 366 * 24 - 1)) 
#else
#    N2=$((N1 + 365 * 24 - 1))
#  fi
  
#  echo N2: $N2 ;
  
#  ncks -d time,$N1,$N2 ${pathin}meansspressure_2016to2020.nc ${pathout}meansspressure_${y}.nc ;
#  ncks --mk_rec_dmn time ${pathout}meansspressure_${y}.nc -o ${pathout}meansspressure_y${y}.nc ;

#  let n++; 
#  if [[ $n -eq 4 ]] ; then
#    n=0 ;
#  fi
#  N1=$((N2+1))
#  echo n: $n ;
