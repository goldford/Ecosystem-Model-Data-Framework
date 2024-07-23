#!/bin/ksh
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;

# GO 2022-01-13 - 
# Processes NEMO output daily 3D files on T grid (e.g. vosaline, votemper) and
# converts to 2D using cdo **for a single selected level (depth) only *
# then it concatenates them together into annual NC files

pathin="/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN202/CDF/"
pathout="../DATA/SS1500-RUN202/NEMO_annualNC/"
pathtemp="../DATA/TEMP/" # for intermed files

NEMOrun="RUN202"

for y in {1979..2018} ; do
  for file in ${pathin}/SalishSea1500-${NEMOrun}_1d_grid_T_${y}*.nc ; do # change if files moved to annual folders
  #for file in ${pathin}${y}/SalishSea1500-{$NEMOrun}_1d_grid_W_*.nc ; do
    echo found file $file ;

    dir="${pathout}" ;

    if [[ ! -e $dir ]]; then
      mkdir $dir
    fi

    # this should be the monthandday
    echo "${file:113:4}" ;
    date1="${file:113:4}" ;
    
    # GO-2021-12-22 two-step process (1/ get only top x layers and 2/ vertmean) 
    # when 'levels' aren't integers then I have to type out the level values exactly, as below 
    # https://code.mpimet.mpg.de/boards/1/topics/786
    # cdo showlevel temp.nc
    #  0.500000298 1.5000031 2.50001144 3.50003052 4.50007057 5.50015068 6.50031042 7.50062323 8.50123596 9.50243282 10.5047655 11.5093117 12.5181665 13.5354118 14.5689821 15.6342878 16.7611732 18.0071354 19.4817848 21.3899784 24.100256 28.2299156 34.6857567 44.5177231 58.484333 76.5855865 98.0629578 121.866516 147.089462 173.114487 199.573044 226.2603 253.066635 279.93454 306.834198 333.750183 360.67453 387.60321 414.534088 441.466095
    
    cdo sellevel,1.5000031 ${file} ${pathtemp}/${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ;
    cdo vertmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ${pathtemp}${NEMOrun}_1d_grid_T_singledepth1m_${y}${date1}.nc ;
    
    #cdo sellevel,4.50007057 ${file} ${pathtemp}/${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ;
    #cdo vertmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ${pathtemp}${NEMOrun}_1d_grid_T_singledepth4m_${y}${date1}.nc ;
    
    #cdo sellevel,9.50243282 ${file} ${pathtemp}/${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ;
    #cdo vertmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth_${y}${date1}.nc ${pathtemp}${NEMOrun}_1d_grid_T_singledepth10m_${y}${date1}.nc ;
  done
  
  # concatenate a year of daily or hourly files into annual files
  echo concatenating y: $y ; 
  
  ncrcat ${pathtemp}${NEMOrun}_1d_grid_T_singledepth1m_${y}* ${pathtemp}${NEMOrun}_1d_grid_T_singledepth1m_${y}.nc ;
  #ncrcat ${pathtemp}${NEMOrun}_1d_grid_T_singledepth4m_${y}* ${pathtemp}${NEMOrun}_1d_grid_T_singledepth4m_${y}.nc ;
  #ncrcat ${pathtemp}${NEMOrun}_1d_grid_T_singledepth10m_${y}* ${pathtemp}${NEMOrun}_1d_grid_T_singledepth10m_${y}.nc ;
  
  echo calculating monmean y: $y ;
  
  cdo monmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth1m_${y}.nc ${pathout}SalishSea1500-${NEMOrun}_MonthlyMean_grid_T_singledepth1m_${y}.nc ;
  #cdo monmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth4m_${y}.nc ${pathout}SalishSea1500-${NEMOrun}_MonthlyMean_grid_T_singledepth4m_${y}.nc ;
  #cdo monmean ${pathtemp}${NEMOrun}_1d_grid_T_singledepth10m_${y}.nc ${pathout}SalishSea1500-${NEMOrun}_MonthlyMean_grid_T_singledepth10m_${y}.nc ;

  cd ${pathtemp} ; 
  rm *.nc
  cd - ; 
  
done


#for y in {1979..1980} ; do
#  echo y: $y ;
#  for m in {1..12} ; do
#    if [[ $m -ge 10 ]] ; then
#      echo m: $m ;
#    else
#      echo m: 0$m ;
#    fi
#  done
#done
  #cdo vertmean ${pathin}ERA5_NEMOgrid_light_${y}.nc ${pathout}ERA5_NEMOgrid_light_daily_${y}.nc;
  #cdo vertmean ~/projects/def-nereusvc/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN200/CDF/1979/SalishSea1500-RUN200_1d_grid_W_19790420-19790420.nc ../DATA/TEMP/SalishSea1500-RUN200_1d_grid_W_19790420-19790420_vertavg.nc

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
#done 

