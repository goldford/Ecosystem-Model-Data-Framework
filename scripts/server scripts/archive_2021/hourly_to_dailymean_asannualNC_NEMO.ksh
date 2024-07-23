#!/bin/ksh
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;
module load nco/4.9.5 ;

# Takes hourly data in daily files and converts to mean monthly data in annual files

pathin="/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN201/CDF/"
pathtemp="/project/6006412/goldford/data_temp/"
pathout="../DATA/NEMOresults_hourlytodailymeans_NC/"

# concatenate a year of daily files into annual files
for y in {1979..1982} ; do 
  echo concatenating y: $y ; 
  ncrcat ${pathin}SalishSea1500-RUN201_1h_grid_T_2D_${y}* ${pathtemp}SalishSea1500-RUN201_1h_grid_T_2D_${y}.nc
done 

# monthly mean 
for y in {1979..1982} ; do
  echo y: $y ;
  cdo daymean ${pathtemp}SalishSea1500-RUN201_1h_grid_T_2D_${y}.nc ${pathout}SalishSea1500-RUN201_DailyMean_grid_T_2D_${y}.nc;
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
done 

