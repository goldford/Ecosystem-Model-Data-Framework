#!/bin/ksh
source $LMOD_PKG/init/ksh
module load cdo/1.9.8 ;

pathin="../DATA/light_daily/"
pathout="../DATA/light_monthly/"
n=0
N1=0

for y in {1979..2020} ; do
  echo y: $y ;
  cdo monmean ${pathin}ERA5_NEMOgrid_light_daily_${y}.nc ${pathout}ERA5_NEMOgrid_light_monthly_${y}.nc;
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

