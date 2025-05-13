#!/bin/bash
#SBATCH --job-name="nemo to daily asc"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=02:00:00 # 
#SBATCH --ntasks=3 # CPUs

python 6_PAR3_VarK_VarZ_day.py