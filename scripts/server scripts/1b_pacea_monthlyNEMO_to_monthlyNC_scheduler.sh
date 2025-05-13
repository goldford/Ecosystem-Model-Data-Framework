#!/bin/bash
#SBATCH --job-name="nemo to mon"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=02:00:00 # 
#SBATCH --ntasks=3 # CPUs

python 1b_pacea_monthlyNEMO_to_monthlyNC.py