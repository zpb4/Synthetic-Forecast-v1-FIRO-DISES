#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=synforc
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --output=synforc.txt
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --hint=nomultithread

module load R/4.1.2
module load gsl/2.6
module load openmpi4/4.0.5

mpirun -quiet Rscript ./src/create_synthetic_forecasts.R




