
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=synforc${SEED}\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --exclusive\n\
#SBATCH --mem-per-cpu=4G\n\
#SBATCH --output=synforc${SEED}.txt\n\
#SBATCH --ntasks-per-node=40\n\
#SBATCH --hint=nomultithread\n\

module load R\n\
module load gsl\n\
module load openmpi4\n\

mpirun Rscript ./src/create_synthetic_forecasts.R"

echo -e $SLURM | sbatch 
sleep 0.5


