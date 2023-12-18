
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=synforc2${SEED}\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --exclusive\n\
#SBATCH --mem-per-cpu=4G\n\
#SBATCH --output=synforc2${SEED}.txt\n\
#SBATCH --ntasks-per-node=80\n\

module load R\n\
module load gsl\n\

Rscript ./src/create_synthetic_forecasts.R"

echo -e $SLURM | sbatch 
sleep 0.5


