
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=datapro\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --nodes=1\n\
#SBATCH --exclusive\n\
#SBATCH --mem-per-cpu=4G\n\
#SBATCH --output=datapro.txt\n\
#SBATCH --ntasks-per-node=80\n\

module load R\n\
module load gsl\n\

Rscript ./src/data_processing.R"

echo -e $SLURM | sbatch 
sleep 0.5


