#!/bin/sh
#SBATCH --job-name=psd-Test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.rutkowski@chem.ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=16
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=800mb
#SBATCH --qos=colina
#SBATCH --account=colina
#SBATCH --time=50:00:00
#SBATCH --output=mpi_test_%j.out

cd $SLURM_SUBMIT_DIR

date

module load ufrc intel

export OMP_NUM_THREADS=16
 
./psd > log.out 2>&1

date
