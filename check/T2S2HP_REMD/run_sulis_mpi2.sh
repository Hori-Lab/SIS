#!/bin/bash
#SBATCH --account=su006-044
#SBATCH --job-name=mpitest
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=devel
#SBATCH --mem-per-cpu=500
#SBATCH --time=0:05:00
#SBATCH -o out.slurm
#SBATCH -e err.slurm

#SBATCH --mail-user=naoto.hori@nottingham.ac.uk
#SBATCH --mail-type=ALL

module purge
module load GCCcore/11.3.0 intel iimpi

srun mpirun2.sh 1> out 2> err
