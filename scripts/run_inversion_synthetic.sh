#!/bin/bash
#SBATCH --account=OD-231390
#SBATCH --job-name=fem_01
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=50g

# Application specific commands:
module load intel-fc
module load intel-cc
module load intel-mkl/2020.1.217
module load openmpi/4.1.1-ofed51-intel20

project=synthetic
cd $project/inversion

mpirun -n 8 /scratch3/sei029/femtic/femtic