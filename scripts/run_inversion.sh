#!/bin/bash
#SBATCH --job-name=inv_01
#SBATCH --time=144:00:00
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=4
#SBATCH --nodes=2
#SBATCH --mem=900g

# Application specific commands:
module load intel-fc
module load intel-cc
module load intel-mkl/2020.1.217
module load openmpi/4.1.1-ofed51-intel20

project=auslamp

cd $project/inversion/inv_01

mpirun -n 28 /scratch3/sei029/femtic/femtic
