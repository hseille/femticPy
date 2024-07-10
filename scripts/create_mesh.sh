#!/bin/bash
#SBATCH --job-name=meshGen
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=2g

project=auslamp

bash meshGen.sh $project

