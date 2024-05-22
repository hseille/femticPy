#!/bin/bash
#SBATCH --account=OD-231390
#SBATCH --job-name=meshGen
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=2g

project=synthetic

bash meshGen.sh $project

