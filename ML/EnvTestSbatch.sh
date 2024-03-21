#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1 -N 1
#SBATCH --mem=1gb
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL
#SBATCH -J EASEnvTest
#SBATCH -o EnvTest.out


module load Anaconda3/2023.03-1
module load TensorFlow/2.11.0-foss-2022a
python EnvTest.py

