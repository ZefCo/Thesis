#!/bin/bash
#SBATCH -J EASEnvTest
#SBATCH -o EnvTest.out
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1 -N 1
#SBATCH --mem=1gb
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL


module load Anaconda3/2023.03-1
module load TensorFlow/2.11.0-foss-2022a
module avail
python EnvTest.py

# For some reason the issue is related to Tensorflow. Somehow it's not properly loading, and is causing a conflict or something to prevent the import of matplotlib and PIL at best, more at worst.
# Going to 2.13 seemed to help but still didn't allow for the import of matplotlib or PIL.
