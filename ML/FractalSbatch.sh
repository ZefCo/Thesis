#!/bin/bash
#SBATCH -J fractal_6x6v3
#SBATCH -o fractal_6x6v3.out
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=10 -N 1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL


module load Anaconda3/2023.03-1
module load TensorFlow/2.11.0-foss-2022a
python FractalModelSlurm.py

