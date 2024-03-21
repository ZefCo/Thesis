#!/bin/bash
#SBATCH -J fractal_GHMSS_v1
#SBATCH -o fractal_GHMSS_v1.out
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=8 -N 1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL


module load Anaconda3/2023.03-1
module load TensorFlow/2.11.0-foss-2022a
python FractalModelSlurm.py

