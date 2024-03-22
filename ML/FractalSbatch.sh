#!/bin/bash
#SBATCH -J fractal_GHMSS_v1
#SBATCH -o fractal_GHMSS_v1.out
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=8 -N 1
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL


module load TensorFlow/2.11.0-foss-2022a
python FractalModelSlurm.py
module unload TensorFlow/2.11.0-foss-2022a

module load Anaconda3
python PrintPlots.py fractal_GHMSS_v1.out
module unload Anaconda3

