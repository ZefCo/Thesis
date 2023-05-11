#!/bin/bash
#SBATCH --job-name=fractal_ml_05032023_1
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=10 -N 1
#SBATCH --mem=5gb
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL

python FractalModelSlurm.py


