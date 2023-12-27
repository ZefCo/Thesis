#!/bin/bash
#SBATCH --job-name=envtest_eas1
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=10 -N 1
#SBATCH --mem=5gb
#SBATCH --mail-user=espeakma@cougarnet.uh.edu
#SBATCH --mail-type=ALL


module load Anaconda3/2023.03-1
module load TensorFlow/2.11.0-foss-2022a
python EnvTest.py

