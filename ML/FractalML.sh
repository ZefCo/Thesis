#!/bin/bash
#SBATCH --job-name=job
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-node=20
#SBATCH --mem=30gb


conda activate thesis

python FractalModel.py


