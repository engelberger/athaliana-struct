#!/bin/bash
#SBATCH --job-name=ESM--test
#SBATCH --output=ESM--test-job-out.%J
#SBATCH --error=example.err.%J
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --gres=gpu:v100:1
#SBATCH --partition=clara

~/micromamba/envs/esmfold/bin/python esmfold.py all.fasta
